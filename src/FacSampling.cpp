#include "FacSampling.h"
#include "utility-function.h"

FacSampling::FacSampling(int N, int K,int Kg, Rcpp::NumericVector cor_fac, bool std, double w0)
  : N(N), K(K),Kg(Kg), cor_fac(cor_fac), std(std), w0(w0) {

  W0 = Eigen::MatrixXd::Constant(K, K, w0);
  W0.diagonal().array() = 1.0;

  C = Eigen::MatrixXd::Identity(K, K);
  phi = inv_phi = C;

  C_gen = Eigen::MatrixXd::Identity(Kg, Kg);
  C_spe = Eigen::MatrixXd::Identity(K-Kg, K-Kg);
  phi_par = inv_phi_par = Eigen::MatrixXd::Identity(K-Kg, K-Kg);

  if(cor_fac.size() == 1) {
    cor_fac = Rcpp::NumericVector(K, cor_fac[0]);
  }

  if (cor_fac.size() != K) {
    Rcpp::stop("cor_fac should be either scalar 0/1, True/False, or an array of 0/1 or True/False with length K.");
  }

  rf = Eigen::VectorXi(K);
  for (int i = 0; i < K; ++i) {
    rf[i] = cor_fac[i];
  }
  K_rf = rf.sum();
  sub_W0 = Eigen::MatrixXd::Identity(K_rf, K_rf);
  int rowIndex = 0;
  for(int i = 0; i < rf.size(); ++i) {
    if(rf[i]) {
      int colIndex = 0;
      for(int j = 0; j < rf.size(); ++j) {
        if(rf[j]) {
          sub_W0(rowIndex, colIndex++) = W0(i, j);
        }
      }
      rowIndex++;
    }
  }

  // initial for ssp
  //Eigen::MatrixXd adj = (abs(as<Eigen::Map<Eigen::MatrixXd>>(phi)) >= 0).cast<double>();
  Eigen::MatrixXd adj = (phi.array().abs() >= 0).cast<double>();
  //Eigen::MatrixXd::Constant(K, K, 1);
  // Calculate pii
  pii = (K >= 4) ? (2.0 / (K - 1)) : 0.3; // Conditional assignment
  pii_par = ((K-Kg) >= 4) ? (2.0 / ((K-Kg) - 1)) : 0.3; // Conditional assignment
  // Create matrices V0 and V1
  lambda = 1;
  double v0tmp = 0.02*0.02;
  double h = 50*50;
  double v1tmp = h * v0tmp;
  V0 = Eigen::MatrixXd::Constant(K, K, v0tmp); // K x K matrix filled with v0
  V1 = Eigen::MatrixXd::Constant(K, K, v1tmp); // K x K matrix filled with v1
  // Initialize taus with V0
  taus = V0;
  // Update taus based on the adjacency matrix
  for (int row = 0; row < K; ++row) {
    for (int col = 0; col < K; ++col) {
      if (adj(row, col) > 0) { // If the adjacency condition is met
        taus(row, col) = v1tmp;   // Set to v1 where adj is TRUE
      }
    }
  }

  taus_par= taus.block(Kg, Kg, K-Kg, K-Kg);
  V0_par= V0.block(Kg, Kg, K-Kg, K-Kg);
  V1_par= V1.block(Kg, Kg, K-Kg, K-Kg);


  //init for horse
  Lambda_sq = Eigen::MatrixXd::Constant(K, K, 1);
  Nu = Eigen::MatrixXd::Constant(K, K, 1);
  xi=1.0;
  regdiag= true; //also for ssp and lasso


  mip=Eigen::MatrixXd::Constant(K, K, 0);
  }

Eigen::MatrixXd FacSampling::gibbs_factor_LD(Eigen::MatrixXd Y, Eigen::MatrixXd A, Eigen::MatrixXd inv_psx) {  // Eigen::VectorXd dV
  Eigen::MatrixXd inv_C = C.inverse(); // 'inv_sympd' can replace 'inv' for symmetric positive definite matrices, which is more efficient and numerically stable.
  Eigen::MatrixXd Sn = (A.transpose() * inv_psx * A + inv_C).inverse(); //the posterior covariance matrix of the latent factors
  Eigen::MatrixXd mean = (Y * inv_psx * A) * Sn;
  Eigen::MatrixXd CSn = Eigen::LLT<Eigen::MatrixXd>(Sn).matrixL(); // Compute the Cholesky decomposition of Sn, Retrieve the lower triangular matrix L

  arma::mat Random_norm = arma::randn(N, K);
  Eigen::MatrixXd Random = Eigen::Map<Eigen::MatrixXd>(Random_norm.memptr(), Random_norm.n_rows, Random_norm.n_cols);
  Eigen::MatrixXd fac = Random * CSn.transpose() + mean;

  return fac;
}

Eigen::MatrixXd FacSampling::cor_cov(Eigen::MatrixXd fac, bool std) {
  Eigen::MatrixXd scaleInv = (fac.transpose() * fac + sub_W0).inverse();
  Eigen::MatrixXd inv_cov=rwish1( K + 2 + N,scaleInv);
  Eigen::MatrixXd new_C = inv_cov.inverse();
  if (std) {
    //////new method to convert correlation
    //Eigen::VectorXd std_dev = new_C.diagonal().array().sqrt(); // Standard deviations
    //Eigen::MatrixXd C1 = Eigen::MatrixXd::Zero(new_C.rows(), new_C.cols());
    //for (int i = 0; i < new_C.rows(); ++i) {
    //  for (int j = 0; j < new_C.cols(); ++j) {
    //    C1(i, j) = new_C(i, j) / (std_dev(i) * std_dev(j));
    //  }
    //}
    /////translate from python
    //Rcpp::Rcout  << "cor_cov--The new_C is: "<< new_C << std::endl;
    Eigen::MatrixXd diag = new_C.diagonal().array().sqrt().matrix().asDiagonal();
    //Eigen::MatrixXd cho = Eigen::LLT<Eigen::MatrixXd>(diag).matrixL();
    //Eigen::MatrixXd tmp = cho.inverse();
    Eigen::MatrixXd tmp = diag.inverse();
    Eigen::MatrixXd C1 = tmp * new_C * tmp;
    double acc = exp((K + 1) / 2 * (log(C1.determinant()) - log(C.determinant())));
    //////////std uniform sampling
    //std::mt19937 gen(seed+iter); // Seed the generator
    //std::uniform_real_distribution<> dis_unif(0.0, 1.0); // Uniform distribution in [0, 1)
    //double uniform_sample = dis_unif(gen);

    double uniform_sample = arma::randu();

    bool acid = acc > uniform_sample;
    new_C = C1 * acid + C * (!acid);

  }
  C = new_C;
  return C;
}
Eigen::MatrixXd FacSampling::cor_cov1(Eigen::MatrixXd fac, bool std) {

  if (K_rf > 1) {
    //Eigen::MatrixXd fac1 = submatrix(fac,Eigen::VectorXi::Ones(1),rf);
    Eigen::MatrixXd fac1(fac.rows(), K_rf);
    int colIndex = 0;
    for(int i = 0; i < rf.size(); ++i) {
      if(rf[i]) {
        fac1.col(colIndex++) = fac.col(i);
      }
    }
    Eigen::MatrixXd scaleInv = (fac1.transpose() * fac1 + sub_W0).inverse();
    Eigen::MatrixXd inv_cov=rwish1( K_rf + 2 + N,scaleInv);
    Eigen::MatrixXd new_C = inv_cov.inverse();

    if (std) {

      //Eigen::MatrixXd C0 = submatrix(C,rf,rf);
      Eigen::MatrixXd C0(K_rf, K_rf);
      int rowIndex = 0;
      for(int i = 0; i < rf.size(); ++i) {
        if(rf[i]) {
          int colIndex = 0;
          for(int j = 0; j < rf.size(); ++j) {
            if(rf[j]) {
              C0(rowIndex, colIndex++) = C(i, j);
            }}
          rowIndex++;}}

      Eigen::MatrixXd diag = new_C.diagonal().array().sqrt().matrix().asDiagonal();
      //Eigen::MatrixXd cho = diag.llt().matrixL();
      //Eigen::MatrixXd tmp = cho.inverse();

      Eigen::MatrixXd tmp = diag.inverse(); //LAWBL
      Eigen::MatrixXd C1 = tmp * new_C * tmp;

      //Eigen::VectorXd std_dev = new_C.diagonal().array().sqrt(); // Standard deviations
      //Eigen::MatrixXd C1 = Eigen::MatrixXd::Zero(new_C.rows(), new_C.cols());
      //for (int i = 0; i < new_C.rows(); ++i) {
      //  for (int j = 0; j < new_C.cols(); ++j) {
      //    C1(i, j) = new_C(i, j) / (std_dev(i) * std_dev(j));
      //  }
      //}

      double acc = exp((K_rf + 1.0) / 2.0 * (log(C1.determinant()) - log(C0.determinant())));
      double uniform_sample = arma::randu();
      bool acid = acc > uniform_sample;
      new_C = C1 * acid + C0 * (!acid);
    }else {
      // Handling the case when std is false and updating self.C with covariance for uncorrelated factors
      Eigen::MatrixXd fac_uc(fac.rows(), fac.cols() - K_rf);
      int colIndex = 0;
      for(int i = 0; i < rf.size(); ++i) {
        if(!rf[i]) {  // Selecting columns based on the negation of rf
          fac_uc.col(colIndex++) = fac.col(i);
        }
      }

      // Calculate covariance for the uncorrelated factors
      // Centering the matrix by subtracting the column means
      Eigen::MatrixXd centered = fac_uc.rowwise() - fac_uc.colwise().mean();
      int n = fac_uc.rows() - 1;
      Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(n);
      colIndex = 0;
      for(int i = 0, updateCol = 0; i < rf.size(); ++i) {
        if(!rf[i]) {
          for(int j = 0, updateRow = 0; j < rf.size(); ++j) {
            if(!rf[j]) {
              C(j, i) = cov(updateRow++, updateCol);
            }
          }
          updateCol++;
        }
      }
    }

    int rowIndex = 0;
    for(int i = 0; i < rf.size(); ++i) {
      if(rf[i]) {
        int colIndex = 0;
        for(int j = 0; j < rf.size(); ++j) {
          if(rf[j]) {
            C(i,j) = new_C(rowIndex, colIndex++);
          }
        }
        rowIndex++;
      }
    }
  }
  return C;
}

Eigen::MatrixXd FacSampling::cor_cov_lasso(Eigen::MatrixXd fac, bool std,Eigen::MatrixXd reg) {

  // Generate ind_nod matrix
  Eigen::MatrixXi ind_nod(K - 1, K);
  for (int k = 0; k < K; ++k) {
    int index = 0;
    for (int i = 0; i < K; ++i) {
      if (i != k) {
        ind_nod(index++, k) = i ; //different with R (i+1), 0-based
      }
    }
  }

  double a_gams = 1.0;
  double b_gams = 0.1;

  // Calculation of S
  Eigen::MatrixXd S = fac.transpose() * fac;

  // Sample gammas
  double apost = a_gams + (K * (K + 1)) / 2.0;
  double bpost = b_gams + (inv_phi.cwiseAbs().sum()) / 2.0;

  std::mt19937 gen;
  std::gamma_distribution<> gamma(apost, 1.0 / bpost);
  double gammas  = gamma(gen);
  // Sample tau off-diagonal

  Eigen::VectorXd Cadj(K * (K - 1) / 2);
  int index = 0;
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < k; ++i) {
      Cadj(index++) = std::max(std::abs(inv_phi(i, k)), 1e-6);
    }
  }
  Eigen::VectorXd mu_p_tmp = gammas / Cadj.array();

  mu_p_tmp = mu_p_tmp.array().min(1e12); // pmin equivalent
  Rcpp::NumericVector mu_p = Rcpp::wrap(mu_p_tmp);

  double gammas_p = std::pow(gammas, 2);
  Rcpp::NumericVector taus_tmp = 1.0/ rinvgauss(mu_p_tmp.size(), mu_p, gammas_p);

  // Convert Rcpp::NumericVector to Eigen::VectorXd for further operations
  //Eigen::Map<Eigen::VectorXd> taus_tmp_eigen(taus_tmp.begin(), taus_tmp.size());

  Eigen::MatrixXd taus = Eigen::MatrixXd::Zero(K, K);
  index = 0;

  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < k; ++i) {
      // Assign taus_tmp values to both upper and lower triangular parts
      taus(i, k) = taus_tmp(index);
      taus(k, i) = taus_tmp(index);
      ++index;
    }
  }
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < k; ++i) {
      if (reg(i, k) == 0) {
        taus(i, k) = 1000; // tau large, inv_tau is small, general part is follow flat normal, not shink to zero
        taus(k, i) = 1000;
      }
    }
  }

  // Set lamdiag based on regdiag flag
  double lamdiag = regdiag ? gammas : 0.0;

  for (int i = 0; i < K; ++i) {
    Eigen::VectorXi ind_noi = ind_nod.col(i);
    // Manually constructing Sig11 as a submatrix of 'phi'
    int size = ind_noi.size();
    Eigen::MatrixXd Sig11(size, size);
    for(int row = 0; row < size; ++row) {
      for(int col = 0; col < size; ++col) {
        Sig11(row, col) = phi(ind_noi(row), ind_noi(col));
      }
    }

    // Extracting Sig12 as a subvector of column 'i' of 'phi'
    Eigen::VectorXd Sig12(size);
    for(int row = 0; row < size; ++row) {
      Sig12(row) = phi(ind_noi(row), i);
    }
    Eigen::MatrixXd invC11 = Sig11 - (Sig12 * Sig12.transpose()) / phi(i, i);
    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Zero(size, size);
    for(int idx = 0; idx < size; ++idx) {
      diagMatrix(idx, idx) = 1.0 / taus(ind_noi(idx), i); // Assuming 'i' is defined elsewhere as the column index
    }

    Eigen::MatrixXd Ci = (S(i, i) + lamdiag) * invC11 + diagMatrix;

    Eigen::MatrixXd Sigma = Ci.llt().solve(Eigen::MatrixXd::Identity(size, size)); // Assuming chol2inv(chol(Ci)) equivalent
    Eigen::VectorXd mu_i(size);
    for (int idx = 0; idx < size; ++idx) {
      mu_i(idx) = -S(ind_noi(idx), i);
    }
    mu_i = Sigma * mu_i;
    arma::vec armaMean(mu_i.data(), mu_i.size(), false, true);
    arma::mat armaCov(Sigma.data(), Sigma.rows(), Sigma.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Rcpp::NumericMatrix betaMat = Rcpp::wrap(armasamples);

    //Rcpp::NumericMatrix betaMat = rmvnorm(1, mu_i, Sigma, seed);
    Eigen::Map<Eigen::MatrixXd> beta(betaMat.begin(), betaMat.nrow(), betaMat.ncol());  //beta is one col here

    for (int idx = 0; idx < size; ++idx) {
      inv_phi(ind_noi(idx), i) = beta(idx,0);
      inv_phi(i, ind_noi(idx)) = beta(idx,0);  // Assuming symmetry
    }

    std::mt19937 gen;
    std::gamma_distribution<double> gamma(N / 2.0 + 1,2.0/(S(i, i) + lamdiag));
    double gam  = gamma(gen);

    Eigen::MatrixXd tmp = beta.transpose() * invC11 * beta;
    inv_phi(i, i) = gam + tmp(0,0);

    // Updating covariance matrix 'phi' according to one-column change in the precision matrix 'inv_phi'
    Eigen::MatrixXd invC11beta = invC11 * beta;
    Sig12 = (-invC11beta / gam).col(0); //1 column
    Eigen::MatrixXd sub_phi =   invC11 + invC11beta * invC11beta.transpose()/gam;
    for (int idx = 0; idx < size; ++idx) {
      // Update off-diagonal elements for column 'i' and symmetric positions
      phi(ind_noi(idx), i) = Sig12(idx);
      phi(i,ind_noi(idx)) = Sig12(idx);
      for (int jdx = 0; jdx < size; ++jdx) {
        // Add invC11beta contribution to both off-diagonal and diagonal elements
        phi( ind_noi(idx), ind_noi(jdx)) = sub_phi(idx,jdx);
      }
    }
    phi(i, i) = 1.0 / gam;
  }

  Eigen::MatrixXd new_C = phi;
  if (std) {
    Eigen::MatrixXd diag = new_C.diagonal().array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd tmp = diag.inverse();
    Eigen::MatrixXd C1 = tmp * new_C * tmp;
    //new_C = C1;
    double acc = exp((K + 1) / 2 * (log(C1.determinant()) - log(C.determinant())));
    double uniform_sample = arma::randu();
    bool acid = acc > uniform_sample;
    new_C = C1 * acid + C * (!acid);
  }
  C=new_C;
  return new_C;

}
Eigen::MatrixXd FacSampling::cor_cov_lasso_par(Eigen::MatrixXd fac, bool std) {

  int Ks = K - Kg;
  Eigen::MatrixXd fac_g = fac.block(0, 0, N, Kg);           // first Kg columns , user need to make the loading in this order, Kg + Ks
  Eigen::MatrixXd fac_s = fac.block(0, Kg, N, Ks);  // remaining Ks columns
  Eigen::MatrixXd sub_W0_par = sub_W0.block(0, 0, Kg, Kg);
  // general factor
  Eigen::MatrixXd scaleInv = (fac_g.transpose() * fac_g + sub_W0_par).inverse();
  Eigen::MatrixXd inv_cov=rwish1( Kg + 2 + N,scaleInv);
  Eigen::MatrixXd new_C_gen = inv_cov.inverse();



  // special factor
  // Generate ind_nod matrix
  Eigen::MatrixXi ind_nod(Ks - 1, Ks);
  for (int k = 0; k < Ks; ++k) {
    int index = 0;
    for (int i = 0; i < Ks; ++i) {
      if (i != k) {
        ind_nod(index++, k) = i ; //different with R (i+1), 0-based
      }
    }
  }

  double a_gams = 1.0;
  double b_gams = 0.1;

  // Calculation of S
  Eigen::MatrixXd S = fac_s.transpose() * fac_s;

  // Sample gammas
  double apost = a_gams + (Ks * (Ks + 1)) / 2.0;
  double bpost = b_gams + (inv_phi_par.cwiseAbs().sum()) / 2.0;

  std::mt19937 gen;
  std::gamma_distribution<> gamma(apost, 1.0 / bpost);
  double gammas  = gamma(gen);
  // Sample tau off-diagonal

  Eigen::VectorXd Cadj(Ks * (Ks - 1) / 2);
  int index = 0;
  for (int k = 0; k < Ks; ++k) {
    for (int i = 0; i < k; ++i) {
      Cadj(index++) = std::max(std::abs(inv_phi_par(i, k)), 1e-6);
    }
  }
  Eigen::VectorXd mu_p_tmp = gammas / Cadj.array();

  mu_p_tmp = mu_p_tmp.array().min(1e12); // pmin equivalent
  Rcpp::NumericVector mu_p = Rcpp::wrap(mu_p_tmp);

  double gammas_p = std::pow(gammas, 2);
  Rcpp::NumericVector taus_tmp = 1.0/ rinvgauss(mu_p_tmp.size(), mu_p, gammas_p);

  // Convert Rcpp::NumericVector to Eigen::VectorXd for further operations
  //Eigen::Map<Eigen::VectorXd> taus_tmp_eigen(taus_tmp.begin(), taus_tmp.size());

  Eigen::MatrixXd taus = Eigen::MatrixXd::Zero(Ks, Ks);
  index = 0;

  for (int k = 0; k < Ks; ++k) {
    for (int i = 0; i < k; ++i) {
      // Assign taus_tmp values to both upper and lower triangular parts
      taus(i, k) = taus_tmp(index);
      taus(k, i) = taus_tmp(index);
      ++index;
    }
  }


  for (int i = 0; i < Ks; ++i) {
    Eigen::VectorXi ind_noi = ind_nod.col(i);
    // Manually constructing Sig11 as a submatrix of 'phi_par'
    int size = ind_noi.size();
    Eigen::MatrixXd Sig11(size, size);
    for(int row = 0; row < size; ++row) {
      for(int col = 0; col < size; ++col) {
        Sig11(row, col) = phi_par(ind_noi(row), ind_noi(col));
      }
    }

    // Extracting Sig12 as a subvector of column 'i' of 'phi_par'
    Eigen::VectorXd Sig12(size);
    for(int row = 0; row < size; ++row) {
      Sig12(row) = phi_par(ind_noi(row), i);
    }
    Eigen::MatrixXd invC11 = Sig11 - (Sig12 * Sig12.transpose()) / phi_par(i, i);
    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Zero(size, size);
    for(int idx = 0; idx < size; ++idx) {
      diagMatrix(idx, idx) = 1.0 / taus(ind_noi(idx), i); // Assuming 'i' is defined elsewhere as the column index
    }

    Eigen::MatrixXd Ci = (S(i, i) + gammas) * invC11 + diagMatrix;

    Eigen::MatrixXd Sigma = Ci.llt().solve(Eigen::MatrixXd::Identity(size, size)); // Assuming chol2inv(chol(Ci)) equivalent
    Eigen::VectorXd mu_i(size);
    for (int idx = 0; idx < size; ++idx) {
      mu_i(idx) = -S(ind_noi(idx), i);
    }
    mu_i = Sigma * mu_i;
    arma::vec armaMean(mu_i.data(), mu_i.size(), false, true);
    arma::mat armaCov(Sigma.data(), Sigma.rows(), Sigma.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Rcpp::NumericMatrix betaMat = Rcpp::wrap(armasamples);

    //Rcpp::NumericMatrix betaMat = rmvnorm(1, mu_i, Sigma, seed);
    Eigen::Map<Eigen::MatrixXd> beta(betaMat.begin(), betaMat.nrow(), betaMat.ncol());  //beta is one col here

    for (int idx = 0; idx < size; ++idx) {
      inv_phi_par(ind_noi(idx), i) = beta(idx,0);
      inv_phi_par(i, ind_noi(idx)) = beta(idx,0);  // Assuming symmetry
    }

    std::mt19937 gen;
    std::gamma_distribution<double> gamma(N / 2.0 + 1,2.0/(S(i, i) + gammas));
    double gam  = gamma(gen);

    Eigen::MatrixXd tmp = beta.transpose() * invC11 * beta;
    inv_phi_par(i, i) = gam + tmp(0,0);

    // Updating covariance matrix 'phi_par' according to one-column change in the precision matrix 'inv_phi_par'
    Eigen::MatrixXd invC11beta = invC11 * beta;
    Sig12 = (-invC11beta / gam).col(0); //1 column
    Eigen::MatrixXd sub_phi =   invC11 + invC11beta * invC11beta.transpose()/gam;
    for (int idx = 0; idx < size; ++idx) {
      // Update off-diagonal elements for column 'i' and symmetric positions
      phi_par(ind_noi(idx), i) = Sig12(idx);
      phi_par(i,ind_noi(idx)) = Sig12(idx);
      for (int jdx = 0; jdx < size; ++jdx) {
        // Add invC11beta contribution to both off-diagonal and diagonal elements
        phi_par( ind_noi(idx), ind_noi(jdx)) = sub_phi(idx,jdx);
      }
    }
    phi_par(i, i) = 1.0 / gam;
  }

  Eigen::MatrixXd new_C_spe = phi_par;

  Eigen::MatrixXd new_C = C;
  new_C.block(0, 0, Kg, Kg) = new_C_gen;
  new_C.block(Kg, Kg, Ks, Ks)=new_C_spe;

  if (std) {
    Eigen::MatrixXd diag_gen = new_C_gen.diagonal().array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd tmp_gen = diag_gen.inverse();
    Eigen::MatrixXd C1_gen = tmp_gen * new_C_gen * tmp_gen;

    Eigen::MatrixXd diag_spe = new_C_spe.diagonal().array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd tmp_spe = diag_spe.inverse();
    Eigen::MatrixXd C1_spe = tmp_spe * new_C_spe * tmp_spe;

    Eigen::MatrixXd C1 = C;
    C1.block(0, 0, Kg, Kg) = C1_gen;
    C1.block(Kg, Kg, Ks, Ks)=C1_spe;

    double acc = exp((K + 1) / 2 * (log(C1.determinant()) - log(C.determinant())));
    double uniform_sample = arma::randu();
    bool acid = acc > uniform_sample;
    new_C = C1 * acid + C * (!acid);

  }
  C = new_C;
  return C;

}

Eigen::MatrixXd FacSampling::cor_cov_horse(Eigen::MatrixXd fac, bool std) {
  // Generate ind_nod matrix
  Eigen::MatrixXi ind_nod(K - 1, K);
  for (int k = 0; k < K; ++k) {
    int index = 0;
    for (int i = 0; i < K; ++i) {
      if (i != k) {
        ind_nod(index++, k) = i ; //different with R (i+1), 0-based
      }
    }
  }

  // Calculation of S
  Eigen::MatrixXd S = fac.transpose() * fac;

  //sample tau_sq and xi

  // Calculate rate
  double rate = 1.0 / xi;
  for (int i = 1; i < K; ++i) {
    for (int j = 0; j < i; ++j) {
      rate += inv_phi(i, j) * inv_phi(i, j) / (2.0 * Lambda_sq(i, j));
    }
  }
  // Gamma distribution for tau_sq and xi
  std::mt19937 gen;
  std::gamma_distribution<double> gamma_tau_sq(K * (K + 1) / 2.0, 1.0/rate);
  double tau_sq = gamma_tau_sq(gen); // Sample from gamma distribution
  std::gamma_distribution<double> gamma_xi(1.0, 1.0 + 1.0 / tau_sq);
  xi = gamma_xi(gen); // Sample from gamma distribution

  // Set lamdiag based on regdiag flag
  double lamdiag = regdiag ? tau_sq : 0.0;
  for (int i = 0; i < K; ++i) {
    Eigen::VectorXi ind_noi = ind_nod.col(i);
    // Manually constructing Sig11 as a submatrix of 'phi'
    int size = ind_noi.size();
    Eigen::MatrixXd Sig11(size, size);
    for(int row = 0; row < size; ++row) {
      for(int col = 0; col < size; ++col) {
        Sig11(row, col) = phi(ind_noi(row), ind_noi(col));
      }
    }

    // Extracting Sig12 as a subvector of column 'i' of 'phi'
    Eigen::VectorXd Sig12(size);
    for(int row = 0; row < size; ++row) {
      Sig12(row) = phi(ind_noi(row), i);
    }
    Eigen::MatrixXd invC11 = Sig11 - (Sig12 * Sig12.transpose()) / phi(i, i);

    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Zero(size, size);
    for(int idx = 0; idx < size; ++idx) {
      diagMatrix(idx, idx) = 1.0 /(Lambda_sq(ind_noi(idx), i)* tau_sq); // Assuming 'i' is defined elsewhere as the column index
    }
    Eigen::MatrixXd Ci = (S(i, i) + lamdiag) * invC11 + diagMatrix;
    Eigen::MatrixXd Sigma = Ci.llt().solve(Eigen::MatrixXd::Identity(size, size)); // Assuming chol2inv(chol(Ci)) equivalent
    Eigen::VectorXd mu_i(size);
    for (int idx = 0; idx < size; ++idx) {
      mu_i(idx) = -S(ind_noi(idx), i);
    }
    mu_i = Sigma * mu_i;
    // bool sys=Sigma.isApprox(Sigma.transpose());
    // if (!sys) {
    //   Rcpp::warning("The covariance matrix is not symmetric.");
    // } else {
    //   Rcpp::Rcout << "The covariance matrix is symmetric." << std::endl;
    // }
    //
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Sigma);
    // bool pd = es.eigenvalues().minCoeff() >= 0;
    // if (!pd) {
    //   Rcpp::stop("The covariance matrix is not positive semi-definite.");
    // } else {
    //   Rcpp::Rcout << "The covariance matrix is positive semi-definite." << std::endl;
    // }
    arma::vec armaMean(mu_i.data(), mu_i.size(), false, true);
    arma::mat armaCov(Sigma.data(), Sigma.rows(), Sigma.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Rcpp::NumericMatrix betaMat = Rcpp::wrap(armasamples);

    //Rcpp::NumericMatrix betaMat = rmvnorm(1, mu_i, Sigma, seed);
    Eigen::Map<Eigen::MatrixXd> beta(betaMat.begin(), betaMat.nrow(), betaMat.ncol());  //beta is one col here
    for (int idx = 0; idx < size; ++idx) {
      inv_phi(ind_noi(idx), i) = beta(idx,0);
      inv_phi(i, ind_noi(idx)) = beta(idx,0);  // Assuming symmetry
    }
    std::mt19937 gen;
    std::gamma_distribution<double> gamma(N / 2.0 + 1,2.0/(S(i, i) + lamdiag));
    double gam  = gamma(gen);
    Eigen::MatrixXd tmp = beta.transpose() * invC11 * beta;
    inv_phi(i, i) = gam + tmp(0,0);

    // Updating covariance matrix 'phi' according to one-column change in the precision matrix 'inv_phi'
    Eigen::MatrixXd invC11beta = invC11 * beta;
    Sig12 = (-invC11beta / gam).col(0); //1 column
    Eigen::MatrixXd sub_phi =   invC11 + invC11beta * invC11beta.transpose()/gam;
    // upadate phi[ind_noi,i]
    for (int idx = 0; idx < size; ++idx) {
      phi(ind_noi(idx), i) = Sig12(idx);
      phi(i, ind_noi(idx)) = Sig12(idx);
      for (int jdx = 0; jdx < size; ++jdx) {
        // Add invC11beta contribution to both off-diagonal and diagonal elements
        phi(ind_noi(idx), ind_noi(jdx)) = sub_phi(idx,jdx);
      }
    }

    // Calculate rate for Lambda_sq
    //Eigen::VectorXd beta_squared = beta.transpose() * beta; // This is a scalar if beta is a vector
    //Eigen::VectorXd rate_lambda(size);
    double beta_product = (beta.transpose() * beta)(0, 0);
    //= (beta_squared / (2.0 * tau_sq)).array() + (1.0 / Nu(row, i)).array();
    // for (int j = 0; j < size; ++j) {
    //   rate_lambda(j) = (beta_product / (2.0 * tau_sq)) + 1.0 / Nu(ind_noi(j), i); // Calculate rate for each index
    // }

    // Sample from gamma distribution for Lambda_sq
    //Eigen::VectorXd sampled_lambda(rate_lambda.size());
    for (int j = 0; j < size; ++j) {
      std::gamma_distribution<double> gamma_sampled_lambda(1.0,(beta_product / (2.0 * tau_sq)) + 1.0 / Nu(ind_noi(j), i)); // Shape = 1, Scale = 1/rate
      //sampled_lambda(j) = gamma_sampled_lambda(gen);
      //double tmp_Lambda_sq= gamma_sampled_lambda(gen);
      Lambda_sq(i, ind_noi(j)) = gamma_sampled_lambda(gen);
      Lambda_sq(ind_noi(j), i) = Lambda_sq(i, ind_noi(j));
      //double rate_nu =1.0 + 1.0 / Lambda_sq(i, ind_noi(j)); // Calculate new rate for Nu
      std::gamma_distribution<double> gamma_Nu(1.0, (1.0 + 1.0 / Lambda_sq(i, ind_noi(j)))); // Shape = 1, Scale = rate_nu
      Nu(i, ind_noi(j)) = gamma_Nu(gen);
      Nu(ind_noi(j), i) = Nu(i, ind_noi(j)); // Update symmetric entry
    }
    phi(i, i) = 1.0 / gam;
  }
  Eigen::MatrixXd new_C = phi;
  if (std) {
    Eigen::MatrixXd diag = new_C.diagonal().array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd tmp = diag.inverse();
    Eigen::MatrixXd C1 = tmp * new_C * tmp;
    new_C = C1;
  }
  C=new_C;
  return new_C;

}
Eigen::MatrixXd FacSampling::cor_cov_horse_par(Eigen::MatrixXd fac, bool std) {

  int Ks = K - Kg;
  Eigen::MatrixXd fac_g = fac.block(0, 0, N, Kg);           // first Kg columns , user need to make the loading in this order, Kg + Ks
  Eigen::MatrixXd fac_s = fac.block(0, Kg, N, Ks);  // remaining Ks columns
  Eigen::MatrixXd sub_W0_par = sub_W0.block(0, 0, Kg, Kg);
  // general factor
  Eigen::MatrixXd scaleInv = (fac_g.transpose() * fac_g + sub_W0_par).inverse();
  Eigen::MatrixXd inv_cov=rwish1( Kg + 2 + N,scaleInv);
  Eigen::MatrixXd new_C_gen = inv_cov.inverse();


  // Generate ind_nod matrix
  Eigen::MatrixXi ind_nod(Ks - 1, Ks);
  for (int k = 0; k < Ks; ++k) {
    int index = 0;
    for (int i = 0; i < Ks; ++i) {
      if (i != k) {
        ind_nod(index++, k) = i ; //different with R (i+1), 0-based
      }
    }
  }

  // Calculation of S
  Eigen::MatrixXd S = fac_s.transpose() * fac_s;

  //sample tau_sq and xi

  // Calculate rate
  double rate = 1.0 / xi;
  for (int i = 1; i < Ks; ++i) {
    for (int j = 0; j < i; ++j) {
      rate += inv_phi_par(i, j) * inv_phi_par(i, j) / (2.0 * Lambda_sq(i, j));
    }
  }
  // Gamma distribution for tau_sq and xi
  std::mt19937 gen;
  std::gamma_distribution<double> gamma_tau_sq(Ks * (Ks + 1) / 2.0, 1.0/rate);
  double tau_sq = gamma_tau_sq(gen); // Sample from gamma distribution
  std::gamma_distribution<double> gamma_xi(1.0, 1.0 + 1.0 / tau_sq);
  xi = gamma_xi(gen); // Sample from gamma distribution

  // Set lamdiag based on regdiag flag
  double lamdiag = regdiag ? tau_sq : 0.0;
  for (int i = 0; i < Ks; ++i) {
    Eigen::VectorXi ind_noi = ind_nod.col(i);
    // Manually constructing Sig11 as a submatrix of 'phi'
    int size = ind_noi.size();
    Eigen::MatrixXd Sig11(size, size);
    for(int row = 0; row < size; ++row) {
      for(int col = 0; col < size; ++col) {
        Sig11(row, col) = phi_par(ind_noi(row), ind_noi(col));
      }
    }

    // Extracting Sig12 as a subvector of column 'i' of 'phi'
    Eigen::VectorXd Sig12(size);
    for(int row = 0; row < size; ++row) {
      Sig12(row) = phi_par(ind_noi(row), i);
    }
    Eigen::MatrixXd invC11 = Sig11 - (Sig12 * Sig12.transpose()) / phi_par(i, i);

    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Zero(size, size);
    for(int idx = 0; idx < size; ++idx) {
      diagMatrix(idx, idx) = 1.0 /(Lambda_sq(ind_noi(idx), i)* tau_sq); // Assuming 'i' is defined elsewhere as the column index
    }
    Eigen::MatrixXd Ci = (S(i, i) + lamdiag) * invC11 + diagMatrix;
    Eigen::MatrixXd Sigma = Ci.llt().solve(Eigen::MatrixXd::Identity(size, size)); // Assuming chol2inv(chol(Ci)) equivalent
    Eigen::VectorXd mu_i(size);
    for (int idx = 0; idx < size; ++idx) {
      mu_i(idx) = -S(ind_noi(idx), i);
    }
    mu_i = Sigma * mu_i;
    // bool sys=Sigma.isApprox(Sigma.transpose());
    // if (!sys) {
    //   Rcpp::warning("The covariance matrix is not symmetric.");
    // } else {
    //   Rcpp::Rcout << "The covariance matrix is symmetric." << std::endl;
    // }
    //
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Sigma);
    // bool pd = es.eigenvalues().minCoeff() >= 0;
    // if (!pd) {
    //   Rcpp::stop("The covariance matrix is not positive semi-definite.");
    // } else {
    //   Rcpp::Rcout << "The covariance matrix is positive semi-definite." << std::endl;
    // }
    arma::vec armaMean(mu_i.data(), mu_i.size(), false, true);
    arma::mat armaCov(Sigma.data(), Sigma.rows(), Sigma.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Rcpp::NumericMatrix betaMat = Rcpp::wrap(armasamples);

    //Rcpp::NumericMatrix betaMat = rmvnorm(1, mu_i, Sigma, seed);
    Eigen::Map<Eigen::MatrixXd> beta(betaMat.begin(), betaMat.nrow(), betaMat.ncol());  //beta is one col here
    for (int idx = 0; idx < size; ++idx) {
      inv_phi_par(ind_noi(idx), i) = beta(idx,0);
      inv_phi_par(i, ind_noi(idx)) = beta(idx,0);  // Assuming symmetry
    }
    std::mt19937 gen;
    std::gamma_distribution<double> gamma(N / 2.0 + 1,2.0/(S(i, i) + lamdiag));
    double gam  = gamma(gen);
    Eigen::MatrixXd tmp = beta.transpose() * invC11 * beta;
    inv_phi_par(i, i) = gam + tmp(0,0);

    // Updating covariance matrix 'phi' according to one-column change in the precision matrix 'inv_phi'
    Eigen::MatrixXd invC11beta = invC11 * beta;
    Sig12 = (-invC11beta / gam).col(0); //1 column
    Eigen::MatrixXd sub_phi =   invC11 + invC11beta * invC11beta.transpose()/gam;
    // upadate phi[ind_noi,i]
    for (int idx = 0; idx < size; ++idx) {
      phi_par(ind_noi(idx), i) = Sig12(idx);
      phi_par(i, ind_noi(idx)) = Sig12(idx);
      for (int jdx = 0; jdx < size; ++jdx) {
        // Add invC11beta contribution to both off-diagonal and diagonal elements
        phi_par(ind_noi(idx), ind_noi(jdx)) = sub_phi(idx,jdx);
      }
    }


    double beta_product = (beta.transpose() * beta)(0, 0);

    for (int j = 0; j < size; ++j) {
      std::gamma_distribution<double> gamma_sampled_lambda(1.0,(beta_product / (2.0 * tau_sq)) + 1.0 / Nu(ind_noi(j), i)); // Shape = 1, Scale = 1/rate
      Lambda_sq(i, ind_noi(j)) = gamma_sampled_lambda(gen);
      Lambda_sq(ind_noi(j), i) = Lambda_sq(i, ind_noi(j));
      //double rate_nu =1.0 + 1.0 / Lambda_sq(i, ind_noi(j)); // Calculate new rate for Nu
      std::gamma_distribution<double> gamma_Nu(1.0, (1.0 + 1.0 / Lambda_sq(i, ind_noi(j)))); // Shape = 1, Scale = rate_nu
      Nu(i, ind_noi(j)) = gamma_Nu(gen);
      Nu(ind_noi(j), i) = Nu(i, ind_noi(j)); // Update symmetric entry
    }
    phi_par(i, i) = 1.0 / gam;
  }
  Eigen::MatrixXd new_C_spe = phi_par;

  Eigen::MatrixXd new_C = C;
  new_C.block(0, 0, Kg, Kg) = new_C_gen;
  new_C.block(Kg, Kg, Ks, Ks)=new_C_spe;

  if (std) {
    Eigen::MatrixXd diag_gen = new_C_gen.diagonal().array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd tmp_gen = diag_gen.inverse();
    Eigen::MatrixXd C1_gen = tmp_gen * new_C_gen * tmp_gen;

    Eigen::MatrixXd diag_spe = new_C_spe.diagonal().array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd tmp_spe = diag_spe.inverse();
    Eigen::MatrixXd C1_spe = tmp_spe * new_C_spe * tmp_spe;

    Eigen::MatrixXd C1 = C;
    C1.block(0, 0, Kg, Kg) = C1_gen;
    C1.block(Kg, Kg, Ks, Ks)=C1_spe;

    double acc = exp((K + 1) / 2 * (log(C1.determinant()) - log(C.determinant())));
    double uniform_sample = arma::randu();
    bool acid = acc > uniform_sample;
    new_C = C1 * acid + C * (!acid);

  }
  C = new_C;
  return C;

}

Eigen::MatrixXd FacSampling::cor_cov_ssp(Eigen::MatrixXd fac, bool std) {
  // Generate ind_nod matrix
  Eigen::MatrixXi ind_nod(K - 1, K);
  for (int k = 0; k < K; ++k) {
    int index = 0;
    for (int i = 0; i < K; ++i) {
      if (i != k) {
        ind_nod(index++, k) = i ; //different with R (i+1), 0-based
      }
    }
  }

  // Calculation of S
  Eigen::MatrixXd S = fac.transpose() * fac;

  double lamdiag = regdiag ? lambda : 0.0;

  for (int i = 0; i < K; ++i) {
    Eigen::VectorXi ind_noi = ind_nod.col(i);
    // Manually constructing Sig11 as a submatrix of 'phi'
    int size = ind_noi.size();
    Eigen::MatrixXd Sig11(size, size);
    for(int row = 0; row < size; ++row) {
      for(int col = 0; col < size; ++col) {
        Sig11(row, col) = phi(ind_noi(row), ind_noi(col));
      }
    }

    // Extracting Sig12 as a subvector of column 'i' of 'phi'
    Eigen::VectorXd Sig12(size);
    for(int row = 0; row < size; ++row) {
      Sig12(row) = phi(ind_noi(row), i);
    }

    Eigen::MatrixXd invC11 = Sig11 - (Sig12 * Sig12.transpose()) / phi(i, i);

    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Zero(size, size);
    for(int idx = 0; idx < size; ++idx) {
      diagMatrix(idx, idx) = 1.0 / taus(ind_noi(idx), i); // Assuming 'i' is defined elsewhere as the column index
    }

    Eigen::MatrixXd Ci = (S(i, i) + lamdiag) * invC11 + diagMatrix;
    Eigen::MatrixXd Sigma = Ci.llt().solve(Eigen::MatrixXd::Identity(size, size)); // Assuming chol2inv(chol(Ci)) equivalent
    Eigen::VectorXd mu_i(size);
    for (int idx = 0; idx < size; ++idx) {
      mu_i(idx) = -S(ind_noi(idx), i);
    }

    mu_i = Sigma * mu_i;
    arma::vec armaMean(mu_i.data(), mu_i.size(), false, true);
    arma::mat armaCov(Sigma.data(), Sigma.rows(), Sigma.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Rcpp::NumericMatrix betaMat = Rcpp::wrap(armasamples);

    //Rcpp::NumericMatrix betaMat = rmvnorm(1, mu_i, Sigma, seed);
    Eigen::Map<Eigen::MatrixXd> beta(betaMat.begin(), betaMat.nrow(), betaMat.ncol());  //beta is one col here
    for (int idx = 0; idx < size; ++idx) {
      inv_phi(ind_noi(idx), i) = beta(idx,0);
      inv_phi(i, ind_noi(idx)) = beta(idx,0);  // Assuming symmetry
    }
    std::mt19937 gen;
    std::gamma_distribution<double> gamma(N / 2.0 + 1,2.0/(S(i, i) + lamdiag));
    double gam  = gamma(gen);
    Eigen::MatrixXd tmp = beta.transpose() * invC11 * beta;
    inv_phi(i, i) = gam + tmp(0,0);

    // Updating covariance matrix 'phi' according to one-column change in the precision matrix 'inv_phi'
    Eigen::MatrixXd invC11beta = invC11 * beta;
    Sig12 = (-invC11beta / gam).col(0); //1 column
    Eigen::MatrixXd sub_phi =   invC11 + invC11beta * invC11beta.transpose()/gam;

    Eigen::VectorXd v0(size);
    Eigen::VectorXd v1(size);

    for (int idx = 0; idx < size; ++idx) {
      // Update off-diagonal elements for column 'i' and symmetric positions
      phi(ind_noi(idx), i) = Sig12(idx);
      phi(i, ind_noi(idx)) = Sig12(idx);
      for (int jdx = 0; jdx < size; ++jdx) {
        // Add invC11beta contribution to both off-diagonal and diagonal elements
        phi(ind_noi(idx), ind_noi(jdx)) = sub_phi(idx,jdx);
      }

      // Extracting v0 and v1
      for (int j = 0; j < size; j++) {
        v0(j) = V0(ind_noi[j], i);
        v1(j) = V1(ind_noi[j], i);
      }

      // Calculate w1 and w2
      Eigen::VectorXd w1(size);
      w1 = -0.5 * v0.array().log() - (0.5 * beta.array().square().col(0) / v0.array()) + log(1 - pii);
      Eigen::VectorXd w2(size);
      w2 = -0.5 * v1.array().log() - (0.5 * beta.array().square().col(0) / v1.array()) + log(pii);
      // Find w_max
      Eigen::VectorXd w_max(size);
      for (int j = 0; j < size; j++) {
        w_max(j) = std::max(w1(j), w2(j));
      }

      // Calculate w
      Eigen::VectorXd exp_w1 = (w1 - w_max).array().exp();
      Eigen::VectorXd exp_w2 = (w2 - w_max).array().exp();
      Eigen::VectorXd rowSums = exp_w1 + exp_w2;
      Eigen::VectorXd w = exp_w2.array() / rowSums.array();

      // Set up random number generator
      std::mt19937 gen; // Seed the generator
      std::uniform_real_distribution<> dis(0.0, 1.0); // Define the range

      // Generate z, taus
      Eigen::VectorXd z(size);
      for (int j = 0; j < size; j++) {
        double u = dis(gen); // Generate a random number
        z(j) = (u < w(j)) ? 1.0 : 0.0; // Compare and assign 1 or 0
        mip(i,ind_noi[j]) =z(j);
        mip(ind_noi[j],i)=z(j);
        if (z(j) == 1) {
          taus(ind_noi[j], i) = v1(j);
          taus(i, ind_noi[j]) = v1(j);

        }else{
          taus(ind_noi[j], i) = v0(j);
          taus(i, ind_noi[j]) = v0(j);
        }
      }
      // // Update v
      // Eigen::VectorXd v = v0;
      // for (int j = 0; j < size; j++) {
      //   if (z(j) == 1) {
      //     v(j) = v1(j);
      //   }
      // }
      // // Update taus
      // for (int j = 0; j < size; j++) {
      //   taus(ind_noi[j], i) = v(j);
      //   taus(i, ind_noi[j]) = v(j);
      // }
    }
    phi(i, i) = 1.0 / gam;
  }
  Eigen::MatrixXd new_C = phi;
  if (std) {
    Eigen::MatrixXd diag = new_C.diagonal().array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd tmp = diag.inverse();
    Eigen::MatrixXd C1 = tmp * new_C * tmp;

    new_C = C1;

    // MH
    // double acc = exp((K + 1) / 2 * (log(C1.determinant()) - log(C.determinant())));
    // double uniform_sample = arma::randu();
    // bool acid = acc > uniform_sample;
    // new_C = C1 * acid + C * (!acid);
  }
  C=new_C;
  return new_C;

}
Eigen::MatrixXd FacSampling::cor_cov_ssp_par(Eigen::MatrixXd fac, bool std) {

  int Ks = K - Kg;
  Eigen::MatrixXd fac_g = fac.block(0, 0, N, Kg);           // first Kg columns , user need to make the loading in this order, Kg + Ks
  Eigen::MatrixXd fac_s = fac.block(0, Kg, N, Ks);  // remaining Ks columns
  Eigen::MatrixXd sub_W0_par = sub_W0.block(0, 0, Kg, Kg);
  // general factor
  Eigen::MatrixXd scaleInv = (fac_g.transpose() * fac_g + sub_W0_par).inverse();
  Eigen::MatrixXd inv_cov=rwish1( Kg + 2 + N,scaleInv);
  Eigen::MatrixXd new_C_gen = inv_cov.inverse();




  // Generate ind_nod matrix
  Eigen::MatrixXi ind_nod(Ks - 1, Ks);
  for (int k = 0; k < Ks; ++k) {
    int index = 0;
    for (int i = 0; i < Ks; ++i) {
      if (i != k) {
        ind_nod(index++, k) = i ; //different with R (i+1), 0-based
      }
    }
  }

  // Calculation of S
  Eigen::MatrixXd S = fac_s.transpose() * fac_s;

  for (int i = 0; i < Ks; ++i) {
    Eigen::VectorXi ind_noi = ind_nod.col(i);
    // Manually constructing Sig11 as a submatrix of 'phi'
    int size = ind_noi.size();
    Eigen::MatrixXd Sig11(size, size);
    for(int row = 0; row < size; ++row) {
      for(int col = 0; col < size; ++col) {
        Sig11(row, col) = phi_par(ind_noi(row), ind_noi(col));
      }
    }

    // Extracting Sig12 as a subvector of column 'i' of 'phi_par'
    Eigen::VectorXd Sig12(size);
    for(int row = 0; row < size; ++row) {
      Sig12(row) = phi_par(ind_noi(row), i);
    }

    Eigen::MatrixXd invC11 = Sig11 - (Sig12 * Sig12.transpose()) / phi_par(i, i);

    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Zero(size, size);
    for(int idx = 0; idx < size; ++idx) {
      diagMatrix(idx, idx) = 1.0 / taus_par(ind_noi(idx), i); // Assuming 'i' is defined elsewhere as the column index
    }

    Eigen::MatrixXd Ci = (S(i, i) + lambda) * invC11 + diagMatrix;
    Eigen::MatrixXd Sigma = Ci.llt().solve(Eigen::MatrixXd::Identity(size, size)); // Assuming chol2inv(chol(Ci)) equivalent
    Eigen::VectorXd mu_i(size);
    for (int idx = 0; idx < size; ++idx) {
      mu_i(idx) = -S(ind_noi(idx), i);
    }

    mu_i = Sigma * mu_i;
    arma::vec armaMean(mu_i.data(), mu_i.size(), false, true);
    arma::mat armaCov(Sigma.data(), Sigma.rows(), Sigma.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Rcpp::NumericMatrix betaMat = Rcpp::wrap(armasamples);

    //Rcpp::NumericMatrix betaMat = rmvnorm(1, mu_i, Sigma, seed);
    Eigen::Map<Eigen::MatrixXd> beta(betaMat.begin(), betaMat.nrow(), betaMat.ncol());  //beta is one col here
    for (int idx = 0; idx < size; ++idx) {
      inv_phi_par(ind_noi(idx), i) = beta(idx,0);
      inv_phi_par(i, ind_noi(idx)) = beta(idx,0);  // Assuming symmetry
    }
    std::mt19937 gen;
    std::gamma_distribution<double> gamma(N / 2.0 + 1,2.0/(S(i, i) + lambda));
    double gam  = gamma(gen);
    Eigen::MatrixXd tmp = beta.transpose() * invC11 * beta;
    inv_phi_par(i, i) = gam + tmp(0,0);

    // Updating covariance matrix 'phi_par' according to one-column change in the precision matrix 'inv_phi'
    Eigen::MatrixXd invC11beta = invC11 * beta;
    Sig12 = (-invC11beta / gam).col(0); //1 column
    Eigen::MatrixXd sub_phi =   invC11 + invC11beta * invC11beta.transpose()/gam;

    Eigen::VectorXd v0(size);
    Eigen::VectorXd v1(size);

    for (int idx = 0; idx < size; ++idx) {
      // Update off-diagonal elements for column 'i' and symmetric positions
      phi_par(ind_noi(idx), i) = Sig12(idx);
      phi_par(i, ind_noi(idx)) = Sig12(idx);
      for (int jdx = 0; jdx < size; ++jdx) {
        // Add invC11beta contribution to both off-diagonal and diagonal elements
        phi_par(ind_noi(idx), ind_noi(jdx)) = sub_phi(idx,jdx);
      }

      // Extracting v0 and v1
      for (int j = 0; j < size; j++) {
        v0(j) = V0_par(ind_noi[j], i);
        v1(j) = V1_par(ind_noi[j], i);
      }

      // Calculate w1 and w2
      Eigen::VectorXd w1(size);
      w1 = -0.5 * v0.array().log() - (0.5 * beta.array().square().col(0) / v0.array()) + log(1 - pii_par);
      Eigen::VectorXd w2(size);
      w2 = -0.5 * v1.array().log() - (0.5 * beta.array().square().col(0) / v1.array()) + log(pii_par);
      // Find w_max
      Eigen::VectorXd w_max(size);
      for (int j = 0; j < size; j++) {
        w_max(j) = std::max(w1(j), w2(j));
      }

      // Calculate w
      Eigen::VectorXd exp_w1 = (w1 - w_max).array().exp();
      Eigen::VectorXd exp_w2 = (w2 - w_max).array().exp();
      Eigen::VectorXd rowSums = exp_w1 + exp_w2;
      Eigen::VectorXd w = exp_w2.array() / rowSums.array();

      // Set up random number generator
      std::mt19937 gen; // Seed the generator
      std::uniform_real_distribution<> dis(0.0, 1.0); // Define the range

      // Generate z, taus_par
      Eigen::VectorXd z(size);
      for (int j = 0; j < size; j++) {
        double u = dis(gen); // Generate a random number
        z(j) = (u < w(j)) ? 1.0 : 0.0; // Compare and assign 1 or 0
        if (z(j) == 1) {
          taus_par(ind_noi[j], i) = v1(j);
          taus_par(i, ind_noi[j]) = v1(j);
        }else{
          taus_par(ind_noi[j], i) = v0(j);
          taus_par(i, ind_noi[j]) = v0(j);
        }
      }
      // // Update v
      // Eigen::VectorXd v = v0;
      // for (int j = 0; j < size; j++) {
      //   if (z(j) == 1) {
      //     v(j) = v1(j);
      //   }
      // }
      // // Update taus
      // for (int j = 0; j < size; j++) {
      //   taus(ind_noi[j], i) = v(j);
      //   taus(i, ind_noi[j]) = v(j);
      // }
    }
    phi_par(i, i) = 1.0 / gam;
  }
  Eigen::MatrixXd new_C_spe = phi_par;



  Eigen::MatrixXd new_C = C;
  new_C.block(0, 0, Kg, Kg) = new_C_gen;
  new_C.block(Kg, Kg, Ks, Ks)=new_C_spe;

  if (std) {
    Eigen::MatrixXd diag_gen = new_C_gen.diagonal().array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd tmp_gen = diag_gen.inverse();
    Eigen::MatrixXd C1_gen = tmp_gen * new_C_gen * tmp_gen;

    Eigen::MatrixXd diag_spe = new_C_spe.diagonal().array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd tmp_spe = diag_spe.inverse();
    Eigen::MatrixXd C1_spe = tmp_spe * new_C_spe * tmp_spe;

    Eigen::MatrixXd C1 = C;
    C1.block(0, 0, Kg, Kg) = C1_gen;
    C1.block(Kg, Kg, Ks, Ks)=C1_spe;

    double acc = exp((K + 1) / 2 * (log(C1.determinant()) - log(C.determinant())));
    double uniform_sample = arma::randu();
    bool acid = acc > uniform_sample;
    new_C = C1 * acid + C * (!acid);

  }
  C = new_C;

  return C;

}
