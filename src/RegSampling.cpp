#include "RegSampling.h"
#include "utility-function.h"

RegSampling::RegSampling(int N,bool LD, double s0, double r0 , double lam_a , double lam_b, double m0, double c0)
  : N(N),LD(LD),s0(s0), r0(r0), lam_a(lam_a), lam_b(lam_b), m0(m0), c0(c0){}

Eigen::MatrixXd RegSampling::initial(int N, Eigen::MatrixXi Q,std::string reglo) {
  regloading=reglo;
  N = N;
  P = Q.cols();
  J = Q.rows();
  sub1 = (Q.array() == -1).cast<int>();
  sub2 = (Q.array() == -2).cast<int>();
  jp1 = sub1.sum();
  jp2 = sub2.sum();

  std::vector<std::pair<int, int>> indices1;
  std::vector<std::pair<int, int>> indices2;
  for (int i = 0; i < J; ++i) {
    for (int j = 0; j < P; ++j) {
      if (sub1(i, j) == 1) {
        indices1.push_back(std::make_pair(i, j));
      }else if (sub2(i, j) == 1) {
        indices2.push_back(std::make_pair(i, j));
      }
    }
  }


  // Convert indices vector to pos1 matrix
  pos1.resize(indices1.size(), 2);
  for (size_t k = 0; k < indices1.size(); ++k) {
    pos1(k, 0) = indices1[k].first;
    pos1(k, 1) = indices1[k].second;
  }
  pos2.resize(indices2.size(), 2);
  for (size_t k = 0; k < indices2.size(); ++k) {
    pos2(k, 0) = indices2[k].first;
    pos2(k, 1) = indices2[k].second;
  }


  // Calculate jnd1 and j1
  std::vector<int> jnd1List;
  for (int i = 0; i < J; ++i) {
    if (sub1.row(i).sum() > 0) {
      jnd1List.push_back(i);
    }
  }
  std::vector<int> jnd2List;
  for (int i = 0; i < J; ++i) {
    if (sub2.row(i).sum() > 0) {
      jnd2List.push_back(i);
    }
  }
  j1 = jnd1List.size();
  j2 = jnd2List.size();
  // Convert jnd1List to Eigen::VectorXi
  jnd1.resize(jnd1List.size());
  for (size_t i = 0; i < jnd1List.size(); ++i) {
    jnd1(i) = jnd1List[i];
  }
  jnd2.resize(jnd2List.size());
  for (size_t i = 0; i < jnd2List.size(); ++i) {
    jnd2(i) = jnd2List[i];
  }

  //Initialize parameter estimates
  tausq = Eigen::MatrixXd::Zero(J, P);
  if (j1>0){
    for (int i = 0; i < jp1; ++i) {
      tausq(pos1(i, 0), pos1(i, 1)) = 1.0; // Example, adjust as needed
    }
  }

  if(regloading=="sspPEFA"){
    if (j2>0){
      for (int i = 0; i < jp2; ++i) {
        tausq(pos2(i, 0), pos2(i, 1)) = 1.0; // Example, adjust as needed
      }
    }
  }
  r1_las = Eigen::VectorXd::Zero(J);
  vtausq = Eigen::VectorXd::Ones(P);
  dsig = Eigen::VectorXd::Constant(J, 0.5);
  psx = dsig.asDiagonal();
  Eigen::LLT<Eigen::MatrixXd> lltOfPSX(psx);
  if(lltOfPSX.info() == Eigen::NumericalIssue) {
    throw std::runtime_error("Matrix is not positive definite.");
  }
  inv_psx = lltOfPSX.solve(Eigen::MatrixXd::Identity(psx.rows(), psx.cols()));

  msig = Eigen::MatrixXd::Identity(J, J);
  lamsq = 1.0;
  Eigen::MatrixXd beta = Q.cast<double>();
  for (int i = 0; i < J; ++i) {
    for (int j = 0; j < P; ++j) {
      if (Q(i, j) == -1) beta(i, j) = 0.1;
      else if (Q(i, j) == -2) beta(i, j) = 0.5;
    }
  }

  //initial for LD and passing to cati
  psxjs= Eigen::VectorXd::Zero(J);
  mst=Eigen::MatrixXd::Zero(N, J);

  // initial for ssp
  Eigen::MatrixXd adj = (psx.array().abs() >= 0).cast<double>();
  //Eigen::MatrixXd::Constant(K, K, 1);
  // Calculate pii
  //pii = (J >= 4) ? (3.0 / (J - 1)) : 0.3; // Conditional assignment
  pii = (J >= 4) ? (2.0 / (J - 1)) : 0.3;
  // Create matrices V0 and V1
  lambda = 1;
  double v0tmp = 0.02*0.02;
  double h = 50*50;
  double v1tmp = h * v0tmp;
  V0 = Eigen::MatrixXd::Constant(J, J, v0tmp); // K x K matrix filled with v0
  V1 = Eigen::MatrixXd::Constant(J, J, v1tmp); // K x K matrix filled with v1
  // Initialize taus with V0
  taus = V0;
  // Update taus based on the adjacency matrix
  for (int row = 0; row < J; ++row) {
    for (int col = 0; col < J; ++col) {
      if (adj(row, col) > 0) { // If the adjacency condition is met
        taus(row, col) = v1tmp;   // Set to v1 where adj is TRUE
      }
    }
  }

  //init for horse
  Lambda_sq = Eigen::MatrixXd::Constant(J, J, 1);
  Nu = Eigen::MatrixXd::Constant(J, J, 1);
  xi=1.0;
  regdiag= true;


  //init for loading horse
  tau=1.0;
  v = Eigen::MatrixXd::Zero(J, P);
  if (j1>0){
    for (int i = 0; i < jp1; ++i) {
      v(pos1(i, 0), pos1(i, 1)) = 0.5; // Example, adjust as needed
    }
  }
  eta = 0.5;

  //init for loading ssp
  sksq = Eigen::VectorXd::Ones(P);
  sksq_t = sksq;
  return beta;
};
Eigen::MatrixXd RegSampling ::bmvreg_plas(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd beta) {
  auto beta_func1= [this,Y,X,&beta](int j) {
    int i = jnd1[j] ;
    Eigen::MatrixXi ones = Eigen::MatrixXi::Ones(sub1.rows(), sub1.cols());
    Eigen::MatrixXi u_sub1 = ones - sub1;
    Eigen::MatrixXd subtausq = selectColumns(tausq,  sub1,  i).row(i).cwiseInverse().asDiagonal(); //*(1.0/tau)
    subtausq = subtausq * (1.0 / tau);

    // Eigen::MatrixXd Sn = ((selectColumns(X,  sub1,  i).transpose() * selectColumns(X,  sub1,  i))/dsig[i] + subtausq).inverse();
    // Eigen::MatrixXd tmp = Y.col(i) - (selectColumns(X,  u_sub1,  i) * selectColumns(beta,  u_sub1,  i).row(i).transpose());
    // Eigen::VectorXd mu = Sn * (selectColumns(X,  sub1,  i).transpose() * tmp/dsig[i]);
    // Eigen::MatrixXd cov = Sn ;

    //FOR SSP
    Eigen::MatrixXd Sn = ((selectColumns(X,  sub1,  i).transpose() * selectColumns(X,  sub1,  i)) + subtausq).inverse();
    Eigen::MatrixXd tmp = Y.col(i) - (selectColumns(X,  u_sub1,  i) * selectColumns(beta,  u_sub1,  i).row(i).transpose());
    Eigen::VectorXd mu = Sn * (selectColumns(X,  sub1,  i).transpose() * tmp);
    Eigen::MatrixXd cov = Sn*dsig[i] ;


    arma::vec armaMean(mu.data(), mu.size(), false, true);
    arma::mat armaCov(cov.data(), cov.rows(), cov.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Rcpp::NumericMatrix sample = Rcpp::wrap(armasamples);
    return sample;
  };

  auto beta_func2 = [this,Y,X,&beta] (int j) {
    int i = jnd2[j] ;
    int size = sub2.row(i).count();
    Eigen::MatrixXd mat_c0 = Eigen::MatrixXd::Identity(size, size) * c0;
    Eigen::MatrixXd mat_m0 = Eigen::MatrixXd::Identity(size, size) * m0;

    Eigen::MatrixXi ones = Eigen::MatrixXi::Ones(sub2.rows(), sub2.cols());
    Eigen::MatrixXi u_sub2 = ones - sub2;
    Eigen::MatrixXd Sn = ((selectColumns(X,  sub2,  i).transpose() * selectColumns(X,  sub2,  i)) + mat_c0).inverse();
    Eigen::MatrixXd tmp = Y.col(i) - (selectColumns(X,  u_sub2,  i) * selectColumns(beta,  u_sub2,  i).row(i).transpose());
    Eigen::VectorXd mu = (Sn * (selectColumns(X, sub2, i).transpose() * tmp )).array() + (mat_m0 * mat_c0).diagonal().array();
    Eigen::MatrixXd cov = Sn* dsig[i] ; //

    arma::vec armaMean(mu.data(), mu.size(), false, true);
    arma::mat armaCov(cov.data(), cov.rows(), cov.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Rcpp::NumericMatrix sample = Rcpp::wrap(armasamples);

    return sample;
  };
  Eigen::MatrixXd beta_list2 = beta;
  Eigen::MatrixXd beta_list1 = beta;

  if (j2 > 0) {
    beta_list2 =beta;
    int index = 0;

    for (int j = 0; j < j2; ++j) {
      Rcpp::NumericMatrix temp = beta_func2(j);
      Eigen::Map<Eigen::MatrixXd> rowVec(temp.begin(), temp.nrow(), temp.ncol());
      for (int k = 0; k < rowVec.rows(); ++k) {  //here k=1, no need loop. but write for later might have two false in one row
        int rowIndex = pos2(index + k, 0);
        int colIndex = pos2(index + k, 1);
        // Check bounds to prevent out-of-bounds access
        if(rowIndex < beta.rows() && colIndex < beta.cols()) {
          beta_list2(rowIndex, colIndex) = rowVec(k, 0);
        }
      }
      index += rowVec.rows();
    }
  }
  beta=beta_list2;

  if (j1 > 0) {
    beta_list1 =beta;
    int index = 0;
    for (int j = 0; j < j1; ++j) {
      Rcpp::NumericMatrix temp = beta_func1(j);
      Eigen::Map<Eigen::MatrixXd> rowVec(temp.begin(), temp.nrow(), temp.ncol());
      for (int k = 0; k < rowVec.rows(); ++k) {
        int rowIndex = pos1(index + k, 0);
        int colIndex = pos1(index + k, 1);
        if(rowIndex < beta.rows() && colIndex < beta.cols()) {
          beta_list1(rowIndex, colIndex) = rowVec(k, 0);
        }
      }
      index += rowVec.rows();
    }
    beta=beta_list1;
    if(regloading=="lasso"){
      std::mt19937 gen;
      Eigen::MatrixXd inv_mu(beta.rows(), beta.cols());
      for (int row = 0; row < beta.rows(); ++row) {
        inv_mu.row(row) = (beta.row(row).array().square() / (dsig[row] * lamsq)).sqrt();
      }
      // Update tausq based on rinvgauss sampling
      Rcpp::NumericVector mu(pos1.rows());
      for (int i = 0; i < pos1.rows(); ++i) {
        int row = pos1(i, 0);
        int col = pos1(i, 1);
        mu[i] = 1.0 / inv_mu(row, col);
      }
      Rcpp::NumericVector sample = rinvgauss(jp1, mu, lamsq);
      for (int k = 0; k < sample.size(); ++k) {
        int rowIndex = pos1(k, 0);
        int colIndex = pos1(k, 1);
        tausq(rowIndex, colIndex) = 1.0/sample[k];

      }
      std::gamma_distribution<> gamma(lam_a + jp1, 1.0 / (lam_b + 0.5 * tausq.sum()));
      lamsq = gamma(gen);
    }else if(regloading=="ssp"){
      std::mt19937 gen;
      double v0 = 0.02 * 0.02;
      double h = 50 * 50;
      double pii = 0.5;
      double v1 = h * v0;
      for (int i = 0; i < jp1; ++i) {
        int row = pos1(i, 0);
        int col = pos1(i, 1);

        double w1 = -0.5 * std::log(v0) - 0.5 *  beta(row, col)*beta(row, col)  / (v0*dsig[row]) + std::log(1.0 - pii);
        double w2 = -0.5 * std::log(v1) - 0.5 *  beta(row, col)*beta(row, col) / (v1*dsig[row]) + std::log(pii);
        double w_max = std::max(w1, w2);
        double exp1 = std::exp(w1 - w_max);
        double exp2 = std::exp(w2 - w_max);
        double prob = exp2 / (exp1 + exp2);

        std::bernoulli_distribution bernoulli(prob);
        bool z = bernoulli(gen);
        tausq(row, col)  = z ? v1 : v0;
      }
    }else if(regloading=="horse"){

        std::mt19937 gen;
        for (int k = 0; k < jp1; ++k) {
          int row = pos1(k, 0);
          int col = pos1(k, 1);
          tausq(row, col) = 1.0/  R::rgamma(1, 1.0/(((beta(row, col) * beta(row, col))/ (2 *tau*dsig[row])) + 1.0 / v(row, col)));
          v(row, col)= 1.0/ R::rgamma(1,1.0/(1 + 1.0/tausq(row, col)));
        }

    }
   for(int i = 0; i < jnd1.size(); ++i) {
     int idx = jnd1[i]; //0-based
     double sum = 0.0;
     sum += (selectColumns(beta, sub1, idx).row(idx).array().square() / selectColumns(tausq, sub1, idx).row(idx).array()).sum();
     r1_las[idx] = sum / tau;
   }
  }


  Eigen::MatrixXd diff = Y - (X * beta.transpose());
  Eigen::VectorXd r1 = diff.array().square().colwise().sum();
  Eigen::VectorXd r1_combined = (r1 + r1_las) / 2.0;

  if(!LD){
    for(int i = 0; i < r1_combined.size(); ++i) {
      dsig[i] =1.0 /R::rgamma(s0 + N / 2.0, 1.0 / (r0 + r1_combined[i]));
      psx(i,i)=dsig[i];
    }

  }


  if(j1 > 0 && regloading=="horse"){
    std::mt19937 gen;
    double total_sum = 0.0;

    for (int j = 0; j < beta.rows(); j++) {
      double sum_k = 0.0;
      for (int k = 0; k < beta.cols(); k++) {
        if (sub1(j, k) == 1) {  // Check if Q[j, k] == -2
          sum_k += (beta(j, k) * beta(j, k))/tausq(j, k);
        }
      }
      total_sum += 0.5 * sum_k / dsig[j];
    }

    tau = 1.0/  R::rgamma((1+jp1)/2.0,  1.0/(1.0/eta+total_sum));
    eta = 1.0/  R::rgamma(1,  1.0/(1 + 1.0/tau));
  }

  return  beta;
}
Eigen::MatrixXd RegSampling::bmvreg_plas_LD(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd beta) {
  auto beta_func1 = [this, Y, X, &beta](int j, Eigen::VectorXd psxjs, Eigen::MatrixXd mst) {
    int i = jnd1[j];
    Eigen::MatrixXi ones = Eigen::MatrixXi::Ones(sub1.rows(), sub1.cols());
    Eigen::MatrixXi u_sub1 = ones - sub1;
    Eigen::MatrixXd subtausq = selectColumns(tausq, sub1, i).row(i).cwiseInverse().asDiagonal();
    subtausq = subtausq * (1.0 / tau);
    Eigen::MatrixXd Sn = ((selectColumns(X, sub1, i).transpose() * selectColumns(X, sub1, i))/psxjs[i] + subtausq).inverse();
    Eigen::MatrixXd tmp = Y.col(i) - mst.col(i) - (selectColumns(X, u_sub1, i) * selectColumns(beta, u_sub1, i).row(i).transpose());
    Eigen::VectorXd mu = Sn * (selectColumns(X, sub1, i).transpose() * tmp/psxjs[i]);
    Eigen::MatrixXd cov = Sn ;
    arma::vec armaMean(mu.data(), mu.size(), false, true);
    arma::mat armaCov(cov.data(), cov.rows(), cov.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Eigen::MatrixXd sample(armasamples.n_rows, armasamples.n_cols);
    for (int i = 0; i < armasamples.n_rows; ++i) {
      for (int j = 0; j < armasamples.n_cols; ++j) {
        sample(i, j) = armasamples(i, j);
      }
    }
    return sample;
  };

  auto beta_func2 = [this, Y, X, &beta](int j, Eigen::VectorXd psxjs, Eigen::MatrixXd mst) {
    int i = jnd2[j];
    int size = sub2.row(i).count();
    Eigen::MatrixXd mat_c0 = Eigen::MatrixXd::Identity(size, size) * c0;//Eigen::MatrixXd::Constant(size, size, c0);
    Eigen::MatrixXd mat_m0 = Eigen::MatrixXd::Identity(size, size) * m0; //Eigen::MatrixXd::Constant(size, size, m0);

    Eigen::MatrixXi ones = Eigen::MatrixXi::Ones(sub2.rows(), sub2.cols());
    Eigen::MatrixXi u_sub2 = ones - sub2;
    Eigen::MatrixXd Sn = ((selectColumns(X, sub2, i).transpose() * selectColumns(X, sub2, i))/psxjs[i] + mat_c0).inverse();
    Eigen::MatrixXd tmp = Y.col(i) - mst.col(i) - (selectColumns(X, u_sub2, i) * selectColumns(beta, u_sub2, i).row(i).transpose());
    Eigen::VectorXd mu = Sn * (selectColumns(X, sub2, i).transpose() * tmp/psxjs[i]) + mat_m0 * mat_c0;
    Eigen::MatrixXd cov = Sn;
    arma::vec armaMean(mu.data(), mu.size(), false, true);
    arma::mat armaCov(cov.data(), cov.rows(), cov.cols(), false, true);
    arma::mat armasamples = arma::mvnrnd(armaMean, armaCov, 1);
    Eigen::MatrixXd sample(armasamples.n_rows, armasamples.n_cols);
    for (int i = 0; i < armasamples.n_rows; ++i) {
      for (int j = 0; j < armasamples.n_cols; ++j) {
        sample(i, j) = armasamples(i, j);
      }
    }
    return sample;
  };

  for (int j = 0; j < J; j++){
    //std::cout << psx << std::endl;
    Eigen::VectorXd psx_j_col = psx.block(0, j, j, 1).eval();
    psx_j_col.conservativeResize(J - 1);
    psx_j_col.tail(J -j - 1) = psx.block(j + 1, j, J - j - 1, 1);

    Eigen::VectorXd psx_j_row = psx.block(j, 0, 1, j).transpose();
    psx_j_row.conservativeResize(J - 1);
    psx_j_row.tail(J - j - 1) = psx.block(j, j + 1, 1, J - j - 1).transpose();

    Eigen::MatrixXd psx_jj = psx.block(0, 0, j, j).eval();
    psx_jj.conservativeResize( J - 1 , J - 1);
    psx_jj.block(0, j, j, J - j - 1) = psx.block(0, j + 1, j, J - j - 1);
    psx_jj.block(j, 0, J - j - 1, j) = psx.block(j+1, 0, J - j - 1, j);
    psx_jj.block(j, j, J - j - 1, J - j - 1) = psx.block(j+1, j + 1, J - j - 1, J - j - 1);

    Eigen::LLT<Eigen::MatrixXd> llt(psx_jj);
    Eigen::VectorXd psx_tmp = psx_j_col.transpose() * llt.solve(Eigen::MatrixXd::Identity(psx_jj.rows(), psx_jj.cols()));

    psxjs(j) = psx(j, j) - (psx_tmp.transpose() * psx_j_row);

    Eigen::MatrixXd Y_j = Y.block(0, 0, N, j);
    Y_j.conservativeResize(N, J - 1);
    Y_j.block(0, j, N, J - j - 1) = Y.block(0, j + 1, N, J - j - 1);

    Eigen::MatrixXd beta_j = beta.block(0, 0, j, P);
    beta_j.conservativeResize(J - 1, P);
    beta_j.block(j, 0, J-j-1, P) = beta.block(j+1, 0, J - j - 1, P);

    mst.col(j) = (Y_j - X * beta_j.transpose()) * psx_tmp;
  }

  Eigen::MatrixXd beta_list2 = beta;
  Eigen::MatrixXd beta_list1 = beta;

  if (j2 > 0) {
    beta_list2 = beta;
    int index = 0;
    for (int j = 0; j < j2; ++j) {
      Eigen::MatrixXd rowVec = beta_func2(j, psxjs,mst);
      for (int k = 0; k < rowVec.rows(); ++k) {  //here k=1, no need loop. but write for later might have two false in one row
        int rowIndex = pos2(index + k, 0);
        int colIndex = pos2(index + k, 1);
        if (rowIndex < beta.rows() && colIndex < beta.cols()) {
          beta_list2(rowIndex, colIndex) = rowVec(k, 0);
        }
      }
      index += rowVec.rows();
    }
  }
  beta = beta_list2;
  if (j1 > 0) {
    beta_list1 = beta;
    int index = 0;
    for (int j = 0; j < j1; ++j) {
      Eigen::MatrixXd rowVec = beta_func1(j, psxjs, mst);
      for (int k = 0; k < rowVec.rows(); ++k) {
        int rowIndex = pos1(index + k, 0);
        int colIndex = pos1(index + k, 1);
        if (rowIndex < beta.rows() && colIndex < beta.cols()) {
          beta_list1(rowIndex, colIndex) = rowVec(k, 0);
        }
      }
      index += rowVec.rows();
    }
    beta = beta_list1;
    if(regloading=="lasso"){

      Eigen::MatrixXd inv_mu(beta.rows(), beta.cols());
      for (int row = 0; row < beta.rows(); ++row) {
        inv_mu.row(row) = (beta.row(row).array().square() / lamsq).sqrt();
     }
      // Update tausq based on rinvgauss sampling
      Rcpp::NumericVector mu(pos1.rows());
      for (int i = 0; i < pos1.rows(); ++i) {
        int row = pos1(i, 0);
        int col = pos1(i, 1);
        mu[i] = 1.0 / inv_mu(row, col);
      }
      Rcpp::NumericVector sample = rinvgauss(jp1, mu, lamsq);
      for (int k = 0; k < sample.size(); ++k) {
        int rowIndex = pos1(k, 0);
        int colIndex = pos1(k, 1);
        tausq(rowIndex, colIndex) = 1.0/sample[k];

      }
      std::mt19937 gen;
      std::gamma_distribution<> gamma(lam_a + jp1, 1.0 / (lam_b + 0.5 * tausq.sum()));
      lamsq = gamma(gen);

    }else if(regloading=="ssp"){
      std::mt19937 gen;
      double v0 = 0.02 * 0.02;
      double h = 50 * 50;
      double pii = 0.5;
      double v1 = h * v0;
      for (int i = 0; i < jp1; ++i) {
        int row = pos1(i, 0);
        int col = pos1(i, 1);

        double w1 = -0.5 * std::log(v0) - 0.5 *  beta(row, col)*beta(row, col)  / v0 + std::log(1.0 - pii);
        double w2 = -0.5 * std::log(v1) - 0.5 *  beta(row, col)*beta(row, col) / v1 + std::log(pii);
        double w_max = std::max(w1, w2);
        double exp1 = std::exp(w1 - w_max);
        double exp2 = std::exp(w2 - w_max);
        double prob = exp2 / (exp1 + exp2);

        std::bernoulli_distribution bernoulli(prob);
        bool z = bernoulli(gen);
        tausq(row, col)  = z ? v1 : v0;

      }
    }else if(regloading=="horse"){

      std::mt19937 gen;
      for (int k = 0; k < jp1; ++k) {
        int row = pos1(k, 0);
        int col = pos1(k, 1);
        tausq(row, col) = 1.0/  R::rgamma(1, 1.0/(((beta(row, col) * beta(row, col))/ (2 *tau)) + 1.0 / v(row, col)));
        v(row, col)= 1.0/ R::rgamma(1,1.0/(1 + 1.0/tausq(row, col)));
      }

    }

  }

  if(j1 > 0 && regloading=="horse"){
    std::mt19937 gen;
    double total_sum = 0.0;
    for (int j = 0; j < beta.rows(); j++) {
      double sum_k = 0.0;
      for (int k = 0; k < beta.cols(); k++) {
        if (sub1(j, k) == 1) {  // Check if Q[j, k] == -2
          sum_k += (beta(j, k) * beta(j, k))/tausq(j, k);
        }
      }
      total_sum += 0.5 * sum_k;
    }
    tau = 1.0/  R::rgamma((1+jp1)/2.0,  1.0/(1.0/eta+total_sum));
    eta = 1.0/  R::rgamma(1,  1.0/(1 + 1.0/tau));
  }

  return  beta;
}

Eigen::MatrixXd RegSampling::bmvreg_PSX(Eigen::MatrixXd Y,   Eigen::MatrixXd X, Eigen::MatrixXd beta,double mu) {
  // Generate ind_nod matrix
  Eigen::MatrixXi ind_nod(J - 1, J);
  for (int j = 0; j < J; ++j) {
    int index = 0;
    for (int i = 0; i < J; ++i) {
      if (i != j) {
        ind_nod(index++, j) = i ; //different with R (i+1), 0-based
      }
    }
  }

  double a_gams = 1.0;
  double b_gams = 0.1;

  // Calculation of S
  Eigen::MatrixXd temp = Y - Eigen::MatrixXd::Ones(Y.rows(), Y.cols()) * mu  - (X*beta.transpose());
  Eigen::MatrixXd S = temp.transpose() * temp;

  // Sample gammas
  double apost = a_gams + (J * (J + 1)) / 2.0;
  double bpost = b_gams + (inv_psx.cwiseAbs().sum()) / 2.0;

  std::mt19937 gen;
  std::gamma_distribution<> gamma(apost, 1.0 / bpost);
  double gammas  = gamma(gen);
  // Sample tau off-diagonal

  Eigen::VectorXd Cadj(J * (J - 1) / 2);
  int index = 0;
  for (int j = 0; j < J; ++j) {
    for (int i = 0; i < j; ++i) {
      Cadj(index++) = std::max(std::abs(inv_psx(i, j)), 1e-6);
    }
  }
  Eigen::VectorXd mu_p_tmp = gammas / Cadj.array();

  mu_p_tmp = mu_p_tmp.array().min(1e12); // pmin equivalent
  Rcpp::NumericVector mu_p = Rcpp::wrap(mu_p_tmp);

  double gammas_p = std::pow(gammas, 2);
  Rcpp::NumericVector taus_tmp = 1.0/ rinvgauss(mu_p_tmp.size(), mu_p, gammas_p);
  Eigen::MatrixXd taus = Eigen::MatrixXd::Zero(J, J);
  index = 0;

  for (int j = 0; j < J; ++j) {
    for (int i = 0; i < j; ++i) {
      // Assign taus_tmp values to both upper and lower triangular parts
      taus(i, j) = taus_tmp(index);
      taus(j, i) = taus_tmp(index);
      ++index;
    }
  }

  for (int i = 0; i < J; ++i) {
    Eigen::VectorXi ind_noi = ind_nod.col(i);
    // Manually constructing Sig11 as a submatrix of 'psx'
    int size = ind_noi.size();
    Eigen::MatrixXd Sig11(size, size);
    for(int row = 0; row < size; ++row) {
      for(int col = 0; col < size; ++col) {
        Sig11(row, col) = psx(ind_noi(row), ind_noi(col));
      }
    }

    // Extracting Sig12 as a subvector of column 'i' of 'psx'
    Eigen::VectorXd Sig12(size);
    for(int row = 0; row < size; ++row) {
      Sig12(row) = psx(ind_noi(row), i);
    }

    Eigen::MatrixXd invC11 = Sig11 - (Sig12 * Sig12.transpose()) / psx(i, i);

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
      inv_psx(ind_noi(idx), i) = beta(idx,0);
      inv_psx(i, ind_noi(idx)) = beta(idx,0);  // Assuming symmetry
    }


    std::mt19937 gen;
    std::gamma_distribution<double> gamma(N / 2.0 + 1,2.0/(S(i, i) + gammas));
    double gam  = gamma(gen);

    Eigen::MatrixXd tmp = beta.transpose() * invC11 * beta;
    inv_psx(i, i) = gam + tmp(0,0);


    // Updating covariance matrix 'psx' according to one-column change in the precision matrix 'inv_psx'
    Eigen::MatrixXd invC11beta = invC11 * beta;
    Sig12 = (-invC11beta / gam).col(0); //1 column
    Eigen::MatrixXd sub_psx =   invC11 + invC11beta * invC11beta.transpose()/gam;

    for (int idx = 0; idx < size; ++idx) {
      int row = ind_noi(idx);

      // Update off-diagonal elements for column 'i' and symmetric positions
      psx(row, i) = Sig12(idx);
      psx(i, row) = Sig12(idx);

      for (int jdx = 0; jdx < size; ++jdx) {
        int col = ind_noi(jdx);
        // Add invC11beta contribution to both off-diagonal and diagonal elements
        psx(row, col) = sub_psx(idx,jdx);
      }
    }

    psx(i, i) = 1.0 / gam;
    dsig[i] = 1.0 / gam;
  }
  Eigen::MatrixXd new_psx = psx;
  return new_psx;
}
Eigen::MatrixXd RegSampling::bmvreg_PSX_ssp(Eigen::MatrixXd Y,   Eigen::MatrixXd X, Eigen::MatrixXd beta,double mu) {
  //Y is N*J
  //X is fac:N*K
  //beta is loading: J*K

  // Generate ind_nod matrix
  Eigen::MatrixXi ind_nod(J - 1, J);
  for (int j = 0; j < J; ++j) {
    int index = 0;
    for (int i = 0; i < J; ++i) {
      if (i != j) {
        ind_nod(index++, j) = i ; //different with R (i+1), 0-based
      }
    }
  }

  // Calculation of S
  Eigen::MatrixXd temp = Y - Eigen::MatrixXd::Ones(Y.rows(), Y.cols()) * mu  - (X*beta.transpose());
  Eigen::MatrixXd S = temp.transpose() * temp;


  for (int i = 0; i < J; ++i) {
    Eigen::VectorXi ind_noi = ind_nod.col(i);
    // Manually constructing Sig11 as a submatrix of 'psx'
    int size = ind_noi.size();
    Eigen::MatrixXd Sig11(size, size);
    for(int row = 0; row < size; ++row) {
      for(int col = 0; col < size; ++col) {
        Sig11(row, col) = psx(ind_noi(row), ind_noi(col));
      }
    }

    // Extracting Sig12 as a subvector of column 'i' of 'psx'
    Eigen::VectorXd Sig12(size);
    for(int row = 0; row < size; ++row) {
      Sig12(row) = psx(ind_noi(row), i);
    }

    Eigen::MatrixXd invC11 = Sig11 - (Sig12 * Sig12.transpose()) / psx(i, i);

    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Zero(size, size);
    for(int idx = 0; idx < size; ++idx) {
      diagMatrix(idx, idx) = 1.0 / taus(ind_noi(idx), i); // Assuming 'i' is defined elsewhere as the column index
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
      inv_psx(ind_noi(idx), i) = beta(idx,0);
      inv_psx(i, ind_noi(idx)) = beta(idx,0);  // Assuming symmetry
    }


    std::mt19937 gen;
    std::gamma_distribution<double> gamma(N / 2.0 + 1,2.0/(S(i, i) + lambda));
    double gam  = gamma(gen);

    Eigen::MatrixXd tmp = beta.transpose() * invC11 * beta;
    inv_psx(i, i) = gam + tmp(0,0);


    // Updating covariance matrix 'psx' according to one-column change in the precision matrix 'inv_psx'
    Eigen::MatrixXd invC11beta = invC11 * beta;
    Sig12 = (-invC11beta / gam).col(0); //1 column
    Eigen::MatrixXd sub_psx =   invC11 + invC11beta * invC11beta.transpose()/gam;

    Eigen::VectorXd v0(size);
    Eigen::VectorXd v1(size);
    for (int idx = 0; idx < size; ++idx) {
      int row = ind_noi(idx);

      // Update off-diagonal elements for column 'i' and symmetric positions
      psx(row, i) = Sig12(idx);
      psx(i, row) = Sig12(idx);

      for (int jdx = 0; jdx < size; ++jdx) {
        int col = ind_noi(jdx);
        // Add invC11beta contribution to both off-diagonal and diagonal elements
        psx(row, col) = sub_psx(idx,jdx);
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
        if (z(j) == 1) {
          taus(ind_noi[j], i) = v1(j);
          taus(i, ind_noi[j]) = v1(j);
        }else{
          taus(ind_noi[j], i) = v0(j);
          taus(i, ind_noi[j]) = v0(j);
        }
      }
    }

    psx(i, i) = 1.0 / gam;
    dsig[i] = 1.0 / gam;
  }

  return psx;
}
Eigen::MatrixXd RegSampling::bmvreg_PSX_horse(Eigen::MatrixXd Y,   Eigen::MatrixXd X, Eigen::MatrixXd beta,double mu) {
  //Y is N*J
  //X is fac:N*K
  //beta is loading: J*K

  // Generate ind_nod matrix
  Eigen::MatrixXi ind_nod(J - 1, J);
  for (int j = 0; j < J; ++j) {
    int index = 0;
    for (int i = 0; i < J; ++i) {
      if (i != j) {
        ind_nod(index++, j) = i ; //different with R (i+1), 0-based
      }
    }
  }

  // Calculation of S
  Eigen::MatrixXd temp = Y - Eigen::MatrixXd::Ones(Y.rows(), Y.cols()) * mu  - (X*beta.transpose());
  Eigen::MatrixXd S = temp.transpose() * temp;

  //sample tau_sq and xi

  // Calculate rate
  double rate = 1.0 / xi;
  for (int i = 1; i < J; ++i) {
    for (int j = 0; j < i; ++j) {
      rate += inv_psx(i, j) * inv_psx(i, j) / (2.0 * Lambda_sq(i, j));
    }
  }
  // Gamma distribution for tau_sq and xi
  std::mt19937 gen;
  std::gamma_distribution<double> gamma_tau_sq(J * (J + 1) / 2.0, 1.0/rate);
  double tau_sq = gamma_tau_sq(gen); // Sample from gamma distribution
  std::gamma_distribution<double> gamma_xi(1.0, 1.0 + 1.0 / tau_sq);
  xi = gamma_xi(gen); // Sample from gamma distribution

  // Set lamdiag based on regdiag flag
  double lamdiag = regdiag ? tau_sq : 0.0;


  for (int i = 0; i < J; ++i) {
    Eigen::VectorXi ind_noi = ind_nod.col(i);
    // Manually constructing Sig11 as a submatrix of 'psx'
    int size = ind_noi.size();
    Eigen::MatrixXd Sig11(size, size);
    for(int row = 0; row < size; ++row) {
      for(int col = 0; col < size; ++col) {
        Sig11(row, col) = psx(ind_noi(row), ind_noi(col));
      }
    }

    // Extracting Sig12 as a subvector of column 'i' of 'psx'
    Eigen::VectorXd Sig12(size);
    for(int row = 0; row < size; ++row) {
      Sig12(row) = psx(ind_noi(row), i);
    }

    Eigen::MatrixXd invC11 = Sig11 - (Sig12 * Sig12.transpose()) / psx(i, i);

    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Zero(size, size);
    for(int idx = 0; idx < size; ++idx) {
      diagMatrix(idx, idx) = 1.0 / (Lambda_sq(ind_noi(idx), i)* tau_sq); // Assuming 'i' is defined elsewhere as the column index
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
      inv_psx(ind_noi(idx), i) = beta(idx,0);
      inv_psx(i, ind_noi(idx)) = beta(idx,0);  // Assuming symmetry
    }


    std::mt19937 gen;
    std::gamma_distribution<double> gamma(N / 2.0 + 1,2.0/(S(i, i) + lamdiag));
    double gam  = gamma(gen);

    Eigen::MatrixXd tmp = beta.transpose() * invC11 * beta;
    inv_psx(i, i) = gam + tmp(0,0);


    // Updating covariance matrix 'psx' according to one-column change in the precision matrix 'inv_psx'
    Eigen::MatrixXd invC11beta = invC11 * beta;
    Sig12 = (-invC11beta / gam).col(0); //1 column
    Eigen::MatrixXd sub_psx =   invC11 + invC11beta * invC11beta.transpose()/gam;

    for (int idx = 0; idx < size; ++idx) {
      int row = ind_noi(idx);

      // Update off-diagonal elements for column 'i' and symmetric positions
      psx(row, i) = Sig12(idx);
      psx(i, row) = Sig12(idx);

      for (int jdx = 0; jdx < size; ++jdx) {
        int col = ind_noi(jdx);
        // Add invC11beta contribution to both off-diagonal and diagonal elements
        psx(row, col) = sub_psx(idx,jdx);
      }
    }

    double beta_product = (beta.transpose() * beta)(0, 0);
    for (int j = 0; j < size; ++j) {
      std::gamma_distribution<double> gamma_sampled_lambda(1.0,(beta_product / (2.0 * tau_sq)) + 1.0 / Nu(ind_noi(j), i)); // Shape = 1, Scale = 1/rate
      Lambda_sq(i, ind_noi(j)) = gamma_sampled_lambda(gen);
      Lambda_sq(ind_noi(j), i) = Lambda_sq(i, ind_noi(j));
      std::gamma_distribution<double> gamma_Nu(1.0, (1.0 + 1.0 / Lambda_sq(i, ind_noi(j)))); // Shape = 1, Scale = rate_nu
      Nu(i, ind_noi(j)) = gamma_Nu(gen);
      Nu(ind_noi(j), i) = Nu(i, ind_noi(j)); // Update symmetric entry
    }
    psx(i, i) = 1.0 / gam;
    dsig[i] = 1.0 / gam;
  }
  return psx;
  //return Rcpp::List::create(Rcpp::Named("obj") = psx_input, Rcpp::Named("inv") = inv_psx_input, Rcpp::Named("gammas") = gammas);
}
