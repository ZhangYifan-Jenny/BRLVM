#include "RegSampling.h"
#include "FacSampling.h"
#include "utility-function.h"
#include "PCPM.h"
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]


void PCPM::pcfa(Eigen::MatrixXd Y,Eigen::MatrixXi Q,bool cati,Eigen::VectorXi cati_v, std::string reglo,bool LD,std::string regphi,int Kg,Eigen::MatrixXd reg,std::string regpsx,
                int iter , int burn,bool std,
                Rcpp::NumericVector cor_fac,bool sign_check, int update ) {
  if (Q.rows() != Y.cols()) throw std::runtime_error("The numbers of items in Y and Q are unequal.");
  if (iter <= 0) throw std::runtime_error("iter must be larger than zero.");

  int N = Y.rows(), J = Y.cols(), K = Q.cols();
  ////missing
  int Nmis = 0; //# of missing indicates
  std::vector<int> mis_rows;
  std::vector<int> mis_cols;
  mis_rows.reserve(J*N);
  mis_cols.reserve(J*N);
  // Iterate through the matrix to find NaN entries
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < J ; ++j) {
      if (std::isnan(Y(i, j))) {
        Nmis++;
        mis_rows.push_back(i);
        mis_cols.push_back(j);
      }
    }
  }

  Rcpp::IntegerMatrix mind(mis_rows.size(), 2);
  for (size_t k = 0; k < mis_rows.size(); ++k) {
    mind(k, 0) = mis_rows[k];
    mind(k, 1) = mis_cols[k];
  }

  int mnoc;
  int Jp;
  double inf = 1e+200;
  Eigen::MatrixXd zind;
  Eigen::MatrixXd THD;
  Eigen::MatrixXd ys;
  if (cati) {
    Jp = cati_v.size();

    if (Jp == 1 && cati_v[0] == -1) {
      cati_v.resize(J);
      for (int i = 0; i < J; ++i) {
        cati_v[i] = i;
      }
      Jp = J;
    }

    if(Jp>0){

      // Identify categories and relabel
      Eigen::VectorXi noc =Eigen::VectorXi::Zero(Jp);
      zind.resize(N,Jp);
      for (int j = 0; j < Jp; j++) { //transfer index start from 1 get 'zind'
        std::vector<double> values;
        values.reserve(N); //reserves capacity
        for (int i = 0; i < N; i++) {
          double val = Y(i, cati_v(j));
          if (!std::isnan(val)) values.push_back(val);
        }

        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end()); //distinct categories
        noc[j] = values.size(); // number of categories
        // relabel categories as 1,2,…,noc[j]
        for (int m = 0; m < noc[j]; m++) {

          for (int i = 0; i < N; i++) {
            if (!std::isnan(Y(i, cati_v(j))) && Y(i, cati_v(j)) == values[m]){
              zind(i, j) = m + 1;
            }else if (std::isnan(Y(i, cati_v(j)))){
              zind(i, j) = NAN;
            }
          }
        }

      }
      mnoc = noc.maxCoeff();

      THD = Eigen::MatrixXd::Zero(Jp, mnoc + 1);

      THD.leftCols(1).array() = -inf;

      for (int j = 0; j < Jp; ++j) {
        THD(j, 0) = -inf; // Set the first column to -inf
        for (int k = 1; k < noc(j); ++k) {
          THD(j, k) = k - 1; // Filling 0 to noc[j]-2
        }
        THD(j, noc(j)) = inf; // Set the last valid category threshold to inf

        if (noc(j) < mnoc) {
          THD.block(j, noc(j) + 1, 1, mnoc - noc(j)) = Eigen::VectorXd::Constant(mnoc - noc(j), inf);
        }
      }
      arma::mat Sigma = arma::eye<arma::mat>(Jp, Jp);
      arma::vec armaMean = arma::zeros(Jp);
      arma::mat ys_sample = arma::mvnrnd(armaMean, Sigma, N);
      ys = Eigen::Map<Eigen::MatrixXd>(ys_sample.memptr(), ys_sample.n_rows, ys_sample.n_cols).transpose();

      for (int i = 0; i < Jp; i++) {
        Y.col(cati_v[i]) = ys.col(i);
      }

      std::vector<bool> isIndexed(Jp, true); // Track whether each column is indexed

      // Mark indexed columns
      for (int i = 0; i < Jp; ++i) {
        isIndexed[cati_v[i]] = false;
      }

      // Allocate matrix for non-indexed columns
      Eigen::MatrixXd ycs(Y.rows(), J-Jp);
      int colIdx = 0;

      // Fill matrix with non-indexed columns
      for (int j = 0; j < J; ++j) {
        if (isIndexed[j]) {
          ycs.col(colIdx++) = Y.col(j);
        }
      }

      //SCALE ycs
      Eigen::VectorXd means(J-Jp);
      Eigen::VectorXd stdDevs(J-Jp);

      for (int j = 0; j < J-Jp; ++j) {
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < N; ++i) {
          if (!std::isnan(ycs(i, j))) {
            sum += ycs(i, j);
            ++count;
          }
        }
        if (count > 0) {
          double mean = sum / count;
          means(j) = mean;
          double sumSquaredDiff = 0.0;
          for (int i = 0; i < N; ++i) {
            if (!std::isnan(ycs(i, j))) {
              double diff = ycs(i, j) - mean;
              sumSquaredDiff += diff * diff;
            }
          }
          stdDevs(j) = sqrt(sumSquaredDiff / (count - 1));

        } else {
          means(j) = NAN;
          stdDevs(j) = NAN;
        }
      }

      // Subtract mean and divide by standard deviation
      for (int j = 0; j < J-Jp; ++j) {
        for (int i = 0; i < N; ++i) {
          if (!std::isnan(ycs(i, j))) {
            ycs(i, j) -= means(j);
            if (std) {
              ycs(i, j) /= stdDevs(j);
            }

          }
        }
      }
      colIdx = 0;
      for (int j = 0; j < J; ++j) {
        if (isIndexed[j]) {
          Y.col(j) =  ycs.col(colIdx++);
        }
      }
    }

  }else{ //  no cati
    Eigen::VectorXd means(J);
    Eigen::VectorXd stdDevs(J);

    for (int j = 0; j < J; ++j) {
      double sum = 0.0;
      int count = 0;
      for (int i = 0; i < N; ++i) {
        if (!std::isnan(Y(i, j))) {
          sum += Y(i, j);
          ++count;
        }
      }
      if (count > 0) {
        double mean = sum / count;
        means(j) = mean;
        double sumSquaredDiff = 0.0;
        for (int i = 0; i < N; ++i) {
          if (!std::isnan(Y(i, j))) {
            double diff = Y(i, j) - mean;
            sumSquaredDiff += diff * diff;
          }
        }
        stdDevs(j) = sqrt(sumSquaredDiff / (count - 1));//N or N-1 //lawbl use N-1, python use N
      } else {
        means(j) = NAN;
        stdDevs(j) = NAN;
      }
    }

    // Subtract mean and divide by standard deviation
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < N; ++i) {
        if (!std::isnan(Y(i, j))) {
          Y(i, j) -= means(j);
          if (std) {
            Y(i, j) /= stdDevs(j);
          }

        }
      }
    }
  }

  if(Nmis>0){
    Rcpp::NumericVector miss_sample =Rcpp::rnorm(Nmis);
    for (size_t k = 0; k < Nmis; ++k) {
      Y(mis_rows[k],mis_cols[k]) = miss_sample[k];
    }
  }

  A_arr = Eigen::MatrixXd::Zero(iter, J * K);
  if(LD){
    V_arr = Eigen::MatrixXd::Zero(iter, J*(J+1)/2); //EPSX in R
  }else{
    V_arr = Eigen::MatrixXd::Zero(iter, J);
  }

  lamA_arr = Eigen::VectorXd::Zero(iter);
  C_arr = Eigen::MatrixXd::Zero(iter, K * (K + 1) / 2);

  mip_arr = Eigen::MatrixXd::Zero(iter, K * (K + 1) / 2);

  // Initialize sampling classes, passing parameters as needed

  FacSampling samF(N, K,Kg, cor_fac, std,w0);
  //RegSampling   (s0=1, r0=0.1, lam_a=1, lam_b=0.1, m0=0, c0=0.1)
  RegSampling samA(N,LD,s0, r0, lam_a, lam_b, m0, c0);

  Eigen::MatrixXd A = samA.initial(N, Q,reglo); //J*K
  Eigen::MatrixXd sfac= Eigen::MatrixXd::Zero(N,  K);

  Eigen::MatrixXd psx = samA.psx;
  Eigen::MatrixXd inv_psx = samA.inv_psx;
  Eigen::VectorXd dV;
  Eigen::MatrixXd fac; //N*K
  Eigen::MatrixXd C;
  //Eigen::MatrixXd mfac;


  auto start = std::chrono::high_resolution_clock::now();
  // Gibbs sampler loop
  for (int i = 0; i < iter + burn; ++i) {
    dV = samA.dsig;
    fac = samF.gibbs_factor_LD(Y, A, inv_psx);
    if (samF.K_rf == K) {
      // No uncorrelated factor
      if(regphi=="lasso"){
        C = samF.cor_cov_lasso(fac, std,reg);
      }else if(regphi=="parlasso"){
        C = samF.cor_cov_lasso_par(fac, std);
      }else if(regphi=="ssp"){
        C = samF.cor_cov_ssp(fac, std);
      }else if(regphi=="parssp"){
        C = samF.cor_cov_ssp_par(fac, std);
      }else if(regphi=="horse"){
        C = samF.cor_cov_horse(fac, std);
      }else if(regphi=="parhorse"){
        C = samF.cor_cov_horse_par(fac, std);
      }else{
        C = samF.cor_cov(fac, std);
      }
    } else {
      C = samF.cor_cov1(fac, std);
    }

    if (LD){
      if(regpsx=="ssp"){
        psx = samA.bmvreg_PSX_ssp(Y,fac,A,0.0);
        inv_psx = samA.inv_psx;
      }else if(regpsx=="horse"){
        psx = samA.bmvreg_PSX_horse(Y,fac,A,0.0);
        inv_psx = samA.inv_psx;
      }else{
        psx = samA.bmvreg_PSX(Y,fac,A,0.0);
        inv_psx = samA.inv_psx;
      }
      A = samA.bmvreg_plas_LD(Y, fac, A);
    }else{
      A = samA.bmvreg_plas(Y, fac, A);

      psx = samA.psx;
      //inv_psx = samA.inv_psx;
      inv_psx = psx.llt().solve(Eigen::MatrixXd::Identity(J, J));
    }

    Eigen::MatrixXd ysa;
    if (Nmis > 0 || cati ){
      Eigen::MatrixXd ysta(N,J);
      Eigen::VectorXd spsxa(J);
      if(LD){
        ysta=fac * A.transpose() + samA.mst;
        spsxa = samA.psxjs.array().sqrt();
      }else{
        ysta =fac * A.transpose();
        spsxa = dV.array().sqrt();
      }
      //Eigen::MatrixXd ysta =fac * A.transpose();
      //Eigen::VectorXd spsxa = dV.array().sqrt();
      arma::mat Random_norm  = arma::randn(N,J);
      ysa = Eigen::Map<Eigen::MatrixXd>(Random_norm.memptr(), Random_norm.n_rows, Random_norm.n_cols);
      for (int j = 0; j < J; ++j) {
        ysa.col(j) += ysta.col(j) / spsxa(j);
      }
      for (int j = 0; j < J; ++j) {
        double sd = sqrt((ysa.col(j).array() - ysa.col(j).mean()).square().sum() / (N-1)); //N or N-1?????
        ysa.col(j) /= sd;
      }

      if(cati){
        Eigen::MatrixXd td1 = THD;
        for (int m = 1; m < mnoc; m++) {
          for (int i = 0; i < Jp; i++) {
            double mean = THD(i, m);
            double sd = 0.2;
            double tmp = mean + sd * as_scalar(arma::randn<arma::mat>(1, 1));
            bool sel = (tmp > td1(i, m - 1)) && (tmp <= THD(i, m + 1));
            td1(i, m) = tmp * sel + THD(i, m) * (1 - sel);
          }
        }
        //Rcpp::Rcout  << "td1 is: "<< td1 << std::endl;
        Eigen::MatrixXd tmp00(N, Jp);
        Eigen::MatrixXd tmp01(N, Jp);
        Eigen::MatrixXd tmp10(N, Jp);
        Eigen::MatrixXd tmp11(N, Jp);
        for (int i = 0; i < N; ++i) {
          for (int j = 0; j < Jp; ++j) {
            if (std::isnan(zind(i, j))) {
              // Apply NaNs to all matrices at once
              tmp00(i, j) = std::numeric_limits<double>::quiet_NaN();
              tmp01(i, j) = std::numeric_limits<double>::quiet_NaN();
              tmp10(i, j) = std::numeric_limits<double>::quiet_NaN();
              tmp11(i, j) = std::numeric_limits<double>::quiet_NaN();
            }else{
              int zind_ij = zind(i, j);
              tmp00(i, j) = THD(j, zind_ij-1);
              tmp10(i, j) = td1(j, zind_ij-1);
              tmp01(i, j) = THD(j, zind_ij);
              tmp11(i, j) = td1(j, zind_ij);
            }

          }
        }

        Eigen::MatrixXd yst(N, Jp);
        Eigen::VectorXd spsx(Jp);
        Eigen::MatrixXd tmp1(N, Jp);
        Eigen::MatrixXd tmp0(N, Jp);


        for (int i = 0; i < Jp; ++i) {
          yst.col(i) = ysta.col(cati_v[i]); // Copy entire row based on index
          spsx(i) = spsxa(cati_v[i]);
        }

        for (int i = 0; i < N; i++) {
          for (int j = 0; j < Jp; j++) {

            std::mt19937 gen;
            double p_upper = 0.5 * std::erfc(-((tmp11(i, j) - yst(i, j)) / spsx(j)) / std::sqrt(2));
            double p_lower = 0.5 * std::erfc(-((tmp10(i, j) - yst(i, j)) / spsx(j)) / std::sqrt(2));

            //double p_upper = R::pnorm((tmp11(i, j) - yst(i, j)) / spsx(j), 0, 1, 1, 0);
            //double p_lower = R::pnorm((tmp10(i, j) - yst(i, j)) / spsx(j), 0, 1, 1, 0);
            tmp1(i, j) = log(p_upper - p_lower);
            if (tmp1(i, j) < -inf) {
              tmp1(i, j) = -inf;
            }

            p_upper = 0.5 * std::erfc(-((tmp01(i, j) - yst(i, j)) / spsx(j)) / std::sqrt(2));
            p_lower = 0.5 * std::erfc(-((tmp00(i, j) - yst(i, j)) / spsx(j)) / std::sqrt(2));
            //p_upper = R::pnorm((tmp01(i, j) - yst(i, j)) / spsx(j), 0, 1, 1, 0);
            //p_lower = R::pnorm((tmp00(i, j) - yst(i, j)) / spsx(j), 0, 1, 1, 0);
            tmp0(i, j) = log(p_upper - p_lower);
            if (tmp0(i, j) < -inf) {
              tmp0(i, j) = -inf;
            }
          }
        }

        Eigen::VectorXd acc(Jp);
        Eigen::VectorXi accind = Eigen::VectorXi::Zero(Jp);
        for (int j = 0; j < Jp; ++j) {
          double sum = 0.0;
          int count = 0;
          for (int i = 0; i < N; ++i) {
            if (std::isnan(tmp1(i, j)) || std::isnan(tmp0(i, j))) {
              continue;
            }
            sum += (tmp1(i, j) - tmp0(i, j));
            count++;
          }
          if (count > 0) {
            acc(j) = exp(sum);
            std::mt19937 gen; // Seed the generator
            // Uniform distribution
            std::uniform_real_distribution<> dis_unif(0.0, 1.0); // Uniform distribution in [0, 1)
            double uniform_sample = dis_unif(gen);
            if(acc(j) > uniform_sample){
              //Rcpp::Rcout  << "true "<< j << std::endl;
              accind[j] = 1;
              THD.row(j) = td1.row(j);
            }
          } else {
            acc(j) = NAN;  // Assign NA if all were NA
          }
        }

        for (int i = 0; i < N; ++i) {
          for (int j = 0; j < Jp; ++j) {
            if (std::isnan(zind(i, j))) {
              // Apply NaNs to all matrices at once
              tmp00(i, j) = std::numeric_limits<double>::quiet_NaN();
              tmp01(i, j) = std::numeric_limits<double>::quiet_NaN();
            }else{
              int zind_ij = zind(i, j);
              tmp00(i, j) = THD(j, zind_ij-1);
              tmp01(i, j) = THD(j, zind_ij);
            }

          }
        }

        for (int i = 0; i < Jp; ++i) {
          ys.col(i) = ysa.col(cati_v[i]);
        }
        Eigen::MatrixXd acc1(N,Jp);
        for (int i = 0; i < N; ++i) {
          for (int j = 0; j < Jp; ++j) {
            if(!std::isnan(tmp00(i, j)) && !std::isnan(tmp01(i, j))){
              bool condition = (ys(i, j) > tmp00(i, j)) && (ys(i, j) <= tmp01(i, j));
              acc1(i, j) = condition ? 1 : 0;
            }else{
              acc1(i, j) = NAN;
            }
          }
        }

        // Applying the condition to update ys
        for (int i = 0; i < N; ++i) {
          for (int j = 0; j < Jp; ++j) {
            if (acc1(i, j) == 0) {
              ys(i, j) = Y(i , cati_v[j]);  // Original value from y
            }
          }
        }


        // Calculating means for accr
        Eigen::VectorXd accr(2);
        accr[0] = accind.cast<double>().mean();  // Mean of acc1, treating true as 1, false as 0

        //Rcpp::Rcout  << "accr[0] is: "<<accr[0] << std::endl;

        int acc1count = 0;
        double acc1sum = 0.0;

        for (int i = 0; i < N; ++i) {
          for (int j = 0; j < Jp; ++j) {
            if (!std::isnan(acc1(i,j))) {  // Check if not NaN
              acc1sum += acc1(i,j);
              ++acc1count;
            }
          }
        }


        accr[1] = acc1sum / acc1count;
        for (int i = 0; i < Jp; i++) {
          Y.col(cati_v[i]) = ys.col(i);
        }
      }



      if (Nmis > 0){
        for (size_t k = 0; k < Nmis; ++k) {
          Y(mis_rows[k],mis_cols[k]) = ysa(mis_rows[k],mis_cols[k]) ;
        }
      }



    }

    if (sign_check) {
    Eigen::VectorXd colSums = A.colwise().sum();
    Eigen::Array<bool, Eigen::Dynamic, 1> chg = (colSums.array() <= -0.5);
    //int sign_sw = 0;
    if (chg.any()) {
      for (int j = 0; j < K; ++j) {
        if (chg(j)) {
          A.col(j) = -A.col(j); // Flip the sign of the column directly
          fac.col(j) = -fac.col(j);
        }
      }
        }
    }

    if (i >= burn) {
      Eigen::Map<Eigen::VectorXd> A_flattened(A.transpose().data(), A.transpose().size());
      A_arr.row(i - burn)= A_flattened;
      if(LD){
        Eigen::VectorXd psx_lower_tri(J*(J+1)/2);
        int index = 0;
        for (int i = 0; i < J; ++i) {
          for (int j = 0; j <= i; ++j) {
            psx_lower_tri(index++) = psx(i, j);
          }
        }
        V_arr.row(i - burn)=psx_lower_tri;
      }else{V_arr.row(i - burn)= psx.diagonal(); }

      lamA_arr(i - burn) = std::sqrt(samA.lamsq);
      sfac += fac;
      int index = 0;
      for (int col = 0; col < C.cols(); ++col) {
        for (int row = col; row < C.rows(); ++row) {
          C_arr(i - burn, index++) = C(row, col);
        }
      }
      int index_mip = 0;
      for (int col = 0; col < C.cols(); ++col) {
        for (int row = col; row < C.rows(); ++row) {
          mip_arr(i - burn, index_mip++) = samF.mip(row, col);
        }
      }


    }

    if ((i + 1) % update == 0) {
      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = finish - start;
      Rcpp::Rcout << "Cumulative iterations: " << i + 1 << ", Time needed: " << elapsed.count() << " s\n";
    }
  }
  mfac = sfac/iter;
  //Rcpp::Rcout  << "mfac is: "<< mfac << std::endl;
}

Eigen::MatrixXd PCPM::stat(std::string parn, bool med , int start, int end , double SL, bool sig ) {
  auto getVariable=[this](std::string parn) {
    Eigen::MatrixXd parn_arr;
    if (parn == "A") {parn_arr=A_arr;}
    else if (parn == "V") {parn_arr=V_arr;}
    else if (parn == "C") {parn_arr= C_arr;}
    else if (parn == "B") {parn_arr=B_arr;}
    else if (parn == "U") {parn_arr= U_arr;}
    else if (parn == "H") {parn_arr=H_arr;}
    else if (parn == "Ap") {parn_arr= Ap_arr;}
    else if (parn == "lamA") {//lamA_arr is Eigen::VectorXd, don't use &
      parn_arr.resize(lamA_arr.size(), 1);
      parn_arr.col(0) = lamA_arr;}
    else if (parn == "lamAp") {
      parn_arr.resize(lamAp_arr.size(), 1);
      parn_arr.col(0) = lamAp_arr;}
    else if (parn == "lamB") {
      parn_arr.resize(lamB_arr.size(), 1);
      parn_arr.col(0) = lamB_arr;}
    else if (parn == "lamH") {
      parn_arr.resize(lamH_arr.size(), 1);
      parn_arr.col(0) = lamH_arr;}
    else {Rcpp::stop("Par " + parn + "_arr does not exist.");}
    return parn_arr;
  };

  Eigen::MatrixXd res;
  try {
    Eigen::MatrixXd dat0 = getVariable(parn);
    Eigen::MatrixXd dat;

    if (end == -1) end = dat0.rows();

    // Single variable case, ensure dat is a column matrix
    if (dat0.cols() == 1) {
      dat = dat0.block(start, 0, end - start, 1);
      double est;

      if (med){
        std::vector<double> data(dat.data(), dat.data() + dat.size());
        std::nth_element(data.begin(), data.begin() + data.size() / 2, data.end());
        est = data[data.size() / 2];
        }
      else{est = dat.mean();}

      Eigen::MatrixXd CI = hpd_interval(dat, 1 - SL);
      double sign = (CI(0, 0) * CI(0, 1) > 0) ? 1 : 0;
      Eigen::VectorXd psrf = compute_psrf(dat);
      res.resize(1, 5);
      res << est, CI(0, 0), CI(0, 1), sign, psrf(0);
      }
    // Multiple variables
    else {

      dat = dat0.block(start, 0, end - start, dat0.cols());
      Eigen::VectorXd est(dat.cols());
      if (med){
        for (int i = 0; i < dat.cols(); ++i) {
        std::vector<double> data(dat.col(i).data(), dat.col(i).data() + dat.col(i).size());
        std::nth_element(data.begin(), data.begin() + data.size() / 2, data.end());
        est(i) = data[data.size() / 2];
        }}
      else{
        for (int i = 0; i < dat.cols(); ++i) {
          est(i) = dat.col(i).mean();
          }}
      Eigen::MatrixXd CI = hpd_interval(dat, 1 - SL); // Adjust hpd_interval to handle matrix input
      Eigen::VectorXd sign = (CI.col(0).array() * CI.col(1).array() > 0).cast<double>();
      Eigen::VectorXd psrf = compute_psrf(dat); // Adjust compute_psrf to handle matrix input
      res.resize(dat.cols(), 6);
      for (int i = 0; i < dat.cols(); ++i) {
        res.row(i) << est(i), CI(i, 0), CI(i, 1), sign(i), psrf(i),i+1;
        }

      //only significant parameter will left
      if (sig) {
        int count = 0;
        for (int i = 0; i < res.rows(); ++i) {
          if (res(i, 3) == 1) {
            ++count;
            }
          }

        // Create a new matrix to hold the filtered results
        Eigen::MatrixXd subres(count, res.cols());

        // Fill the new matrix with rows that meet the condition
        int rowIndex = 0;
        for (int i = 0; i < res.rows(); ++i) {
          if (res(i, 3) == 1) {
            subres.row(rowIndex++) = res.row(i);
            }
          }
        res= subres;
        }
      }
    }
  catch (const std::exception& e) {
    Rcpp::Rcout << "Error: " << e.what() << std::endl;
    }
  return res;
}
Eigen::MatrixXd PCPM::getVariable (std::string parn){
    Eigen::MatrixXd parn_arr;
    if (parn == "A") {parn_arr=A_arr;}
    else if (parn == "V") {parn_arr=V_arr;}
    else if (parn == "C") {parn_arr= C_arr;}
    else if (parn == "B") {parn_arr=B_arr;}
    else if (parn == "U") {parn_arr= U_arr;}
    else if (parn == "H") {parn_arr=H_arr;}
    else if (parn == "Ap") {parn_arr= Ap_arr;}
    else if (parn == "mip") {parn_arr= mip_arr;}
    else if (parn == "lamA") {//lamA_arr is Eigen::VectorXd, don't use &
      parn_arr.resize(lamA_arr.size(), 1);
      parn_arr.col(0) = lamA_arr;}
    else if (parn == "lamAp") {
      parn_arr.resize(lamAp_arr.size(), 1);
      parn_arr.col(0) = lamAp_arr;}
    else if (parn == "lamB") {
      parn_arr.resize(lamB_arr.size(), 1);
      parn_arr.col(0) = lamB_arr;}
    else if (parn == "lamH") {
      parn_arr.resize(lamH_arr.size(), 1);
      parn_arr.col(0) = lamH_arr;}
    else if (parn == "mfac") {parn_arr= mfac;}
    else {Rcpp::stop("Par " + parn + "_arr does not exist.");}

    return parn_arr;
}



