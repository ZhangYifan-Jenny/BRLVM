#ifndef FAC_SAMPLING_H
#define FAC_SAMPLING_H



#include <RcppArmadillo.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]




class FacSampling {
public:

  double w0;
  int N,K;
  Eigen::MatrixXd W0, sub_W0,C;
  Eigen::MatrixXd phi,inv_phi;  //covariance matrix for reg phi
  Rcpp::NumericVector cor_fac;
  bool std;
  Eigen::VectorXi rf; // Assuming cor_fac is a vector of integers 0 or 1, or boolean values
  int K_rf;
  // init for ssp
  double pii,lambda,pii_par;
  Eigen::MatrixXd V0, V1;
  Eigen::MatrixXd taus;
  // init for horse
  Eigen::MatrixXd Lambda_sq, Nu;
  double xi;
  bool regdiag;

  // init for reg matrix - try partial phi lasso
  Eigen::MatrixXd C_gen, C_spe;
  int Kg;
  Eigen::MatrixXd phi_par,inv_phi_par;  //covariance matrix for partial reg phi


  // init for partial ssp
  Eigen::MatrixXd V0_par, V1_par;
  Eigen::MatrixXd taus_par;

  Eigen::MatrixXd  mip;

  FacSampling(int N, int K,int Kg, Rcpp::NumericVector cor_fac=1, bool std=TRUE, double w0=0.1);
  Eigen::MatrixXd gibbs_factor_LD(Eigen::MatrixXd Y, Eigen::MatrixXd A, Eigen::MatrixXd inv_psx);
  Eigen::MatrixXd cor_cov(Eigen::MatrixXd fac,bool std);
  Eigen::MatrixXd cor_cov1(Eigen::MatrixXd fac,bool std) ;


  Eigen::MatrixXd cor_cov_lasso(Eigen::MatrixXd fac, bool std,Eigen::MatrixXd reg);
  Eigen::MatrixXd cor_cov_ssp(Eigen::MatrixXd fac, bool std);
  Eigen::MatrixXd cor_cov_horse(Eigen::MatrixXd fac, bool std);
  Eigen::MatrixXd cor_cov_horse_par(Eigen::MatrixXd fac, bool std);
  Eigen::MatrixXd cor_cov_lasso_par(Eigen::MatrixXd fac, bool std);
  Eigen::MatrixXd cor_cov_ssp_par(Eigen::MatrixXd fac, bool std);
};

#endif // FAC_SAMPLING_H
