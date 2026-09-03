#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include<random>
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
Eigen::MatrixXd rwish1(double v, Eigen::MatrixXd S);
Rcpp::NumericVector rinvgauss(int n, Rcpp::NumericVector mu, double lambda);
Eigen::MatrixXd  selectColumns (Eigen::MatrixXd X, Eigen::MatrixXi sub1, int i);
Eigen::MatrixXd  submatrix (Eigen::MatrixXd C,Eigen::VectorXi row, Eigen::VectorXi col);
Eigen::MatrixXd hpd_interval(Eigen::MatrixXd mc_samples, double prob = 0.95 );
Eigen::VectorXd compute_psrf(Eigen::MatrixXd chain) ;