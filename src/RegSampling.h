#ifndef REG_SAMPLING_H
#define REG_SAMPLING_H


#include <RcppArmadillo.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]



class RegSampling {
  public:
    double s0,r0,lam_a,lam_b,m0,c0;
    int N, P, J, jp1, j1, jp2, j2;
    Eigen::MatrixXi pos1, pos2, sub1, sub2;
    Eigen::VectorXi jnd1, jnd2;
    Eigen::MatrixXd tausq,  msig,psx,inv_psx;
    Eigen::VectorXd r1_las,vtausq,dsig;
    double lamsq;
    bool LD;
    // init for ssp psx
    double pii,lambda;
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXd taus;
    // init for horse psx
    Eigen::MatrixXd Lambda_sq, Nu;
    double xi;
    bool regdiag;

    std::string regloading;
    //init horse loading
    Eigen::MatrixXd v;
    double tau,eta;

    //init ssp loading
    Eigen::VectorXd sksq;
    Eigen::VectorXd sksq_t;

    //global for passing to cati
    Eigen::VectorXd psxjs;
    Eigen::MatrixXd mst;

    RegSampling(int N,bool LD=false, double s0=1, double r0=0.1 , double lam_a=1 , double lam_b=0.1, double m0=0, double c0=0.1);
    Eigen::MatrixXd initial(int N, Eigen::MatrixXi Q,std::string reglo);
    Eigen::MatrixXd bmvreg_plas(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd beta);
    Eigen::MatrixXd bmvreg_plas_LD(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd beta);
    Eigen::MatrixXd bmvreg_PSX(Eigen::MatrixXd Y,  Eigen::MatrixXd X, Eigen::MatrixXd beta, double mu=0);
    Eigen::MatrixXd bmvreg_PSX_ssp(Eigen::MatrixXd Y,  Eigen::MatrixXd X, Eigen::MatrixXd beta,double mu=0);
    Eigen::MatrixXd bmvreg_PSX_horse(Eigen::MatrixXd Y,  Eigen::MatrixXd X, Eigen::MatrixXd beta,double mu=0);


  };

#endif // REG_SAMPLING_H
