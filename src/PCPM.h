
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::plugins("cpp11")]]
class PCPM{
public:
  Eigen::MatrixXd A_arr,V_arr,C_arr,Y,res, Ap_arr,B_arr,U_arr,H_arr,H,mfac,mip_arr;// mip_arr try mip for ssp
  Eigen::VectorXd lamA_arr,lamAp_arr,lamB_arr,lamH_arr ;
  int iter ,burn ,update,start,end,SL;
  double w0 = 0.1;
  double s0 = 1.0,  r0 = 0.1,  lam_a = 1.0,  lam_b = 0.1,  m0 = 0.0,  c0 =0.1;
  double s0B = 1.0,  r0B = 0.1,  lam_aB = 1.0,  lam_bB = 0.1,  m0B = 0.0,  c0B =0.1;
  //double w0; //priorF
  //double s0, r0, lam_a, lam_b, m0, c0 ; //prior A
  //double s0B, r0B, lam_aB, lam_bB, m0B, c0B ; //prior B
  bool std,med,sig,LD;
  std::vector<double> priorF,priorA,priorB,priorC;
  Rcpp::NumericVector cor_fac;
  std::string parn, regphi,regpsx;

  int distance = 1000, consecutive=100; //cati
  double threshold =0.0001;

  PCPM(std::vector<double> priorF,std::vector<double> priorA,std::vector<double> priorB)
    :priorF({0.1}),priorA({1.0, 0.1, 1.0, 0.1, 0.0, 0.1}),priorB({1.0, 0.1, 1.0, 0.1, 0.0, 0.1}){
    w0=priorF[0];
    s0=priorA[0];r0=priorA[1];lam_a =priorA[2];lam_b =priorA[3];m0=priorA[4];c0=priorA[5];
    s0B=priorB[0];r0B=priorB[1];lam_aB=priorB[2];lam_bB=priorB[3];m0B=priorB[4];c0B=priorB[5];
  }
  //PCPM(){}  //comment in package

  void pcfa(Eigen::MatrixXd Y,Eigen::MatrixXi Q,bool cati,Eigen::VectorXi cati_v, std::string reglo,bool LD,std::string regphi,int Kg,Eigen::MatrixXd reg,std::string regpsx,int iter, int burn ,bool std , Rcpp::NumericVector cor_fac,bool sign_check, int update=1000);
  //Eigen::MatrixXd stat(std::string parn, bool med = true, int start = 0, int end = -1, double SL = 0.05, bool sig = true);
  Eigen::MatrixXd stat(std::string parn, bool med, int start , int end , double SL , bool sig);
  Eigen::MatrixXd getVariable (std::string parn);
private:
  std::map<std::string, Rcpp::List> priors;

};
