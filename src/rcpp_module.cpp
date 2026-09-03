
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "PCPM.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// Rcpp Module definition
//RCPP_EXPOSED_CLASS(PCPM)

RCPP_MODULE(PCPM) {
  using namespace Rcpp;
  class_<PCPM>("PCPM")
    .constructor<std::vector<double>,std::vector<double>,std::vector<double>>()
    //.constructor()
    .field("mfac", &PCPM::mfac) //parameter public to user
    .method("pcfa", &PCPM::pcfa)
    .method("stat", &PCPM::stat)
    .method("getVariable", &PCPM::getVariable);
}


