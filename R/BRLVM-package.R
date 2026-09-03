#' BRLVM: Bayesian Regularized Latent Variable Modeling
#'
#' @description
#' BRLVM provides tools for Bayesian regularization in latent variable modeling,
#' including regularization of factor loadings, factor correlation structures,
#' and residual covariance structures.
#'
#' @details
#' BRLVM was developed as the software implementation accompanying Yifan Zhang's
#' PhD thesis, \emph{Bayesian Regularization for Latent Variable Modeling:
#' A Unified Framework and Prior Comparison}.
#'
#' The package builds on the partially confirmatory factor analysis framework
#' implemented in the LAWBL package and extends it by supporting regularization
#' of factor correlation structures, including full and partial regularization.
#' It also provides multiple shrinkage priors, including lasso, horseshoe, and
#' spike-and-slab priors, which can be selected and combined across different
#' model components.
#'
#' @useDynLib BRLVM, .registration = TRUE
#'
#' @import Rcpp
#' @import methods
#' @import coda
#' @import clue
#' @import MASS
#'
#' @importFrom graphics hist par
#' @importFrom stats median rbinom sd var
#'
#' @keywords internal
#'
"_PACKAGE"
