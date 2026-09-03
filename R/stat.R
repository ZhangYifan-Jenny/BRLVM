#' Posterior Summary Statistics for Model Parameters
#'
#' @description
#' This function computes posterior summary statistics for model parameters
#' stored in an object containing MCMC samples. It supports different
#' parameter types (e.g., loadings, factor correlations, residual covariances) and
#' automatically constructs index labels based on the parameter structure.
#'
#' For each parameter element, the function computes:
#'   - Posterior mean or median
#'   - Posterior standard deviation
#'   - Highest Posterior Density (HPD) interval
#'   - Significance indicator (HPD interval excluding zero)
#'   - Potential scale reduction factor (PSRF)
#'
#' @param classname An object containing MCMC samples and a
#'   \code{getVariable()} method to extract parameter chains.
#' @param varName Character string indicating the parameter name
#' @param J Number of items.
#' @param K Number of latent factors.
#' @param P Number of predictors.
#' @param Q Loading pattern matrix for factor loadings.
#' @param Qb Loading pattern matrix for structure coefficients.
#' @param LD Logical. Whether local dependence structure is included.
#' @param med Logical. If TRUE, compute posterior median instead of mean.
#' @param start Starting iteration for posterior summary.
#' @param end Ending iteration for posterior summary.
#' @param SL Significance level for HPD interval (default 0.05).
#' @param sig Logical. If TRUE, return only significant parameters.
#'
#' @details
#' Let theta denote a model parameter. Posterior summaries are computed
#' from MCMC samples {theta^(s)} as:
#'
#'   Mean:        E(theta | data)
#'   SD:          sd(theta | data)
#'   HPD interval: 100(1 - SL)% credible interval
#'
#' A parameter is considered significant if its HPD interval does not
#' contain zero.
#'
#' Index positions are automatically constructed depending on the
#' parameter type:
#'   A: factor loadings
#'   C: factor correlation elements
#'   V: residual covariance elements
#'
#' @return
#' A matrix containing posterior summaries for each parameter element,
#' including estimates, standard deviations, HPD interval bounds,
#' significance indicator, PSRF, and corresponding index labels.
#'
#' @export
stat <- function(classname, varName,J=NULL,K=NULL, P=NULL,  Q=NULL,Qb=NULL, LD=F, med =FALSE,  start = 00,  end = -1,  SL = 0.05,  sig =F){
  dat0 <- classname$getVariable(varName)
  if (end == -1 | end > nrow(dat0)){
    end <- nrow(dat0)
  }

  ind<-1:ncol(dat0)
  if(!is.null(Q)){
    J <- nrow(Q)
    K <- ncol(Q)
    if (varName=='A'){
      ind <- which(!is.na(Q), arr.ind = TRUE)
      colnames(ind) <- c("Item", "Factor")
    }
  }

  if(!is.null(Qb)){
    K <- nrow(Qb)
    P <- ncol(Qb)
    if(varName=='B'){
      ind <- which(!is.na(t(Qb)), arr.ind = TRUE)
      colnames(ind) <- c("Predictor", "Factor")
    }
  }

  if(varName=='C'&& !is.null(K)){
    pos <- lower.tri(matrix(0, K, K),diag = T)
    ind <- which(pos, arr.ind = TRUE)
    colnames(ind) <- c("F", "F")
  }else if(varName=='U'&& !is.null(K)){
    pos <- lower.tri(matrix(0, K, K),diag = T)
    ind <- which(pos, arr.ind = TRUE)
    colnames(ind) <- c("F", "F")
  } else if(varName=='V' && LD && !is.null(J)){
    pos <- upper.tri(matrix(0, J, J), diag = TRUE)
    ind <- which(pos, arr.ind = TRUE)
    colnames(ind) <- c("Item", "Item")
  }else if(varName=='H'&& !is.null(J)&& !is.null(P)){
    pos <- matrix(0, J, P)
    ind <- which(!is.na(pos), arr.ind = TRUE)
    colnames(ind) <- c("Item", "Predictor")
  }



  if (ncol(dat0) == 1) {
    dat <- dat0[start:end]
    if (med) {
      est=median(dat)
    } else {
      est=mean(dat)
    }
    sd = sd(dat)
    CI <- coda::HPDinterval(mcmc(dat), prob = 1 - SL)
    sign <- as.numeric(CI[1] * CI[2] > 0)
    psrf <- compute_psrf(dat)
    res <- c(est, sd, CI, sign, psrf)
  } else if (ncol(dat0)>1) {
    dat <- dat0[start:end, ]
    if (med) {
      est=apply(dat, 2, median)
    } else {
      est=apply(dat, 2, mean)
    }
    sd = apply(dat, 2, sd)
    CI <- coda::HPDinterval(mcmc(dat), prob = 1 - SL)
    sign <- as.numeric(CI[, 1] * CI[, 2] > 0)
    psrf <- compute_psrf(dat)
    res <- cbind(est, sd, CI, sign, psrf,ind)
    if (sig) {
      sub <- (res[, 5] == 1)
      res <- res[sub, ]
    }
    return(res)
  }
}
