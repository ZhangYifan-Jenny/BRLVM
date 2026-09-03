#' Generate Sample Data
#'
#' @description
#' This function generates observed sample data from a specified factor model.
#' It supports continuous or categorical indicators, missingness, correlated latent factors,
#' optional cross-loadings, local dependence, and structured residual covariance.
#' The function can simulate:
#' - Standard CFA models
#' - Bifactor models
#' - MTMM models
#' - Models with cross-loading perturbations
#' - Models with local residual correlations
#' - Categorical response data
#' - Missing data mechanisms
#'
#' @param N Integer. Sample size.
#' @param lam Matrix. True loading
#' @param phi Numeric/Matrix. Correlation among latent factors.
#' @param ecm Matrix. Residual covariance matrix. If NULL, generated internally.
#' @param ecr Numeric. Residual correlation strength.
#' @param necw Integer. Number of residual covariance elements within blocks.
#' @param necb Integer. Number of residual covariance elements between blocks.
#' @param fixloading Logical or numeric indicator for identification constraints.
#' @param Kg Integer. Number of general factors (0 corresponds to simple CFA).
#' @param std Logical. Whether to standardize latent variables.
#' @param cati Numeric vector. The set of polytomous items in sequence number (i.e., can be any number set
#'            in between 1 and {J}); \code{NULL} for no and -1 for all .
#' @param noc Numeric vector. Number of levels for polytomous items. (if categorical).
#' @param misp Numeric. Proportion of missingness..
#' @param fac_score Logical. Whether to return true latent factor scores.
#' @param rseed Integer. Random seed for reproducibility.
#' @param digits Integer. Rounding precision of output.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{lam}{Loading matrix used for data generation.}
#'   \item{PHI}{Latent factor correlation matrix.}
#'   \item{Eigen}{common variance contributed by each factor.}
#'   \item{scale}{Empirical item standard deviations.}
#'   \item{var_ie}{residual (unique) variance proportion.}
#'   \item{fac}{True latent factor scores (returned if fac_score = TRUE).}
#'   \item{dat}{Simulated observed data matrix.}
#'   \item{fixloading}{Identification constraint indicator for fixed loadings.}
#'   \item{std}{Logical flag indicating whether latent variables were standardized.}
#' }
#'
#' @details
#' Data are generated according to the factor model:
#'
#'   X = Lambda F + epsilon
#'
#'
#' @export
sim_data <- function(N = 1000, lam, phi = .3,ecm=NULL,ecr = .0,
                     necw=0, necb=0, fixloading = 0,Kg=0,std = T, cati = NULL, noc = c(4), misp = 0,fac_score = FALSE, rseed = 333,digits = 4) {
  if (is.scalar(lam))
    stop("lam should be a J*K loading matrix.", call. = FALSE)

  if (exists(".Random.seed", .GlobalEnv))
    oldseed <- .GlobalEnv$.Random.seed else oldseed <- NULL
    set.seed(rseed)
    set.seed <- rseed

    oo <- options()       # code line i
    on.exit(options(oo))  # code line i+1
    # old_digits <- getOption("digits")
    options(digits = digits)


    lam<-as.matrix(lam)
    K <- ncol(lam)
    J <- nrow(lam)


    PHI <- cm_check(phi,K)
    if(any(PHI<=-1))
      stop("phi should be a correlation scalar or matrix (K*K & PD).", call. = FALSE)


    # ## set fix loading matrix for non-standard
    # ## if Kg>0, only fix specific factor with 1 loading, keep general factor standardize
    # if(!is.scalar(fixloading)){
    #   if (!all(dim(fixloading) == dim(lam)))
    #     stop("fixloading should have same dimension with loading (J*K).", call. = FALSE)
    #   fixloading<-as.matrix(fixloading)
    # }else{
    #   lam_par<-lam[,(Kg+1):K]
    #   fixloading<-matrix(0,nrow=nrow(lam_par),ncol=ncol(lam_par))
    #   for (j in 1:ncol(lam_par)) {
    #     idx <- which(lam_par[, j] != 0 & lam_par[, j]!=lac)[1]  # first non-zero row index
    #     if (!is.na(idx)) {
    #       fixloading[idx, j] <- 1
    #     }
    #   }
    #   fixloading<-cbind(matrix(0,J,Kg), fixloading)
    # }

    # Convert to non-standardized model
    if(std == FALSE){
      # 1. Sub columns for non-standardized
      ## for bifactor use, only specific factor change standardization
      cols_non_standardized <- which(colSums(fixloading == 1) == 1)
      if (length(cols_non_standardized) > 1) {
        cat("Factor", paste(cols_non_standardized, collapse = ", "), "are set to non-standardized.\n")
      }else{
        cat("Less than 2 factor are set to non-standardized.\n")
      }

      # 2. Extract the corresponding columns from lam
      lam_par <- lam[, cols_non_standardized, drop = FALSE]
      fixloading_par<-fixloading[, cols_non_standardized, drop = FALSE]

      D <- diag(lam_par[which(fixloading_par==1,arr.ind = T)])

      Lambda_star <- lam_par %*% solve(D)
      Phi_star <- D %*% PHI[cols_non_standardized,cols_non_standardized] %*% D

      lam[, cols_non_standardized]<-Lambda_star
      PHI[cols_non_standardized,cols_non_standardized]<-Phi_star
    }

    evr<-1 - diag(lam %*% PHI %*% t(lam))
    # if(is.null(ecm)){
    #   ## residual matrix
    #   ipf <- round(J / K)
    #   ecm <- matrix(0, J, J)
    #   if (ecr > 0) {
    #     # iecb <- iecw <- NULL
    #     li <- c(1:K) * ipf - cpf
    #     for (i in 1:necw) {
    #       if(necw>0) ecm[li[i], li[i] - 1] <- ecm[li[i] - 1,li[i]] <- ecr
    #     }
    #     for (i in 1:necb) {
    #       i1 <- i%%K + 1
    #       if(necb>0) ecm[li[i] - 2, li[i1] - 3] <- ecm[li[i1] - 3,li[i] - 2]<- ecr
    #     }
    #   }
    #   diag(ecm) <- evr
    # }

    err <- mvrnorm(N, rep(0, J), ecm)

    fac <- mvrnorm(N, rep(0, K), PHI)
    y <- t(lam %*% t(fac)) + err

    scale <- apply(y, 2, sd)
    eigen <- diag(crossprod(lam))


    out <- list(lam=lam, PHI = PHI, Eigen = eigen, scale = scale, var_ie=evr)
    if (fac_score)
      out$fac <- fac

    pos <- lower.tri(ecm)
    ind <- which(pos, arr.ind = T)
    rind <- which(ecm[pos] > 0)
    out$loc_dep <- ind[rind, ]

    UL <- 3;LL <- -3
    Jp <- length(cati)
    if (Jp > 0)
    {
      if (Jp == 1 && cati == -1) {
        cati <- c(1:J)
        Jp <- J
      }
      yc <- y[, cati]
      len <- length(noc)
      if (Jp%%len != 0)
        stop("cati is not divisable by length(noc).", call. = F)
      for (i in 1:len) {
        M <- noc[i]  #categories
        stp <- LL + c(1:(M - 1)) * (UL - LL)/M  #M-1 step points
        i0 <- (i - 1) * Jp/len + 1
        i1 <- i * Jp/len
        for (j in i0:i1) {
          tmp <- yc[, j] > matrix(stp, N, (M - 1), byrow = T)
          yc[, j] <- rowSums(tmp) + 1  # value starting from 1
        }
      }  #end len
      y[, cati] <- yc
      out$noc <- noc
      out$cati <- cati
    }  #end Jp

    if (misp > 0){
      mind <- matrix(rbinom(N * J, 1, misp), N, J)
      y0<-y[,1:J]
      y0[mind == 1] <- NA
      y[,1:J] <- y0
    }

    colnames(y)<-paste0("y",c(1:ncol(y)))
    out$dat <- y
    out$fixloading<-fixloading
    out$std <- std
    if (!is.null(oldseed))
      .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)

    return(out)
}
