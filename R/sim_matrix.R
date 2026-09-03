#' Generate Population Parameter Matrices for Simulation
#'
#' @description
#' This function generates population-level parameter matrices for simulation
#' studies in partially confirmatory factor analysis (PCFA), bifactor models,
#' and MTMM-type structures. This function constructs the true loading matrix,
#' latent factor correlation matrix, and residual covariance matrix.
#'
#' The function supports:
#' - Multiple general factors Multiple specific (method) factors
#' - Flexible loading magnitudes
#' - Correlated general and/or specific factors
#' - Optional sparse specific factor structure
#' - Residual covariance structures
#' - Identification constraints for CFA-type estimation
#'
#' @param n_generals Integer. Number of general factors.
#' @param items_per_general Integer. Number of items loading on each general factor.
#' @param n_specifics Integer. Number of specific (group/testlet/method) factors.
#' @param items_per_specific Integer. Number of items loading on each specific factor.
#' @param blocks_spe Integer. Number of blocks used to evenly distribute specific-factor loadings across item clusters.
#'        When greater than 1, each specific factor is partitioned into multiple non-zero blocks. Useful for MTMM designs.
#' @param loadings_g Numeric vector. Population loading values for general factors.
#' @param loadings_s Numeric vector. Population loading values for specific factors.
#' @param lac_spe Numeric. Cross-loading value.
#' @param cpf_spe Integer. Cross-loading per factor for specific factors.
#' @param generals_rho Numeric. Correlation among general factors.
#' @param specifics_rho Numeric. Correlation among specific factors.
#' @param sparse_spe Integer. number of nonzero factor correlation for special factor.
#' Controls sparsity of specific factor correlation (large value means dense). default= 999: full correlated
#' @param seed Integer. Random seed for reproducibility.
#' @param ecr Numeric. Residual correlation strength.
#' @param necw Integer. Number of residual covariance elements within factor
#' @param necb Integer. Number of residual covariance elements between factor
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{A}{Loading matrix.}
#'   \item{R}{Implied population correlation matrix.}
#'   \item{Phi}{Latent factor covariance matrix.}
#'   \item{uniquenesses}{Diagonal residual variances (uniqueness terms).}
#'   \item{fixloading}{Indicator matrix specifying which loadings are fixed to 1
#'     for model identification (especially useful for non-standard CFA or MTMM estimation).}
#'   \item{ecm}{Residual covariance matrix (error covariance matrix).}
#' }
#'
#' @details
#' If no general factor is desired (i.e., simple CFA), set the
#' general factor component to zero. n_generals =0. The resulting structure
#' corresponds to a standard multi-factor CFA model.
#'
#' @examples
#' params <- sim_matrix(
#'   n_generals = 1,
#'   items_per_general = 6,
#'   n_specifics = 2,
#'   items_per_specific = 3
#' )
#'
#' str(params)
#'
#' @export
sim_matrix <- function(n_generals=0, items_per_general=0, n_specifics, items_per_specific,blocks_spe=1,
                       loadings_g =  c(.7, .7), loadings_s = c(.3, .3),lac_spe=0, cpf_spe=0,
                       generals_rho = 0, specifics_rho = 0,sparse_spe=999, seed=123,
                       ecr = .0,necw=0, necb=0) {

  ### this function is to generate population loading and phi
  ### if n_generals = 0, it is a simple CFA
  ### general factor now cannot have cross-loading
  if (!is.null(seed)) set.seed(seed)

  if(n_specifics==0 ){
    stop("n_specifics cannot be 0. If generate CFA n_specifics = n_factor")
  }

  # Numero items:
  if(n_generals!=0  &  n_generals * items_per_general != n_specifics*items_per_specific){
    stop("items number are wrong")
  }

  n_items <- n_specifics * items_per_specific

  # number total factors:
  n_factors <- n_generals + n_specifics

  # initial loading matrix
  A <- matrix(NA, nrow = n_items, ncol = n_factors)
  fixloading<- matrix(0, nrow = n_items, ncol = n_factors)


  ##### Specific factor loadings####
  # generate special loading:
  sequen_spe <- seq(loadings_s[1], loadings_s[2], length.out = items_per_specific)
  # for(i in 0:(n_specifics-1)) {
  #   start_row <- 1 + i*items_per_specific
  #   end_row <- start_row + items_per_specific - 1
  #   A[start_row:end_row , 1+i+n_generals] <- sequen_spe
  # }
  if(n_specifics == 1){
    mla <- matrix(sequen_spe, n_items, n_specifics)
  }else{
    lam0 <- c(sequen_spe, rep(0, items_per_specific - cpf_spe), rep(lac_spe, cpf_spe), rep(0, (n_specifics - 2) * items_per_specific))
    mla1 <- matrix(lam0, items_per_specific, n_specifics)
    mla <- c()
    for (i in n_specifics:1) {
      # i <- 1
      ind <- (c(1:n_specifics) + i)%%n_specifics
      ind[ind == 0] <- n_specifics
      mla <- rbind(mla, mla1[, ind])
    }
  }

  if(blocks_spe>1){  ### for separate in specifal factor, chap 1
    items_in_one_block<-items_per_specific/blocks_spe  #2
    items_per_block <- items_in_one_block * n_specifics   #6

    chunks <- split(sequen_spe, rep(seq_len(blocks_spe), each = items_in_one_block)) #separte loading value

    M <- matrix(0L, nrow = n_items, ncol = n_specifics)
    for (b in 0:(blocks_spe - 1)) {
      base <- b * items_per_block
      for (k in 1:n_specifics) {
        start <- base + (k - 1) * items_in_one_block + 1
        end   <- base + k * items_in_one_block
        idx <- seq.int(start, end)   # safer than ":" with +/()
        M[idx, k] <- chunks[[b + 1]]
      }
    }

    mla<-M
  }

  A[,( 1+n_generals):n_factors]<-mla

  ##### General factor loadings####
  if(n_generals!=0){
    # generate general loading:
    sequen_gen <- seq(loadings_g[1], loadings_g[2], length.out = items_per_general)
    for(i in 0:(n_generals-1)) {

      start_row <- 1 + i*items_per_general
      end_row <- start_row + items_per_general - 1
      A[start_row:end_row , i+1] <- sequen_gen # runif(items_per_general, loadings_g[1], loadings_g[2])

    }

    colnames(A) <- c(paste("G", 1:n_generals, sep = ""), paste("S", 1:n_specifics, sep = ""))

  }else{
    colnames(A) <- paste("K", 1:n_specifics, sep = "")
  }


  A[is.na(A)] <- 0
  rownames(A) <- paste("item", 1:nrow(A), sep = "")


  ##### fixloading for standarization####
  if(blocks_spe==1){  ## for blocks>1, don't use fixloading, keep std =T
    for(i in 1:n_specifics) {
      fix_col <- n_generals + i
      fix_row <- 1 + (i - 1) * items_per_specific
      if(A[fix_row,fix_col]==0){
        stop('the fixloading setting is problematic')
      }
      fixloading[fix_row,fix_col]<-1
    }
  }

  # Fix one loading per general factor to 1, avoiding already fixed items
  if(n_generals != 0) {
    for(i in 1:n_generals) {
      fix_col <- i
      fix_row <- 1 + (i - 1) * items_per_general

      if(any(fixloading[fix_row,] == TRUE)){
        fix_row<-fix_row+1
      }
      if(A[fix_row,fix_col]==0){
        stop('the fixloading setting is problematic')
      }
      fixloading[fix_row,fix_col]<-1
    }
  }

  # Correlaciones entre factores:
  Phi <- matrix(0, n_factors, n_factors)

  if(sparse_spe>=1 && sparse_spe<999){

    Phi_spe <- matrix(0,n_specifics,n_specifics)
    upper_indices <- which(upper.tri(Phi_spe), arr.ind = TRUE)
    selected <- upper_indices[sample(nrow(upper_indices), sparse_spe), , drop = FALSE]

    for (i in seq_len(sparse_spe)) {
      r <- selected[i, "row"]
      c <- selected[i, "col"]
      Phi_spe[r, c] <- Phi_spe[c, r] <- specifics_rho
    }

  }else if(sparse_spe==999){
    Phi_spe <-  matrix(specifics_rho,n_specifics,n_specifics)
  }else{
    stop("sparse: number of non-zero correlation for specific factor, 999: fully correlated")
  }

  Phi[1:n_generals, 1:n_generals] <- generals_rho
  Phi[(1+n_generals):n_factors, (1+n_generals):n_factors] <- Phi_spe
  diag(Phi) <- 1

  # Matriz de correlaciones entre items poblacional:

  R <- A %*% Phi %*% t(A)
  uniquenesses <- 1 - diag(R)
  diag(R) <- 1

  ecm <- matrix(0, n_items, n_items)
  if(ecr!=0){
    ## residual matrix
    ipf <- round(n_items / n_specifics)

    if (ecr > 0) {
      # iecb <- iecw <- NULL
      li <- c(1:n_specifics) * ipf - cpf_spe
      for (i in 1:necw) {
        if(necw>0) ecm[li[i], li[i] - 1] <- ecm[li[i] - 1,li[i]] <- ecr
      }
      for (i in 1:necb) {
        i1 <- i%%n_specifics + 1
        if(necb>0) ecm[li[i] - 2, li[i1] - 3] <- ecm[li[i1] - 3,li[i] - 2]<- ecr
      }
    }
  }
  diag(ecm)<-uniquenesses

  # if( any(uniquenesses < 0) ) {
  #
  #   warning("At least a communality greater than 1 found \n Resampling...")
  #
  #   sim <- sim_bifactor(n_generals, n_specifics, items_per_specific,
  #                       loadings, generals_rho, specifics_rho,
  #                       crossloadings, pure)
  #
  #   A = sim$A; R = sim$R; Phi = sim$Phi; uniquenesses = sim$uniquenesses;
  #
  # }

  return( list(A = A, R = R, Phi = Phi, uniquenesses = uniquenesses, fixloading=fixloading,ecm=ecm) )

}
