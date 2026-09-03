#library(MASS)
is.scalar <- function(x) is.atomic(x) && length(x) == 1L
cm_check<-function(cm,K,diag=1,type=1){
  if (is.scalar(cm)){
    cm<-matrix(cm, K, K)
    diag(cm) <- diag
  }
  c0 <- any(dim(cm)!=K) #not K x K
  c1<- any(cm!=t(cm)) #not symmetric
  if (type ==1){
    c2<- any(diag(cm)!=1)
  }else{
    c2 <- any(diag(cm)<=0)
  }
  if (c0 || c1 || c2) return(-1)

  if (any(eigen(cm)$values<=0)) return(-2) #not positive definite

  return(cm)
}


plot_pcpm <- function(classname, varName, ind) {
  data <- classname$getVariable(varName)
  par(mfrow = c(2, 1), mar = c(2, 2, 2, 2))
  # Trace plot
  plot(data[, ind], type = 'l', main = paste("Trace Plot_",varName,ind), xlab = "Iteration", ylab = "Value")
  # Density plot
  hist(data[, ind], probability = TRUE, main = "Density Plot", xlab = "Value")
}

compute_psrf <- function(chain) {
  chain <- as.matrix(chain)
  num_vars <- ncol(chain)
  r_hats <- c()
  for (i in 1:num_vars) {
    variable_chain <- chain[, i]
    N <- length(variable_chain)
    chain1 <- variable_chain[1:(N/2)]
    chain2 <- variable_chain[(N/2 + 1):N]
    W <- 0.5 * (var(chain1) + var(chain2))
    if (W != 0) {
      B <- (N / 2) * (mean(chain1) - mean(chain2))^2
      V <- (1 - 1/N) * W + (1/N) * B
      r_hat <- sqrt(V / W)
    } else {
      r_hat <- 0
    }
    r_hats <- c(r_hats, r_hat)
  }

  if (num_vars > 1) {
    return(r_hats)
  } else {
    return(r_hats[1])
  }
}

SignSwitch<-function(LA, LB=NULL, Phi=NULL,reorder=F){
  reord<-function(matrix){
    uni<-unique(apply(matrix, 1, which.max))
    missing_columns <- setdiff(1:ncol(matrix), uni)
    matrix <- matrix[,c(uni,missing_columns)]
    return(matrix)
  }
  chg <- (colSums(LA)<= -0.1)
  if (any(chg)) {
    sign <- diag(1 - 2 * chg)
    LA <- LA %*% sign
    if(!is.null(LB)){
      LB <- t(t(LB) %*% sign)
    }
    if(!is.null(Phi)){
      Phi <- sign  %*% Phi %*% sign
      if( !all(Phi>=0)){
        print("illegal swich")
      }
    }

    print(chg)
  }

  if(reorder){
    LA <- reord(LA)
    if(!is.null(LB)){
      if(reorder){LB <- t(reord(t(LB)))}
    }
  }
  list(LA = LA, LB = LB,Phi = Phi,chg=chg)
}

#library(clue)
ColumnSwitch<-function(F1, F2, Phi2 = NULL){ # F1 and F2 should already do sign switch
  Nfac <- ncol(F1)
  UniqueMatch <- TRUE
  #cost = 1 - cor(F1, F2) #plus case cannot use this
  cost<-matrix(0,Nfac,Nfac)
  for(i in 1:Nfac){
    for(j in 1:Nfac){
      cost[i,j]<-sum((F1[,i]-F2[,j])^2)
    }
  }
  Qmatch <- apply(cost, 1, which.min)
  if (length(unique(Qmatch)) != Nfac) {
    UniqueMatch <- FALSE
  }
  if (UniqueMatch == FALSE) {
    Qmatch <- clue::solve_LSAP(cost, maximum = FALSE)
  }
  F2 <- F2[, Qmatch]
  FactorMap <- rbind(1:Nfac, Qmatch)
  rownames(FactorMap) <- c("Original Order", "Sorted Order")
  if (!is.null(Phi2)) {
    Phi2 <- Phi2[Qmatch, Qmatch]
  }
  list(F2 = F2, Phi2 = Phi2, FactorMap = FactorMap, UniqueMatch = UniqueMatch,Qmatch=Qmatch)
}

#ctrl+shift+c
# matrix_data <- matrix(c(
#   0.4869690, 0.2632395, 0.0000000,
#   0.4185412, 0.3540344, 0.0000000,
#   0.5444653, 0.1604072, 0.0000000,
#   0.5371706, 0.0000000, 0.0000000,
#   0.4646691, 0.2200288, 0.2221012,
#   0.3868229, 0.2228869, 0.1895790,
#   0.4016635, 0.3043689, 0.1805089,
#   0.6276564, 0.0000000, 0.0000000,
#   0.3726699, 0.0000000, 0.3734036,
#   0.0000000, 0.6795011, 0.0000000,
#   0.0000000, 0.5876358, 0.0000000,
#   0.0000000, 0.5672942, 0.1833906,
#   0.0000000, 0.5249648, 0.0000000,
#   0.0000000, 0.4858463, 0.0000000,
#   0.0000000, 0.6617170, 0.0000000,
#   0.1649421, 0.2498064, 0.0000000,
#   0.0000000, 0.5888468, 0.0000000,
#   0.1244326, 0.5528141, 0.0000000,
#   0.3391217, 0.0000000, 0.2407297,
#   0.1994041, 0.0000000, 0.3633352,
#   0.0000000, 0.0000000, 0.6577342,
#   0.0000000, -0.1558295, 0.8133125,
#   0.0000000, 0.0000000, 0.6548891,
#   0.0000000, 0.0000000, 0.5738575,
#   0.0000000, 0.0000000, 0.6212890,
#   0.0000000, 0.0000000, 0.5090497,
#   0.0000000, 0.0000000, 0.5470675), ncol = 3, byrow = TRUE)
