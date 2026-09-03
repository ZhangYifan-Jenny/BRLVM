#' @title PCPMr R6 Class
#'
#' @importFrom R6 R6Class
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp Module
#'
#' @field N Number of observations.
#' @field J Number of items.
#' @field K Number of latent factors.
#' @field P Number of predictors.
#' @field Q Loading pattern matrix.
#' @field Qb Regression pattern matrix.
#' @field LD Logical indicating local dependence.
#' @field DIF Logical indicating differential item functioning.
#' @field time Computation time.
#' @field MODEL Model type identifier.
#' @field FacSco Stored factor scores.
#'
#' @description This class wraps the PCPM C++ class using R6.
#' @examples
#' # example code
#'
#' \donttest{
#' truemodel<-sim_matrix(n_specifics=3, items_per_specific=6,loadings_s = c(.7,.7),lac_spe=0.2, cpf_spe=2,specifics_rho = 0.3, seed=123)
#' dat <- sim_data( N = 100,lam = truemodel[["A"]],phi = truemodel[["Phi"]],ecm =truemodel[["ecm"]])
#' J<-nrow(dat$lam);K<-ncol(dat$lam)
#' Q <- matrix(-1, nrow = J, ncol = K)
#' Q[dat$lam!=0]=-2
#' Y<-dat$dat
#' a<-PCPMr$new()
#' a$pcfa(Y, Q, LD=F,iter=100, burn=0,std=T, cor_fac=1, update = 100)
#' out<-a$pcpmsummary(sig=F)
#' }
#' @export
PCPMr<- R6::R6Class(
  classname = "PCPMr",
  public = list(
    N='int',
    J='int',
    K='int',
    P='int',
    Q='double',
    Qb='double',
    LD='bool',
    DIF='bool',
    time='double',
    MODEL='character',
    FacSco='double',

    #' @description Initializes a new PCPM object.
    #' @param priorF A numeric vector for priorF (default: c(0.1, 1.0, 1.0)). w0
    #' @param priorA A numeric vector for priorA (default: c(1.0, 1.0, 1.0, 1.0, 0.0, 1.0)).s0;r0;lam_a;lam_b;m0;c0;
    #' @param priorB A numeric vector for priorB (default: c(1.0, 1.0, 1.0, 1.0, 0.0, 1.0)).s0B;r0B;lam_aB;lam_bB;m0B;c0B;
    initialize=function(priorF = c(0.1),
             priorA = c(1.0, 0.1, 1.0, 0.1, 0.0, 0.1),
             priorB = c(1.0, 0.1, 1.0, 0.1, 0.0, 0.1)) {

      if(is.null(priorF)){
        priorF=c(w0=0.1)
      }else if(length(priorF)!=1){
        stop("The numbers of priorF should be 1", call. = FALSE)
      }

      if(is.null(priorA)){
        priorA=c(1.0, 0.1, 1.0, 0.1, 0.0, 0.1)
      }else if(length(priorA)!=6){
        stop("The numbers of priorA should be 6", call. = FALSE)
      }

      if(is.null(priorB)){
        priorB=c(1.0, 0.1, 1.0, 0.1, 0.0, 0.1)
      }else if(length(priorB)!=6){
        stop("The numbers of priorB should be 6", call. = FALSE)
      }


      #loadModule("PCPM", TRUE)
      #RCPP<- new(PCPM,priorF,priorA,priorB)
      private$PCPM_Module <- new(PCPM,priorF,priorA,priorB)
    },

    #' @description Performs PCFA on the provided data.
    #' @param Y A numeric matrix.
    #' @param Q An integer matrix.
    #' @param cati A logical value, wheter item is categorical.
    #' @param cati_v A vector or -1. which item is categorical, -1 means all.
    #' @param reglo A string. which prior for loading. 'lasso','horse','ssp'
    #' @param LD A logical value.
    #' @param regphi A string. which prior for factor correlation. 'none', lasso','horse','ssp','parlasso','parhorse','parssp'
    #' @param Kg A numeric. how many general factor, if none, 0
    #' @param reg An integer matrix. set for par phi, not used here. default as matrix full with 1.
    #' @param regpsx A string. which prior for loading. 'lasso','horse','ssp'
    #' @param iter An integer value for the number of iterations.
    #' @param burn An integer value for the burn-in period.
    #' @param std A logical value indicating whether to standardize.
    #' @param cor_fac A numeric vector for correlation factors.
    #' @param sign_check A logical value.
    #' @param update An integer value for the update frequency (default: 1000).
    #'
    #'
    #'
    pcfa=function(Y, Q, cati=F,cati_v=-1,reglo = 'lasso',LD=F,regphi = 'lasso',Kg = 0,
                  reg,regpsx = 'lasso',iter =1000,  burn = 0000, std = T, cor_fac =1,  sign_check = F, update = 1000) {
                self$N <- nrow(Y)
                self$J <- ncol(Y)
                self$K <- ncol(Q)
                self$Q<-Q
                reg<-matrix(1,ncol(Q),ncol(Q))
                self$LD<-LD
                self$MODEL<-'pcfa'
                private$PCPM_Module$pcfa(Y, Q, cati,cati_v,reglo,LD,regphi,Kg,
                                         reg,regpsx,iter,burn, std, cor_fac, sign_check, update)
                self$FacSco = private$PCPM_Module$mfac
    },

    #' @description get posterial summary
    #' @param varName aaa.
    getVariable = function(varName) {
      private$PCPM_Module$getVariable(varName)
    },

    #' @description Performs pcpmsummary.
    #' @param med aaa.
    #' @param start aaa.
    #' @param end aaaa
    #' @param SL aaa
    #' @param sig aaa
    #' @param detail aaa
    #' @param SignSwitch aaa
    #' @param reorder aaa
    pcpmsummary= function(med =FALSE,  start = 00,  end = -1,  SL = 0.05,  sig =T,
                          detail=F,SignSwitch=F,reorder=F) {
      J<-self$J
      K<-self$K
      P<-self$P
      LD<-self$LD
      DIF<-self$DIF
      out <- list(model=self$MODEL)
      fullout <- list(model=self$MODEL)
      #### factorial correlation (pcfa)
      if(self$MODEL=='pcfa'|self$MODEL=='pcfa_cati'){
        A<-stat(classname=private$PCPM_Module, varName='A',Q=self$Q, med =med,  start = start,  end = end,  SL = SL,  sig =sig)
        if(sig){
          loadA<-matrix(0,nrow=J,ncol=K)
          for(i in 1:nrow(A)){
            loadA[A[i,7],A[i,8]]<-A[i,1]
          }
        }else{
          loadA<-matrix(A[,1],nrow=J,ncol=K)
        }

        C<-stat(classname=private$PCPM_Module, varName='C',K=K,med =med,  start = start,  end = end,  SL = SL,  sig =sig)
        phi<-matrix(0,nrow=K,ncol=K)
        for(i in 1:nrow(C)){
          phi[C[i,7],C[i,8]]<-C[i,1]
        }

        #### sign switch (pcfa)
        if(SignSwitch){
          SIGN<-SignSwitch(loadA,Phi=phi,reorder=reorder)
          loadA<-SIGN$LA
          phi<-SIGN$Phi
        }


        #### error covariance
        V<-stat(classname=private$PCPM_Module, varName='V',J=J,LD=LD,med =med,  start = start,  end = end,  SL = SL,  sig =sig)

        if(LD){
          PSX<-matrix(0,nrow=J,ncol=J)
          for(i in 1:nrow(V)){
            PSX[V[i,7],V[i,8]]<-V[i,1]
          }
        }else{
          PSX<-V[,1]
        }


        out$loadA<-loadA
        out$PHI  <-phi
        out$PSX  <- PSX

        fullout$loadA<-A
        fullout$PHI  <-C
        fullout$PSX  <-V
      }

      if(detail==F){
        return(out)
      }else{
        return(fullout)
      }
      }


  ),
  private = list(
    PCPM_Module=NULL

  )
)
