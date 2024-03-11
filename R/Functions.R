

#==========================================================================#
# AuxCE-PH using CV: Fit the Cox PH Model with Auxiliary SP using Control Variate method
#==========================================================================#

#==== The main function for transfering information from external sources to target data ====#
#' @title Cox proportional hazards model with adaptive information transfer
#'
#' @description Fit the Cox proportional hazards model with aggregate estimate from multiple external sources.
#'
#' @aliases PH.Transfer.CV.fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param aux indicates the aggregate statistics extracted from external sources (based on transitional models). It is a list of lists, and each sub-list represents information from a study.
#' @param uncertain a list that contains two elements:
#'   \code{do} is a logical to indicate whether the uncertainty issue should be addressed;
#'   \code{haveV} is another logical value to show whether we have external estimated variance-covariance matrices.
#' @param hetero a list that contains two elements:
#'   \code{do} is a logical to indicate whether the heterogeneity issue should be addressed;
#'   \code{group} is another logical value to show whether the group Lasso penalty should be used.
#' @param Lams a list that contains \code{tm}, which denotes a vector time points where it will evaluate the baseline cumulative hazard function on.
#' @param eps the error of convergence.
#' @param maxit the maximum number of iteration.
#' @param trace a logical to indicate whether the traces of model fitting should be printed.
#'
#' @examples
#' # Here is an example (to appear).
#'
#' @export PH.Transfer.CV.fit
PH.Transfer.CV.fit <- function(
    yobs,delta,X, # the internal dataset
    aux,uncertain=list(do=TRUE,haveV=TRUE),# information about auxiliary information
    hetero=list(do=TRUE,group=TRUE), # information about heterogeneity
    Lams=list(tm=NULL), # whether to output the estimation of lambda
    eps=1e-6,maxit=5e4,
    trace=TRUE
){

  # Tips:
  #   hetero = FALSE | uncertain = FALSE               : Homogeneous Auxiliary Information Synthesis without Uncertainty
  #   hetero = FALSE | uncertain = TRUE (haveV = FALSE): Homogeneous Auxiliary Information Synthesis with Uncertainty
  #   hetero = FALSE | uncertain = TRUE (haveV = TRUE) : Federated Transfer Learning (Homogeneous)
  #   hetero = TRUE  | uncertain = FALSE               : Heterogeneous Auxiliary Information Synthesis without Uncertainty
  #   hetero = TRUE  | uncertain = TRUE (haveV = FALSE): [Intractable Case]
  #   hetero = TRUE  | uncertain = TRUE (haveV = TRUE) : Federated Transfer Learning (Heterogeneous)

  ### ---- Specify several basic elements ---- ###
  pbet <- ncol(X)   # number of parameters
  yobs[yobs<=0] <- 1e-6
  bet.names <- colnames(X)
  N <- length(yobs)   # sample size in internal data
  K <- length(aux) # number of external studies
  source.index <- do.call(c,lapply(1:K,function(ksource){rep(ksource,length(aux[[ksource]]$ce))})) # a vector to indicate the study number of ce
  methods.suffix <- ifelse(hetero$do==TRUE,"hetero","homo")

  ### ---- various quantities (internal data only) ---- ###
  orimodel.fit <- PH.fit(yobs,delta,X,eps=eps,maxit=maxit)
  bet.ori <- orimodel.fit$res[,1]
  Influs.ori <- PH.Influence(yobs,delta,X,bet=bet.ori)
  Sigma.ori <- t(Influs.ori)%*%Influs.ori/N # orimodel.fit$VarCov #
  if(!is.null(Lams$tm)){
    Lam.ori <- PH.Lam(Lams$tm,yobs,delta,X,bet.ori)
    Influs.Lam.ori <- PH.Influence.Lam(Lams$tm,yobs,delta,X,bet=bet.ori)
    Sigma.Lam <- t(Influs.Lam.ori)%*%Influs.Lam.ori/N
  }

  ### ---- expand the aux: ce.inn + Influs.inn ---- ###
  for(ksource in 1:K){ # ksource <- 1;
    aux.c <- aux[[ksource]]
    # corresponding estimates of reduced coefficients using internal data
    ceall.inn.c <- PH.fit(yobs,delta,X[,aux.c$idxall,drop=F],
                          Var=FALSE,maxit=maxit,eps=eps)$res[,1]
    ce.inn.c <- ceall.inn.c[aux.c$idxall %in% aux.c$idxobs]
    aux[[ksource]]$ceall.inn <- ceall.inn.c
    aux[[ksource]]$ce.inn <- ce.inn.c
    # the influence functions of the above reduced coefficients
    Influs.c <- PH.Influence(yobs,delta,X[,aux.c$idxall,drop=F],bet=ceall.inn.c)
    Influs.inn.c <- Influs.c[,aux.c$idxall %in% aux.c$idxobs]
    aux[[ksource]]$Influs.inn <- Influs.inn.c
  }

  ### ---- get desired elements and matrices ---- ###
  # - for control variates
  ce.ext <- do.call(c,lapply(1:K,function(ksource){aux[[ksource]]$ce}))
  ce.inn <- do.call(c,lapply(1:K,function(ksource){aux[[ksource]]$ce.inn}))
  cva <- ce.inn-ce.ext
  # - for variance of the control variates
  Influs.inn <- do.call(cbind,lapply(1:K,function(ksource){aux[[ksource]]$Influs.inn}))
  V.inn <- t(Influs.inn)%*%Influs.inn/N #
  if(uncertain$do==TRUE){
    if(uncertain$haveV==TRUE){
      V.ext.rho <- as.matrix(do.call(Matrix::bdiag,lapply(1:K,function(ksource){aux[[ksource]]$ceV*N/aux[[ksource]]$M})))
    }else{
      V.ext.rho <- array(0,dim=dim(V.inn))
      for(ksource in 1:K){
        idxk <- (source.index==ksource)
        V.ext.rho[idxk,idxk] <- V.inn[idxk,idxk]*(N/aux[[ksource]]$M)
      }
    }
  }else{
    V.ext.rho <- 0
  }
  V <- V.inn + V.ext.rho # indeed here is Vinn+rho*Vext
  V.inv <- solve(V) # eigen(V.new)$values solve(V)  V %*% V.inv %*% V
  # - for covariance between auxiliary counterpart and original estimator
  Gam <- t(Influs.inn)%*%Influs.ori/N
  if(!is.null(Lams$tm)){
    Gam.Lam <- t(Influs.inn)%*%Influs.Lam.ori/N
  }
  # - for combined quantities
  Gamt.V.inv <- t(V.inv%*%Gam)
  if(!is.null(Lams$tm)){
    Gamt.V.inv.Lam <- t(V.inv%*%Gam.Lam)
  }

  ### ---- estimators with auxiliary information - homo ---- ###
  if("homo" %in% methods.suffix){

    # obtain estimates
    bet.homo <- as.vector(bet.ori-Gamt.V.inv%*%cva)
    # calculate the corresponding variances
    A <- Sigma.ori-t(Gam)%*%V.inv%*%Gam
    bet.homo.se  <- sqrt( diag(A)/N )
    # obtain estimates (for lam, using plug-in principle)
    if(!is.null(Lams$tm)){
      # calculate estimates at different tms (plug-in)
      Lam.homo <- as.vector(Lam.ori-Gamt.V.inv.Lam%*%cva)
      # calculate standard errors
      A.Lam <- Sigma.Lam-t(Gam.Lam)%*%V.inv%*%Gam.Lam
      Lam.homo.se  <- sqrt( diag((A.Lam)/N) )
    }
    # other elements
    info.homo <- list()

  }

  ### ---- estimators with auxiliary information - hetero ---- ###
  if("hetero" %in% methods.suffix){

    # transform necessary matrices
    V.inv.SVD <- svd(V.inv)
    V.inv.root <- V.inv.SVD$u%*%diag(sqrt(V.inv.SVD$d))%*%t(V.inv.SVD$v)
    # solve adaptive lasso using the well-known R package 'lars'
    sol.penfit <- PH.Transfer.CV.Pen.Quadratic.fit(
      bet=bet.ori,cva=cva,V.inv.root=V.inv.root,
      Gamt.V.inv=Gamt.V.inv,N=N,group=hetero$group,group.id=source.index)
    # obtain final estimates
    bet.hetero <- sol.penfit$bet; bet.hetero
    tau.hetero <- sol.penfit$tau; tau.hetero
    # estimator of variance (using oracle properties directly)
    homo.idx <- (tau.hetero==0)
    V.inv.homo <- solve(V[homo.idx,homo.idx,drop=F])
    if(any(homo.idx)){
      Gam.homo <- Gam[homo.idx,,drop=F]
      A.hetero <- Sigma.ori-t(Gam.homo)%*%V.inv.homo%*%Gam.homo
      bet.hetero.se  <- sqrt( diag(A.hetero)/N )
    }else{
      bet.hetero.se <-  sqrt(diag(Sigma.ori)/N)
    }
    # obtain estimates (for lam, using plug in principle)
    if(!is.null(Lams$tm)){
      # calculate estimates at different tms
      Lam.hetero <- as.vector(Lam.ori-Gamt.V.inv.Lam%*%(cva-tau.hetero))
      # calculate standard errors
      if(any(homo.idx)){
        Gam.Lam.homo <- Gam.Lam[homo.idx,,drop=F]
        A.Lam <- Sigma.Lam-t(Gam.Lam.homo)%*%V.inv.homo%*%Gam.Lam.homo
        Lam.hetero.se  <- sqrt( diag(A.Lam)/N )
      }else{
        Lam.hetero.se  <- sqrt( diag(Sigma.Lam)/N )
      }
    }
    # other elements
    tau.hetero.list <- rep(list(list()),K)
    names(tau.hetero.list) <- names(aux)
    for(ksource in 1:K){
      tau.hetero.list[[ksource]] <- tau.hetero[source.index==ksource]
    }
    info.hetero <- list(tau=tau.hetero.list)

  }

  ### ---- do inference and combine results into pre-specified style ---- ###
  suffix <- methods.suffix
  # for estimates (bet)
  bet.c <- get(paste(c("bet"),suffix,sep="."))
  bet.c.se <- get(paste(c("bet"),suffix,"se",sep="."))
  zvalue.bet.c <- bet.c/bet.c.se
  pvalue.bet.c <- 2*(1-pnorm(abs(zvalue.bet.c)))
  res.c <- data.frame(Est=bet.c,SE=bet.c.se,zvalue=zvalue.bet.c,pvalue=pvalue.bet.c,row.names=bet.names)
  # for estimates (Lam)
  if(!is.null(Lams$tm)){
    Lam.c <- get(paste(c("Lam"),suffix,sep="."))
    Lam.c.se <- get(paste(c("Lam"),suffix,"se",sep="."))
    zvalue.Lam.c <- Lam.c/Lam.c.se
    pvalue.Lam.c <- 2*(1-pnorm(abs(zvalue.Lam.c)))
    res.Lam.c <- data.frame(Est=Lam.c,SE=Lam.c.se,zvalue=zvalue.Lam.c,pvalue=pvalue.Lam.c,
                            row.names=paste("t=",Lams$tm,sep=""))
  }else{
    res.Lam.c <- NULL
  }
  # combine inference results
  info.c <- c(get(paste(c("info"),suffix,sep=".")),Lam=list(res.Lam.c))
  out <- c(list(res = res.c),info.c)

  ### ---- extract output values ---- ###
  return(out)

}


#==== Sub-function: for fitting penalized estimator ====#
PH.Transfer.CV.Pen.Quadratic.fit <- function(
    bet,cva,V.inv.root,Gamt.V.inv,N,
    group=FALSE,group.id=NULL
){

  y.tilde <- as.vector( V.inv.root%*%cva )
  X.tilde <- V.inv.root

  # solve the penalized problem
  if(group==FALSE){

    # solve adaptive lasso using lars
    w <- 1/abs(cva)
    X.tilde.star <- t(t(X.tilde)/w)
    sol.lars <- lars::lars(X.tilde.star,y.tilde,trace=FALSE,normalize=FALSE,intercept=FALSE)
    tau.path <- t(as.matrix(sol.lars$beta))/w # each
    tau.path.RSS <- apply(X.tilde %*% tau.path - y.tilde,2,function(x){sum(x^2)})
    tau.path.Card <- apply(tau.path,2,function(x){sum(x!=0)})
    IC.all <- as.numeric( tau.path.RSS + log(N)/N * tau.path.Card ) # BIC type criterion
    min_IC.idx <- which.min( IC.all  )
    tau <- tau.path[,min_IC.idx]

  }else{

    # solve adaptive group lasso using gglasso
    tau.nopen <- cva
    group.num.each <- as.numeric(table(group.id))
    group.num <- length(group.num.each)
    tau.nopen.group.norm <- sapply(1:group.num,function(igroup){sqrt(sum(tau.nopen[group.id==igroup]^2))})
    w <- 1/tau.nopen.group.norm # sqrt(group.num.each)/tau.nopen.group.norm #
    sol.gglasso <- gglasso::gglasso(x=X.tilde,y=y.tilde,group=group.id,pf=w,intercept=FALSE)
    tau.path <- sol.gglasso$beta
    tau.path.RSS <- apply(X.tilde %*% tau.path - y.tilde,2,function(x){sum(x^2)})
    tau.path.group.norm <- apply(tau.path,2,function(tau.pen.i){
      sapply(1:group.num,function(igroup){sqrt(sum(tau.pen.i[group.id==igroup]^2))})
    })
    tau.path.group.Card <- apply(tau.path.group.norm,2,function(x){sum(x!=0)})
    degrees.freedom <- tau.path.group.Card +
      apply(tau.path.group.norm*(group.num.each-1)/tau.nopen.group.norm,2,sum)
    IC.all <- as.numeric( tau.path.RSS + log(N)/N * degrees.freedom ) # BIC type criterion
    min_IC.idx <- which.min( IC.all  )
    tau <- tau.path[,min_IC.idx]

  }

  # output: return final estimates
  out <- list(
    bet=as.vector(bet-Gamt.V.inv%*%(cva-tau)),
    tau=tau
  )
  return(out)

}





#==========================================================================#
# Cox proportional hazards model (PH) -- based on: Profile likelihood
#==========================================================================#

#==== Function for fitting Cox proportional hazards model with partial likelihood ====#
#' @title Cox proportional hazards model based on the profiled likelihood estimation procedure
#'
#' @description Fit Cox proportional hazards model based on the profiled likelihood estimation procedure.
#'
#' @aliases PH.fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates
#' @param Var a logical value. The default setup is \code{TRUE}, indicating that the standard errors of the estimated regression coefficients will be calculated.
#' @param Var.Robust a logical value. The default setup is \code{TRUE}, indicating that we assume that the model is correctly specified, so that there is no need to use the robust estimator of standard errors.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}.
#'
#' @return The fitted results are returned (a list). The estimates of regression coefficients is \code{res}.
#'
#' @export PH.fit
PH.fit <- function(yobs,delta,X,Var=TRUE,Var.Robust=FALSE,eps=1e-6,maxit=5e4){
  # need: survival

  ### Preparations
  N <- length(yobs)
  yobs[yobs<=0] <- 1e-6
  pbet <- ncol(X)

  ### calculate initial values for bet
  bet.init <- rep(0,pbet)

  ### calculate MLE for beta
  numit <- 1
  bet.old <- bet.init
  repeat{
    InfoMScore <- PH.ScoreInfoM(bet=bet.old,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
    Score <- InfoMScore$Score
    InfoM <- InfoMScore$InfoM
    dev <- MASS::ginv(InfoM)%*%Score
    bet <- bet.old + as.vector(dev)
    if( max(abs(bet-bet.old))>eps & numit<maxit ){
      bet.old <- bet
      numit <- numit + 1
    }else{
      break
    }
  }
  convergence <- (numit<maxit)

  ### calculate SEs or not (and tidy them)
  if(Var==TRUE){

    ## calculate SEs for bet using explicit formula !!!
    InfoMScore <- PH.ScoreInfoM(bet=bet,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
    InfoM <- InfoMScore$InfoM
    Score <- InfoMScore$Score
    InfoM.inv <- MASS::ginv(InfoM)
    if(Var.Robust==TRUE){
      Score.Influs <- PH.Influence.EE(yobs=yobs,delta=delta,X=X,bet=bet)
      Score.Influs.Cov <- t(Score.Influs)%*%Score.Influs/N
      VarCov <- InfoM.inv %*% Score.Influs.Cov %*% InfoM.inv
    }else{
      VarCov <- InfoM.inv
    }
    bet.se <- sqrt(diag(VarCov)/N)

    ### tidy the results: inference
    zvalue.bet <- bet/bet.se
    pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
    res <- data.frame(Est=bet,SE=bet.se,zvalue=zvalue.bet,pvalue=pvalue.bet,
                      row.names=colnames(X))

  }else{

    # tidy the results directly
    VarCov <- NULL
    res <- data.frame(Est=bet,row.names=colnames(X))

  }

  ### output
  out <- list(
    info = list(
      convergence = convergence,
      bet.init = bet.init
    ),
    res=res,
    Score=Score,
    InfoM=InfoM,
    VarCov=VarCov
  )
  return(out)

}


#==== Profiled log-likelihood function in Cox proportional hazards model ====#
#' @title Profiled log-likelihood function in Cox proportional hazards model
#'
#' @description Calculate profiled log-likelihood function in Cox proportional hazards model.
#'
#' @aliases PH.Beta.LogLik
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates
#'
#' @export PH.Beta.LogLik
PH.Beta.LogLik <- function(bet,yobs,delta,X){
  # for maximization

  ## prepare
  N <- length(yobs)
  Xbet <- as.vector(X %*% bet)
  log.SS0 <- sapply(yobs,function(Yi){log(mean(exp(Xbet)*(yobs>=Yi)))})

  ## calculate the partial likelihood function for beta (log form)
  val.log <- sum((Xbet-log.SS0)*delta)

  # output
  return(val.log)

}


#==== Obtain score vector and information matrix ====#
#' @title Score vector and information matrix in Cox proportional hazards model
#'
#' @description Calculate score vector and information matrix in Cox proportional hazards model.
#'
#' @aliases PH.ScoreInfoM
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param IsScore whether the score vector is calculated or not.
#' @param IsInfoM whether the information matrix is calculated or not.
#'
#' @export PH.ScoreInfoM
PH.ScoreInfoM <-  function(bet,yobs,delta,X,IsScore=FALSE,IsInfoM=TRUE){
  # Score is the first  derivative of [positive] log-likelihood
  # InfoM is the second derivative of [negative] log-likelihood

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])})/N # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)}))/N # [S(1)(t,bet)] at pre-specified yobs and beta
  out <- list()

  ## prepare information matrix
  if(IsInfoM==TRUE){
    SS2.vec <- do.call(rbind,lapply(yobs,function(Yi){
      yobsGYi <- yobs>=Yi
      as.vector( t(X[yobsGYi,,drop=F])%*%(XexpXbet[yobsGYi,,drop=F]) )
    }))/N
    I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=pbet,ncol=pbet)
    I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
    InfoM <- I1-I2
    out <- c(out,list(InfoM=InfoM))
  }

  ## prepare score vector
  if(IsScore==TRUE){
    U <- (X - SS1/SS0)*delta
    Score <- apply(U,2,mean)
    out <- c(out,list(Score=Score))
  }


  ## output
  return(out)

}


#==== Obtain individual level score vector ====#
#' @title Individual level score vector in Cox proportional hazards model
#'
#' @description Calculate the individual level score vector in Cox proportional hazards model.
#'
#' @aliases PH.Score.Individual
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param Iscov whether convert the individual level form into a covariance form.
#'
#' @export PH.Score.Individual
PH.Score.Individual <-  function(bet,yobs,delta,X,Iscov=TRUE){
  # Score is the first  derivative of [positive] log-likelihood

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])})/N # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)}))/N # [S(1)(t,bet)] at pre-specified yobs and beta
  U <- (X - SS1/SS0)*delta
  if(Iscov==TRUE){
    out <- t(U)%*%U/N
  }else{
    out <- U
  }

  ## output
  return(out)

}


#==== Nonparametric component ====#
#' @title Nonparametric component in Cox proportional hazards model
#'
#' @description  Nonparametric component for the Cox proportional hazards model at specified time points.
#'
#' @aliases PH.Lam
#'
#' @param tm The time points that the nonparametric baseline function will be estimated at.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet unknown parameters corresponding to the model.
#' @param type if \code{type="right"}, right-continuous results will be calculated, and if \code{type="left"}, left-continuous results will be calculated.
#'
#' @export PH.Lam
PH.Lam <- function(tm,yobs,delta,X,bet=NULL,type="right"){

  # prepare
  if(is.null(bet)){
    bet <- PH.fit(yobs,delta,X)$res[,1]
  }
  N <- length(yobs)
  expXbet <- as.vector(exp(X %*% bet))
  SS0 <- sapply(yobs,function(Yi){mean(expXbet*(yobs>=Yi))}) # [S(0)(t,bet)] at pre-specified yobs and beta

  # calculate the Lam(t) at specified time points tm
  if(type=="right"){
    Lam <- sapply(tm,function(tmj){
      sum( (yobs<=tmj)[delta==1]/SS0[delta==1] ) / N
    })
  }else if(type=="left"){
    Lam <- sapply(tm,function(tmj){
      sum( (yobs< tmj)[delta==1]/SS0[delta==1] ) / N
    })
  }


  # output
  return(Lam)

}


#==== The influences for proportional hazards model ====#
#' @title  Influence functions for the estimator of regression coefficients in the Cox proportional hazards model
#'
#' @description Calculate the influence functions for the estimator of regression coefficients in the Cox proportional hazards model.
#'
#' @aliases PH.Influence
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet unknown parameters corresponding to the model.
#'
#' @export PH.Influence
PH.Influence <- function(yobs,delta,X,bet=NULL){

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  if(is.null(bet)){
    bet <- PH.fit(yobs,delta,X)$res[,1]
  }
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta

  ## calculate main part (except an inverse of information matrix)
  UU1 <- (X - SS1/SS0)*delta
  UU2 <- X * sapply(yobs,function(Yi){sum(delta*(yobs<=Yi)/SS0)})
  UU3 <- do.call(rbind,lapply(yobs,function(Yi){apply(delta*(yobs<=Yi)*SS1/SS0^2,2,sum)}))
  U <- UU1-(UU2-UU3)*expXbet

  ## prepare information matrix
  SS2.vec <- do.call(rbind,lapply(yobs,function(Yi){
    yobsGYi <- yobs>=Yi
    as.vector( t(X[yobsGYi,,drop=F])%*%(XexpXbet[yobsGYi,,drop=F]) )
  }))
  I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=pbet,ncol=pbet)
  I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
  InfoM <- I1-I2

  ## final influence functions (individual level) and output
  Influs <- U%*%MASS::ginv(InfoM)
  return(Influs)

}


#==== The influences for score of proportional hazards model ====#
#' @title  Influence functions for the estimator of the score equation in the Cox proportional hazards model
#'
#' @description Calculate the influence functions for the estimator of the score equation in the Cox proportional hazards model.
#'
#' @aliases PH.Influence.EE
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet unknown parameters corresponding to the model.
#'
#' @export PH.Influence.EE
PH.Influence.EE <- function(yobs,delta,X,bet=NULL){
  # The incluence function for estimating equations

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  if(is.null(bet)){
    bet <- PH.fit(yobs,delta,X)$res[,1]
  }
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta

  ## calculate main part
  UU1 <- (X - SS1/SS0)*delta
  UU2 <- X * sapply(yobs,function(Yi){sum(delta*(yobs<=Yi)/SS0)})
  UU3 <- do.call(rbind,lapply(yobs,function(Yi){apply(delta*(yobs<=Yi)*SS1/SS0^2,2,sum)}))
  U <- UU1-(UU2-UU3)*expXbet

  ## final influence functions (individual level) and output
  return(U)

}


#==== The influences for nonparametric part of proportional hazards model ====#
#' @title  Influence functions for the estimator of the nonparametric part in the Cox proportional hazards model
#'
#' @description Calculate the influence functions for the estimator of the nonparametric part in the Cox proportional hazards model.
#'
#' @aliases PH.Influence.Lam
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet unknown parameters corresponding to the model.
#'
#' @export PH.Influence.Lam
PH.Influence.Lam <- function(tm,yobs,delta,X,bet=NULL){

  # tm <- c(0.2,0.4,0.6,0.8,1)

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  if(is.null(bet)){
    bet <- PH.fit(yobs,delta,X)$res[,1]
  }
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta

  ## prepare information matrix
  SS2.vec <- do.call(rbind,lapply(yobs,function(Yi){
    yobsGYi <- yobs>=Yi
    as.vector( t(X[yobsGYi,,drop=F])%*%(XexpXbet[yobsGYi,,drop=F]) )
  }))
  I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=pbet,ncol=pbet)
  I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
  InfoM <- I1-I2

  ## calculate main part 1 (same as bet)
  UU1 <- (X - SS1/SS0)*delta
  UU2 <- X * sapply(yobs,function(Yi){sum(delta*(yobs<=Yi)/SS0)})
  UU3 <- do.call(rbind,lapply(yobs,function(Yi){apply(delta*(yobs<=Yi)*SS1/SS0^2,2,sum)}))
  U <- UU1-(UU2-UU3)*expXbet
  Influs.bet <- U%*%MASS::ginv(InfoM)

  ## calculate main part 2 (with tm vary)
  LL1.tm <- do.call(cbind,lapply(tm,function(tmi){(yobs<=tmi)*delta/(SS0/N)}))
  LL2.tm <- do.call(cbind,lapply(tm,function(tmi){
    sapply(yobs,function(Yj){sum(delta*(yobs<=min(tmi,Yj))/(SS0^2/N))})
  }))*expXbet
  LL3.pre <- do.call(cbind,lapply(tm,function(tmi){apply(SS1*delta*(yobs<=tmi)/SS0^2,2,sum)}))
  LL3.tm <- Influs.bet %*% LL3.pre

  ## influence function for Lam (individual level) and output
  Influs.Lam <- LL1.tm - LL2.tm - LL3.tm

  ## final influence functions
  return(Influs.Lam)

}




#==========================================================================#
# Calculate Subgroup survival rates using KM
#==========================================================================#
#' @title Subgroup survival rates
#'
#' @description Calculate Subgroup survival rates basedon the Kaplan-Meier estimation procedure
#'
#' @aliases St.Sub.KM
#'
#' @param tstar time points that the survival rates will be estimated at.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param G a matrix used to indicate which subgroups he/she belongs to for each of these subjects.
#'
#' @export St.Sub.KM
St.Sub.KM <- function(tstar,yobs,delta,G){

  K <- nrow(G)
  tSP <- rep(0, K)
  for(k in 1:K){
    idx <- (G[k,]==1)
    fit.k <- summary(survival::survfit(survival::Surv(yobs, delta) ~ 1, subset=idx))
    tSP[k] <- min( c(1,fit.k$surv)[ c(0,fit.k$time) <= tstar ] )
  }
  return( tSP )
}


#==== Sub-function: disturbing function used in the bootstrap procedure ====# (insert the elements of this function into the main body)
#' @title Disturbing function used in the bootstrap procedure
#'
#' @description Disturb the auxiliary in the bootstrap procedure
#'
#' @aliases AuxSP.disturb.KM
#'
#' @param nboot specifies the number of bootstrap sampling. The default \code{nboot = 100}.
#' @param V.KM the estimated variance-covariance matrix of the auxiliary information based on Kaplan-Meier estimator
#' @param auxinfo a matrix that collects the auxiliary information in a pre-specified way.
#'
#' @export AuxSP.disturb.KM
AuxSP.disturb.KM <- function(nboot,V.KM,auxinfo){

  ScaleP <- diag(1/sqrt(auxinfo[,'M']))
  disturbs <- mvtnorm::rmvnorm(nboot,mean=rep(0,nrow(auxinfo)),
                               sigma=ScaleP %*% V.KM %*% ScaleP)
  return(disturbs)

}








