#' Estimate
#'
#' This function does the following
#'
#' @param Y Matrix of the output variable rows=length of time series, cols=number of individuals
#' @param AX Matrix of the input variable rows=length of time series, cols=number of individuals x number of explanatory variables
#' @param LAB vector of group membership of each individual
#' @param k number of?
#' @param ks vector of number of group specific factors
#' @param lambda regularization parameter
#' @return A list containing PredL, PredXB, PredG and Bes
#' @export
Estimation_Function<-function(Y,AX,LAB,k,ks,lambda=0.1){
  num_labels <- length(ks)
  N <- nrow(Y)
  Nc <- ncol(Y)
  p <- ncol(AX)/Nc
  AY <- Y
  
  PredXB <- matrix(0,nrow=N,ncol=Nc)
  
  Bes <- 0*matrix(runif((p+1)*Nc,-1,1),ncol=Nc)
  
  for(j in 1:Nc){
    X <- AX[,(p*(j-1)+1):(p*j)]; y <- Y[,j]
    fit <- ncvreg::ncvreg(X, y, family="gaussian", penalty="SCAD",lambda=c(1,lambda))
    Bes[,j] <- (fit$beta)[,2]
    PredXB[,j] <- cbind(1,X)%*%as.vector(Bes[,j])
  }
  
  Y <- AY-PredXB
  
  VEC <- eigen(Y%*%t(Y))$vectors
  F <- sqrt(N)*(VEC)[,1:k]
  L <- t(t(F)%*%Y/N)
  FG <- F; LG <- L
  
  PredG <- matrix(0,nrow=N,ncol=Nc)
  
  for(j in 1:Nc){PredG[,j] <- FG%*%LG[j,]}
  
  Y <- AY-PredXB-PredG
  FS <- matrix(0,nrow=N,ncol=sum(ks))
  LS <- matrix(0,nrow=Nc,ncol=max(ks))
  
  for(i in 1:length(ks)){
    index <- subset(1:(Nc),LAB==i)
    Z <- Y[,index]
    VEC <- eigen(Z%*%t(Z))$vectors
    F <- sqrt(N)*(VEC)[,1:ks[i]]
    L <- t(t(F)%*%Z/N)
    
    LS[index,1:ks[i]] <- L
    if(i==1){FS[,1:ks[1]] <- F}
    if(i!=1){FS[,(sum(ks[1:(i-1)])+1):(sum(ks[1:i]))] <- F}
    
  }
  
  PredL <- matrix(0,nrow=N,ncol=Nc)
  
  for(i in 1:length(ks)){
    index <- subset(1:(Nc),LAB==i)
    L <- LS[index,1:ks[i]]
    if(i==1){F <- FS[,1:ks[1]]}
    if(i!=1){F <- FS[,(sum(ks[1:(i-1)])+1):(sum(ks[1:i]))]}
    PredL[,index] <- F%*%t(L)
  }
  
  
  #########################Estimation
  
  
  for(ite in 1:100){
    
    Bes.old <- Bes
    
    ###----------Beta
    
    Y <- AY-PredG-PredL
    
    for(j in 1:Nc){
      X <- AX[,(p*(j-1)+1):(p*j)]; y <- Y[,j]
      fit <- ncvreg::ncvreg(X, y, family="gaussian", penalty="SCAD",lambda=c(1,lambda))
      Bes[,j] <- (fit$beta)[,2]
      PredXB[,j] <- cbind(1,X)%*%as.vector(Bes[,j])
    }
    
    ###----------Global Factors
    
    Y <- AY-PredXB-PredL
    
    VEC <- eigen(Y%*%t(Y))$vectors
    F <- sqrt(N)*(VEC)[,1:k]
    L <- t(t(F)%*%Y/N)
    
    FG <- F; LG <- L
    
    PredG <- matrix(0,nrow=N,ncol=Nc)
    
    for(j in 1:Nc){PredG[,j] <- FG%*%LG[j,]}
    
    ###----------Local Factors
    
    Y <- AY-PredXB-PredG
    
    FS <- matrix(0,nrow=N,ncol=sum(ks))
    LS <- matrix(0,nrow=Nc,ncol=max(ks))
    
    for(i in 1:length(ks)){
      
      index <- subset(1:(Nc),LAB==i)
      
      Z <- Y[,index]
      VEC <- eigen(Z%*%t(Z))$vectors
      F <- sqrt(N)*(VEC)[,1:ks[i]]
      L <- t(t(F)%*%Z/N)
      
      LS[index,1:ks[i]] <- L
      if(i==1){FS[,1:ks[1]] <- F}
      if(i!=1){FS[,(sum(ks[1:(i-1)])+1):(sum(ks[1:i]))] <- F}
      
    }
    
    PredL <- matrix(0,nrow=N,ncol=Nc)
    
    for(i in 1:length(ks)){
      index <- subset(1:(Nc),LAB==i)
      L <- LS[index,1:ks[i]]
      if(i==1){F <- FS[,1:ks[1]]}
      if(i!=1){F <- FS[,(sum(ks[1:(i-1)])+1):(sum(ks[1:i]))]}
      PredL[,index] <- F%*%t(L)
    }
    
    if(sum( abs(Bes.old-Bes) )<=10^-2){break}
    
  }
  return(list("PredL"=PredL,"PredXB"=PredXB,"PredG"=PredG,"Bes"=Bes))
}
