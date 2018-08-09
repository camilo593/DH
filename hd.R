VARmakexy <- function(DATA,lags,c_case){
  
  nobs <- nrow(DATA)
  
  #Y matrix 
  Y <- DATA[(lags+1):nrow(DATA),]
  Y <- DATA[-c(1:lags),]
  
  #X-matrix 
  if (c_case==0){
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    } 
  } else if(c_case==1){ #constante
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    X <- cbind(matrix(1,(nobs-lags),1), X) 
  } else if(c_case==2){ # tendencia y constante
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    trend <- c(1:nrow(X))
    X <-cbind(matrix(1,(nobs-lags),1), t(trend))
  }
  A <- (t(X) %*% as.matrix(X)) 
  B <- (as.matrix(t(X)) %*% as.matrix(Y))
  
  Ft <- ginv(A) %*% B
  
  retu <- list(X=X,Y=Y, Ft=Ft)
  return(retu)
}


#funcion para crear la matriz de compañia
companionmatrix <- function (x) 
{
  if (!(class(x) == "varest")) {
    stop("\nPoner un VAR generado por el paquete VARS.\n")
  }
  K <- x$K #dimension del var estimado
  p <- x$p #numero de rezagos
  A <- unlist(Acoef(x)) # sacamos los coeficinetes estimados del VAR
  companion <- matrix(0, nrow = K * p, ncol = K * p) # se crea una matriz vacia con las dimensiones necesarias
  companion[1:K, 1:(K * p)] <- A # se llena con los coeficinetes sacados y se empieza a llenar la matriz de compañia
  if (p > 1) {
    j <- 0
    for (i in (K + 1):(K * p)) {
      j <- j + 1
      companion[i, j] <- 1
    }
  }
  return(companion)
}




VARhd <- function(Estimation){
  
  ## make X and Y
  nlag    <- Estimation$p   # numero de rezagos
  DATA    <- Estimation$y   # datos
  QQ      <- VARmakexy(DATA,nlag,1) #dejamos en caso uno porque vamos a correr con constante
  
  
  
  invA    <- t(chol(as.matrix(summary(Estimation)$covres)))   # inversa de la matriz A (covarianza residuos) se hace con  la descomposición de cholesky (triangular baja)
  Fcomp   <- companionmatrix(Estimation)                      # Matriz de compañia con la funcion que creamos arriba para el VAR
  
  #det     <- c_case                                          # constante o tendencia
  F1      <- t(QQ$Ft)                                         
  eps     <- ginv(invA) %*% t(residuals(Estimation))          # Los errores estructurales multiplicamos la matriz ortogonalizada por los residuos  
  nvar    <- Estimation$K                                     # variables endoges
  nvarXeq <- nvar * nlag                                      # numero de endogenas rezagadas por ecuacion
  nvar_ex <- 0                                                # exogenas
  Y       <- QQ$Y                                             
  #X       <- QQ$X[,(1+det):(nvarXeq+det)]                    
  nobs    <- nrow(Y)                                          # observaciones
  
  
  ## Calculo de la descomposición histórica
  
  # contribución de cada shock
  invA_big <- matrix(0,nvarXeq,nvar)
  invA_big[1:nvar,] <- invA
  Icomp <- cbind(diag(nvar), matrix(0,nvar,(nlag-1)*nvar))
  HDshock_big <- array(0, dim=c(nlag*nvar,nobs+1,nvar))
  HDshock <- array(0, dim=c(nvar,(nobs+1),nvar))
  
  for (j in 1:nvar){  # para cada variable
    eps_big <- matrix(0,nvar,(nobs+1)) 
    eps_big[j,2:ncol(eps_big)] <- eps[j,]
    for (i in 2:(nobs+1)){
      HDshock_big[,i,j] <- invA_big %*% eps_big[,i] + Fcomp %*% HDshock_big[,(i-1),j]
      HDshock[,i,j] <-  Icomp %*% HDshock_big[,i,j]
    } 
    
  } 
  
  HD.shock <- array(0, dim=c((nobs+nlag),nvar,nvar))   # Es necesario crear un objeto con tres dimensiones observaciones, shock y variable
  
  for (i in 1:nvar){
    
    for (j in 1:nvar){
      HD.shock[,j,i] <- c(rep(NA,nlag), HDshock[i,(2:dim(HDshock)[2]),j])
    }
  }
  
  return(HD.shock)
  
}