iGS.eps <- function(X=NULL,K=NULL,Y=NULL,Var=NULL,alg="ai"){
  # Objects: prediction using AI-REML, may use multiple random effects
  # Design marix for random effects is identy matrix
  
  # Input: 
  #     K: list, Kinship with inference and reference
  #     mm: list, generated by calcAIvaf
  #     Y: matrix, n by 2, inference individuals values are NA
  #     X: n by p, p is the number of covariates
  
  # Output:
  #     BLUP: blup values for all individuals; use the same function as GAPIT
  #     phe: predicted phenotype for all individuals
  
  y = Y[,2]
  index.y = which(!is.na(y))
  n.random = length(K)
  if(n.random==0) stop("This no random effect! Kinship Matrix is needed.")
  for( i in 1:n.random){
    A = list(A=K[[i]][index.y,index.y])
    if(i==1){
      Knew= A
    }else{
      Knew = append(Knew,A)
    }
  }
  if(!is.null(X)){
    Xnew = as.matrix(X[index.y,])
  }else{
    Xnew = NULL
  }
  ynew = as.matrix(y[index.y])
  mm = calcAIvar(y = ynew,B = Knew, X = Xnew,Var=Var,alg=alg)
  gblup = gBLUP.eps(y = y, K = K, X = X, mm = mm)
  BLUP = gblup$BLUP
  phe = gblup$phe
  Var = mm$Var
  LL = mm$llik
  return(list(BLUP = BLUP, phe = phe, Var = Var, LL = LL))
}
