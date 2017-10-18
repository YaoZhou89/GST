# general Mixed model code for variance component analysis

# y named(!) phenotype vector (n)
# B named list of covariance structures (m x m matrices)
# X a incidence matrix for fixed effects (n x p fixed effects)
# Z a incidence matrix for random effects (n x m) 
# Ze a incidence matrix for random effects for error term only
# reps  (T or F, are there replicates)
# alg derivative based algorithm ('ai' = average information, 'fs' = fisher scoring, or 'nr' = newton-rhapson (fastest to slowest), sqrt  )
# conv.val value for convergence
# Var vector of initialization values for variance components

library(Matrix)
# gpu only works for smaller datasets, otherwise vid card doesn't have enough memory
# advise not using
#library(gputools)

calcAIvar = function(y, B=NULL,X=NULL, Z=NULL, Ze=NULL, reps=FALSE,
                  alg='ai', conv.val=1e-4, Var=NULL){
  doGPU = FALSE
  # y = as.matrix(y)
  n = length(y)
  #for constructing Strain Variance component
  Ind = Matrix(diag(n), sparse=T)
  
  # If B is NULL then add a covariance term that is the identity matix - will then calculate effect of strain (broad-sense heritability)
  if(is.null(B)) {
    B=list(Ind=Ind);
  }# If Z is null and there are no replicates this will make Z a diagonal incidence matrix, otherwise this constructs an incidence matrix based on strain names
  if(is.null(Z)) {   Z = Ind }
  # If X is null assume one fixed effect of population mean
  if(is.null(X) ) {  X = model.matrix(y~1)}
  # If Ze is null assume no error structure
  if(is.null(Ze)) {  Ze = Ind }
  
  #number of terms in the structured covariance
  alg.test = TRUE
  while(alg.test){
    alg.test = F
    VC.names=paste('sigma', c(names(B), 'E'), sep='')
    N.s.c=length(B)
    Vcmp.cnt=N.s.c+1
    # starting values for VC estimates as 1 / (#of VCs including residual error term)
    if(is.null(Var) ) { Var=rep(1/Vcmp.cnt, Vcmp.cnt) }
    I = matrix(0, ncol= Vcmp.cnt, nrow= Vcmp.cnt)
    s = matrix(0, ncol=1, nrow= Vcmp.cnt)
    
    diffs=rep(10,  Vcmp.cnt)
    # second derivatives of V with respect to the variance components (Lynch and Walsh 27.15)
    
    VV = list()
    for(i in 1:N.s.c) {
      VV[[i]]=Z %*% tcrossprod(B[[i]],Z)
    }
    VV[[ Vcmp.cnt ]]=Ze 
    
    i = 0
    theLoop = 100
    # while the differences haven't converged 
    var.old = 0
    alg.all = c("fs","ai","nr")
    while ( sum(ifelse(diffs<conv.val, TRUE,FALSE)) <  Vcmp.cnt & i < theLoop) { 
      i = i + 1
      V=matrix(0,n, n)
      for( vcs in 1:length(VV)) {  V=V+(VV[[vcs]]*Var[vcs]) }
      print('Inverting V')
      Vinv = solve(V)
      print('Done inverting V')
      tXVinvX=t(X) %*% Vinv %*% X
      # print(tXVinvX)
      inv.tXVinvX =  ginv(as.matrix(tXVinvX))
      # ginvXVX = inv.tXVinvX # ginvXVX = inv.tXVinvX
      VinvX = Vinv %*% X
      itv = inv.tXVinvX %*% t(X)%*%Vinv
      P = Vinv - Vinv %*% X %*% itv
      # P = Vinv - (VinvX %*% (ginvXVX %*% t(VinvX))) 
      if(alg=='fs') {print("Fisher scoring algorithm: calculating expected VC Hessian") }
      if(alg=='nr') {print("Netwon rhapson algorithm: calculating observed VC Hessian") }
      if(alg=='ai') {print("Average information algorithm: calculating avg of expected and observed VC Hessians") }
      Py = P %*% y
      for(ii in 1:Vcmp.cnt) {
        for(jj in ii:Vcmp.cnt) {
          if (alg=='fs') {    I[ii,jj]= 0.5*sum(diag( ((P%*%VV[[ii]]) %*%P )%*%VV[[jj]])) }
          if (alg=='nr') {    I[ii,jj]=-0.5*sum(diag(P%*%VV[[ii]]%*%P%*%VV[[jj]])) + (t(y)%*%P%*%VV[[ii]]%*%P%*%VV[[jj]]%*%P%*%y)[1,1] }
          if (alg=='ai') {
            if(doGPU)  {    I[ii,jj]= 0.5*gpuMatMult(gpuMatMult(gpuMatMult(gpuMatMult(gpuMatMult(gpuMatMult(t(y),P),VV[[ii]]),P),VV[[jj]]),P),y) }
            else{
              I[ii,jj]= 0.5*( t(Py)%*%VV[[ii]]%*%P%*%VV[[jj]]%*%Py)[1,1]
            }
          }
          print(paste(ii, jj))
          I[jj,ii]=I[ii,jj]
        }
      }
      for(ii in 1:Vcmp.cnt){
        s[ii,1]= -0.5*sum(diag(P%*%VV[[ii]])) + .5*(t(Py)%*%VV[[ii]]%*%Py )[1,1]
        # s[ii,1]<- - 0.5 * c(sum(apply(P * VV[[ii]], 1,sum)) - t(Py) %*% (VV[[ii]] %*% Py), sum(diag(P)) -  sum(Py^2))
      }
      print(s)
      # s[ii,1]<- - 0.5 * c(sum(apply(P * VV[[ii]], 1,sum)) - t(Py) %*% (VV[[ii]] %*% Py), sum(diag(P)) -  sum(Py^2))
      # print("Average information algorithm: calculating avg of expected and observed VC Hessians") 
      # temp0 = P%*%y
      # for(ii in 1:Vcmp.cnt) {
      #   temp2  = P%*%VV[[ii]]
      #   temp1 = t(y)%*% temp2
      #   for(jj in ii:Vcmp.cnt) {
      #     I[ii,jj]= 0.5*(temp1%*%P%*%VV[[jj]]%*%temp0)[1,1] 
      #     print(paste(ii, jj))
      #     I[jj,ii]=I[ii,jj]
      #   }
      #   s[ii,1]= -0.5*sum(diag(temp2)) + .5*(temp1%*%temp0 )[1,1]
      # }
      cat("I is:","\n",I,"\n")
      invI = ginv(I)
      newVar = Var + invI%*%s
      #newVar[newVar<0]=2.6e-9
      newVar[newVar<0] = 2.6e-9
      for(d in 1:length(diffs)) { diffs[d]=abs(Var[d] - newVar[d]) }
      Var = newVar
      
      cat('\n')
      cat("iteration ", i, '\n')
      cat(VC.names, '\n')
      cat(Var, '\n')
      
      Bhat= itv %*% y
      cat("Fixed Effects, Bhat = ", as.matrix(Bhat), '\n')
      det.tXVinvX=determinant(tXVinvX, logarithm=TRUE)
      det.tXVinvX=det.tXVinvX$modulus*det.tXVinvX$sign
      det.V =determinant(V, logarithm=TRUE)
      det.V=det.V$modulus*det.V$sign
      LL = -.5 * (det.tXVinvX + det.V + log(t(y) %*% P %*% y) )
      cat("Log Likelihood = " , as.matrix(LL), '\n')
      cat("VC convergence vals", '\n')
      cat(diffs, '\n')
      if ( i ==50 ){
        d1 = diffs
        alg.all = setdiff(alg.all,alg)
      }
      if (i == 51){
        d2 = diffs
      }
      if(!is.na(alg.all)){
        if (i == 52 & (sum(d2>d1)>1)){
          alg = alg.all[1]
          alg.test = TRUE
        }
      }
    }
  }
  return(list(Var=Var, invI=invI, W=Vinv, Bhat=Bhat, llik=LL))
}



