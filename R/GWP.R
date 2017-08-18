`GWP` <- function(Y = NULL, LD= 0.7,GM = NULL, SNPnum = NULL, K= NULL, nrep = 100, nfold = 5, GD = NULL,
                  GDcali = NULL, GDvali = NULL,GWAS.method = NULL, GS.method = "gBLUP",CV = NULL,h2 = NULL,maxLoop = 10,
                  nIter= 20000,burnIn = 15000,thin = 10, saveAt=''){

  # Objects: Genome Wide Prediction, do GS analysis based on GWAS results
  # Input:
  ##   GWAS.method: GWAS methods,GLM, FarmCPU, Blink
  ##   GS.method: GS methods, gBLUP, rrBLUP, Bayes A, Bayes B, Bayes Cpi, Bayes LASSO
  ## 	 GD: Genotype, n by m, all genotype, will do cross validation
  ## 	 GDcali: Genotype for calibration in GS; if provided, GDvali must be provided, and will ignore GD
  ##			 Y must be the same order as GDcali
  ## 	 GDvali: Genotype for validation in GS; if provided, GDcali must be provided, and will ignore GD
  ## 	 Y: phenotype data, n x (p+1) matrix; n individuals, p phenotypes
  ## 	 nfold: only valid when GD exist, for cross validation;
  ##   nrep: only valid when gd exist, for cross validation;
  ##   LD: threshold used for filtering SNPs;
  ##   SNPnum: # of SNPs will be kept for GS analysis, gBLUP and rrBLUP will use all markers and ignore this
  ##   K: n by n, the same order as Y, kinship used for gBLUP and rrBLUP
  ##   h2: heritability, should be the same length as phenotypes
  ##	 other values are the same for GWAS methods
  # Output:
  ##   EBV: Estimated BLUP value, only provided when GDcali used
  ##	 accuracy: accuracy of GS
  ##   fit: model fitness
  ##   accuracy_se: standard error of accuracy, only provided for cross validation
  ##   fit_se: standard error of fitness, only provided for cross validation
  # Authors: Yao Zhou
  # Last update: 01/25/2017

  ## Detection
  try( if(is.null(GD) & is.null(GDcali))  stop("GD and GDcali must be provided at least one"))
  if(!is.null(GD) & !is.null(GDcali)) warning("GD and GDcali both provided, only GD will be used")
  try(if(is.null(GD) & !is.null(GDcali) & is.null(GDvali)) stop("GDcali and GDvali must be provided at the same time"))
  if(!is.null(CV)) try(if(nrow(Y)!=nrow(CV)) stop("CV and Y must be provided at the same individuals"))
  ## Before analysis
  PheNum = ncol(Y) - 1
  PheName = colnames(Y)[-1]
  if(is.null(h2)|length(h2)!=PheNum) h2 = rep(0.5, PheNum)
  # diagnose
  if(!is.null(GD)) try(if(ncol(GD)!=nrow(Y)) stop("Number of individuals in Y and GD do not match"))
  if(is.null(GD)) try(if(ncol(GDcali)!=nrow(Y)) stop("Number of individuals in Y and GDcali do not match"))
  if(!is.null(K)) try(if(ncol(K)!=nrow(Y)) stop("Number of individuals in K and Y do not match"))
  ## do all phenotypes provided
  for( phe in 1:PheNum){
    cat("Phenotype: ",PheName[phe],"\n")
    ## prepare phenotype data
    myY = Y[,c(1,(1+phe))]
    n = sum(!is.na(myY[,2]))
    index.y = !is.na(myY[,2])
    myY = myY[index.y,]
    myCV = CV[index.y,]

    if(!is.null(GD) & is.null(K) & (GS.method == "gBLUP")) K = A.mat(t(as.matrix(GD)))
    if (!is.null(GWAS.method)){
      if(is.null(K) & GWAS.method =="MLM"){
        if (!is.null(GD)) {
          K = A.mat(t(as.matrix(GD)))
        }else{
          K = A.mat(t(as.matrix(GDcali)))
        }
      }
    }
    if(is.null(GD) & !is.null(GDcali) & is.null(K) & GS.method == "gBLUP") K = A.mat(t(as.matrix(GDcali)))
    myK <- K[index.y,index.y]
    if(!is.null(K)) myK = as.matrix(myK)
    if(!is.null(GD)){
      myGD = deepcopy(GD, cols = index.y)
      acc = matrix(NA,nrep*nfold,4)
      for ( rep in 1:nrep){
        sets=sample(cut(1:n,nfold,labels=FALSE),n)
        for ( i in 1:nfold){
          ## GWAS
          index.ref <- (sets != i)
          Zcali = deepcopy(myGD, cols = index.ref)
          Zvali = deepcopy(myGD, cols = !index.ref)
          Ycali = myY[index.ref,]
          Yvali = myY[!index.ref,]
          CVcali = myCV[index.ref,]
          # change GWAS methods
          if(!is.null(GWAS.method)& GS.method!="gBLUP" & GS.method != "rrBLUP"){
            if(GWAS.method == "GLM"){
              maxLoop = 1
              GWAS = Blink(CV = CVcali, Y = Ycali,GD = Zcali,GM = GM,maxLoop = maxLoop,model = "A",time.cal = T,BIC.method = "naive",file.output = F)
              GWAS = GWAS$GWAS
            } else if(GWAS.method == "Blink"){
              GWAS = Blink(CV = CVcali, Y = Ycali,GD = Zcali,GM = GM,maxLoop = maxLoop,model = "A",time.cal = T,BIC.method = "naive",file.output = F)
              GWAS = GWAS$GWAS
            } else if(GWAS.method == "MLM"){
              pheno = myY
              pheno[!index.ref,2] = NA
              geno = as.matrix(myGD)
              colnames(geno) = as.character(pheno[,1])
              geno = cbind(myGM,geno)
              GWAS = GWAS(pheno = pheno, geno = geno, K = K,plot=F)
            } else{
              print("GWAS method not provided")
            }

            ##	order
            GWAS[is.na(GWAS[,4]),4] = 1
            if(is.null(SNPnum)){
              lenSNP = ceiling(nrow(myY)/log(nrow(myY)))
            }else{
              lenSNP = SNPnum
            }
            index.col = order(GWAS[,4],decreasing =F, na.last = T)
            ## LD remove
            add = TRUE
            j = 1
            while (add == TRUE ){
              lend = 5*j*lenSNP
              if(lend > length(index.col)){
                lend = length(index.col)
                add = FALSE
              }else{
                j = j + 1
              }
              Porder = index.col[1:lend]
              GDneo = deepcopy(Zcali,rows = Porder)
              Psort=Blink.LDRemove(Porder=Porder,GDneo = GDneo,LD = LD,orientation="row")
              if(length(Psort) > lenSNP){
                Psort = Psort[1:lenSNP]
                break
              }
            }
            Acali = deepcopy(Zcali, rows = Psort)
            Zcali = as.matrix(Acali)
            Zcali = t(Zcali)
            GDvali = deepcopy(Zvali, rows= Psort)
            GDvali = as.matrix(GDvali)
            Zvali = t(GDvali)
          }else{
            Zcali = as.matrix(Zcali)
            Zcali = t(Zcali)
            Zvali = as.matrix(Zvali)
            Zvali = t(Zvali)
          }
          ## GS
          if(GS.method == "gBLUP"){
            ynew <- myY
            ynew[!index.ref,2] = NA
            rrmod <- mixed.solve(y=ynew[,2], K = myK,X = myCV)
            g_hat <- rrmod$u
            prediction.vali <- g_hat[!index.ref]
            prediction.cali <- g_hat[index.ref]
            accuracy <- cor(Yvali[,2],prediction.vali)
            fit <- cor(Ycali[,2],prediction.cali)
            acc[((rep-1)*nfold+i),1] = accuracy
            acc[((rep-1)*nfold+i),2] = fit
          }else if(GS.method == "rrBLUP"){
            ynew <- myY
            ynew[!index.ref,2] = NA
            rrmod <- mixed.solve(y=ynew[,2], Z = t(as.matrix(myGD)),X = myCV)
            g_hat <- t(as.matrix(myGD)) %*% rrmod$u
            prediction.vali <- g_hat[!index.ref]
            prediction.cali <- g_hat[index.ref]
            accuracy <- cor(Yvali[,2],prediction.vali)
            fit <- cor(Ycali[,2],prediction.cali)
            acc[((rep-1)*nfold+i),1] = accuracy
            acc[((rep-1)*nfold+i),2] = fit
          }else{
            R2 = h2[phe]
            if(GS.method == "BayesA"){
              ETA<-list(list(X = Zcali,model='BayesA'))
              if(!is.null(CVcali)) ETA<-list(list(X = CVcali,model="FIXED"),list(X = Zcali,model='BayesA'))
            }else if(GS.method == "BayesB"){
              ETA<-list(list(X = Zcali,model='BayesB',probIn=0.8))
              if(!is.null(CVcali)) ETA<-list(list(X = CVcali,model="FIXED"),list(X = Zcali,model='BayesB',probIn=0.8))
            }else if(GS.method == "BayesCpi"){
              ETA<-list(list(X = Zcali,model='BayesCpi'))
              if(!is.null(CVcali)) ETA<-list(list(X = CVcali,model="FIXED"),list(X = Zcali,model='BayesCpi'))
            }else if(GS.method =="BayesLASSO"){
              ETA<-list(list(X = Zcali,model='BL'))
              if(!is.null(CVcali)) ETA<-list(list(X = CVcali,model="FIXED"),list(X = Zcali,model='BL'))
            }else{
              stop("GS method only support limited methods: gBLUP, rrBLUP, BayesA, BayesB, BayesCpi, BayesLASSO")
            }
            fit = BGLR(y=Ycali[,2],ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=NULL,weights=NULL,R2=R2)
            prediction.vali <- Zvali%*%fit$ETA[[1]]$b
            prediction.cali <- Zcali%*%fit$ETA[[1]]$b
            accuracy <- cor(Yvali[,2],prediction.vali)
            fit <- cor(Ycali[,2],prediction.cali)
            acc[((rep-1)*nfold+i),1] = accuracy
            acc[((rep-1)*nfold+i),2] = fit
          }
          cat( "Repetition:", rep, "\n")
          cat("Iteration:", i, "\n")
        }

      }
      phe.name = list(BLUP = NULL, accuracy = acc)
    }else{
      Zcali = GDcali
      Zvali = GDvali
      if(!is.null(GWAS.method)& GS.method!="gBLUP" & GS.method != "rrBLUP"){
        if(GWAS.method == "GLM"){
          maxLoop = 1
          GWAS = Blink(CV = CV, Y = Y,GD = Zcali, GM = GM,maxLoop = maxLoop,model = "A",time.cal = T,BIC.method = "naive",file.output = F)
          GWAS = GWAS$GWAS
        }else if(GWAS.method == "Blink"){
          GWAS = Blink(CV = CV, Y = Y,GD = Zcali, GM = GM,maxLoop = maxLoop,model = "A",time.cal = T,BIC.method = "naive",file.output = F)
          GWAS = GWAS$GWAS
        }else if(GWAS.method == "MLM"){
          pheno = Y
          geno = as.matrix(Zcali)
          colnames(geno) = as.character(pheno[,1])
          geno = cbind(myGM,geno)
          GWAS = GWAS(pheno = pheno, geno = geno, K = K,plot=F)
        } else{
          print("GWAS method not provided")
        }

        ##	order
        GWAS[is.na(GWAS[,4]),4] = 1
        if(is.null(SNPnum)){
          lenSNP = ceiling(nrow(myY)/log(nrow(myY)))
        }else{
          lenSNP = SNPnum
        }
        index.col = order(GWAS[,4],decreasing =F, na.last = T)
        ## LD remove
        add = TRUE
        j = 1
        while (add == TRUE ){
          lend = 5*j*lenSNP
          if(lend > length(index.col)){
            lend = length(index.col)
            add = FALSE
          }else{
            j = j + 1
          }
          Porder = index.col[1:lend]
          GDneo = deepcopy(Zcali,rows = Porder)
          Psort=Blink.LDRemove(Porder=Porder,GDneo = GDneo,LD = LD,orientation="row")
          if(length(Psort) > lenSNP){
            Psort = Psort[1:lenSNP]
            break
          }
        }
        Acali = deepcopy(Zcali, rows = Psort)
        Zcali = as.matrix(Acali)
        Zcali = t(Zcali)
        GDvali = deepcopy(Zvali, rows= Psort)
        GDvali = as.matrix(GDvali)
        Zvali = t(GDvali)
      }else{
        Zcali = as.matrix(Zcali)
        Zcali = t(Zcali)
        Zvali = as.matrix(Zvali)
        Zvali = t(Zvali)
      }
      ## GS
      print(dim(Zcali))
      if(GS.method == "gBLUP"){
        ynew = matrix(NA,(nrow(Zvali)+nrow(Zcali)),2)
        ynew[1:nrow(Zcali),2] = myY[,2]
        myK = A.mat(cbind(Zcali,Zvali))
        rrmod <- mixed.solve(y=ynew[,2], K = myK,X = myCV)
        g_hat <- GDvali %*% rrmod$u
        mean.pop <- rrmod$beta
        prediction.vali <- g_hat[!index.ref]
        BLUP = prediction.vali
      }else if(GS.method == "rrBLUP"){
        ynew <- myY
        rrmod <- mixed.solve(y=ynew[,2], Z = t(as.matrix(myGD)),X = myCV)
        g_hat <- GDvali %*% rrmod$u
        mean.pop <- rrmod$beta
        prediction.vali <- g_hat
        BLUP = prediction.vali
      }else{
        R2 = h2[phe]
        if(GS.method == "BayesA"){
          ETA<-list(list(X = Zcali,model='BayesA'))
          if(!is.null(CV)) ETA<-list(list(X = CVcali,model="FIXED"),list(X = Zcali,model='BayesA'))
        }else if(GS.method == "BayesB"){
          ETA<-list(list(X = Zcali,model='BayesB',probIn=0.8))
          if(!is.null(CV)) ETA<-list(list(X = CVcali,model="FIXED"),list(X = Zcali,model='BayesB',probIn=0.8))
        }else if(GS.method == "BayesCpi"){
          ETA<-list(list(X = Zcali,model='BayesCpi'))
          if(!is.null(CV)) ETA<-list(list(X = CVcali,model="FIXED"),list(X = Zcali,model='BayesCpi'))
        }else if(GS.method =="BayesLASSO"){
          ETA<-list(list(X = Zcali,model='BL'))
          if(!is.null(CV)) ETA<-list(list(X = CVcali,model="FIXED"),list(X = Zcali,model='BL'))
        }else{
          stop("GS method only support limited methods: gBLUP, rrBLUP, BayesA, BayesB, BayesCpi, BayesLASSO")
        }
        fit = BGLR(y=myY[,2],ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=NULL,weights=NULL,R2=R2)
        BLUP <- Zvali%*%fit$ETA[[1]]$b
      }
      phe.name = list(BLUP=BLUP,accuracy = NULL)
    }
    if(phe == 1){
      result = phe.name
    }else{
      result = c(result,phe.name)
    }
  }
  return(result)
}
