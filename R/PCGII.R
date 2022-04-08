require(tidyverse)
require(corpcor)
require(GeneNet)
require(FastGGM)
require(glmnet)

unMat=function(X_est, X_p){
  # X is a p by p matrix
  p=dim(X_est)[2]
  out=matrix(NA, nrow=p*(p-1)/2, ncol=4)
  ind=1
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      out[ind,1]=i
      out[ind,2]=j
      out[ind,3]=X_est[i,j]
      out[ind,4]=X_p[i,j]
      ind=ind+1
    }
  }
  colnames(out)=c("node1","node2","pcor","pval")
  out
}

sigs2vec=function(sigs, P){
  require(corpcor)
  m=matrix(0,P,P)
  for (h in 1: dim(sigs)[1]){
    m[sigs[h,1],sigs[h,2]]=1
  }
  sm2vec(m)
}


PCGII_estimating=function(df, IIC=TRUE, prior){
  n = dim(df)[1]; p = dim(df)[2]
  t0=2
  IndMatrix = matrix(1, p, p) - diag(rep(1, p))
  Eresidual = matrix(0, n, p) # regression residuals matrix n*p
  CoefMatrix = matrix(0, p, p - 1) # regression coefficient matrix p*p-1
  meanX = colMeans(df)
  X = t(t(df) - meanX)
  XS = matrix(0, n, p)
  # XS: Standardized X
  for (i in 1 : p){
    XS[, i] = X[, i] / sd(X[, i])
  }
  
  colnames(X)=colnames(df)
  colnames(XS)=colnames(df)
  
  shat=sqrt(n/(log(p)^3))
  lambda=sqrt(2*(2+0.01)*log(p/shat)/n)
  # lambda=sqrt(2*log(p)/n)
  default_penalty=rep(1,p-1)
  
  if (IIC==T){
    for (i in 1 : p){
      penalty_fac=default_penalty
      temp.node=prior[with(prior,row==i),'col']
      for(nds in temp.node){
        if (nds < i) penalty_fac[nds]=0 else penalty_fac[nds-1]=0
      }
      out = glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda, penalty.factor=penalty_fac)
      Coef = out$beta
      Predict = predict(out, XS[, -i], type = "link")
      CoefMatrix[i, ] = as.numeric(Coef) / apply(X[, -i], 2, sd)
      Eresidual[, i] = X[, i] - Predict
    }
  } else {
    for (i in 1 : p){
      out = glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda)
      Coef = out$beta
      Predict = predict(out, XS[, -i], type = "link")
      CoefMatrix[i, ] = as.numeric(Coef) / apply(X[, -i], 2, sd)
      Eresidual[, i] = X[, i] - Predict
    }
  }
  
  CovRes = t(Eresidual) %*% Eresidual / n # residuals covariance
  Est = matrix(1, p, p) # estimated partial correlation (rho hat in the paper )
  
  for (i in 1 : (p - 1)){
    for (j in (i + 1) : p){
      temp = Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1]
      Est[i, j] = mean(temp) / sqrt(diag(CovRes)[i] * diag(CovRes)[j])
      Est[j, i] = Est[i, j]
    }
  }
  
  # sparse partial correlation estimation with threshold (rho ~ in the paper)
  EstThresh = Est * ( abs(Est) >= (t0 * sqrt(log(p) / n) * IndMatrix) ) 

  kappa = (n / 3) * mean( colSums(Eresidual^4) / (colSums(Eresidual^2))^2 )  # forth moment, a number 
  
  SE=sqrt((kappa*(1-EstThresh^2))^2/n)
  
  tscore=Est/SE
  
  return(list(Est=Est,
              tscore=tscore,
              kappa=kappa,
              EstThresh=EstThresh,
              n=n, p=p))
}  

PCGII_estimating_cv=function(df, IIC=TRUE, prior, degree_freedom){
  require(glmnet)
  n = dim(df)[1]; p = dim(df)[2]
  t0=2
  IndMatrix = matrix(1, p, p) - diag(rep(1, p))
  Eresidual = matrix(0, n, p) # regression residuals matrix n*p
  CoefMatrix = matrix(0, p, p - 1) # regression coefficient matrix p*p-1
  meanX = colMeans(df)
  X = t(t(df) - meanX)
  XS = matrix(0, n, p) # XS: Standardized X
  for (i in 1 : p){
    XS[, i] = X[, i] / sd(X[, i])
  }
  
  colnames(X)=colnames(df)
  colnames(XS)=colnames(df)
  
  default_penalty=rep(1,p-1)
  
  if (IIC==T){
    for (i in 1 : p){
      penalty_fac=default_penalty
      temp.node=prior[with(prior,row==i),'col']
      for(nds in temp.node){
        if (nds < i) {penalty_fac[nds]=0} else penalty_fac[nds-1]=0
      }
      out.cv = cv.glmnet(XS[, -i], X[, i], family = "gaussian", nfolds = n, penalty.factor=penalty_fac)
      Coef = coef(out.cv,s= "lambda.min")[-1]
      Predict = as.matrix(predict(out.cv,XS[, -i], s= "lambda.min"))
      
      # if (length(which(coef(out.cv, s='lambda.min')[-1]!=0))<=degree_freedom){ # cv
      #   Coef = coef(out.cv,s= "lambda.min")[-1]
      #   Predict = as.matrix(predict(out.cv,XS[, -i], s= "lambda.min"))
      # } else{ # dfmax
      #   out = glmnet(XS[, -i], X[, i], family = "gaussian", dfmax = degree_freedom, penalty.factor=penalty_fac)
      #   Coef=coef(out,s=out$lambda[which.max(out$df[out$df<=degree_freedom])])[-1]
      #   Predict = as.matrix(predict(out,XS[, -i],s=out$lambda[which.max(out$df[out$df<=degree_freedom])]))
      # }
      CoefMatrix[i,] = Coef / apply(X[, -i], 2, sd)
      Eresidual[,i] = X[, i] - Predict 
    }
  } else {
    for (i in 1 : p){
      out.cv = cv.glmnet(XS[, -i], X[, i], family = "gaussian", nfolds = n)
      Coef = coef(out.cv,s= "lambda.1se")[-1]
      Predict = as.matrix(predict(out.cv,XS[, -i], s= "lambda.1se"))
      
      # if (length(which(coef(out.cv, s='lambda.min')[-1]!=0))<=degree_freedom){ # cv
      #   Coef = coef(out.cv,s= "lambda.min")[-1]
      #   Predict = as.matrix(predict(out.cv,XS[, -i], s= "lambda.min"))
      # } else{ # dfmax
      #   out = glmnet(XS[, -i], X[, i], family = "gaussian", dfmax = degree_freedom)
      #   Coef=coef(out,s=out$lambda[which.max(out$df[out$df<=degree_freedom])])[-1]
      #   Predict = as.matrix(predict(out,XS[, -i],s=out$lambda[which.max(out$df[out$df<=degree_freedom])]))
      # }
      CoefMatrix[i,] = Coef / apply(X[, -i], 2, sd)
      Eresidual[,i] = X[, i] - Predict 
    }  
  }

  CovRes = t(Eresidual) %*% Eresidual / n # residuals covariance
  Est = matrix(1, p, p) # estimated partial correlation (rho hat in the paper )
  
  for (i in 1 : (p - 1)){
    for (j in (i + 1) : p){
      temp = Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1]
      Est[i, j] = mean(temp) / sqrt(diag(CovRes)[i] * diag(CovRes)[j])
      Est[j, i] = Est[i, j]
    }
  }
  
  EstThresh = Est * ( abs(Est) >= (t0 * sqrt(log(p) / n) * IndMatrix) ) 

  kappa = (n / 3) * mean( colSums(Eresidual^4) / (colSums(Eresidual^2))^2 )  # forth moment, a number 
  SE=sqrt((kappa*(1-EstThresh^2))^2/n)
  tscore=Est/SE
  
  return(list(Est=Est,
              tscore=tscore,
              kappa=kappa,
              EstThresh=EstThresh,
              n=n, p=p))
}  


model.eval=function(X, omega, rep=5, degfree=20, prior_prop=0.05){

  p=dim(omega)[1] 
  truth<-sm2vec(omega) 
  TP<- which(truth!=0) 
  CN=p*(p-1)/2-length(TP) 
  CP=length(TP)
  
  prior=which(omega != 0, arr.ind = TRUE) %>% as.data.frame() %>% 
    subset(row != col) %>% sample_n(size=round(CP*prior_prop,1)) %>%
    transform(row = pmin(row, col), col = pmax(row, col)) %>% 
    arrange(row, col) 
  
  al=seq(0.001, 0.2005, 0.0035) # inference sig levels
  cutt=seq(0.001,1,0.001) # desired FPR
  
  {
    ### ROC 
    ENF.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    ENF.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    ENF.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    ENF.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    MLE.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    MLE.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    MLE.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    MLE.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    cLevel_theo.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_theo.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_theo.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_theo.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    cLevel_cv.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_cv.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_cv.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_cv.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_theo.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_theo.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_theo.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_theo.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_cv.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_cv.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_cv.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_cv.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    F_GGM.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    F_GGM.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    F_GGM.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    F_GGM.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    tempstore.ENF=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.MLE=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points 
    tempstore.cLevel_theo=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.cLevel_cv=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_theo=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_cv=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.F_GGM=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    
    ### Inference
    
    # Number of total correctly selected edges
    all.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    all.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.mle.p=matrix(Inf, nrow=length(al), ncol=rep) 
    all.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.cLevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    all.cLevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_theo=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_cv=matrix(Inf, nrow=length(al), ncol=rep)
    all.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # true positives
    tp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.cLevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    tp.cLevel_cv=matrix(Inf, nrow=length(al), ncol=rep)    
    tp.PCGII_theo=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_cv=matrix(Inf, nrow=length(al), ncol=rep)
    tp.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false positives
    fp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.cLevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    fp.cLevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_theo=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_cv=matrix(Inf, nrow=length(al), ncol=rep)
    fp.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false negatives
    tn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.cLevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    tn.cLevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_theo=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_cv=matrix(Inf, nrow=length(al), ncol=rep)
    tn.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false negatives
    fn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.cLevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    fn.cLevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_theo=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_cv=matrix(Inf, nrow=length(al), ncol=rep)
    fn.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
  }
  
  
  for (k in 1:rep){
    print(paste0(k,"th rep starts!!!!!"))
    sim.data=X[[k]]
    
    ### Fast GGM
    out=FastGGM(sim.data)
    F_GGM.test=unMat(X_est=out$partialCor, X_p=out$p_partialCor) # Est, pvals
    F_GGM.test=cbind.data.frame(F_GGM.test,truth=truth, qval=p.adjust(F_GGM.test[,4], method="BH"))
    print("Fast GGM done")
    
    ### cLevel
    cLevel_theo=CLEVEL_estimating(sim.data)
    cLevel_cv=CLEVEL_estimating_cv(sim.data,degfree)
    # estimates by cLevel
    Est_cLevel_theo=sm2vec(cLevel_theo$Est)
    Est_cLevel_cv=sm2vec(cLevel_cv$Est)
    # test statistics of cLevel
    tscore_cLevel_theo=sm2vec(cLevel_theo$tscore)
    tscore_cLevel_cv=sm2vec(cLevel_cv$tscore)
    print("cLevel done")
    
    ### PCGII
    PCGII_theo=PCGII_estimating(df=sim.data, prior = prior)
    PCGII_cv=PCGII_estimating_cv(df=sim.data, prior = prior, degree_freedom =  degfree)
    # estimates by PCGII
    Est_PCGII_theo=sm2vec(PCGII_theo$Est)
    Est_PCGII_cv=sm2vec(PCGII_cv$Est)
    # test statistics of PCGII
    tscore_PCGII_theo=sm2vec(PCGII_theo$tscore)
    tscore_PCGII_cv=sm2vec(PCGII_cv$tscore)
    print("PCGII done")
    
    ### Shrunk Partial Cor
    GGM=pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
    lambda=attr(GGM, "lambda") 
    while (lambda == 1){lambda=0.99999}
    shrunk_p=sm2vec(GGM) # off diagonal elements of estimated shrunk partial corr matrix
    # P values by Empirical null fitting (ENF)
    ENF.test=network.test.edges(GGM, fdr=TRUE, plot=FALSE,verbose=FALSE) 
    ENF.test=ENF.test[order(ENF.test$node1,ENF.test$node2),] # node1 is col; WATCH OUT:  the test must be in node's order. Sort by node 2 and by node 1
    # P values by Shrunk MLE
    pval.shrunk=p.shrunk(shrunk_p,p,n,lambda)
    print("Shrunk MLE done")
    
    
    results=cbind.data.frame(truth, Est_cLevel_theo, Est_cLevel_cv, 
                                    Est_PCGII_theo, Est_PCGII_cv, 
                                    Est_F_GGM=F_GGM.test$pcor,
                             shrunk_p, # Shrunk_pcor
                             tscore_cLevel_theo, tscore_cLevel_cv,
                             tscore_PCGII_theo, tscore_PCGII_cv,
                             F_GGM_pval=F_GGM.test$pval,
                             F_GGM_qval=F_GGM.test$qval,
                             ENF_p=ENF.test$pval, ENF_q=p.adjust(ENF.test$pval, method="BH"),
                             ShrunkMLE_p=pval.shrunk, ShrunkMLE_q=p.adjust(pval.shrunk, method="BH"))
    
    
    ## Ranking
    F_GGM.ordered<-results[order(results$F_GGM_qval, results$F_GGM_pval, decreasing = F),c("truth","F_GGM_qval","F_GGM_pval")]
    cLevel_theo.ordered<-results[order(abs(results$tscore_cLevel_theo), decreasing = T),c("truth","Est_cLevel_theo","tscore_cLevel_theo")]
    cLevel_cv.ordered<-results[order(abs(results$tscore_cLevel_cv), decreasing = T),c("truth","Est_cLevel_cv","tscore_cLevel_cv")]
    PCGII_theo.ordered<-results[order(abs(results$tscore_PCGII_theo), decreasing = T),c("truth","Est_PCGII_theo","tscore_PCGII_theo")]
    PCGII_cv.ordered<-results[order(abs(results$tscore_PCGII_cv), decreasing = T),c("truth","Est_PCGII_cv","tscore_PCGII_cv")]
    ENF.ordered<-results[order(results$ENF_q, results$ENF_p, decreasing = F),c("truth","ENF_q","ENF_p")]
    shrunk.MLE.ordered<-results[order(results$ShrunkMLE_q, results$ShrunkMLE_p, decreasing = F),c("truth","ShrunkMLE_q","ShrunkMLE_p")]
    
    
    ##### ROC
    print("ROC starts")
    for (loop in seq(1,p*(p-1)/2,1)){
      ENF.ordered.TPR[loop, k]=sum(ENF.ordered[1:loop,]$truth!=0)/CP 
      ENF.ordered.FPR[loop, k]=sum(ENF.ordered[1:loop,]$truth==0)/CN 
      ENF.ordered.PPV[loop, k]=sum(ENF.ordered[1:loop,]$truth!=0)/loop 
      ENF.ordered.FDR[loop, k]=sum(ENF.ordered[1:loop,]$truth==0)/loop   
      
      MLE.ordered.TPR[loop, k]=sum(shrunk.MLE.ordered[1:loop,]$truth!=0)/CP 
      MLE.ordered.FPR[loop, k]=sum(shrunk.MLE.ordered[1:loop,]$truth==0)/CN 
      MLE.ordered.PPV[loop, k]=sum(shrunk.MLE.ordered[1:loop,]$truth!=0)/loop
      MLE.ordered.FDR[loop, k]=sum(shrunk.MLE.ordered[1:loop,]$truth==0)/loop
      
      cLevel_theo.ordered.TPR[loop, k]=sum(cLevel_theo.ordered[1:loop,]$truth!=0)/CP 
      cLevel_theo.ordered.FPR[loop, k]=sum(cLevel_theo.ordered[1:loop,]$truth==0)/CN 
      cLevel_theo.ordered.PPV[loop, k]=sum(cLevel_theo.ordered[1:loop,]$truth!=0)/loop
      cLevel_theo.ordered.FDR[loop, k]=sum(cLevel_theo.ordered[1:loop,]$truth==0)/loop
      
      cLevel_cv.ordered.TPR[loop, k]=sum(cLevel_cv.ordered[1:loop,]$truth!=0)/CP 
      cLevel_cv.ordered.FPR[loop, k]=sum(cLevel_cv.ordered[1:loop,]$truth==0)/CN 
      cLevel_cv.ordered.PPV[loop, k]=sum(cLevel_cv.ordered[1:loop,]$truth!=0)/loop
      cLevel_cv.ordered.FDR[loop, k]=sum(cLevel_cv.ordered[1:loop,]$truth==0)/loop
      
      PCGII_theo.ordered.TPR[loop, k]=sum(PCGII_theo.ordered[1:loop,]$truth!=0)/CP 
      PCGII_theo.ordered.FPR[loop, k]=sum(PCGII_theo.ordered[1:loop,]$truth==0)/CN 
      PCGII_theo.ordered.PPV[loop, k]=sum(PCGII_theo.ordered[1:loop,]$truth!=0)/loop
      PCGII_theo.ordered.FDR[loop, k]=sum(PCGII_theo.ordered[1:loop,]$truth==0)/loop
      
      PCGII_cv.ordered.TPR[loop, k]=sum(PCGII_cv.ordered[1:loop,]$truth!=0)/CP 
      PCGII_cv.ordered.FPR[loop, k]=sum(PCGII_cv.ordered[1:loop,]$truth==0)/CN 
      PCGII_cv.ordered.PPV[loop, k]=sum(PCGII_cv.ordered[1:loop,]$truth!=0)/loop
      PCGII_cv.ordered.FDR[loop, k]=sum(PCGII_cv.ordered[1:loop,]$truth==0)/loop
      
      F_GGM.ordered.TPR[loop, k]=sum(F_GGM.ordered[1:loop,]$truth!=0)/CP 
      F_GGM.ordered.FPR[loop, k]=sum(F_GGM.ordered[1:loop,]$truth==0)/CN 
      F_GGM.ordered.PPV[loop, k]=sum(F_GGM.ordered[1:loop,]$truth!=0)/loop
      F_GGM.ordered.FDR[loop, k]=sum(F_GGM.ordered[1:loop,]$truth==0)/loop
    }
    
    for (c in 1:1000){
      tempstore.ENF[c,k]=max(ENF.ordered.TPR[ENF.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.MLE[c,k]=max(MLE.ordered.TPR[MLE.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.cLevel_theo[c,k]=max(cLevel_theo.ordered.TPR[cLevel_theo.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.cLevel_cv[c,k]=max(cLevel_cv.ordered.TPR[cLevel_cv.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_theo[c,k]=max(PCGII_theo.ordered.TPR[PCGII_theo.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_cv[c,k]=max(PCGII_cv.ordered.TPR[PCGII_cv.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.F_GGM[c,k]=max(F_GGM.ordered.TPR[F_GGM.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
    }
    print("ROC done")
    
    #### FDR 
    print("FDR starts")
    for (a in 1:length(al)){
      temp=inference(cLevel_theo, alpha = al[a])$sigs # c=0
      sigs_cLevel_theo=sigs2vec(temp, p) # significant edges
      
      temp=inference(cLevel_cv, alpha = al[a])$sigs # c=0
      sigs_cLevel_cv=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_theo, alpha = al[a])$sigs # c=0
      sigs_PCGII_theo=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_cv, alpha = al[a])$sigs # c=0
      sigs_PCGII_cv=sigs2vec(temp, p) # significant edges
      
      # Number of total selected edges
      all.enf.p[a,k]=sum(results$ENF_p<=al[a])
      all.enf.q[a,k]=sum(results$ENF_q<=al[a]) 
      all.mle.p[a,k]=sum(results$ShrunkMLE_p<=al[a]) 
      all.mle.q[a,k]=sum(results$ShrunkMLE_q<=al[a])
      all.cLevel_theo[a,k]=sum(sigs_cLevel_theo==1)
      all.cLevel_cv[a,k]=sum(sigs_cLevel_cv==1)
      all.PCGII_theo[a,k]=sum(sigs_PCGII_theo==1)
      all.PCGII_cv[a,k]=sum(sigs_PCGII_cv==1)
      all.F_GGM[a,k]=sum(results$F_GGM_qval<=al[a])
      
      # true positives
      tp.enf.p[a,k]=sum(which(results$ENF_p<=al[a]) %in% TP)
      tp.enf.q[a,k]=sum(which(results$ENF_q<=al[a]) %in% TP)
      tp.mle.p[a,k]=sum(which(results$ShrunkMLE_p<=al[a]) %in% TP)
      tp.mle.q[a,k]=sum(which(results$ShrunkMLE_q<=al[a]) %in% TP)
      tp.cLevel_theo[a,k]=sum(which(sigs_cLevel_theo==1) %in% TP)
      tp.cLevel_cv[a,k]=sum(which(sigs_cLevel_cv==1) %in% TP)
      tp.PCGII_theo[a,k]=sum(which(sigs_PCGII_theo==1) %in% TP)
      tp.PCGII_cv[a,k]=sum(which(sigs_PCGII_cv==1) %in% TP)
      tp.F_GGM[a,k]=sum(which(results$F_GGM_qval<=al[a]) %in% TP)
      
      # false positives
      fp.enf.p[a,k]=sum(!which(results$ENF_p<=al[a]) %in% TP)
      fp.enf.q[a,k]=sum(!which(results$ENF_q<=al[a]) %in% TP)
      fp.mle.p[a,k]=sum(!which(results$ShrunkMLE_p<=al[a]) %in% TP)
      fp.mle.q[a,k]=sum(!which(results$ShrunkMLE_q<=al[a]) %in% TP)
      fp.cLevel_theo[a,k]=sum(!which(sigs_cLevel_theo==1) %in% TP)
      fp.cLevel_cv[a,k]=sum(!which(sigs_cLevel_cv==1) %in% TP)
      fp.PCGII_theo[a,k]=sum(!which(sigs_PCGII_theo==1) %in% TP)
      fp.PCGII_cv[a,k]=sum(!which(sigs_PCGII_cv==1) %in% TP)
      fp.F_GGM[a,k]=sum(!which(results$F_GGM_qval<=al[a]) %in% TP)
      
      # true negatives
      tn.enf.p[a,k]=sum(!which(results$ENF_p>al[a]) %in% TP)
      tn.enf.q[a,k]=sum(!which(results$ENF.qval>al[a]) %in% TP)
      tn.mle.p[a,k]=sum(!which(results$ShrunkMLE_p>al[a]) %in% TP)
      tn.mle.q[a,k]=sum(!which(results$ShrunkMLE_q>al[a]) %in% TP)
      tn.cLevel_theo[a,k]=sum(!which(sigs_cLevel_theo!=1) %in% TP)
      tn.cLevel_cv[a,k]=sum(!which(sigs_cLevel_cv!=1) %in% TP)
      tn.PCGII_theo[a,k]=sum(!which(sigs_PCGII_theo!=1) %in% TP)
      tn.PCGII_cv[a,k]=sum(!which(sigs_PCGII_cv!=1) %in% TP)
      tn.F_GGM[a,k]=sum(!which(results$F_GGM_qval>al[a]) %in% TP)
      
      # false negatives
      fn.enf.p[a,k]=sum(which(results$ENF_p>al[a]) %in% TP)
      fn.enf.q[a,k]=sum(which(results$ENF_q>al[a]) %in% TP)
      fn.mle.p[a,k]=sum(which(results$ShrunkMLE_p>al[a]) %in% TP)
      fn.mle.q[a,k]=sum(which(results$ShrunkMLE_q>al[a]) %in% TP)
      fn.cLevel_theo[a,k]=sum(which(sigs_cLevel_theo!=1) %in% TP)
      fn.cLevel_cv[a,k]=sum(which(sigs_cLevel_cv!=1) %in% TP)
      fn.PCGII_theo[a,k]=sum(which(sigs_PCGII_theo!=1) %in% TP)
      fn.PCGII_cv[a,k]=sum(which(sigs_PCGII_cv!=1) %in% TP)
      fn.F_GGM[a,k]=sum(which(results$F_GGM_qval>al[a]) %in% TP)
    }
    print("FDR done")
    
  }
  print("Summary starts")
  {
    rank.ENF=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                              TPR=rowMeans(ENF.ordered.TPR),
                              FPR=rowMeans(ENF.ordered.FPR),
                              PPV=rowMeans(ENF.ordered.PPV),
                              FDR=rowMeans(ENF.ordered.FDR))
    rank.MLE=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                              TPR=rowMeans(MLE.ordered.TPR),
                              FPR=rowMeans(MLE.ordered.FPR),
                              PPV=rowMeans(MLE.ordered.PPV),
                              FDR=rowMeans(MLE.ordered.FDR))
    rank.cLevel_theo=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                      TPR=rowMeans(cLevel_theo.ordered.TPR),
                                      FPR=rowMeans(cLevel_theo.ordered.FPR),
                                      PPV=rowMeans(cLevel_theo.ordered.PPV),
                                      FDR=rowMeans(cLevel_theo.ordered.FDR))
    rank.cLevel_cv=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                    TPR=rowMeans(cLevel_cv.ordered.TPR),
                                    FPR=rowMeans(cLevel_cv.ordered.FPR),
                                    PPV=rowMeans(cLevel_cv.ordered.PPV),
                                    FDR=rowMeans(cLevel_cv.ordered.FDR))
    rank.PCGII_theo=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                      TPR=rowMeans(PCGII_theo.ordered.TPR),
                                      FPR=rowMeans(PCGII_theo.ordered.FPR),
                                      PPV=rowMeans(PCGII_theo.ordered.PPV),
                                      FDR=rowMeans(PCGII_theo.ordered.FDR))
    rank.PCGII_cv=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                    TPR=rowMeans(PCGII_cv.ordered.TPR),
                                    FPR=rowMeans(PCGII_cv.ordered.FPR),
                                    PPV=rowMeans(PCGII_cv.ordered.PPV),
                                    FDR=rowMeans(PCGII_cv.ordered.FDR))
    rank.F_GGM=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                TPR=rowMeans(F_GGM.ordered.TPR),
                                FPR=rowMeans(F_GGM.ordered.FPR),
                                PPV=rowMeans(F_GGM.ordered.PPV),
                                FDR=rowMeans(F_GGM.ordered.FDR))
    
    rank.table=cbind.data.frame(methods=c(rep("ENF",p*(p-1)/2),
                                          rep("MLE",p*(p-1)/2),
                                          rep("cLevel_theo",p*(p-1)/2),
                                          rep("cLevel_cv",p*(p-1)/2),
                                          rep("PCGII_theo",p*(p-1)/2),
                                          rep("PCGII_cv",p*(p-1)/2),
                                          rep("F_GGM",p*(p-1)/2)),
                                rbind.data.frame(rank.ENF, rank.MLE, rank.cLevel_theo, rank.cLevel_cv, rank.PCGII_theo, rank.PCGII_cv, rank.F_GGM))
    
    
    ROC.out=cbind.data.frame(FPR=seq(0.001,1,0.001), # Sample FPR cut points
                             ENF=rowMeans(tempstore.ENF), # averaged max TPR at cutted FPR level
                             MLE=rowMeans(tempstore.MLE),
                             cLevel_theo=rowMeans(tempstore.cLevel_theo),
                             cLevel_cv=rowMeans(tempstore.cLevel_cv),
                             PCGII_theo=rowMeans(tempstore.PCGII_theo),
                             PCGII_cv=rowMeans(tempstore.PCGII_cv),
                             F_GGM=rowMeans(tempstore.F_GGM))
    #ROC.table=tidyr::gather(ROC.out, methods, TPR, ENF:cPCG_df)
    
    
    # Number of total correctly selected edges
    All.enf.p=rowMeans(all.enf.p)
    All.enf.q=rowMeans(all.enf.q)
    All.mle.p=rowMeans(all.mle.p)
    All.mle.q=rowMeans(all.mle.q)
    All.cLevel.theo=rowMeans(all.cLevel_theo)
    All.cLevel.cv=rowMeans(all.cLevel_cv)
    All.PCGII.theo=rowMeans(all.PCGII_theo)
    All.PCGII.cv=rowMeans(all.PCGII_cv)
    All.f_ggm=rowMeans(all.F_GGM)
    
    # true positives
    TP.enf.p=rowMeans(tp.enf.p)
    TP.enf.q=rowMeans(tp.enf.q)
    TP.mle.p=rowMeans(tp.mle.p)
    TP.mle.q=rowMeans(tp.mle.q)
    TP.cLevel.theo=rowMeans(tp.cLevel_theo)
    TP.cLevel.cv=rowMeans(tp.cLevel_cv)
    TP.PCGII.theo=rowMeans(tp.PCGII_theo)
    TP.PCGII.cv=rowMeans(tp.PCGII_cv)
    TP.f_ggm=rowMeans(tp.F_GGM)
    
    # false positives
    FP.enf.p=rowMeans(fp.enf.p)
    FP.enf.q=rowMeans(fp.enf.q)
    FP.mle.p=rowMeans(fp.mle.p)
    FP.mle.q=rowMeans(fp.mle.q)
    FP.cLevel.theo=rowMeans(fp.cLevel_theo)
    FP.cLevel.cv=rowMeans(fp.cLevel_cv)
    FP.PCGII.theo=rowMeans(fp.PCGII_theo)
    FP.PCGII.cv=rowMeans(fp.PCGII_cv)
    FP.f_ggm=rowMeans(fp.F_GGM)
    
    # true negatives
    TN.enf.p=rowMeans(tn.enf.p)
    TN.enf.q=rowMeans(tn.enf.q)
    TN.mle.p=rowMeans(tn.mle.p)
    TN.mle.q=rowMeans(tn.mle.q)
    TN.cLevel.theo=rowMeans(tn.cLevel_theo)
    TN.cLevel.cv=rowMeans(tn.cLevel_cv)
    TN.PCGII.theo=rowMeans(tn.PCGII_theo)
    TN.PCGII.cv=rowMeans(tn.PCGII_cv)
    TN.f_ggm=rowMeans(tn.F_GGM)
    
    # false negatives
    FN.enf.p=rowMeans(fn.enf.p)
    FN.enf.q=rowMeans(fn.enf.q)
    FN.mle.p=rowMeans(fn.mle.p)
    FN.mle.q=rowMeans(fn.mle.q)
    FN.cLevel.theo=rowMeans(fn.cLevel_theo)
    FN.cLevel.cv=rowMeans(fn.cLevel_cv)
    FN.PCGII.theo=rowMeans(fn.PCGII_theo)
    FN.PCGII.cv=rowMeans(fn.PCGII_cv)
    FN.f_ggm=rowMeans(fn.F_GGM)
    
    # empirical fdr
    fdr.enf.p=FP.enf.p/All.enf.p
    fdr.enf.p[is.na(fdr.enf.p)]=0
    fdr.enf.q=FP.enf.q/All.enf.q
    fdr.enf.q[is.na(fdr.enf.q)]=0
    fdr.mle.p=FP.mle.p/All.mle.p
    fdr.mle.p[is.na(fdr.mle.p)]=0
    fdr.mle.q=FP.mle.q/All.mle.q
    fdr.mle.q[is.na(fdr.mle.q)]=0
    fdr.cLevel.theo=FP.cLevel.theo/All.cLevel.theo
    fdr.cLevel.theo[is.na(fdr.cLevel.theo)]=0
    fdr.cLevel.cv=FP.cLevel.cv/All.cLevel.cv
    fdr.cLevel.cv[is.na(fdr.cLevel.cv)]=0
    fdr.PCGII.theo=FP.PCGII.theo/All.PCGII.theo
    fdr.PCGII.theo[is.na(fdr.PCGII.theo)]=0
    fdr.PCGII.cv=FP.PCGII.cv/All.PCGII.cv
    fdr.PCGII.cv[is.na(fdr.PCGII.cv)]=0
    fdr.f_ggm=FP.f_ggm/All.f_ggm
    fdr.f_ggm[is.na(fdr.f_ggm)]=0
    
    
    fdr.table=cbind.data.frame(
      FDR=rep(al,9), # nominal FDR
      methods=c(rep("ENF.p",length(al)),
                rep("ENF.q",length(al)),
                rep("MLE.p",length(al)),
                rep("MLE.q",length(al)),
                rep("cLevel.theo",length(al)),
                rep("cLevel.cv",length(al)),
                rep("PCGII.theo",length(al)),
                rep("PCGII.cv",length(al)),
                rep("F_GGM",length(al))),
      fdr=c(fdr.enf.p, fdr.enf.q, fdr.mle.p, fdr.mle.q, fdr.cLevel.theo, fdr.cLevel.cv, fdr.PCGII.theo, fdr.PCGII.cv, fdr.f_ggm))
  } # summary
  print("Summary done")
  {
    return(
      tepp=list(
        ROC=list(
          # ROC
          ENF.ordered.TPR=ENF.ordered.TPR, ENF.ordered.FPR=ENF.ordered.FPR, 
          ENF.ordered.PPV=ENF.ordered.PPV, ENF.ordered.FDR=ENF.ordered.FDR,
          MLE.ordered.TPR=MLE.ordered.TPR, MLE.ordered.FPR=MLE.ordered.FPR, 
          MLE.ordered.PPV=MLE.ordered.PPV, MLE.ordered.FDR=MLE.ordered.FDR,
          cLevel_theo.ordered.TPR=cLevel_theo.ordered.TPR, cLevel_theo.ordered.FPR=cLevel_theo.ordered.FPR, 
          cLevel_theo.ordered.PPV=cLevel_theo.ordered.PPV, cLevel_theo.ordered.FDR=cLevel_theo.ordered.FDR,
          cLevel_cv.ordered.TPR=cLevel_cv.ordered.TPR, cLevel_cv.ordered.FPR=cLevel_cv.ordered.FPR, 
          cLevel_cv.ordered.PPV=cLevel_cv.ordered.PPV, cLevel_cv.ordered.FDR=cLevel_cv.ordered.FDR,
          PCGII_theo.ordered.TPR=PCGII_theo.ordered.TPR, PCGII_theo.ordered.FPR=PCGII_theo.ordered.FPR, 
          PCGII_theo.ordered.PPV=PCGII_theo.ordered.PPV, PCGII_theo.ordered.FDR=PCGII_theo.ordered.FDR,
          PCGII_cv.ordered.TPR=PCGII_cv.ordered.TPR, PCGII_cv.ordered.FPR=PCGII_cv.ordered.FPR, 
          PCGII_cv.ordered.PPV=PCGII_cv.ordered.PPV, PCGII_cv.ordered.FDR=PCGII_cv.ordered.FDR,
          F_GGM.ordered.TPR=F_GGM.ordered.TPR, F_GGM.ordered.FPR=F_GGM.ordered.FPR, 
          F_GGM.ordered.PPV=F_GGM.ordered.PPV, F_GGM.ordered.FDR=F_GGM.ordered.FDR,
          rank.table=rank.table,
          
          # empirical TPR
          tempstore.ENF=tempstore.ENF, tempstore.ENF=tempstore.ENF, tempstore.MLE=tempstore.MLE, 
          tempstore.cLevel_theo=tempstore.cLevel_theo, tempstore.cLevel_cv=tempstore.cLevel_cv, 
          tempstore.PCGII_theo=tempstore.PCGII_theo, tempstore.PCGII_cv=tempstore.PCGII_cv, 
          tempstore.F_GGM=tempstore.F_GGM,
          ROC.out=ROC.out),
        
        # FDR
        metrics=list( 
          # Number of total correctly selected edges
          all.enf.p=all.enf.p, all.enf.q=all.enf.q,  
          all.mle.p=all.mle.p,  all.mle.q=all.mle.q,
          all.cLevel_theo=all.cLevel_theo, all.cLevel_cv=all.cLevel_cv, 
          all.PCGII_theo=all.PCGII_theo, all.PCGII_cv=all.PCGII_cv, 
          all.F_GGM=all.F_GGM,
          
          # true positives
          tp.enf.p=tp.enf.p, tp.enf.q=tp.enf.q, 
          tp.mle.p=tp.mle.p, tp.mle.q=tp.mle.q,
          tp.cLevel_theo=tp.cLevel_theo, tp.cLevel_cv=tp.cLevel_cv, 
          tp.PCGII_theo=tp.PCGII_theo, tp.PCGII_cv=tp.PCGII_cv, 
          tp.F_GGM=tp.F_GGM,
          
          # false positives
          fp.enf.p=fp.enf.p, fp.enf.q=fp.enf.q, 
          fp.mle.p=fp.mle.p, fp.mle.q=fp.mle.q,
          fp.cLevel_theo=fp.cLevel_theo, fp.cLevel_cv=fp.cLevel_cv, 
          fp.PCGII_theo=fp.PCGII_theo, fp.PCGII_cv=fp.PCGII_cv, 
          fp.F_GGM=fp.F_GGM,
          
          # false negatives
          tn.enf.p=tn.enf.p, tn.enf.q=tn.enf.q, 
          tn.mle.p=tn.mle.p, tn.mle.q=tn.mle.q,
          tn.cLevel_theo=tn.cLevel_theo, tn.cLevel_cv=tn.cLevel_cv, 
          tn.PCGII_theo=tn.PCGII_theo, tn.PCGII_cv=tn.PCGII_cv, 
          tn.F_GGM=tn.F_GGM,
          
          # false negatives
          fn.enf.p=fn.enf.p, fn.enf.q=fn.enf.q,
          fn.mle.p=fn.mle.p, fn.mle.q=fn.mle.q,
          fn.cLevel_theo=fn.cLevel_theo, fn.cLevel_cv=fn.cLevel_cv, 
          fn.PCGII_theo=fn.PCGII_theo, fn.PCGII_cv=fn.PCGII_cv, 
          fn.F_GGM=fn.F_GGM,
          
          All.enf.p=All.enf.p, All.enf.q=All.enf.q, 
          All.mle.p=All.mle.p, All.mle.q=All.mle.q,
          All.cLevel.theo=All.cLevel.theo, All.cLevel.cv=All.cLevel.cv, 
          All.PCGII.theo=All.PCGII.theo, All.PCGII.cv=All.PCGII.cv, 
          All.f_ggm=All.f_ggm,
          
          # true positives
          TP.enf.p=TP.enf.p,TP.enf.q=TP.enf.q, 
          TP.mle.p=TP.mle.p,TP.mle.q=TP.mle.q, 
          TP.cLevel.theo=TP.cLevel.theo,TP.cLevel.cv=TP.cLevel.cv,
          TP.PCGII.theo=TP.PCGII.theo,TP.PCGII.cv=TP.PCGII.cv,
          TP.f_ggm=TP.f_ggm,
          
          # false positives
          FP.enf.p=FP.enf.p,FP.enf.q=FP.enf.q, 
          FP.mle.p=FP.mle.p, FP.mle.q=FP.mle.q,
          FP.cLevel.theo=FP.cLevel.theo, FP.cLevel.cv=FP.cLevel.cv,
          FP.PCGII.theo=FP.PCGII.theo, FP.PCGII.cv=FP.PCGII.cv,
          FP.f_ggm=FP.f_ggm,
          
          # true negatives
          TN.enf.p=TN.enf.p, TN.enf.q=TN.enf.q, 
          TN.mle.p=TN.mle.p, TN.mle.q=TN.mle.q,
          TN.cLevel.theo=TN.cLevel.theo, TN.cLevel.cv=TN.cLevel.cv, 
          TN.PCGII.theo=TN.PCGII.theo, TN.PCGII.cv=TN.PCGII.cv, 
          TN.f_ggm=TN.f_ggm,
          
          # false negatives
          FN.enf.p=FN.enf.p, FN.enf.q=FN.enf.q, 
          FN.mle.p=FN.mle.p,FN.mle.q=FN.mle.q,
          FN.cLevel.theo=FN.cLevel.theo,FN.cLevel.cv=FN.cLevel.cv, 
          FN.PCGII.theo=FN.PCGII.theo,FN.PCGII.cv=FN.PCGII.cv, 
          FN.f_ggm=FN.f_ggm
        ),
        # empirical fdr
        fdr=list(
          fdr.enf.p=fdr.enf.p, fdr.enf.q=fdr.enf.q, 
          fdr.mle.p=fdr.mle.p, fdr.mle.q=fdr.mle.q,
          fdr.cLevel.theo=fdr.cLevel.theo, fdr.cLevel.cv=fdr.cLevel.cv, 
          fdr.PCGII.theo=fdr.PCGII.theo, fdr.PCGII.cv=fdr.PCGII.cv, 
          fdr.f_ggm=fdr.f_ggm,
      
          fdr.table=fdr.table
        )
      )
    )
  } # return
}
