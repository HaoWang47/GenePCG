library(tidyverse)
library(GeneNet)
library(FastGGM)
library(corpcor)
library(glmnet)

source("R/Qiu.R")
source("R/shrinkagefunctions.R") # https://github.com/V-Bernal/GGM-Shrinkage


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

model.eval=function(X, pcor, rep=5, degfree=5){
  require(GeneNet)
  
  p=dim(pcor)[1]
  truth<-sm2vec(pcor) 
  TP<- which(truth!=0) 
  CN=p*(p-1)/2-length(TP) 
  CP=length(TP)
  
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
    
    F_GGM.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    F_GGM.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    F_GGM.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    F_GGM.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    tempstore.ENF=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.MLE=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points 
    tempstore.cLevel_theo=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.cLevel_cv=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.F_GGM=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    
    ### Inference
    
    # Number of total correctly selected edges
    all.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    all.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.mle.p=matrix(Inf, nrow=length(al), ncol=rep) 
    all.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    all.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    all.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # true positives
    tp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    tp.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    tp.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false positives
    fp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    fp.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    fp.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false negatives
    tn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    tn.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    tn.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false negatives
    fn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    fn.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    fn.F_GGM=matrix(Inf, nrow=length(al), ncol=rep)
  }
  
  
  for (k in 1:rep){
    print(paste0(k,"th rep starts!!!!!"))
    sim.data=X[[k]]
    
    ### Fast GGM
    out=FastGGM(sim.data)
    F_GGM.test=unMat(X_est=out$partialCor, X_p=out$p_partialCor) # Est, pvals
    F_GGM.test=cbind.data.frame(F_GGM.test,truth=truth, qval=p.adjust(F_GGM.test[,4], method="BH"))
    
    
    ### cPCG
    cPCG_theo=estimating(sim.data)
    cPCG_cv=estimating_cv(sim.data,degfree)
    # estimates by cPCG
    Est_theo=sm2vec(cPCG_theo$Est)
    Est_cv=sm2vec(cPCG_cv$Est)
    # test statistics of cPCG
    tscore_theo=sm2vec(cPCG_theo$tscore)
    tscore_cv=sm2vec(cPCG_cv$tscore)
    
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
    
    
    results=cbind.data.frame(truth, Est_theo, Est_cv, Est_F_GGM=F_GGM.test$pcor,
                             shrunk_p, # Shrunk_pcor
                             tscore_theo, tscore_cv,
                             F_GGM_pval=F_GGM.test$pval,
                             F_GGM_qval=F_GGM.test$qval,
                             ENF_p=ENF.test$pval, ENF_q=p.adjust(ENF.test$pval, method="BH"),
                             ShrunkMLE_p=pval.shrunk, ShrunkMLE_q=p.adjust(pval.shrunk, method="BH"))
    
    
    
    
    ## Ranking
    F_GGM.ordered<-results[order(results$F_GGM_qval, results$F_GGM_pval, decreasing = F),c("truth","F_GGM_qval","F_GGM_pval")]
    cLevel_theo.ordered<-results[order(abs(results$tscore_theo), decreasing = T),c("truth","Est_theo","tscore_theo")]
    cLevel_cv.ordered<-results[order(abs(results$tscore_cv), decreasing = T),c("truth","Est_theo","tscore_theo")]
    ENF.ordered<-results[order(results$ENF_q, results$ENF_p, decreasing = F),c("truth","ENF_q","ENF_p")]
    shrunk.MLE.ordered<-results[order(results$ShrunkMLE_q, results$ShrunkMLE_p, decreasing = F),c("truth","ShrunkMLE_q","ShrunkMLE_p")]

    
    ##### ROC
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
      tempstore.F_GGM[c,k]=max(F_GGM.ordered.TPR[F_GGM.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
    }
    
    #### FDR 
    for (a in 1:length(al)){
      temp=inference(cPCG_theo, alpha = al[a])$sigs # c=0
      sigs_theo=sigs2vec(temp, p) # significant edges
      
      temp=inference(cPCG_cv, alpha = al[a])$sigs # c=0
      sigs_cv=sigs2vec(temp, p) # significant edges

      # Number of total selected edges
      all.enf.p[a,k]=sum(results$ENF_p<=al[a])
      all.enf.q[a,k]=sum(results$ENF_q<=al[a]) 
      all.mle.p[a,k]=sum(results$ShrunkMLE_p<=al[a]) 
      all.mle.q[a,k]=sum(results$ShrunkMLE_q<=al[a])
      all.clevel_theo[a,k]=sum(sigs_theo==1)
      all.clevel_cv[a,k]=sum(sigs_cv==1)
      all.F_GGM[a,k]=sum(results$F_GGM_qval<=al[a])
      
      # true positives
      tp.enf.p[a,k]=sum(which(results$ENF_p<=al[a]) %in% TP)
      tp.enf.q[a,k]=sum(which(results$ENF_q<=al[a]) %in% TP)
      tp.mle.p[a,k]=sum(which(results$ShrunkMLE_p<=al[a]) %in% TP)
      tp.mle.q[a,k]=sum(which(results$ShrunkMLE_q<=al[a]) %in% TP)
      tp.clevel_theo[a,k]=sum(which(sigs_theo==1) %in% TP)
      tp.clevel_cv[a,k]=sum(which(sigs_cv==1) %in% TP)
      tp.F_GGM[a,k]=sum(which(results$F_GGM_qval<=al[a]) %in% TP)
      
      # false positives
      fp.enf.p[a,k]=sum(!which(results$ENF_p<=al[a]) %in% TP)
      fp.enf.q[a,k]=sum(!which(results$ENF_q<=al[a]) %in% TP)
      fp.mle.p[a,k]=sum(!which(results$ShrunkMLE_p<=al[a]) %in% TP)
      fp.mle.q[a,k]=sum(!which(results$ShrunkMLE_q<=al[a]) %in% TP)
      fp.clevel_theo[a,k]=sum(!which(sigs_theo==1) %in% TP)
      fp.clevel_cv[a,k]=sum(!which(sigs_cv==1) %in% TP)
      fp.F_GGM[a,k]=sum(!which(results$F_GGM_qval<=al[a]) %in% TP)
      
      # true negatives
      tn.enf.p[a,k]=sum(!which(results$ENF_p>al[a]) %in% TP)
      tn.enf.q[a,k]=sum(!which(results$ENF.qval>al[a]) %in% TP)
      tn.mle.p[a,k]=sum(!which(results$ShrunkMLE_p>al[a]) %in% TP)
      tn.mle.q[a,k]=sum(!which(results$ShrunkMLE_q>al[a]) %in% TP)
      tn.clevel_theo[a,k]=sum(!which(sigs_theo!=1) %in% TP)
      tn.clevel_cv[a,k]=sum(!which(sigs_cv!=1) %in% TP)
      tn.F_GGM[a,k]=sum(!which(results$F_GGM_qval>al[a]) %in% TP)
      
      # false negatives
      fn.enf.p[a,k]=sum(which(results$ENF_p>al[a]) %in% TP)
      fn.enf.q[a,k]=sum(which(results$ENF_q>al[a]) %in% TP)
      fn.mle.p[a,k]=sum(which(results$ShrunkMLE_p>al[a]) %in% TP)
      fn.mle.q[a,k]=sum(which(results$ShrunkMLE_q>al[a]) %in% TP)
      fn.clevel_theo[a,k]=sum(which(sigs_theo!=1) %in% TP)
      fn.clevel_cv[a,k]=sum(which(sigs_cv!=1) %in% TP)
      fn.F_GGM[a,k]=sum(which(results$F_GGM_qval>al[a]) %in% TP)
    }
    
  }
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
    rank.F_GGM=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                    TPR=rowMeans(F_GGM.ordered.TPR),
                                    FPR=rowMeans(F_GGM.ordered.FPR),
                                    PPV=rowMeans(F_GGM.ordered.PPV),
                                    FDR=rowMeans(F_GGM.ordered.FDR))
    rank.table=cbind.data.frame(methods=c(rep("ENF",p*(p-1)/2),
                                          rep("MLE",p*(p-1)/2),
                                          rep("cPCG_theo",p*(p-1)/2),
                                          rep("cPCG_cv",p*(p-1)/2),
                                          rep("F_GGM",p*(p-1)/2)),
                                rbind.data.frame(rank.ENF, rank.MLE, rank.cLevel_theo, rank.cLevel_cv, rank.F_GGM))
    
    
    ROC.out=cbind.data.frame(FPR=seq(0.001,1,0.001), # Sample FPR cut points
                             ENF=rowMeans(tempstore.ENF), # averaged max TPR at cutted FPR level
                             MLE=rowMeans(tempstore.MLE),
                             cPCG_theo=rowMeans(tempstore.cLevel_theo),
                             cPCG_cv=rowMeans(tempstore.cLevel_cv),
                             F_GGM=rowMeans(tempstore.F_GGM))
    #ROC.table=tidyr::gather(ROC.out, methods, TPR, ENF:cPCG_df)
    
    
    # Number of total correctly selected edges
    All.enf.p=rowMeans(all.enf.p)
    All.enf.q=rowMeans(all.enf.q)
    All.mle.p=rowMeans(all.mle.p)
    All.mle.q=rowMeans(all.mle.q)
    All.clevel.theo=rowMeans(all.clevel_theo)
    All.clevel.cv=rowMeans(all.clevel_cv)
    All.f_ggm=rowMeans(all.F_GGM)
    
    # true positives
    TP.enf.p=rowMeans(tp.enf.p)
    TP.enf.q=rowMeans(tp.enf.q)
    TP.mle.p=rowMeans(tp.mle.p)
    TP.mle.q=rowMeans(tp.mle.q)
    TP.clevel.theo=rowMeans(tp.clevel_theo)
    TP.clevel.cv=rowMeans(tp.clevel_cv)
    TP.f_ggm=rowMeans(tp.F_GGM)
    
    # false positives
    FP.enf.p=rowMeans(fp.enf.p)
    FP.enf.q=rowMeans(fp.enf.q)
    FP.mle.p=rowMeans(fp.mle.p)
    FP.mle.q=rowMeans(fp.mle.q)
    FP.clevel.theo=rowMeans(fp.clevel_theo)
    FP.clevel.cv=rowMeans(fp.clevel_cv)
    FP.f_ggm=rowMeans(fp.F_GGM)
    
    # true negatives
    TN.enf.p=rowMeans(tn.enf.p)
    TN.enf.q=rowMeans(tn.enf.q)
    TN.mle.p=rowMeans(tn.mle.p)
    TN.mle.q=rowMeans(tn.mle.q)
    TN.clevel.theo=rowMeans(tn.clevel_theo)
    TN.clevel.cv=rowMeans(tn.clevel_cv)
    TN.f_ggm=rowMeans(tn.F_GGM)
    
    # false negatives
    FN.enf.p=rowMeans(fn.enf.p)
    FN.enf.q=rowMeans(fn.enf.q)
    FN.mle.p=rowMeans(fn.mle.p)
    FN.mle.q=rowMeans(fn.mle.q)
    FN.clevel.theo=rowMeans(fn.clevel_theo)
    FN.clevel.cv=rowMeans(fn.clevel_cv)
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
    fdr.clevel.theo=FP.clevel.theo/All.clevel.theo
    fdr.clevel.theo[is.na(fdr.clevel.theo)]=0
    fdr.clevel.cv=FP.clevel.cv/All.clevel.cv
    fdr.clevel.cv[is.na(fdr.clevel.cv)]=0
    fdr.f_ggm=FP.f_ggm/All.f_ggm
    fdr.f_ggm[is.na(fdr.f_ggm)]=0
    
    
    fdr.table=cbind.data.frame(
      FDR=rep(al,7), # nominal FDR
      methods=c(rep("ENF.p",length(al)),
                rep("ENF.q",length(al)),
                rep("MLE.p",length(al)),
                rep("MLE.q",length(al)),
                rep("cPCG.theo",length(al)),
                rep("cPCG.cv",length(al)),
                rep("F_GGM",length(al))),
      fdr=c(fdr.enf.p, fdr.enf.q, fdr.mle.p, fdr.mle.q, fdr.clevel.theo, fdr.clevel.cv, fdr.f_ggm))
  } # summary
  
  {
    return(
      list(
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
          F_GGM.ordered.TPR=F_GGM.ordered.TPR, F_GGM.ordered.FPR=F_GGM.ordered.FPR, 
          F_GGM.ordered.PPV=F_GGM.ordered.PPV, F_GGM.ordered.FDR=F_GGM.ordered.FDR,
          rank.table=rank.table,
          
          # empirical TPR
          tempstore.ENF=tempstore.ENF, tempstore.ENF=tempstore.ENF, tempstore.MLE=tempstore.MLE, 
          tempstore.cLevel_theo=tempstore.cLevel_theo, tempstore.cLevel_cv=tempstore.cLevel_cv, 
          tempstore.F_GGM=tempstore.F_GGM,
          ROC.out=ROC.out),
        
        # FDR
        metrics=list( 
          # Number of total correctly selected edges
          all.enf.p=all.enf.p, all.enf.q=all.enf.q,  
          all.mle.p=all.mle.p,  all.mle.q=all.mle.q,
          all.clevel_theo=all.clevel_theo, all.clevel_cv=all.clevel_cv, 
          all.F_GGM=all.F_GGM,
          
          # true positives
          tp.enf.p=tp.enf.p, tp.enf.q=tp.enf.q, 
          tp.mle.p=tp.mle.p, tp.mle.q=tp.mle.q,
          tp.clevel_theo=tp.clevel_theo, tp.clevel_cv=tp.clevel_cv, 
          tp.F_GGM=tp.F_GGM,
          
          # false positives
          fp.enf.p=fp.enf.p, fp.enf.q=fp.enf.q, 
          fp.mle.p=fp.mle.p, fp.mle.q=fp.mle.q,
          fp.clevel_theo=fp.clevel_theo, fp.clevel_cv=fp.clevel_cv, 
          fp.F_GGM=fp.F_GGM,
          
          # false negatives
          tn.enf.p=tn.enf.p, tn.enf.q=tn.enf.q, 
          tn.mle.p=tn.mle.p, tn.mle.q=tn.mle.q,
          tn.clevel_theo=tn.clevel_theo, tn.clevel_cv=tn.clevel_cv, 
          tn.F_GGM=tn.F_GGM,
          
          # false negatives
          fn.enf.p=fn.enf.p, fn.enf.q=fn.enf.q,
          fn.mle.p=fn.mle.p, fn.mle.q=fn.mle.q,
          fn.clevel_theo=fn.clevel_theo, fn.clevel_cv=fn.clevel_cv, 
          fn.F_GGM=fn.F_GGM,
          
          All.enf.p=All.enf.p, All.enf.q=All.enf.q, 
          All.mle.p=All.mle.p, All.mle.q=All.mle.q,
          All.clevel.theo=All.clevel.theo, All.clevel.cv=All.clevel.cv, 
          All.f_ggm=All.f_ggm,
          
          # true positives
          TP.enf.p=TP.enf.p,TP.enf.q=TP.enf.q, 
          TP.mle.p=TP.mle.p,TP.mle.q=TP.mle.q, 
          TP.clevel.theo=TP.clevel.theo,TP.clevel.cv=TP.clevel.cv,
          TP.f_ggm=TP.f_ggm,
          
          # false positives
          FP.enf.p=FP.enf.p,FP.enf.q=FP.enf.q, 
          FP.mle.p=FP.mle.p, FP.mle.q=FP.mle.q,
          FP.clevel.theo=FP.clevel.theo, FP.clevel.cv=FP.clevel.cv,
          FP.f_ggm=FP.f_ggm,
          
          # true negatives
          TN.enf.p=TN.enf.p, TN.enf.q=TN.enf.q, 
          TN.mle.p=TN.mle.p, TN.mle.q=TN.mle.q,
          TN.clevel.theo=TN.clevel.theo, TN.clevel.cv=TN.clevel.cv, 
          TN.f_ggm=TN.f_ggm,
          
          # false negatives
          FN.enf.p=FN.enf.p, FN.enf.q=FN.enf.q, 
          FN.mle.p=FN.mle.p,FN.mle.q=FN.mle.q,
          FN.clevel.theo=FN.clevel.theo,FN.clevel.cv=FN.clevel.cv, 
          FN.f_ggm=FN.f_ggm
        ),
        # empirical fdr
        fdr=list(
          fdr.enf.p=fdr.enf.p, fdr.enf.q=fdr.enf.q, 
          fdr.mle.p=fdr.mle.p, fdr.mle.q=fdr.mle.q,
          fdr.clevel.theo=fdr.clevel.theo, fdr.clevel.cv=fdr.clevel.cv, 
          fdr.f_ggm=fdr.f_ggm,
          
          fdr.table=fdr.table
        )
      )
    )
  } # return
}


# Block Diagonal ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
#min.beta = 0.3 # v2
#max.beta = 1
#blocksize = c(4,10)
min.beta = 0.3 # v3, v4 
max.beta = 0.9
blocksize = c(4,8,10)
for (n in nl) {
  for (p in pl) {
    for (e in blocksize) {
      load( file=sprintf("simu_data/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta) )
      degf=5
      RESULT=model.eval(X, pcor, rep = 50, degfree = degf)
      save(RESULT, file=sprintf("results/BlockDiag_FGGM_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta) )
      
    }
  }
}


# Scale Free ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta = 0.3
edge = c(1, 2)
for (n in nl) {
  for (p in pl) {
    for (e in edge) {
      load( file=sprintf("simu_data/ScaleFree_simu_n%d_p%d_e%d_min_beta%g.RData", n, p, e, min.beta) )
      degf=10
      RESULT=model.eval(X, pcor, rep = 50, degfree = degf)
      save(RESULT, file=sprintf("results/ScaleFree_FGGM_n%d_p%d_e%d_min_beta%g.RData", n, p, e, min.beta) )
    }
  }
}


# Random Structure ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
etal = c(.01, .02, .03)
min.beta = 0.3
for (n in nl) {
  for (p in pl) {
    for (eta in etal) {
      load( file=sprintf("simu_data/Random_simu_n%d_p%d_eta%g_min_beta%g.RData", n, p, eta, min.beta) )
      degf=10
      RESULT=model.eval(X, pcor, rep = 50, degfree = degf)
      save(RESULT, file=sprintf("results/Random_FGGM_n%d_p%d_eta%g_min_beta%g.RData", n, p, eta, min.beta) )
      
    }
  }
}




