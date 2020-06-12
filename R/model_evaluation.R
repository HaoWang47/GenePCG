library(ggplot2)
library(tidyverse)
library(GeneNet)

source("/Users/Wanghao/Documents/ISU/RA/Fall-2019/Qiu.R")
#source("R/Qiu.R")
source("R/shrinkagefunctions.R") # https://github.com/V-Bernal/GGM-Shrinkage

sigs2vec=function(sigs, P){
  require(corpcor)
  m=matrix(0,P,P)
  for (h in 1: dim(sigs)[1]){
    m[sigs[h,1],sigs[h,2]]=1
  }
  sm2vec(m)
}


model.eval=function(X, pcor, rep=10, degfree=5){
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
    
    cLevel_df.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_df.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_df.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_df.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    tempstore.ENF=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.MLE=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points 
    tempstore.cLevel_theo=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.cLevel_cv=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.cLevel_df=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    
    ### Inference
    
    # Number of total correctly selected edges
    all.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    all.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.mle.p=matrix(Inf, nrow=length(al), ncol=rep) 
    all.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    all.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    all.clevel_df=matrix(Inf, nrow=length(al), ncol=rep)
    
    # true positives
    tp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    tp.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    tp.clevel_df=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false positives
    fp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    fp.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    fp.clevel_df=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false negatives
    tn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    tn.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    tn.clevel_df=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false negatives
    fn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.clevel_theo=matrix(Inf, nrow=length(al), ncol=rep)
    fn.clevel_cv=matrix(Inf, nrow=length(al), ncol=rep)
    fn.clevel_df=matrix(Inf, nrow=length(al), ncol=rep)
  }
  
  
  for (k in 1:rep){
    print(paste0(k,"th rep starts!!!!!"))
    sim.data=X[[k]]
    
    #cPCG_theo=cPCG_theo_lambda(sim.data)
    cPCG_theo=estimating(sim.data)
    #cPCG_cv=cPCG_cv_dfmax(sim.data)
    cPCG_cv=estimating_cv(sim.data)
    #cPCG_df=cPCG_cv_dfmax(sim.data, degree_freedom = degfree)
    cPCG_df=estimating_q(sim.data,degfree)
    
    # estimates by cPCG
    Est_theo=sm2vec(cPCG_theo$Est)
    Est_cv=sm2vec(cPCG_cv$Est)
    Est_df=sm2vec(cPCG_df$Est)
    
    # test statistics of cPCG
    tscore_theo=sm2vec(cPCG_theo$tscore)
    tscore_cv=sm2vec(cPCG_cv$tscore)
    tscore_df=sm2vec(cPCG_df$tscore)
    
    ### Shrunk Partial Cor
    GGM=pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
    lambda=attr(GGM, "lambda") 
    while (lambda == 1){lambda=0.99999}
    shrunk_p=sm2vec(GGM) # off diagonal elements of estimated shrunk partial corr matrix
    
    ### P values by Empirical null fitting (ENF)
    ENF.test=network.test.edges(GGM, fdr=TRUE, plot=FALSE,verbose=FALSE) 
    ENF.test=ENF.test[order(ENF.test$node1,ENF.test$node2),] # node1 is col; WATCH OUT:  the test must be in node's order. Sort by node 2 and by node 1
    
    ### P values by Shrunk MLE
    pval.shrunk=p.shrunk(shrunk_p,p,n,lambda)
    
    results=cbind.data.frame(truth, Est_theo, Est_cv, Est_df, 
                             shrunk_p, 
                             tscore_theo, tscore_cv, tscore_df,
                             ENF_p=ENF.test$pval, ENF_q=p.adjust(ENF.test$pval, method="BH"),
                             ShrunkMLE_p=pval.shrunk, ShrunkMLE_q=p.adjust(pval.shrunk, method="BH"))
    
    ## Ranking
    ENF.ordered<-results[order(results$ENF_q, results$ENF_p, decreasing = F),c("truth","ENF_p","ENF_q")]
    shrunk.MLE.ordered<-results[order(results$ShrunkMLE_q, results$ShrunkMLE_p, decreasing = F),c("truth","ShrunkMLE_p","ShrunkMLE_q")]
    
    cLevel_theo.ordered<-results[order(abs(results$tscore_theo), decreasing = T),]
    cLevel_cv.ordered<-results[order(abs(results$tscore_cv), decreasing = T),]
    cLevel_df.ordered<-results[order(abs(results$tscore_df), decreasing = T),]
    
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
      
      cLevel_df.ordered.TPR[loop, k]=sum(cLevel_df.ordered[1:loop,]$truth!=0)/CP 
      cLevel_df.ordered.FPR[loop, k]=sum(cLevel_df.ordered[1:loop,]$truth==0)/CN 
      cLevel_df.ordered.PPV[loop, k]=sum(cLevel_df.ordered[1:loop,]$truth!=0)/loop
      cLevel_df.ordered.FDR[loop, k]=sum(cLevel_df.ordered[1:loop,]$truth==0)/loop
    }
    
    for (c in 1:1000){
      tempstore.ENF[c,k]=max(ENF.ordered.TPR[ENF.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.MLE[c,k]=max(MLE.ordered.TPR[MLE.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.cLevel_theo[c,k]=max(cLevel_theo.ordered.TPR[cLevel_theo.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.cLevel_cv[c,k]=max(cLevel_cv.ordered.TPR[cLevel_cv.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.cLevel_df[c,k]=max(cLevel_df.ordered.TPR[cLevel_df.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
    }
    
    #### FDR 
    for (a in 1:length(al)){
      temp=inference(cPCG_theo, alpha = al[a])$sigs # c=0
      sigs_theo=sigs2vec(temp, p) # significant edges
      
      temp=inference(cPCG_cv, alpha = al[a])$sigs # c=0
      sigs_cv=sigs2vec(temp, p) # significant edges
      
      temp=inference(cPCG_df, alpha = al[a])$sigs # c=0
      sigs_df=sigs2vec(temp, p) # significant edges
      
      # Number of total correctly selected edges
      all.enf.p[a,k]=sum(results$ENF_p<=al[a])
      all.enf.q[a,k]=sum(results$ENF_q<=al[a]) 
      all.mle.p[a,k]=sum(results$ShrunkMLE_p<=al[a]) 
      all.mle.q[a,k]=sum(results$ShrunkMLE_q<=al[a])
      all.clevel_theo[a,k]=sum(sigs_theo==1)
      all.clevel_cv[a,k]=sum(sigs_cv==1)
      all.clevel_df[a,k]=sum(sigs_df==1)
      
      # true positives
      tp.enf.p[a,k]=sum(which(results$ENF_p<=al[a]) %in% TP)
      tp.enf.q[a,k]=sum(which(results$ENF_q<=al[a]) %in% TP)
      tp.mle.p[a,k]=sum(which(results$ShrunkMLE_p<=al[a]) %in% TP)
      tp.mle.q[a,k]=sum(which(results$ShrunkMLE_q<=al[a]) %in% TP)
      tp.clevel_theo[a,k]=sum(which(sigs_theo==1) %in% TP)
      tp.clevel_cv[a,k]=sum(which(sigs_cv==1) %in% TP)
      tp.clevel_df[a,k]=sum(which(sigs_df==1) %in% TP)
      
      # false positives
      fp.enf.p[a,k]=sum(!which(results$ENF_p<=al[a]) %in% TP)
      fp.enf.q[a,k]=sum(!which(results$ENF_q<=al[a]) %in% TP)
      fp.mle.p[a,k]=sum(!which(results$ShrunkMLE_p<=al[a]) %in% TP)
      fp.mle.q[a,k]=sum(!which(results$ShrunkMLE_q<=al[a]) %in% TP)
      fp.clevel_theo[a,k]=sum(!which(sigs_theo==1) %in% TP)
      fp.clevel_cv[a,k]=sum(!which(sigs_cv==1) %in% TP)
      fp.clevel_df[a,k]=sum(!which(sigs_df==1) %in% TP)
      
      # true negatives
      tn.enf.p[a,k]=sum(!which(results$ENF_p>al[a]) %in% TP)
      tn.enf.q[a,k]=sum(!which(results$ENF.qval>al[a]) %in% TP)
      tn.mle.p[a,k]=sum(!which(results$ShrunkMLE_p>al[a]) %in% TP)
      tn.mle.q[a,k]=sum(!which(results$ShrunkMLE_q>al[a]) %in% TP)
      tn.clevel_theo[a,k]=sum(!which(sigs_theo!=1) %in% TP)
      tn.clevel_cv[a,k]=sum(!which(sigs_cv!=1) %in% TP)
      tn.clevel_df[a,k]=sum(!which(sigs_df!=1) %in% TP)
      
      # false negatives
      fn.enf.p[a,k]=sum(which(results$ENF_p>al[a]) %in% TP)
      fn.enf.q[a,k]=sum(which(results$ENF_q>al[a]) %in% TP)
      fn.mle.p[a,k]=sum(which(results$ShrunkMLE_p>al[a]) %in% TP)
      fn.mle.q[a,k]=sum(which(results$ShrunkMLE_q>al[a]) %in% TP)
      fn.clevel_theo[a,k]=sum(which(sigs_theo!=1) %in% TP)
      fn.clevel_cv[a,k]=sum(which(sigs_cv!=1) %in% TP)
      fn.clevel_df[a,k]=sum(which(sigs_df!=1) %in% TP)
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
    rank.cLevel_df=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                    TPR=rowMeans(cLevel_df.ordered.TPR),
                                    FPR=rowMeans(cLevel_df.ordered.FPR),
                                    PPV=rowMeans(cLevel_df.ordered.PPV),
                                    FDR=rowMeans(cLevel_df.ordered.FDR))
    rank.table=cbind.data.frame(methods=c(rep("ENF",p*(p-1)/2),
                                          rep("MLE",p*(p-1)/2),
                                          rep("cPCG_theo",p*(p-1)/2),
                                          rep("cPCG_cv",p*(p-1)/2),
                                          rep("cPCG_df",p*(p-1)/2)),
                                rbind.data.frame(rank.ENF, rank.MLE, rank.cLevel_theo, rank.cLevel_cv, rank.cLevel_df))
    
    
    ROC.out=cbind.data.frame(FPR=seq(0.001,1,0.001), # Sample FPR cut points
                             ENF=rowMeans(tempstore.ENF), # averaged max TPR at cutted FPR level
                             MLE=rowMeans(tempstore.MLE),
                             cPCG_theo=rowMeans(tempstore.cLevel_theo),
                             cPCG_cv=rowMeans(tempstore.cLevel_cv),
                             cPCG_df=rowMeans(tempstore.cLevel_df))
    #ROC.table=tidyr::gather(ROC.out, methods, TPR, ENF:cPCG_df)
    
    
    # Number of total correctly selected edges
    All.enf.p=rowMeans(all.enf.p)
    All.enf.q=rowMeans(all.enf.q)
    All.mle.p=rowMeans(all.mle.p)
    All.mle.q=rowMeans(all.mle.q)
    All.clevel.theo=rowMeans(all.clevel_theo)
    All.clevel.cv=rowMeans(all.clevel_cv)
    All.clevel.df=rowMeans(all.clevel_df)
    
    # true positives
    TP.enf.p=rowMeans(tp.enf.p)
    TP.enf.q=rowMeans(tp.enf.q)
    TP.mle.p=rowMeans(tp.mle.p)
    TP.mle.q=rowMeans(tp.mle.q)
    TP.clevel.theo=rowMeans(tp.clevel_theo)
    TP.clevel.cv=rowMeans(tp.clevel_cv)
    TP.clevel.df=rowMeans(tp.clevel_df)
    
    # false positives
    FP.enf.p=rowMeans(fp.enf.p)
    FP.enf.q=rowMeans(fp.enf.q)
    FP.mle.p=rowMeans(fp.mle.p)
    FP.mle.q=rowMeans(fp.mle.q)
    FP.clevel.theo=rowMeans(fp.clevel_theo)
    FP.clevel.cv=rowMeans(fp.clevel_cv)
    FP.clevel.df=rowMeans(fp.clevel_df)
    
    # true negatives
    TN.enf.p=rowMeans(tn.enf.p)
    TN.enf.q=rowMeans(tn.enf.q)
    TN.mle.p=rowMeans(tn.mle.p)
    TN.mle.q=rowMeans(tn.mle.q)
    TN.clevel.theo=rowMeans(tn.clevel_theo)
    TN.clevel.cv=rowMeans(tn.clevel_cv)
    TN.clevel.df=rowMeans(tn.clevel_df)
    
    # false negatives
    FN.enf.p=rowMeans(fn.enf.p)
    FN.enf.q=rowMeans(fn.enf.q)
    FN.mle.p=rowMeans(fn.mle.p)
    FN.mle.q=rowMeans(fn.mle.q)
    FN.clevel.theo=rowMeans(fn.clevel_theo)
    FN.clevel.cv=rowMeans(fn.clevel_cv)
    FN.clevel.df=rowMeans(fn.clevel_df)
    
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
    fdr.clevel.df=FP.clevel.df/All.clevel.df
    fdr.clevel.df[is.na(fdr.clevel.df)]=0
    
    
    fdr.table=cbind.data.frame(
      FDR=rep(al,7), # nominal FDR
      methods=c(rep("ENF.p",length(al)),
                rep("ENF.q",length(al)),
                rep("MLE.p",length(al)),
                rep("MLE.q",length(al)),
                rep("cPCG.theo",length(al)),
                rep("cPCG.cv",length(al)),
                rep("cPCG.df",length(al))),
      fdr=c(fdr.enf.p, fdr.enf.q, fdr.mle.p, fdr.mle.q, fdr.clevel.theo, fdr.clevel.cv, fdr.clevel.df))
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
        cLevel_df.ordered.TPR=cLevel_df.ordered.TPR, cLevel_df.ordered.FPR=cLevel_df.ordered.FPR, 
        cLevel_df.ordered.PPV=cLevel_df.ordered.PPV, cLevel_df.ordered.FDR=cLevel_df.ordered.FDR,
        rank.table=rank.table,
        
        # empirical TPR
        tempstore.ENF=tempstore.ENF, tempstore.ENF=tempstore.ENF, tempstore.MLE=tempstore.MLE, 
        tempstore.cLevel_theo=tempstore.cLevel_theo, tempstore.cLevel_cv=tempstore.cLevel_cv, 
        tempstore.cLevel_df=tempstore.cLevel_df,
        ROC.out=ROC.out),
      
        # FDR
      metrics=list( 
        # Number of total correctly selected edges
        all.enf.p=all.enf.p, all.enf.q=all.enf.q,  
        all.mle.p=all.mle.p,  all.mle.q=all.mle.q,
        all.clevel_theo=all.clevel_theo, all.clevel_cv=all.clevel_cv, all.clevel_df=all.clevel_df,
        
        # true positives
        tp.enf.p=tp.enf.p, tp.enf.q=tp.enf.q, 
        tp.mle.p=tp.mle.p, tp.mle.q=tp.mle.q,
        tp.clevel_theo=tp.clevel_theo, tp.clevel_cv=tp.clevel_cv, tp.clevel_df=tp.clevel_df,
        
        # false positives
        fp.enf.p=fp.enf.p, fp.enf.q=fp.enf.q, 
        fp.mle.p=fp.mle.p, fp.mle.q=fp.mle.q,
        fp.clevel_theo=fp.clevel_theo, fp.clevel_cv=fp.clevel_cv, fp.clevel_df=fp.clevel_df,
        
        # false negatives
        tn.enf.p=tn.enf.p, tn.enf.q=tn.enf.q, 
        tn.mle.p=tn.mle.p, tn.mle.q=tn.mle.q,
        tn.clevel_theo=tn.clevel_theo, tn.clevel_cv=tn.clevel_cv, tn.clevel_df=tn.clevel_df,
        
        # false negatives
        fn.enf.p=fn.enf.p, fn.enf.q=fn.enf.q,
        fn.mle.p=fn.mle.p, fn.mle.q=fn.mle.q,
        fn.clevel_theo=fn.clevel_theo, fn.clevel_cv=fn.clevel_cv, fn.clevel_df=fn.clevel_df,
        
        All.enf.p=All.enf.p, All.enf.q=All.enf.q, 
        All.mle.p=All.mle.p, All.mle.q=All.mle.q,
        All.clevel.theo=All.clevel.theo, All.clevel.cv=All.clevel.cv, All.clevel.df=All.clevel.df,
        
        # true positives
        TP.enf.p=TP.enf.p,TP.enf.q=TP.enf.q, TP.mle.p=TP.mle.p,TP.mle.q=TP.mle.q, 
        TP.clevel.theo=TP.clevel.theo,TP.clevel.cv=TP.clevel.cv,TP.clevel.df=TP.clevel.df,
        
        # false positives
        FP.enf.p=FP.enf.p,FP.enf.q=FP.enf.q, 
        FP.mle.p=FP.mle.p, FP.mle.q=FP.mle.q,
        FP.clevel.theo=FP.clevel.theo, FP.clevel.cv=FP.clevel.cv,FP.clevel.df=FP.clevel.df,
        
        # true negatives
        TN.enf.p=TN.enf.p, TN.enf.q=TN.enf.q, 
        TN.mle.p=TN.mle.p, TN.mle.q=TN.mle.q,
        TN.clevel.theo=TN.clevel.theo, TN.clevel.cv=TN.clevel.cv, TN.clevel.df=TN.clevel.df,
        
        # false negatives
        FN.enf.p=FN.enf.p, FN.enf.q=FN.enf.q, 
        FN.mle.p=FN.mle.p,FN.mle.q=FN.mle.q,
        FN.clevel.theo=FN.clevel.theo,FN.clevel.cv=FN.clevel.cv, FN.clevel.df=FN.clevel.df
        ),
      # empirical fdr
      fdr=list(
        fdr.enf.p=fdr.enf.p, fdr.enf.q=fdr.enf.q, 
        fdr.mle.p=fdr.mle.p, fdr.mle.q=fdr.mle.q,
        fdr.clevel.theo=fdr.clevel.theo, fdr.clevel.cv=fdr.clevel.cv, fdr.clevel.df=fdr.clevel.df,
        fdr.table=fdr.table
        )
      )
    )
  } # return
}



# Block Diagonal ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
#pl = c(100)
#min.beta = 0.3
#max.beta = 1
#blocksize = c(4,10)
min.beta = 0.05
max.beta = .5
blocksize = c(4,8,10)
for (n in nl) {
  for (p in pl) {
    for (e in blocksize) {
      load( file=sprintf("simu_data_v2/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta) )
      
      RESULT=model.eval(X, pcor, rep = 10)
      save(RESULT, file=sprintf("temp/6_12/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta) )
      
    }
  }
}


# Scale Free ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
#pl = c(100)
min.beta = 0.3
edge = c(1)
for (n in nl) {
  for (p in pl) {
    for (e in edge) {
      load( file=sprintf("simu_data_v2/ScaleFree_simu_n%d_p%d_e%d_min_beta%g.RData", n, p, e, min.beta) )
      
      RESULT=model.eval(X, pcor, rep = 10)
      save(RESULT, file=sprintf("temp/6_12/ScaleFree_simu_n%d_p%d_e%d_min_beta%g.RData", n, p, e, min.beta) )
    }
  }
}


# Random Structure ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
#pl = c(100)
etal = c(.01, .03, .05)
min.beta = 0.3
for (n in nl) {
  for (p in pl) {
    for (eta in etal) {
      load( file=sprintf("simu_data_v2/Random_simu_n%d_p%d_eta%g_min_beta%g.RData", n, p, eta, min.beta) )
      
      RESULT=model.eval(X, pcor, rep = 5)
      save(RESULT, file=sprintf("temp/6_12/Random_simu_n%d_p%d_eta%g_min_beta%g.RData", n, p, eta, min.beta) )
      
    }
  }
}


