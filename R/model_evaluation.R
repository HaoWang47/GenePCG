library(ggplot2)
library(tidyverse)
library(GeneNet)

source("R/Qiu.R")
source("R/shrinkagefunctions.R") # https://github.com/V-Bernal/GGM-Shrinkage

sigs2vec=function(sigs){
  require(corpcor)
  m=matrix(0,P,P)
  for (h in 1: dim(sigs)[1]){
    m[sigs[h,1],sigs[h,2]]=1
  }
  sm2vec(m)
}

model_results=function(sim.data, pcor) {
  require(GeneNet)
  P=dim(sim.data)[2]; N=dim(sim.data)[1];
  truth<-sm2vec(pcor) 
  TP<- which(truth!=0) 
  
  cPCG_theo=cPCG_theo(sim.data)
  cPCG_cv=cPCG_cv_dfmax(sim.data)
  cPCG_df=cPCG_cv_dfmax(sim.data, degree_freedom = 5)
  
  # estimates by cPCG
  Est_theo=sm2vec(cPCG_theo$Est)
  Est_cv=sm2vec(cPCG_cv$Est)
  Est_df=sm2vec(cPCG_df$Est)
  
  # test statistics of cPCG
  tscore_theo=sm2vec(cPCG_theo$tscore)
  tscore_cv=sm2vec(cPCG_cv$tscore)
  tscore_df=sm2vec(cPCG_df$tscore)
  
  ### Shrunk Partial Cor
  { GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
    lambda<-attr(GGM, "lambda") 
    while (lambda == 1) { # exclude the case of complete shrinkage (lambda=1)
      sim.data=mvrnorm(n=N,mu=rep(0,P),Sigma=Cov)
      GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
      lambda <- attr(GGM, "lambda") 
      print("TRY AGAIN")
    }  
  }
  shrunk_p=sm2vec(GGM) # off diagonal elements of estimated shrunk partial corr matrix
  
  ### P values by Empirical null fitting (ENF)
  ENF.test=network.test.edges(GGM, fdr=TRUE, plot=FALSE,verbose=FALSE) 
  ENF.test=ENF.test[order(ENF.test$node1,ENF.test$node2),] # node1 is col; WATCH OUT:  the test must be in node's order. Sort by node 2 and by node 1
  
  ### P values by Shrunk MLE
  pval.shrunk=p.shrunk(shrunk_p,P,N,lambda)
  
  ### cPCG inference
  PCG_theo_sigs=inference(cPCG_theo)
  PCG_cv_sigs=inference(cPCG_cv)
  PCG_df_sigs=inference(cPCG_df)
  
  results=cbind.data.frame(truth, Est_theo, Est_cv, Est_df, shrunk_p, 
                           tscore_theo, tscore_cv, tscore_df,
                           sigs_theo=sigs2vec(PCG_theo_sigs$sigs),
                           sigs_cv=sigs2vec(PCG_cv_sigs$sigs),
                           sigs_df=sigs2vec(PCG_df_sigs$sigs),
                           ENF_p=ENF.test$pval, ENF_q=p.adjust(ENF.test$pval, method="BH"),
                           ShrunkMLE_p=pval.shrunk, ShrunkMLE_q=p.adjust(pval.shrunk, method="BH"))
  
  return(results)
}

# Random Structure ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
etal = c(0.01,0.02,0.03)
min.beta = 0.3
for (n in nl) {
  for (p in pl) {
    for (eta in etal) {
      load( file=sprintf("simu_data/Random_simu_n%d_p%d_eta%g_min_beta%g.RData", n, p, eta, min.beta) )
    }
  }
}


# Block Diagonal ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta = 0.3
max.beta = 1
blocksize = c(4,10)
for (n in nl) {
  for (p in pl) {
    for (e in blocksize) {
      load( file=sprintf("simu_data/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta) )
      RESULT=list()
      for (i in 1:2){
        RESULT[[i]]=model_results(sim.data = X[[i]], pcor = PCOR[[i]])
      }
      save(RESULT, file=sprintf("results/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta))
    }
  }
}


# Scale Free ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta = 0.1
edge = c(1,2)
for (n in nl) {
  for (p in pl) {
    for (e in edge) {
      load( file=sprintf("simu_data/ScaleFree_simu_n%d_p%d_e%d_min_beta%g.RData", n, p, e, min.beta) )
    }
  }
}
