library(tidyverse)
library(GeneNet)
library(FastGGM)
library(corpcor)
library(glmnet)

source("R/Qiu.R")
source("R/PCGII.R")
source("R/shrinkagefunctions.R") # https://github.com/V-Bernal/GGM-Shrinkage

# Scale Free ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
el = c(1, 2, 3)
for (n in nl) {
  for (p in pl) {
    for (e in el) {
      load( file=sprintf("simu_data/ScaleFree_simu_n%d_p%d_e%d.RData", n, p, e) )
      degf=n/5
      RESULT=model.eval(X, omega, rep = 20, degfree = degf, prior_prop=0.1)
      save(RESULT, file=sprintf("results/ScaleFree_n%d_p%d_e%d.RData", n, p, e) )
    }
  }
}


# Random ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
etal = c(.01, .02, .03)
for (n in nl) {
  for (p in pl) {
    for (eta in etal) {
      load( file=sprintf("simu_data/Random_simu_n%d_p%d_eta%g.RData", n, p, eta) )
      degf=n/5
      RESULT=model.eval(X, omega, rep = 20, degfree = degf, prior_prop=0.1)
      save(RESULT, file=sprintf("results/Random_n%d_p%d_eta%g.RData", n, p, eta) )
    }
  }
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
      degf=n/5
      RESULT=model.eval(X, omega, rep = 20, degfree = degf, prior_prop=0.1)
      save(RESULT, file=sprintf("results/BlockDiag_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta) )
    }
  }
}



# plots ----
ScaleFree=list()
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
el = c(1, 2, 3)
i=1
for (n in nl) {
  for (p in pl) {
    for (e in el) {
      load(file=sprintf("results/ScaleFree_n%d_p%d_e%d.RData", n, p, e))
      ScaleFree[[i]]=list(RESULT=RESULT, n=n, p=p, e=e, model="scale free")
      i=i+1
    }
  }
}

Random=list()
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
el = c(.01, .02, .03)
i=1
for (n in nl) {
  for (p in pl) {
    for (e in el) {
      load(file=sprintf("results/Random_n%d_p%d_eta%g.RData", n, p, e))
      Random[[i]]=list(RESULT=RESULT, n=n, p=p, e=e, model="random")
      i=i+1
    }
  }
}

Summary_ROC=function(mylist){
  tempres=mylist$RESULT
  tempn=mylist$n
  tempp=mylist$p
  tempe=mylist$e
  model=mylist$model
  tempROC=tempres$ROC$ROC.out %>% 
    tidyr::gather("methods", "TPR", ENF:F_GGM) %>%
    filter(!methods %in% c("ENF")) %>% 
    mutate(n=tempn,p=tempp,e=tempe, model=model)
  tempROC
}

Summary_FDR=function(mylist){
  tempres=mylist$RESULT
  tempn=mylist$n
  tempp=mylist$p
  tempe=mylist$e
  model=mylist$model
  tempfdr=tempres$fdr$fdr.table %>% 
    filter(!methods %in% c("ENF.p","MLE.p")) %>% 
    mutate(n=tempn,p=tempp,e=tempe, model=model)
  tempfdr
}


lapply(ScaleFree, Summary_ROC) %>% bind_rows() %>% 
  ggplot(aes(x=FPR, y=TPR,col=methods)) + 
  geom_line(aes(linetype=methods))+
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +
  scale_color_manual(values=c(6,6,4,3,2,2)) +
  ggtitle("Scale Free, ROC") +
  geom_abline(slope = 1, intercept = 0, size=0.3) + facet_wrap(n~p*e, ncol = 3)

lapply(Random, Summary_ROC) %>% bind_rows() %>% 
  ggplot(aes(x=FPR, y=TPR,col=methods)) + 
  geom_line(aes(linetype=methods))+
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +
  scale_color_manual(values=c(6,6,4,3,2,2)) +
  ggtitle("Random, ROC") +
  geom_abline(slope = 1, intercept = 0, size=0.3) + facet_wrap(n~p*e, ncol = 3)

lapply(ScaleFree, Summary_FDR) %>% bind_rows() %>% 
  ggplot(aes(x=FDR, y=fdr, col=methods))+
  geom_line(aes(linetype=methods))+
  scale_linetype_manual(values=c(1,2,7,3,4,5,6)) +
  scale_color_manual(values=c(6,6,5,3,4,2,2)) +
  ggtitle("Scale Free, FDR") +
  geom_abline(slope = 1, intercept = 0)+ facet_wrap(n~p*e, ncol = 3)

lapply(Random, Summary_FDR) %>% bind_rows() %>% 
  ggplot(aes(x=FDR, y=fdr, col=methods))+
  geom_line(aes(linetype=methods))+
  scale_linetype_manual(values=c(1,2,7,3,4,5,6)) +
  scale_color_manual(values=c(6,6,5,3,4,2,2)) +
  ggtitle("Random, FDR") +
  geom_abline(slope = 1, intercept = 0)+ facet_wrap(n~p*e, ncol = 3)

my_power=function(RESULT, alpha=0.05){
  ref=seq(0.001, 0.2005, 0.0035)
  ind=which(ref==alpha)
  recall=cbind(RESULT$metrics$tp.enf.p[ind,]/(RESULT$metrics$tp.enf.p[ind,]+RESULT$metrics$fn.enf.p[ind,]), # TPR
               RESULT$metrics$tp.enf.q[ind,]/(RESULT$metrics$tp.enf.q[ind,]+RESULT$metrics$fn.enf.q[ind,]),
               RESULT$metrics$tp.mle.p[ind,]/(RESULT$metrics$tp.mle.p[ind,]+RESULT$metrics$fn.mle.p[ind,]),
               RESULT$metrics$tp.mle.q[ind,]/(RESULT$metrics$tp.mle.q[ind,]+RESULT$metrics$fn.mle.q[ind,]),
               RESULT$metrics$tp.cLevel_theo[ind,]/(RESULT$metrics$tp.cLevel_theo[ind,]+RESULT$metrics$fn.cLevel_theo[ind,]),
               RESULT$metrics$tp.cLevel_cv[ind,]/(RESULT$metrics$tp.cLevel_cv[ind,]+RESULT$metrics$fn.cLevel_cv[ind,]),
               RESULT$metrics$tp.PCGII_theo[ind,]/(RESULT$metrics$tp.PCGII_theo[ind,]+RESULT$metrics$fn.PCGII_theo[ind,]),
               RESULT$metrics$tp.PCGII_cv[ind,]/(RESULT$metrics$tp.PCGII_cv[ind,]+RESULT$metrics$fn.PCGII_cv[ind,]),
               RESULT$metrics$tp.F_GGM[ind,]/(RESULT$metrics$tp.F_GGM[ind,]+RESULT$metrics$fn.F_GGM[ind,])
  )
  
  
  power=as.data.frame(recall)
  colnames(power)=c("ENF.p","ENF.q","MLE.p","MLE.q","cLevel.theo","cLevel.cv","PCGII.theo","PCGII.cv","FGGM")
  power
}
Summary_power=function(mylist){
  tempres=mylist$RESULT
  tempn=mylist$n
  tempp=mylist$p
  tempe=mylist$e
  model=mylist$model
  temppower=my_power(tempres) %>% 
    rename(F_GGM=FGGM) %>% 
    tidyr::gather("methods", "power", ENF.p:F_GGM) %>%
    filter(!methods %in% c("ENF.p","MLE.p")) %>% 
    mutate(n=tempn,p=tempp,e=tempe, model=model)
  temppower
}
lapply(ScaleFree, Summary_power) %>% bind_rows() %>% 
  ggplot(aes(x=methods, y=power, col=methods))+
  geom_boxplot() + 
  scale_color_manual(values=c(6,6,5,3,4,2,2)) +
  ggtitle("Scale Free, Power") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) ) +
  facet_wrap(n~p*e, ncol = 3)

lapply(Random, Summary_power) %>% bind_rows() %>% 
  ggplot(aes(x=methods, y=power, col=methods))+
  geom_boxplot() + 
  scale_color_manual(values=c(6,6,5,3,4,2,2)) +
  ggtitle("Random, Power") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) ) +
  facet_wrap(n~p*e, ncol = 3)

# PPV=1-FDR=TP/(TP+FP)
my_PPV=function(RESULT, alpha=0.05){
  ref=seq(0.001, 0.2005, 0.0035)
  ind=which(ref==alpha)
  is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
  PPV=cbind.data.frame(RESULT$metrics$tp.enf.p[ind,]/(RESULT$metrics$tp.enf.p[ind,]+RESULT$metrics$fp.enf.p[ind,]), # TPR
                       RESULT$metrics$tp.enf.q[ind,]/(RESULT$metrics$tp.enf.q[ind,]+RESULT$metrics$fp.enf.q[ind,]),
                       RESULT$metrics$tp.mle.p[ind,]/(RESULT$metrics$tp.mle.p[ind,]+RESULT$metrics$fp.mle.p[ind,]),
                       RESULT$metrics$tp.mle.q[ind,]/(RESULT$metrics$tp.mle.q[ind,]+RESULT$metrics$fp.mle.q[ind,]),
                       RESULT$metrics$tp.cLevel_theo[ind,]/(RESULT$metrics$tp.cLevel_theo[ind,]+RESULT$metrics$fp.cLevel_theo[ind,]),
                       RESULT$metrics$tp.cLevel_cv[ind,]/(RESULT$metrics$tp.cLevel_cv[ind,]+RESULT$metrics$fp.cLevel_cv[ind,]),
                       RESULT$metrics$tp.PCGII_theo[ind,]/(RESULT$metrics$tp.PCGII_theo[ind,]+RESULT$metrics$fp.PCGII_theo[ind,]),
                       RESULT$metrics$tp.PCGII_cv[ind,]/(RESULT$metrics$tp.PCGII_cv[ind,]+RESULT$metrics$fp.PCGII_cv[ind,]),
                       RESULT$metrics$tp.F_GGM[ind,]/(RESULT$metrics$tp.F_GGM[ind,]+RESULT$metrics$fp.F_GGM[ind,])
  )
  
  colnames(PPV)=c("ENF.p","ENF.q","MLE.p","MLE.q","cLevel.theo","cLevel.cv","PCGII.theo","PCGII.cv","FGGM")
  PPV[is.nan(PPV)] <- 0
  PPV
}
my_fdr=function(RESULT, alpha=0.05){
  ref=seq(0.001, 0.2005, 0.0035)
  ind=which(ref==alpha)
  is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
  fdr=cbind.data.frame(RESULT$metrics$fp.enf.p[ind,]/(RESULT$metrics$tp.enf.p[ind,]+RESULT$metrics$fp.enf.p[ind,]), 
                       RESULT$metrics$fp.enf.q[ind,]/(RESULT$metrics$tp.enf.q[ind,]+RESULT$metrics$fp.enf.q[ind,]),
                       RESULT$metrics$fp.mle.p[ind,]/(RESULT$metrics$tp.mle.p[ind,]+RESULT$metrics$fp.mle.p[ind,]),
                       RESULT$metrics$fp.mle.q[ind,]/(RESULT$metrics$tp.mle.q[ind,]+RESULT$metrics$fp.mle.q[ind,]),
                       RESULT$metrics$fp.cLevel_theo[ind,]/(RESULT$metrics$tp.cLevel_theo[ind,]+RESULT$metrics$fp.cLevel_theo[ind,]),
                       RESULT$metrics$fp.cLevel_cv[ind,]/(RESULT$metrics$tp.cLevel_cv[ind,]+RESULT$metrics$fp.cLevel_cv[ind,]),
                       RESULT$metrics$fp.PCGII_theo[ind,]/(RESULT$metrics$tp.PCGII_theo[ind,]+RESULT$metrics$fp.PCGII_theo[ind,]),
                       RESULT$metrics$fp.PCGII_cv[ind,]/(RESULT$metrics$tp.PCGII_cv[ind,]+RESULT$metrics$fp.PCGII_cv[ind,]),
                       RESULT$metrics$fp.F_GGM[ind,]/(RESULT$metrics$tp.F_GGM[ind,]+RESULT$metrics$fp.F_GGM[ind,])
  )
  
  colnames(fdr)=c("ENF.p","ENF.q","MLE.p","MLE.q","cLevel.theo","cLevel.cv","PCGII.theo","PCGII.cv","FGGM")
  fdr[is.nan(fdr)] <- 0
  fdr
}
Summary_PPV=function(mylist){
  tempres=mylist$RESULT
  tempn=mylist$n
  tempp=mylist$p
  tempe=mylist$e
  model=mylist$model
  tempppv=my_PPV(tempres) %>% 
    rename(F_GGM=FGGM) %>% 
    tidyr::gather("methods", "PPV", ENF.p:F_GGM) %>%
    filter(!methods %in% c("ENF.p","MLE.p")) %>% 
    mutate(n=tempn,p=tempp,e=tempe, model=model)
  tempppv
}

Summary_fdr=function(mylist){
  tempres=mylist$RESULT
  tempn=mylist$n
  tempp=mylist$p
  tempe=mylist$e
  model=mylist$model
  tempfdr=my_fdr(tempres) %>% 
    rename(F_GGM=FGGM) %>% 
    tidyr::gather("methods", "fdr", ENF.p:F_GGM) %>%
    filter(!methods %in% c("ENF.p","MLE.p")) %>% 
    mutate(n=tempn,p=tempp,e=tempe, model=model)
  tempfdr
}

lapply(ScaleFree, Summary_fdr) %>% bind_rows() %>% 
  ggplot(aes(x=methods, y=fdr, col=methods))+
  geom_boxplot() + 
  geom_abline(intercept = 0.05, slope=0) +
  scale_color_manual(values=c(6,6,5,3,4,2,2)) +
  ggtitle("Scale Free, fdr") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) ) +
  facet_wrap(n~p*e, ncol = 3)

lapply(Random, Summary_fdr) %>% bind_rows() %>% 
  ggplot(aes(x=methods, y=power, col=methods))+
  geom_boxplot() + 
  geom_abline(intercept = 0.05, slope=0) +
  scale_color_manual(values=c(6,6,5,3,4,2,2)) +
  ggtitle("Random, fdr") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) ) +
  facet_wrap(n~p*e, ncol = 3)



