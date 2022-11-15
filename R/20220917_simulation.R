############# Last updated on 20220919 #############;

library(igraph)
library(tidyverse)
library(GeneNet)
library(FastGGM)
library(corpcor)
library(glmnet)
library(mvtnorm)

source("./R/shrinkagefunctions.R") # https://github.com/V-Bernal/GGM-Shrinkage
source("./R/PCGII.R")
source("./R/Utility.R")
############# Simulate Data, 2022/05/12 #############;

makeBlockDiag=function(blocksize=4, p=20, min.beta=0.2, max.beta=0.5){ # blocksize has to be a factor of p
  reps=p/blocksize
  S=list()
  for (i in 1:reps) {
    bd=matrix(0,blocksize,blocksize)
    for(j1 in 1:(blocksize-1)){
      for (j2 in (j1+1):blocksize) {
        bd[j1,j2]=runif(1,min.beta,max.beta)
        bd[j2,j1]=bd[j1,j2]
      }
    }
    diag(bd)=runif(1,1,1.25)
    while(!is.positive.definite(bd)){diag(bd)=diag(bd)+0.01}
    S[[i]]=bd
  }
  as.matrix(Matrix::bdiag(S))
}

#### Strong Signal ####

## Scale Free 1
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220520)

for(p in pl){
  for(e in 1:3){
    g <- sample_pa(p, power=1, m=e, directed = FALSE)
    omega=as_adjacency_matrix(g) %>% as.matrix()
    for(h1 in 1:(p-1)){
      for(h2 in (h1+1):p){
        if(omega[h1,h2]!=0){
          temp=runif(1, 0.2, 0.5)*sample(c(-1,1),size=1)
          omega[h1,h2]=temp
          omega[h2,h1]=temp
        }
      }
    }
    
    diag(omega)=rowSums(abs(omega))
    if(e==1){diag(omega)=diag(omega)+0.10}
    if(e==2){diag(omega)=diag(omega)+0.1}
    if(e==3){diag(omega)=diag(omega)+0.1}
    
    temp=eigen(omega)$values
    print(paste0("p=",p,",e=",e,",ratio=",max(temp)/min(temp)))
    # diag(omega)=4
    # increment=.05
    # temp=eigen(omega)$values
    # while (min(temp)<0 | max(temp)/min(temp)>50) {
    #   diag(omega)=4+increment
    #   increment=increment+0.05
    #   temp=eigen(omega)$values
    # }
    ppc=-sm2vec(cov2cor(omega))
    ppc=ppc[which(ppc!=0)]
    print(summary(abs(ppc)))
    sparsity=round(length(ppc)/(p*(p-1)/2),3)
    print(paste0("number of true edges: ",length(ppc)))
    Sigma=solve(omega)
    #is.positive.definite(Sigma)
    for(n in nl){  
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/ScaleFree1_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}


## Scale Free.5
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220520)

for(p in pl){
  for(e in 1:3){
    g <- sample_pa(p, power=.5, m=e, directed = FALSE)
    omega=as_adjacency_matrix(g) %>% as.matrix()
    for(h1 in 1:(p-1)){
      for(h2 in (h1+1):p){
        if(omega[h1,h2]!=0){
          temp=runif(1, 0.2, 0.5)*sample(c(-1,1),size=1)
          omega[h1,h2]=temp
          omega[h2,h1]=temp
        }
      }
    }
    
    diag(omega)=rowSums(abs(omega))
    if(e==1){diag(omega)=diag(omega)+0.1}
    if(e==2){diag(omega)=diag(omega)+0.1}
    if(e==3){diag(omega)=diag(omega)+0.1}
    
    temp=eigen(omega)$values
    print(paste0("p=",p,",e=",e,",ratio=",max(temp)/min(temp)))
    ppc=-sm2vec(cov2cor(omega))
    ppc=ppc[which(ppc!=0)]
    print(summary(abs(ppc)))
    sparsity=round(length(ppc)/(p*(p-1)/2),3)
    print(paste0("number of true edges: ",length(ppc)))
    Sigma=solve(omega)
    #is.positive.definite(Sigma)
    for(n in nl){    
      
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/ScaleFree.5_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Scale Free.1
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220520)

for(p in pl){
  for(e in 1:3){
    g <- sample_pa(p, power=.1, m=e, directed = FALSE)
    omega=as_adjacency_matrix(g) %>% as.matrix()
    for(h1 in 1:(p-1)){
      for(h2 in (h1+1):p){
        if(omega[h1,h2]!=0){
          temp=runif(1, 0.2, 0.5)*sample(c(-1,1),size=1)
          omega[h1,h2]=temp
          omega[h2,h1]=temp
        }
      }
    }
    
    diag(omega)=rowSums(abs(omega))
    if(e==1){diag(omega)=diag(omega)+0.1}
    if(e==2){diag(omega)=diag(omega)+0.1}
    if(e==3){diag(omega)=diag(omega)+0.1}
    
    temp=eigen(omega)$values
    print(paste0("p=",p,",e=",e,",ratio=",max(temp)/min(temp)))
    ppc=-sm2vec(cov2cor(omega))
    ppc=ppc[which(ppc!=0)]
    print(summary(abs(ppc)))
    sparsity=round(length(ppc)/(p*(p-1)/2),3)
    print(paste0("number of true edges: ",length(ppc)))
    Sigma=solve(omega)
    #is.positive.definite(Sigma)
    for(n in nl){    
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/ScaleFree.1_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}


## Random 20220521
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220520)

for(p in pl){
  for(eta in c(0.01,0.02,0.03)){
    g <- sample_gnp(p, eta, directed = FALSE)
    omega=as_adjacency_matrix(g) %>% as.matrix()
    for(h1 in 1:(p-1)){
      for(h2 in (h1+1):p){
        if(omega[h1,h2]!=0){
          temp=runif(1, 0.2, 0.5)*sample(c(-1,1),size=1)
          omega[h1,h2]=temp
          omega[h2,h1]=temp
        }
      }
    }
    
    diag(omega)=rowSums(abs(omega))
    if(e==1){diag(omega)=diag(omega)+0.1}
    if(e==2){diag(omega)=diag(omega)+0.1}
    if(e==3){diag(omega)=diag(omega)+0.1}
    
    temp=eigen(omega)$values
    print(paste0("p=",p,",eta=",eta,",ratio=",max(temp)/min(temp)))
    ppc=-sm2vec(cov2cor(omega))
    ppc=ppc[which(ppc!=0)]
    print(summary(abs(ppc)))
    sparsity=round(length(ppc)/(p*(p-1)/2),3)
    print(paste0("number of true edges: ",length(ppc)))
    Sigma=solve(omega)
    #is.positive.definite(Sigma)
    for(n in nl){
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta)) 
    }
  }
}

## Random 20220616
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220616)

for(p in pl){
  for(eta in c(0.01,0.02,0.03)){
    g <- sample_gnp(p, eta, directed = FALSE)
    omega=as_adjacency_matrix(g) %>% as.matrix()
    for(h1 in 1:(p-1)){
      for(h2 in (h1+1):p){
        if(omega[h1,h2]!=0){
          temp=runif(1, 0.2, 0.5)*sample(c(-1,1),size=1)
          omega[h1,h2]=temp
          omega[h2,h1]=temp
        }
      }
    }
    
    diag(omega)=rowSums(abs(omega))
    diag(omega)=diag(omega)+0.1

    
    temp=eigen(omega)$values
    print(paste0("p=",p,",eta=",eta,",ratio=",max(temp)/min(temp)))
    ppc=-sm2vec(cov2cor(omega))
    ppc=ppc[which(ppc!=0)]
    print(summary(abs(ppc)))
    sparsity=round(length(ppc)/(p*(p-1)/2),3)
    print(paste0("number of true edges: ",length(ppc)))
    Sigma=solve(omega)
    #is.positive.definite(Sigma)
    for(n in nl){
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta)) 
    }
  }
}


## Block Diagonal
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220520)

for(p in pl){
  for(e in c(4,5,8,10)){
    Sigma=makeBlockDiag(blocksize=e, p=p, min.beta=.3, max.beta=.5)
    omega=solve(Sigma)
    temp=eigen(omega)$values
    print(paste0("p=",p,",e=",e,",ratio=",max(temp)/min(temp)))
    ppc=-sm2vec(cov2cor(omega))
    ppc=ppc[which(ppc!=0)]
    print(summary(abs(ppc)))
    sparsity=round(length(ppc)/(p*(p-1)/2),3)
    print(paste0("number of true edges: ",length(ppc)))

    for(n in nl){  
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/BlockDiag_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}


#### Functions ####

eval_models=function(X, omega, rep=3, Seed=1234){
  set.seed(Seed)
  # here, we want to check how incorrect prior affects our approach
  
  p=dim(omega)[1] 
  n=dim(X[[1]])[1]
  truth<-sm2vec(omega, diag = F) 
  TP<-which(truth!=0) 
  CN=p*(p-1)/2-length(TP) 
  CP=length(TP)
  
  # lambda smaller, less penalty, less correlation in residuals
  lam=sqrt(2*log(p/sqrt(n))/n) ## FastGGM default lambda
    
  true.par=-cov2cor(omega)
  diag(true.par)=0
  
  prior.all=cbind.data.frame(which(true.par!=0, arr.ind = TRUE), signal=round(true.par[which(true.par!=0)],6)) %>% 
    transform(row = pmin(row, col), col = pmax(row, col)) %>% 
    arrange(row, col) %>% 
    unique() %>%
    arrange(desc(abs(signal))) ## all off-diag !=0
  prior.wrong=cbind.data.frame(which(true.par==0, arr.ind = TRUE), signal=0) %>%
    transform(row = pmin(row, col), col = pmax(row, col)) %>% 
    arrange(row, col) %>% 
    unique() %>% 
    filter(row!=col)
  prior.wrong=prior.wrong[sample(nrow(prior.wrong)),] # shuffle
  
  if(dim(prior.all)[1]!=CP) {stop("Prior information error!")}

  #create ID
  IDs=lapply(c(1,.7,.3), function(x){ # all from first CP's edges
    sample(1:CP, round(CP*x,0))
  })

  prior_all=lapply(IDs, function(x){ #100-0, 70-30, 30-70
    temp=prior.all[x,]
    if(nrow(temp)<CP){
      temp=rbind(temp,prior.wrong[1:(CP-nrow(temp)),])
    }
    temp
  })
  names(prior_all)=c("100%","70%","30%")
  
  prior_pct=lapply(list(IDs[[2]],IDs[[3]]), function(ID){ # "70%","30%"
    prior=lapply(c(1,0.7,0.3), function(x){ 
      temp=prior.all[ID,] %>% filter(signal!=0)
      if(x!=1){
        wrong=nrow(temp)*(1-x)
        temp[1:wrong,]=prior.wrong[1:wrong,]  
      }
      temp
    })
    names(prior)=c("100%right","70%right","30%right")
    prior
  })
  names(prior_pct)=c("70%ofALL","30%ofALL")
  
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
    
    cLevel.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_all.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_all.30.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all.30.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all.30.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all.30.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_all.70.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all.70.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all.70.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all.70.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_30perct.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_30perct.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_30perct.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_30perct.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_30perct.30.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_30perct.30.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_30perct.30.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_30perct.30.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_30perct.70.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_30perct.70.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_30perct.70.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_30perct.70.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_70perct.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_70perct.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_70perct.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_70perct.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_70perct.30.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_70perct.30.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_70perct.30.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_70perct.30.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_70perct.70.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_70perct.70.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_70perct.70.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_70perct.70.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    FGGM.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep) 
    
    tempstore.ENF=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.MLE=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points 
    tempstore.cLevel=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_all=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_all.30=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_all.70=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points

    tempstore.PCGII_30perct=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_30perct.30=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_30perct.70=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_70perct=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_70perct.30=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_70perct.70=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.FGGM=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    
    ### Inference
    
    # Number of total correctly selected edges
    all.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    all.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.mle.p=matrix(Inf, nrow=length(al), ncol=rep) 
    all.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.cLevel=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_all=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_all.30=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_all.70=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_30perct=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_30perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_30perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_70perct=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_70perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_70perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    all.FGGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # true positives
    tp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.cLevel=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_all=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_all.30=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_all.70=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_30perct=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_30perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_30perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_70perct=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_70perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_70perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    tp.FGGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false positives
    fp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.cLevel=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_all=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_all.30=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_all.70=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_30perct=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_30perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_30perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_70perct=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_70perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_70perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    fp.FGGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # true negatives
    tn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.cLevel=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_all=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_all.30=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_all.70=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_30perct=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_30perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_30perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_70perct=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_70perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_70perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    tn.FGGM=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false negatives
    fn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.cLevel=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_all=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_all.30=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_all.70=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_30perct=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_30perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_30perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_70perct=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_70perct.30=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_70perct.70=matrix(Inf, nrow=length(al), ncol=rep)
    fn.FGGM=matrix(Inf, nrow=length(al), ncol=rep)
  }
  
  
  for (k in 1:rep){
    print(paste0(k,"th rep starts!!!!!"))
    sim.data=X[[k]]
    
    ### Fast GGM
    out=FastGGM(sim.data, lambda=lam)
    FGGM.test=unMat(X_est=out$partialCor, X_p=out$p_partialCor) # Est, pvals
    FGGM.test=cbind.data.frame(FGGM.test, truth=truth, qval=p.adjust(FGGM.test[,4], method="BH"))
    
    print("Fast GGM done")
    
    ### cLevel
    cLevel=clevel(sim.data, lambda = lam)
    # estimates by cLevel
    Est_cLevel=sm2vec(cLevel$Est, diag = F)
    # test statistics of cLevel
    tscore_cLevel=sm2vec(cLevel$tscore, diag = F)
    print("cLevel done")
    
    ### PCGII-prior.all
    PCGII_all=lapply(prior_all, function(pri){
      temp=PCGII(df=sim.data, prior=double_prior(pri), lambda=lam)
      list(out=temp, Est=sm2vec(temp$Est, diag=F), tscore=sm2vec(temp$tscore, diag = F))
    })

    ### PCGII-70/30 of true edges
    PCGII_pct=lapply(prior_pct, function(pri){
      lapply(pri, function(prii){
        temp=PCGII(df=sim.data, prior=double_prior(prii), lambda=lam)
        list(out=temp, Est=sm2vec(temp$Est, diag=F), tscore=sm2vec(temp$tscore, diag = F))
      })
    })
    
    print("PCGII done")
    
    ### Shrunk Partial Cor
    GGM=pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
    lambda=attr(GGM, "lambda") 
    while (lambda == 1){lambda=0.99999}
    shrunk_p=sm2vec(GGM, diag = F) # off diagonal elements of estimated shrunk partial corr matrix
    # P values by Empirical null fitting (ENF)
    ENF.test=network.test.edges(GGM, fdr=TRUE, plot=FALSE,verbose=FALSE) 
    ENF.test=ENF.test[order(ENF.test$node1,ENF.test$node2),] # node1 is col; WATCH OUT:  the test must be in node's order. Sort by node 2 and by node 1
    # P values by Shrunk MLE
    pval.shrunk=p.shrunk(shrunk_p,p,n,lambda)
    print("Shrunk MLE done")
    
    
    results=cbind.data.frame(truth,
                             Est_cLevel, 
                             Est_PCGII_all=PCGII_all$`100%`$Est, 
                             Est_PCGII_all.70=PCGII_all$`70%`$Est, 
                             Est_PCGII_all.30=PCGII_all$`30%`$Est, 
                             Est_PCGII_70perct=PCGII_pct$`70%ofALL`$`100%right`$Est, 
                             Est_PCGII_70perct.70=PCGII_pct$`70%ofALL`$`70%right`$Est, 
                             Est_PCGII_70perct.30=PCGII_pct$`70%ofALL`$`30%right`$Est, 
                             Est_PCGII_30perct=PCGII_pct$`30%ofALL`$`100%right`$Est, 
                             Est_PCGII_30perct.70=PCGII_pct$`30%ofALL`$`70%right`$Est, 
                             Est_PCGII_30perct.30=PCGII_pct$`30%ofALL`$`30%right`$Est, 

                             Est_FGGM=FGGM.test$pcor, 
                             shrunk_p, # Shrunk_pcor
                             
                             tscore_cLevel, 
                             tscore_PCGII_all=PCGII_all$`100%`$tscore, 
                             tscore_PCGII_all.70=PCGII_all$`70%`$tscore, 
                             tscore_PCGII_all.30=PCGII_all$`30%`$tscore, 
                             tscore_PCGII_70perct=PCGII_pct$`70%ofALL`$`100%right`$tscore, 
                             tscore_PCGII_70perct.70=PCGII_pct$`70%ofALL`$`70%right`$tscore,
                             tscore_PCGII_70perct.30=PCGII_pct$`70%ofALL`$`30%right`$tscore,
                             tscore_PCGII_30perct=PCGII_pct$`30%ofALL`$`100%right`$tscore, 
                             tscore_PCGII_30perct.70=PCGII_pct$`30%ofALL`$`70%right`$tscore,
                             tscore_PCGII_30perct.30=PCGII_pct$`30%ofALL`$`30%right`$tscore,
                             
                             FGGM_pval=FGGM.test$pval,
                             FGGM_qval=FGGM.test$qval,
                             
                             ENF_p=ENF.test$pval, ENF_q=p.adjust(ENF.test$pval, method="BH"),
                             ShrunkMLE_p=pval.shrunk, ShrunkMLE_q=p.adjust(pval.shrunk, method="BH"))
    
    
    ## Ranking
    FGGM.ordered<-results[order(results$FGGM_qval, results$FGGM_pval, decreasing = F),c("truth","FGGM_qval","FGGM_pval")]
    
    
    cLevel.ordered<-results[order(abs(results$tscore_cLevel), decreasing = T),c("truth","Est_cLevel","tscore_cLevel")]
    
    PCGII_all.ordered<-results[order(abs(results$tscore_PCGII_all), decreasing = T),c("truth","Est_PCGII_all","tscore_PCGII_all")]
    PCGII_all.70.ordered<-results[order(abs(results$tscore_PCGII_all.70), decreasing = T),c("truth","Est_PCGII_all.70","tscore_PCGII_all.70")]
    PCGII_all.30.ordered<-results[order(abs(results$tscore_PCGII_all.30), decreasing = T),c("truth","Est_PCGII_all.30","tscore_PCGII_all.30")]

    
    PCGII_30perct.ordered<-results[order(abs(results$tscore_PCGII_30perct), decreasing = T),c("truth","Est_PCGII_30perct","tscore_PCGII_30perct")]
    PCGII_30perct.30.ordered<-results[order(abs(results$tscore_PCGII_30perct.30), decreasing = T),c("truth","Est_PCGII_30perct.30","tscore_PCGII_30perct.30")]
    PCGII_30perct.70.ordered<-results[order(abs(results$tscore_PCGII_30perct.70), decreasing = T),c("truth","Est_PCGII_30perct.70","tscore_PCGII_30perct.70")]
    
    PCGII_70perct.ordered<-results[order(abs(results$tscore_PCGII_70perct), decreasing = T),c("truth","Est_PCGII_70perct","tscore_PCGII_70perct")]
    PCGII_70perct.30.ordered<-results[order(abs(results$tscore_PCGII_70perct.30), decreasing = T),c("truth","Est_PCGII_70perct.30","tscore_PCGII_70perct.30")]
    PCGII_70perct.70.ordered<-results[order(abs(results$tscore_PCGII_70perct.70), decreasing = T),c("truth","Est_PCGII_70perct.70","tscore_PCGII_70perct.70")]
    
    
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
      
      cLevel.ordered.TPR[loop, k]=sum(cLevel.ordered[1:loop,]$truth!=0)/CP 
      cLevel.ordered.FPR[loop, k]=sum(cLevel.ordered[1:loop,]$truth==0)/CN 
      cLevel.ordered.PPV[loop, k]=sum(cLevel.ordered[1:loop,]$truth!=0)/loop
      cLevel.ordered.FDR[loop, k]=sum(cLevel.ordered[1:loop,]$truth==0)/loop
      
      PCGII_all.ordered.TPR[loop, k]=sum(PCGII_all.ordered[1:loop,]$truth!=0)/CP 
      PCGII_all.ordered.FPR[loop, k]=sum(PCGII_all.ordered[1:loop,]$truth==0)/CN 
      PCGII_all.ordered.PPV[loop, k]=sum(PCGII_all.ordered[1:loop,]$truth!=0)/loop
      PCGII_all.ordered.FDR[loop, k]=sum(PCGII_all.ordered[1:loop,]$truth==0)/loop
      
      
      PCGII_all.30.ordered.TPR[loop, k]=sum(PCGII_all.30.ordered[1:loop,]$truth!=0)/CP 
      PCGII_all.30.ordered.FPR[loop, k]=sum(PCGII_all.30.ordered[1:loop,]$truth==0)/CN 
      PCGII_all.30.ordered.PPV[loop, k]=sum(PCGII_all.30.ordered[1:loop,]$truth!=0)/loop
      PCGII_all.30.ordered.FDR[loop, k]=sum(PCGII_all.30.ordered[1:loop,]$truth==0)/loop
      
      PCGII_all.70.ordered.TPR[loop, k]=sum(PCGII_all.70.ordered[1:loop,]$truth!=0)/CP 
      PCGII_all.70.ordered.FPR[loop, k]=sum(PCGII_all.70.ordered[1:loop,]$truth==0)/CN 
      PCGII_all.70.ordered.PPV[loop, k]=sum(PCGII_all.70.ordered[1:loop,]$truth!=0)/loop
      PCGII_all.70.ordered.FDR[loop, k]=sum(PCGII_all.70.ordered[1:loop,]$truth==0)/loop
      
      PCGII_30perct.ordered.TPR[loop, k]=sum(PCGII_30perct.ordered[1:loop,]$truth!=0)/CP 
      PCGII_30perct.ordered.FPR[loop, k]=sum(PCGII_30perct.ordered[1:loop,]$truth==0)/CN 
      PCGII_30perct.ordered.PPV[loop, k]=sum(PCGII_30perct.ordered[1:loop,]$truth!=0)/loop
      PCGII_30perct.ordered.FDR[loop, k]=sum(PCGII_30perct.ordered[1:loop,]$truth==0)/loop
      
      PCGII_30perct.30.ordered.TPR[loop, k]=sum(PCGII_30perct.30.ordered[1:loop,]$truth!=0)/CP 
      PCGII_30perct.30.ordered.FPR[loop, k]=sum(PCGII_30perct.30.ordered[1:loop,]$truth==0)/CN 
      PCGII_30perct.30.ordered.PPV[loop, k]=sum(PCGII_30perct.30.ordered[1:loop,]$truth!=0)/loop
      PCGII_30perct.30.ordered.FDR[loop, k]=sum(PCGII_30perct.30.ordered[1:loop,]$truth==0)/loop
      
      PCGII_30perct.70.ordered.TPR[loop, k]=sum(PCGII_30perct.70.ordered[1:loop,]$truth!=0)/CP 
      PCGII_30perct.70.ordered.FPR[loop, k]=sum(PCGII_30perct.70.ordered[1:loop,]$truth==0)/CN 
      PCGII_30perct.70.ordered.PPV[loop, k]=sum(PCGII_30perct.70.ordered[1:loop,]$truth!=0)/loop
      PCGII_30perct.70.ordered.FDR[loop, k]=sum(PCGII_30perct.70.ordered[1:loop,]$truth==0)/loop
      
      PCGII_70perct.ordered.TPR[loop, k]=sum(PCGII_70perct.ordered[1:loop,]$truth!=0)/CP 
      PCGII_70perct.ordered.FPR[loop, k]=sum(PCGII_70perct.ordered[1:loop,]$truth==0)/CN 
      PCGII_70perct.ordered.PPV[loop, k]=sum(PCGII_70perct.ordered[1:loop,]$truth!=0)/loop
      PCGII_70perct.ordered.FDR[loop, k]=sum(PCGII_70perct.ordered[1:loop,]$truth==0)/loop
      
      PCGII_70perct.30.ordered.TPR[loop, k]=sum(PCGII_70perct.30.ordered[1:loop,]$truth!=0)/CP 
      PCGII_70perct.30.ordered.FPR[loop, k]=sum(PCGII_70perct.30.ordered[1:loop,]$truth==0)/CN 
      PCGII_70perct.30.ordered.PPV[loop, k]=sum(PCGII_70perct.30.ordered[1:loop,]$truth!=0)/loop
      PCGII_70perct.30.ordered.FDR[loop, k]=sum(PCGII_70perct.30.ordered[1:loop,]$truth==0)/loop
      
      PCGII_70perct.70.ordered.TPR[loop, k]=sum(PCGII_70perct.70.ordered[1:loop,]$truth!=0)/CP 
      PCGII_70perct.70.ordered.FPR[loop, k]=sum(PCGII_70perct.70.ordered[1:loop,]$truth==0)/CN 
      PCGII_70perct.70.ordered.PPV[loop, k]=sum(PCGII_70perct.70.ordered[1:loop,]$truth!=0)/loop
      PCGII_70perct.70.ordered.FDR[loop, k]=sum(PCGII_70perct.70.ordered[1:loop,]$truth==0)/loop
      
      
      FGGM.ordered.TPR[loop, k]=sum(FGGM.ordered[1:loop,]$truth!=0)/CP 
      FGGM.ordered.FPR[loop, k]=sum(FGGM.ordered[1:loop,]$truth==0)/CN 
      FGGM.ordered.PPV[loop, k]=sum(FGGM.ordered[1:loop,]$truth!=0)/loop
      FGGM.ordered.FDR[loop, k]=sum(FGGM.ordered[1:loop,]$truth==0)/loop
    }
    
    for (c in 1:1000){
      tempstore.ENF[c,k]=max(ENF.ordered.TPR[ENF.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.MLE[c,k]=max(MLE.ordered.TPR[MLE.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.cLevel[c,k]=max(cLevel.ordered.TPR[cLevel.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR

      tempstore.PCGII_all[c,k]=max(PCGII_all.ordered.TPR[PCGII_all.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_all.30[c,k]=max(PCGII_all.30.ordered.TPR[PCGII_all.30.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_all.70[c,k]=max(PCGII_all.70.ordered.TPR[PCGII_all.70.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR

      tempstore.PCGII_30perct[c,k]=max(PCGII_30perct.ordered.TPR[PCGII_30perct.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_30perct.30[c,k]=max(PCGII_30perct.30.ordered.TPR[PCGII_30perct.30.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_30perct.70[c,k]=max(PCGII_30perct.70.ordered.TPR[PCGII_30perct.70.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_70perct[c,k]=max(PCGII_70perct.ordered.TPR[PCGII_70perct.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_70perct.30[c,k]=max(PCGII_70perct.30.ordered.TPR[PCGII_70perct.30.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_70perct.70[c,k]=max(PCGII_70perct.70.ordered.TPR[PCGII_70perct.70.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      
      tempstore.FGGM[c,k]=max(FGGM.ordered.TPR[FGGM.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
    }
    print("ROC done")
    
    #### FDR 
    print("FDR starts")
    for (a in 1:length(al)){
      temp=inference(cLevel, alpha = al[a])$sigs # c=0
      sigs_cLevel=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_all$`100%`$out, alpha = al[a])$sigs # c=0
      sigs_PCGII_all=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_all$`30%`$out, alpha = al[a])$sigs # c=0
      sigs_PCGII_all.30=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_all$`70%`$out, alpha = al[a])$sigs # c=0
      sigs_PCGII_all.70=sigs2vec(temp, p) # significant edges
    
      
      temp=inference(PCGII_pct$`30%ofALL`$`100%right`$out, alpha = al[a])$sigs # c=0
      sigs_PCGII_30perct=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_pct$`30%ofALL`$`30%right`$out, alpha = al[a])$sigs # c=0
      sigs_PCGII_30perct.30=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_pct$`30%ofALL`$`70%right`$out, alpha = al[a])$sigs # c=0
      sigs_PCGII_30perct.70=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_pct$`70%ofALL`$`100%right`$out, alpha = al[a])$sigs # c=0
      sigs_PCGII_70perct=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_pct$`70%ofALL`$`30%right`$out, alpha = al[a])$sigs # c=0
      sigs_PCGII_70perct.30=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_pct$`70%ofALL`$`70%right`$out, alpha = al[a])$sigs # c=0
      sigs_PCGII_70perct.70=sigs2vec(temp, p) # significant edges
    
      # Number of total selected edges
      all.enf.p[a,k]=sum(results$ENF_p<=al[a])
      all.enf.q[a,k]=sum(results$ENF_q<=al[a]) 
      all.mle.p[a,k]=sum(results$ShrunkMLE_p<=al[a]) 
      all.mle.q[a,k]=sum(results$ShrunkMLE_q<=al[a])
      all.cLevel[a,k]=sum(sigs_cLevel==1)
      
      all.PCGII_all[a,k]=sum(sigs_PCGII_all==1)
      all.PCGII_all.30[a,k]=sum(sigs_PCGII_all.30==1)
      all.PCGII_all.70[a,k]=sum(sigs_PCGII_all.70==1)
      
      all.PCGII_30perct[a,k]=sum(sigs_PCGII_30perct==1)
      all.PCGII_30perct.30[a,k]=sum(sigs_PCGII_30perct.30==1)
      all.PCGII_30perct.70[a,k]=sum(sigs_PCGII_30perct.70==1)
      
      all.PCGII_70perct[a,k]=sum(sigs_PCGII_70perct==1)
      all.PCGII_70perct.30[a,k]=sum(sigs_PCGII_70perct.30==1)
      all.PCGII_70perct.70[a,k]=sum(sigs_PCGII_70perct.70==1)
      
      all.FGGM[a,k]=sum(results$FGGM_qval<=al[a])
      
      # true positives
      tp.enf.p[a,k]=sum(which(results$ENF_p<=al[a]) %in% TP)
      tp.enf.q[a,k]=sum(which(results$ENF_q<=al[a]) %in% TP)
      tp.mle.p[a,k]=sum(which(results$ShrunkMLE_p<=al[a]) %in% TP)
      tp.mle.q[a,k]=sum(which(results$ShrunkMLE_q<=al[a]) %in% TP)
      tp.cLevel[a,k]=sum(which(sigs_cLevel==1) %in% TP)
      
      tp.PCGII_all[a,k]=sum(which(sigs_PCGII_all==1) %in% TP)
      tp.PCGII_all.30[a,k]=sum(which(sigs_PCGII_all.30==1) %in% TP)
      tp.PCGII_all.70[a,k]=sum(which(sigs_PCGII_all.70==1) %in% TP)
      
      tp.PCGII_30perct[a,k]=sum(which(sigs_PCGII_30perct==1) %in% TP)
      tp.PCGII_30perct.30[a,k]=sum(which(sigs_PCGII_30perct.30==1) %in% TP)
      tp.PCGII_30perct.70[a,k]=sum(which(sigs_PCGII_30perct.70==1) %in% TP)
      tp.PCGII_70perct[a,k]=sum(which(sigs_PCGII_70perct==1) %in% TP)
      tp.PCGII_70perct.30[a,k]=sum(which(sigs_PCGII_70perct.30==1) %in% TP)
      tp.PCGII_70perct.70[a,k]=sum(which(sigs_PCGII_70perct.70==1) %in% TP)
      
      tp.FGGM[a,k]=sum(which(results$FGGM_qval<=al[a]) %in% TP)
      
      # false positives
      fp.enf.p[a,k]=sum(!which(results$ENF_p<=al[a]) %in% TP)
      fp.enf.q[a,k]=sum(!which(results$ENF_q<=al[a]) %in% TP)
      fp.mle.p[a,k]=sum(!which(results$ShrunkMLE_p<=al[a]) %in% TP)
      fp.mle.q[a,k]=sum(!which(results$ShrunkMLE_q<=al[a]) %in% TP)
      fp.cLevel[a,k]=sum(!which(sigs_cLevel==1) %in% TP)

      fp.PCGII_all[a,k]=sum(!which(sigs_PCGII_all==1) %in% TP)
      fp.PCGII_all.30[a,k]=sum(!which(sigs_PCGII_all.30==1) %in% TP)
      fp.PCGII_all.70[a,k]=sum(!which(sigs_PCGII_all.70==1) %in% TP)

      fp.PCGII_30perct[a,k]=sum(!which(sigs_PCGII_30perct==1) %in% TP)
      fp.PCGII_30perct.30[a,k]=sum(!which(sigs_PCGII_30perct.30==1) %in% TP)
      fp.PCGII_30perct.70[a,k]=sum(!which(sigs_PCGII_30perct.70==1) %in% TP)
      fp.PCGII_70perct[a,k]=sum(!which(sigs_PCGII_70perct==1) %in% TP)
      fp.PCGII_70perct.30[a,k]=sum(!which(sigs_PCGII_70perct.30==1) %in% TP)
      fp.PCGII_70perct.70[a,k]=sum(!which(sigs_PCGII_70perct.70==1) %in% TP)

      fp.FGGM[a,k]=sum(!which(results$FGGM_qval<=al[a]) %in% TP)
      
      # true negatives
      tn.enf.p[a,k]=sum(!which(results$ENF_p>al[a]) %in% TP)
      tn.enf.q[a,k]=sum(!which(results$ENF.qval>al[a]) %in% TP)
      tn.mle.p[a,k]=sum(!which(results$ShrunkMLE_p>al[a]) %in% TP)
      tn.mle.q[a,k]=sum(!which(results$ShrunkMLE_q>al[a]) %in% TP)
      tn.cLevel[a,k]=sum(!which(sigs_cLevel!=1) %in% TP)
      
      tn.PCGII_all[a,k]=sum(!which(sigs_PCGII_all!=1) %in% TP)
      tn.PCGII_all.30[a,k]=sum(!which(sigs_PCGII_all.30!=1) %in% TP)
      tn.PCGII_all.70[a,k]=sum(!which(sigs_PCGII_all.70!=1) %in% TP)

      tn.PCGII_30perct[a,k]=sum(!which(sigs_PCGII_30perct!=1) %in% TP)
      tn.PCGII_30perct.30[a,k]=sum(!which(sigs_PCGII_30perct.30!=1) %in% TP)
      tn.PCGII_30perct.70[a,k]=sum(!which(sigs_PCGII_30perct.70!=1) %in% TP)
      tn.PCGII_70perct[a,k]=sum(!which(sigs_PCGII_70perct!=1) %in% TP)
      tn.PCGII_70perct.30[a,k]=sum(!which(sigs_PCGII_70perct.30!=1) %in% TP)
      tn.PCGII_70perct.70[a,k]=sum(!which(sigs_PCGII_70perct.70!=1) %in% TP)
      
      tn.FGGM[a,k]=sum(!which(results$FGGM_qval>al[a]) %in% TP)
      
      # false negatives
      fn.enf.p[a,k]=sum(which(results$ENF_p>al[a]) %in% TP)
      fn.enf.q[a,k]=sum(which(results$ENF_q>al[a]) %in% TP)
      fn.mle.p[a,k]=sum(which(results$ShrunkMLE_p>al[a]) %in% TP)
      fn.mle.q[a,k]=sum(which(results$ShrunkMLE_q>al[a]) %in% TP)
      fn.cLevel[a,k]=sum(which(sigs_cLevel!=1) %in% TP)
      
      fn.PCGII_all[a,k]=sum(which(sigs_PCGII_all!=1) %in% TP)
      fn.PCGII_all.30[a,k]=sum(which(sigs_PCGII_all.30!=1) %in% TP)
      fn.PCGII_all.70[a,k]=sum(which(sigs_PCGII_all.70!=1) %in% TP)

      fn.PCGII_30perct[a,k]=sum(which(sigs_PCGII_30perct!=1) %in% TP)
      fn.PCGII_30perct.30[a,k]=sum(which(sigs_PCGII_30perct.30!=1) %in% TP)
      fn.PCGII_30perct.70[a,k]=sum(which(sigs_PCGII_30perct.70!=1) %in% TP)
      fn.PCGII_70perct[a,k]=sum(which(sigs_PCGII_70perct!=1) %in% TP)
      fn.PCGII_70perct.30[a,k]=sum(which(sigs_PCGII_70perct.30!=1) %in% TP)
      fn.PCGII_70perct.70[a,k]=sum(which(sigs_PCGII_70perct.70!=1) %in% TP)
      
      fn.FGGM[a,k]=sum(which(results$FGGM_qval>al[a]) %in% TP)
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
    rank.cLevel=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                      TPR=rowMeans(cLevel.ordered.TPR),
                                      FPR=rowMeans(cLevel.ordered.FPR),
                                      PPV=rowMeans(cLevel.ordered.PPV),
                                      FDR=rowMeans(cLevel.ordered.FDR))
    rank.PCGII_all=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                     TPR=rowMeans(PCGII_all.ordered.TPR),
                                     FPR=rowMeans(PCGII_all.ordered.FPR),
                                     PPV=rowMeans(PCGII_all.ordered.PPV),
                                     FDR=rowMeans(PCGII_all.ordered.FDR))
    rank.PCGII_all.30=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                       TPR=rowMeans(PCGII_all.30.ordered.TPR),
                                       FPR=rowMeans(PCGII_all.30.ordered.FPR),
                                       PPV=rowMeans(PCGII_all.30.ordered.PPV),
                                       FDR=rowMeans(PCGII_all.30.ordered.FDR))
    rank.PCGII_all.70=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_all.70.ordered.TPR),
                                            FPR=rowMeans(PCGII_all.70.ordered.FPR),
                                            PPV=rowMeans(PCGII_all.70.ordered.PPV),
                                            FDR=rowMeans(PCGII_all.70.ordered.FDR))

    rank.PCGII_30perct=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_30perct.ordered.TPR),
                                            FPR=rowMeans(PCGII_30perct.ordered.FPR),
                                            PPV=rowMeans(PCGII_30perct.ordered.PPV),
                                            FDR=rowMeans(PCGII_30perct.ordered.FDR))
    rank.PCGII_30perct.30=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_30perct.30.ordered.TPR),
                                            FPR=rowMeans(PCGII_30perct.30.ordered.FPR),
                                            PPV=rowMeans(PCGII_30perct.30.ordered.PPV),
                                            FDR=rowMeans(PCGII_30perct.30.ordered.FDR))
    rank.PCGII_30perct.70=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_30perct.70.ordered.TPR),
                                            FPR=rowMeans(PCGII_30perct.70.ordered.FPR),
                                            PPV=rowMeans(PCGII_30perct.70.ordered.PPV),
                                            FDR=rowMeans(PCGII_30perct.70.ordered.FDR))
    rank.PCGII_70perct=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_70perct.ordered.TPR),
                                            FPR=rowMeans(PCGII_70perct.ordered.FPR),
                                            PPV=rowMeans(PCGII_70perct.ordered.PPV),
                                            FDR=rowMeans(PCGII_70perct.ordered.FDR))
    rank.PCGII_70perct.30=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_70perct.30.ordered.TPR),
                                            FPR=rowMeans(PCGII_70perct.30.ordered.FPR),
                                            PPV=rowMeans(PCGII_70perct.30.ordered.PPV),
                                            FDR=rowMeans(PCGII_70perct.30.ordered.FDR))
    rank.PCGII_70perct.70=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_70perct.70.ordered.TPR),
                                            FPR=rowMeans(PCGII_70perct.70.ordered.FPR),
                                            PPV=rowMeans(PCGII_70perct.70.ordered.PPV),
                                            FDR=rowMeans(PCGII_70perct.70.ordered.FDR))

    rank.FGGM=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                TPR=rowMeans(FGGM.ordered.TPR),
                                FPR=rowMeans(FGGM.ordered.FPR),
                                PPV=rowMeans(FGGM.ordered.PPV),
                                FDR=rowMeans(FGGM.ordered.FDR))
    
    rank.table=cbind.data.frame(methods=c(rep("ENF",p*(p-1)/2),
                                          rep("MLE",p*(p-1)/2),
                                          rep("cLevel",p*(p-1)/2),
                                          rep("PCGII_all",p*(p-1)/2),
                                          rep("PCGII_all.70",p*(p-1)/2),
                                          rep("PCGII_all.30",p*(p-1)/2),
                                          rep("PCGII_70perct",p*(p-1)/2),
                                          rep("PCGII_70perct.70",p*(p-1)/2),
                                          rep("PCGII_70perct.30",p*(p-1)/2),
                                          rep("PCGII_30perct",p*(p-1)/2),
                                          rep("PCGII_30perct.70",p*(p-1)/2),
                                          rep("PCGII_30perct.30",p*(p-1)/2),
                                          rep("FGGM",p*(p-1)/2)),
                                rbind.data.frame(rank.ENF, rank.MLE, 
                                                 rank.cLevel,
                                                 rank.PCGII_all, 
                                                 rank.PCGII_all.70, 
                                                 rank.PCGII_all.30, 
                                                 rank.PCGII_70perct, rank.PCGII_70perct.70, rank.PCGII_70perct.30,
                                                 rank.PCGII_30perct, rank.PCGII_30perct.70, rank.PCGII_30perct.30,
                                                 rank.FGGM))
    
    
    ROC.out=cbind.data.frame(FPR=seq(0.001,1,0.001), # Sample FPR cut points
                             ENF=rowMeans(tempstore.ENF), # averaged max TPR at cutted FPR level
                             MLE=rowMeans(tempstore.MLE),
                             cLevel=rowMeans(tempstore.cLevel),
                             PCGII_all=rowMeans(tempstore.PCGII_all),
                             PCGII_all.70=rowMeans(tempstore.PCGII_all.70),
                             PCGII_all.30=rowMeans(tempstore.PCGII_all.30),
                             PCGII_70perct=rowMeans(tempstore.PCGII_70perct),
                             PCGII_70perct.70=rowMeans(tempstore.PCGII_70perct.70),
                             PCGII_70perct.30=rowMeans(tempstore.PCGII_70perct.30),
                             PCGII_30perct=rowMeans(tempstore.PCGII_30perct),
                             PCGII_30perct.70=rowMeans(tempstore.PCGII_30perct.70),
                             PCGII_30perct.30=rowMeans(tempstore.PCGII_30perct.30),

                             FGGM=rowMeans(tempstore.FGGM))
    #ROC.table=tidyr::gather(ROC.out, methods, TPR, ENF:cPCG_df)
    
    
    # Number of total correctly selected edges
    All.enf.p=rowMeans(all.enf.p)
    All.enf.q=rowMeans(all.enf.q)
    All.mle.p=rowMeans(all.mle.p)
    All.mle.q=rowMeans(all.mle.q)
    All.cLevel=rowMeans(all.cLevel)
    
    All.PCGII_all=rowMeans(all.PCGII_all)
    All.PCGII_all.30=rowMeans(all.PCGII_all.30)
    All.PCGII_all.70=rowMeans(all.PCGII_all.70)
    
    All.PCGII_30perct=rowMeans(all.PCGII_30perct)
    All.PCGII_30perct.30=rowMeans(all.PCGII_30perct.30)
    All.PCGII_30perct.70=rowMeans(all.PCGII_30perct.70)
    All.PCGII_70perct=rowMeans(all.PCGII_70perct)
    All.PCGII_70perct.30=rowMeans(all.PCGII_70perct.30)
    All.PCGII_70perct.70=rowMeans(all.PCGII_70perct.70)
    
    All.FGGM=rowMeans(all.FGGM)
    
    # true positives
    TP.enf.p=rowMeans(tp.enf.p)
    TP.enf.q=rowMeans(tp.enf.q)
    TP.mle.p=rowMeans(tp.mle.p)
    TP.mle.q=rowMeans(tp.mle.q)
    TP.cLevel=rowMeans(tp.cLevel)
    
    TP.PCGII_all=rowMeans(tp.PCGII_all)
    TP.PCGII_all.30=rowMeans(tp.PCGII_all.30)
    TP.PCGII_all.70=rowMeans(tp.PCGII_all.70)

    TP.PCGII_30perct=rowMeans(tp.PCGII_30perct)
    TP.PCGII_30perct.30=rowMeans(tp.PCGII_30perct.30)
    TP.PCGII_30perct.70=rowMeans(tp.PCGII_30perct.70)
    TP.PCGII_70perct=rowMeans(tp.PCGII_70perct)
    TP.PCGII_70perct.30=rowMeans(tp.PCGII_70perct.30)
    TP.PCGII_70perct.70=rowMeans(tp.PCGII_70perct.70)
    
    TP.FGGM=rowMeans(tp.FGGM)
    
    # false positives
    FP.enf.p=rowMeans(fp.enf.p)
    FP.enf.q=rowMeans(fp.enf.q)
    FP.mle.p=rowMeans(fp.mle.p)
    FP.mle.q=rowMeans(fp.mle.q)
    FP.cLevel=rowMeans(fp.cLevel)
    
    FP.PCGII_all=rowMeans(fp.PCGII_all)
    FP.PCGII_all.30=rowMeans(fp.PCGII_all.30)
    FP.PCGII_all.70=rowMeans(fp.PCGII_all.70)

    FP.PCGII_30perct=rowMeans(fp.PCGII_30perct)
    FP.PCGII_30perct.30=rowMeans(fp.PCGII_30perct.30)
    FP.PCGII_30perct.70=rowMeans(fp.PCGII_30perct.70)
    FP.PCGII_70perct=rowMeans(fp.PCGII_70perct)
    FP.PCGII_70perct.30=rowMeans(fp.PCGII_70perct.30)
    FP.PCGII_70perct.70=rowMeans(fp.PCGII_70perct.70)
    
    FP.FGGM=rowMeans(fp.FGGM)
    
    # true negatives
    TN.enf.p=rowMeans(tn.enf.p)
    TN.enf.q=rowMeans(tn.enf.q)
    TN.mle.p=rowMeans(tn.mle.p)
    TN.mle.q=rowMeans(tn.mle.q)
    TN.cLevel=rowMeans(tn.cLevel)
    
    TN.PCGII_all=rowMeans(tn.PCGII_all)
    TN.PCGII_all.30=rowMeans(tn.PCGII_all.30)
    TN.PCGII_all.70=rowMeans(tn.PCGII_all.70)

    TN.PCGII_30perct=rowMeans(tn.PCGII_30perct)
    TN.PCGII_30perct.30=rowMeans(tn.PCGII_30perct.30)
    TN.PCGII_30perct.70=rowMeans(tn.PCGII_30perct.70)
    TN.PCGII_70perct=rowMeans(tn.PCGII_70perct)
    TN.PCGII_70perct.30=rowMeans(tn.PCGII_70perct.30)
    TN.PCGII_70perct.70=rowMeans(tn.PCGII_70perct.70)

    TN.FGGM=rowMeans(tn.FGGM)
    
    # false negatives
    FN.enf.p=rowMeans(fn.enf.p)
    FN.enf.q=rowMeans(fn.enf.q)
    FN.mle.p=rowMeans(fn.mle.p)
    FN.mle.q=rowMeans(fn.mle.q)
    FN.cLevel=rowMeans(fn.cLevel)
    
    FN.PCGII_all=rowMeans(fn.PCGII_all)
    FN.PCGII_all.30=rowMeans(fn.PCGII_all.30)
    FN.PCGII_all.70=rowMeans(fn.PCGII_all.70)

    FN.PCGII_30perct=rowMeans(fn.PCGII_30perct)
    FN.PCGII_30perct.30=rowMeans(fn.PCGII_30perct.30)
    FN.PCGII_30perct.70=rowMeans(fn.PCGII_30perct.70)
    FN.PCGII_70perct=rowMeans(fn.PCGII_70perct)
    FN.PCGII_70perct.30=rowMeans(fn.PCGII_70perct.30)
    FN.PCGII_70perct.70=rowMeans(fn.PCGII_70perct.70)
    
    FN.FGGM=rowMeans(fn.FGGM)
    
    # empirical fdr
    fdr.enf.p=FP.enf.p/All.enf.p
    fdr.enf.p[is.na(fdr.enf.p)]=0
    fdr.enf.q=FP.enf.q/All.enf.q
    fdr.enf.q[is.na(fdr.enf.q)]=0
    fdr.mle.p=FP.mle.p/All.mle.p
    fdr.mle.p[is.na(fdr.mle.p)]=0
    fdr.mle.q=FP.mle.q/All.mle.q
    fdr.mle.q[is.na(fdr.mle.q)]=0
    fdr.cLevel=FP.cLevel/All.cLevel
    fdr.cLevel[is.na(fdr.cLevel)]=0

    
    fdr.PCGII_all=FP.PCGII_all/All.PCGII_all
    fdr.PCGII_all[is.na(fdr.PCGII_all)]=0
    fdr.PCGII_all.30=FP.PCGII_all.30/All.PCGII_all.30
    fdr.PCGII_all.30[is.na(fdr.PCGII_all.30)]=0
    fdr.PCGII_all.70=FP.PCGII_all.70/All.PCGII_all.70
    fdr.PCGII_all.70[is.na(fdr.PCGII_all.70)]=0

    fdr.PCGII_30perct=FP.PCGII_30perct/All.PCGII_30perct
    fdr.PCGII_30perct[is.na(fdr.PCGII_30perct)]=0
    fdr.PCGII_30perct.30=FP.PCGII_30perct.30/All.PCGII_30perct.30
    fdr.PCGII_30perct.30[is.na(fdr.PCGII_30perct.30)]=0
    fdr.PCGII_30perct.70=FP.PCGII_30perct.70/All.PCGII_30perct.70
    fdr.PCGII_30perct.70[is.na(fdr.PCGII_30perct.70)]=0
    
    fdr.PCGII_70perct=FP.PCGII_70perct/All.PCGII_70perct
    fdr.PCGII_70perct[is.na(fdr.PCGII_70perct)]=0
    fdr.PCGII_70perct.30=FP.PCGII_70perct.30/All.PCGII_70perct.30
    fdr.PCGII_70perct.30[is.na(fdr.PCGII_70perct.30)]=0
    fdr.PCGII_70perct.70=FP.PCGII_70perct.70/All.PCGII_70perct.70
    fdr.PCGII_70perct.70[is.na(fdr.PCGII_70perct.70)]=0
    
    fdr.FGGM=FP.FGGM/All.FGGM
    fdr.FGGM[is.na(fdr.FGGM)]=0
    
    
    fdr.table=cbind.data.frame(
      FDR=rep(al,15), # nominal FDR
      methods=c(rep("ENF.p",length(al)),
                rep("ENF.q",length(al)),
                rep("MLE.p",length(al)),
                rep("MLE.q",length(al)),
                rep("cLevel",length(al)), #5
                
                rep("PCGII_all",length(al)),
                rep("PCGII_all.70",length(al)),
                rep("PCGII_all.30",length(al)),
                rep("PCGII_70perct",length(al)),
                rep("PCGII_70perct.30",length(al)),
                rep("PCGII_70perct.70",length(al)), 
                rep("PCGII_30perct",length(al)),
                rep("PCGII_30perct.30",length(al)),
                rep("PCGII_30perct.70",length(al)),# 3*3 
                
                rep("FGGM",length(al))),
      fdr=c(fdr.enf.p, fdr.enf.q, fdr.mle.p, fdr.mle.q, 
            fdr.cLevel, 
            
            fdr.PCGII_all, fdr.PCGII_all.70, fdr.PCGII_all.30,  
            fdr.PCGII_70perct, fdr.PCGII_70perct.70,  fdr.PCGII_70perct.30, 
            fdr.PCGII_30perct, fdr.PCGII_30perct.70,  fdr.PCGII_30perct.30, 
            
            fdr.FGGM))
  } # summary
  print("Summary done")
  {
    return(
      list(
        ROC=list(
          # ROC
          ENF.ordered.TPR=ENF.ordered.TPR, ENF.ordered.FPR=ENF.ordered.FPR, 
          ENF.ordered.PPV=ENF.ordered.PPV, ENF.ordered.FDR=ENF.ordered.FDR,
          MLE.ordered.TPR=MLE.ordered.TPR, MLE.ordered.FPR=MLE.ordered.FPR, 
          MLE.ordered.PPV=MLE.ordered.PPV, MLE.ordered.FDR=MLE.ordered.FDR,
          
          cLevel.ordered.TPR=cLevel.ordered.TPR, cLevel.ordered.FPR=cLevel.ordered.FPR, 
          cLevel.ordered.PPV=cLevel.ordered.PPV, cLevel.ordered.FDR=cLevel.ordered.FDR,
          
          PCGII_all.ordered.TPR=PCGII_all.ordered.TPR, PCGII_all.ordered.FPR=PCGII_all.ordered.FPR, 
          PCGII_all.ordered.PPV=PCGII_all.ordered.PPV, PCGII_all.ordered.FDR=PCGII_all.ordered.FDR,
          PCGII_all.30.ordered.TPR=PCGII_all.30.ordered.TPR, PCGII_all.30.ordered.FPR=PCGII_all.30.ordered.FPR, 
          PCGII_all.30.ordered.PPV=PCGII_all.30.ordered.PPV, PCGII_all.30.ordered.FDR=PCGII_all.30.ordered.FDR,
          PCGII_all.70.ordered.TPR=PCGII_all.70.ordered.TPR, PCGII_all.70.ordered.FPR=PCGII_all.70.ordered.FPR, 
          PCGII_all.70.ordered.PPV=PCGII_all.70.ordered.PPV, PCGII_all.70.ordered.FDR=PCGII_all.70.ordered.FDR,

          PCGII_30perct.ordered.TPR=PCGII_30perct.ordered.TPR, PCGII_30perct.ordered.FPR=PCGII_30perct.ordered.FPR, 
          PCGII_30perct.ordered.PPV=PCGII_30perct.ordered.PPV, PCGII_30perct.ordered.FDR=PCGII_30perct.ordered.FDR,
          PCGII_30perct.30.ordered.TPR=PCGII_30perct.30.ordered.TPR, PCGII_30perct.30.ordered.FPR=PCGII_30perct.30.ordered.FPR, 
          PCGII_30perct.30.ordered.PPV=PCGII_30perct.30.ordered.PPV, PCGII_30perct.30.ordered.FDR=PCGII_30perct.30.ordered.FDR,
          PCGII_30perct.70.ordered.TPR=PCGII_30perct.70.ordered.TPR, PCGII_30perct.70.ordered.FPR=PCGII_30perct.70.ordered.FPR, 
          PCGII_30perct.70.ordered.PPV=PCGII_30perct.70.ordered.PPV, PCGII_30perct.70.ordered.FDR=PCGII_30perct.70.ordered.FDR,
          
          PCGII_70perct.ordered.TPR=PCGII_70perct.ordered.TPR, PCGII_70perct.ordered.FPR=PCGII_70perct.ordered.FPR, 
          PCGII_70perct.ordered.PPV=PCGII_70perct.ordered.PPV, PCGII_70perct.ordered.FDR=PCGII_70perct.ordered.FDR,
          PCGII_70perct.30.ordered.TPR=PCGII_70perct.30.ordered.TPR, PCGII_70perct.30.ordered.FPR=PCGII_70perct.30.ordered.FPR, 
          PCGII_70perct.30.ordered.PPV=PCGII_70perct.30.ordered.PPV, PCGII_70perct.30.ordered.FDR=PCGII_70perct.30.ordered.FDR,
          PCGII_70perct.70.ordered.TPR=PCGII_70perct.70.ordered.TPR, PCGII_70perct.70.ordered.FPR=PCGII_70perct.70.ordered.FPR, 
          PCGII_70perct.70.ordered.PPV=PCGII_70perct.70.ordered.PPV, PCGII_70perct.70.ordered.FDR=PCGII_70perct.70.ordered.FDR,

          FGGM.ordered.TPR=FGGM.ordered.TPR, FGGM.ordered.FPR=FGGM.ordered.FPR, 
          FGGM.ordered.PPV=FGGM.ordered.PPV, FGGM.ordered.FDR=FGGM.ordered.FDR,
          rank.table=rank.table,
          
          # empirical TPR
          tempstore.ENF=tempstore.ENF, tempstore.MLE=tempstore.MLE, 
          tempstore.cLevel=tempstore.cLevel, 
          
          tempstore.PCGII_all=tempstore.PCGII_all, 
          tempstore.PCGII_all.30=tempstore.PCGII_all.30, 
          tempstore.PCGII_all.70=tempstore.PCGII_all.70, 

          tempstore.PCGII_30perct=tempstore.PCGII_30perct, 
          tempstore.PCGII_30perct.30=tempstore.PCGII_30perct.30, 
          tempstore.PCGII_30perct.70=tempstore.PCGII_30perct.70, 
          tempstore.PCGII_70perct=tempstore.PCGII_70perct, 
          tempstore.PCGII_70perct.30=tempstore.PCGII_70perct.30, 
          tempstore.PCGII_70perct.70=tempstore.PCGII_70perct.70, 
          
          tempstore.FGGM=tempstore.FGGM,
          
          ROC.out=ROC.out),
        
        # FDR
        metrics=list( 
          # Number of total correctly selected edges
          all.enf.p=all.enf.p, all.enf.q=all.enf.q,  
          all.mle.p=all.mle.p,  all.mle.q=all.mle.q,
          all.cLevel=all.cLevel, 
          
          all.PCGII_all=all.PCGII_all, 
          all.PCGII_all.30=all.PCGII_all.30, 
          all.PCGII_all.70=all.PCGII_all.70,

          all.PCGII_30perct=all.PCGII_30perct, 
          all.PCGII_30perct.30=all.PCGII_30perct.30, 
          all.PCGII_30perct.70=all.PCGII_30perct.70,
          all.PCGII_70perct=all.PCGII_70perct, 
          all.PCGII_70perct.30=all.PCGII_70perct.30, 
          all.PCGII_70perct.70=all.PCGII_70perct.70,
          
          all.FGGM=all.FGGM,
          
          # true positives
          tp.enf.p=tp.enf.p, tp.enf.q=tp.enf.q, 
          tp.mle.p=tp.mle.p, tp.mle.q=tp.mle.q,
          tp.cLevel=tp.cLevel, 

          tp.PCGII_all=tp.PCGII_all, 
          tp.PCGII_all.30=tp.PCGII_all.30, 
          tp.PCGII_all.70=tp.PCGII_all.70,
          
          tp.PCGII_30perct=tp.PCGII_30perct, 
          tp.PCGII_30perct.30=tp.PCGII_30perct.30, 
          tp.PCGII_30perct.70=tp.PCGII_30perct.70,
          tp.PCGII_70perct=tp.PCGII_70perct, 
          tp.PCGII_70perct.30=tp.PCGII_70perct.30, 
          tp.PCGII_70perct.70=tp.PCGII_70perct.70,
          
          tp.FGGM=tp.FGGM,

          # false positives
          fp.enf.p=fp.enf.p, fp.enf.q=fp.enf.q, 
          fp.mle.p=fp.mle.p, fp.mle.q=fp.mle.q,
          fp.cLevel=fp.cLevel, 
          
          fp.PCGII_all=fp.PCGII_all, 
          fp.PCGII_all.30=fp.PCGII_all.30, 
          fp.PCGII_all.70=fp.PCGII_all.70,
          
          fp.PCGII_30perct=fp.PCGII_30perct, 
          fp.PCGII_30perct.30=fp.PCGII_30perct.30, 
          fp.PCGII_30perct.70=fp.PCGII_30perct.70,
          fp.PCGII_70perct=fp.PCGII_70perct, 
          fp.PCGII_70perct.30=fp.PCGII_70perct.30, 
          fp.PCGII_70perct.70=fp.PCGII_70perct.70,
          
          fp.FGGM=fp.FGGM,
          
          # false negatives
          tn.enf.p=tn.enf.p, tn.enf.q=tn.enf.q, 
          tn.mle.p=tn.mle.p, tn.mle.q=tn.mle.q,
          tn.cLevel=tn.cLevel, 
          
          tn.PCGII_all=tn.PCGII_all, 
          tn.PCGII_all.30=tn.PCGII_all.30, 
          tn.PCGII_all.70=tn.PCGII_all.70,

          tn.PCGII_30perct=tn.PCGII_30perct, 
          tn.PCGII_30perct.30=tn.PCGII_30perct.30, 
          tn.PCGII_30perct.70=tn.PCGII_30perct.70,
          tn.PCGII_70perct=tn.PCGII_70perct, 
          tn.PCGII_70perct.30=tn.PCGII_70perct.30, 
          tn.PCGII_70perct.70=tn.PCGII_70perct.70,
          
          tn.FGGM=tn.FGGM,
          
          # false negatives
          fn.enf.p=fn.enf.p, fn.enf.q=fn.enf.q,
          fn.mle.p=fn.mle.p, fn.mle.q=fn.mle.q,
          fn.cLevel=fn.cLevel, 
          
          fn.PCGII_all=fn.PCGII_all, 
          fn.PCGII_all.30=fn.PCGII_all.30, 
          fn.PCGII_all.70=fn.PCGII_all.70,

          fn.PCGII_30perct=fn.PCGII_30perct, 
          fn.PCGII_30perct.30=fn.PCGII_30perct.30, 
          fn.PCGII_30perct.70=fn.PCGII_30perct.70,
          fn.PCGII_70perct=fn.PCGII_70perct, 
          fn.PCGII_70perct.30=fn.PCGII_70perct.30, 
          fn.PCGII_70perct.70=fn.PCGII_70perct.70,

          fn.FGGM=fn.FGGM,
          
          All.enf.p=All.enf.p, All.enf.q=All.enf.q, 
          All.mle.p=All.mle.p, All.mle.q=All.mle.q,
          All.cLevel=All.cLevel, 
          
          All.PCGII_all=All.PCGII_all, 
          All.PCGII_all.30=All.PCGII_all.30, 
          All.PCGII_all.70=All.PCGII_all.70,

          All.PCGII_30perct=All.PCGII_30perct, 
          All.PCGII_30perct.30=All.PCGII_30perct.30, 
          All.PCGII_30perct.70=All.PCGII_30perct.70,
          All.PCGII_70perct=All.PCGII_70perct, 
          All.PCGII_70perct.30=All.PCGII_70perct.30, 
          All.PCGII_70perct.70=All.PCGII_70perct.70,
          
          All.FGGM=All.FGGM,
          
          # true positives
          TP.enf.p=TP.enf.p,TP.enf.q=TP.enf.q, 
          TP.mle.p=TP.mle.p,TP.mle.q=TP.mle.q, 
          TP.cLevel=TP.cLevel,
          
          TP.PCGII_all=TP.PCGII_all, 
          TP.PCGII_all.30=TP.PCGII_all.30, 
          TP.PCGII_all.70=TP.PCGII_all.70,

          TP.PCGII_30perct=TP.PCGII_30perct, 
          TP.PCGII_30perct.30=TP.PCGII_30perct.30, 
          TP.PCGII_30perct.70=TP.PCGII_30perct.70,
          TP.PCGII_70perct=TP.PCGII_70perct, 
          TP.PCGII_70perct.30=TP.PCGII_70perct.30, 
          TP.PCGII_70perct.70=TP.PCGII_70perct.70,
          
          TP.FGGM=TP.FGGM,
          
          # false positives
          FP.enf.p=FP.enf.p,FP.enf.q=FP.enf.q, 
          FP.mle.p=FP.mle.p, FP.mle.q=FP.mle.q,
          FP.cLevel=FP.cLevel, 
          
          FP.PCGII_all=FP.PCGII_all, 
          FP.PCGII_all.30=FP.PCGII_all.30, 
          FP.PCGII_all.70=FP.PCGII_all.70,

          FP.PCGII_30perct=FP.PCGII_30perct, 
          FP.PCGII_30perct.30=FP.PCGII_30perct.30, 
          FP.PCGII_30perct.70=FP.PCGII_30perct.70,
          FP.PCGII_70perct=FP.PCGII_70perct, 
          FP.PCGII_70perct.30=FP.PCGII_70perct.30, 
          FP.PCGII_70perct.70=FP.PCGII_70perct.70,
          
          FP.FGGM=FP.FGGM,
          
          # true negatives
          TN.enf.p=TN.enf.p, TN.enf.q=TN.enf.q, 
          TN.mle.p=TN.mle.p, TN.mle.q=TN.mle.q,
          TN.cLevel=TN.cLevel, 
          
          TN.PCGII_all=TN.PCGII_all, 
          TN.PCGII_all.30=TN.PCGII_all.30, 
          TN.PCGII_all.70=TN.PCGII_all.70,

          TN.PCGII_30perct=TN.PCGII_30perct, 
          TN.PCGII_30perct.30=TN.PCGII_30perct.30, 
          TN.PCGII_30perct.70=TN.PCGII_30perct.70,
          TN.PCGII_70perct=TN.PCGII_70perct, 
          TN.PCGII_70perct.30=TN.PCGII_70perct.30, 
          TN.PCGII_70perct.70=TN.PCGII_70perct.70,
          
          TN.FGGM=TN.FGGM,

          # false negatives
          FN.enf.p=FN.enf.p, FN.enf.q=FN.enf.q, 
          FN.mle.p=FN.mle.p,FN.mle.q=FN.mle.q,
          FN.cLevel=FN.cLevel,
          
          FN.PCGII_all=FN.PCGII_all, 
          FN.PCGII_all.30=FN.PCGII_all.30, 
          FN.PCGII_all.70=FN.PCGII_all.70,

          FN.PCGII_30perct=FN.PCGII_30perct, 
          FN.PCGII_30perct.30=FN.PCGII_30perct.30, 
          FN.PCGII_30perct.70=FN.PCGII_30perct.70,
          FN.PCGII_70perct=FN.PCGII_70perct, 
          FN.PCGII_70perct.30=FN.PCGII_70perct.30, 
          FN.PCGII_70perct.70=FN.PCGII_70perct.70,
          
          FN.FGGM=FN.FGGM
        ),
        # empirical fdr
        fdr=list(
          fdr.enf.p=fdr.enf.p, fdr.enf.q=fdr.enf.q, 
          fdr.mle.p=fdr.mle.p, fdr.mle.q=fdr.mle.q,
          fdr.cLevel=fdr.cLevel, 
          
          fdr.PCGII_all=fdr.PCGII_all, 
          fdr.PCGII_all.30=fdr.PCGII_all.30, 
          fdr.PCGII_all.70=fdr.PCGII_all.70,

          fdr.PCGII_30perct=fdr.PCGII_30perct, 
          fdr.PCGII_30perct.30=fdr.PCGII_30perct.30, 
          fdr.PCGII_30perct.70=fdr.PCGII_30perct.70,
          fdr.PCGII_70perct=fdr.PCGII_70perct, 
          fdr.PCGII_70perct.30=fdr.PCGII_70perct.30, 
          fdr.PCGII_70perct.70=fdr.PCGII_70perct.70,
          
          fdr.FGGM=fdr.FGGM,
          
          fdr.table=fdr.table
        )
      )
    )
  } # return
}


#### Eval_Strong_Signal #### 

#### Scale Free 1
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      load(file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/ScaleFree1_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=50)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220917/strong_signal/Results/ScaleFree1_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}


#### Scale Free .5
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      load(file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/ScaleFree.5_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=25)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220917/strong_signal/Results/ScaleFree.5_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}

#### Scale Free .1
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      load(file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/ScaleFree.1_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=25)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220917/strong_signal/Results/ScaleFree.1_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}

#### Random

nl = c(60) # Sample Size
pl = c(100, 200)  # Number of Genes

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.02,0.03)){
      load(file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta))
      res=eval_models(X, omega, rep=50)
      save(res,file=sprintf("~/Desktop/GenePCG/R/20220917/strong_signal/Results/Random_n%d_p%d_eta%g.RData", n, p, eta))
    }
  }
}

nl = c(80) # Sample Size
pl = c(100, 200)  # Number of Genes

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.02,0.03)){
      load(file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta))
      res=eval_models(X, omega, rep=50)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220917/strong_signal/Results/Random_n%d_p%d_eta%g.RData", n, p, eta))
    }
  }
}

#### BlockDiag
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in c(4,8,10)){
      load(file=sprintf("~/Desktop/GenePCG/R/20220521/strong_signal/BlockDiag_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=50)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220917/strong_signal/Results/BlockDiag_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}


#### Eval_Weak_Signal #### 

#### Scale Free 1
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      load(file=sprintf("~/Desktop/GenePCG/R/weak_signal/ScaleFree1_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=20)
      save(res, file=sprintf("~/Desktop/GenePCG/R/weak_signal/Results/ScaleFree1_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}


#### Scale Free .5
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      load(file=sprintf("~/Desktop/GenePCG/R/weak_signal/ScaleFree.5_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=20)
      save(res, file=sprintf("~/Desktop/GenePCG/R/weak_signal/Results/ScaleFree.5_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}

#### Scale Free .1
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      load(file=sprintf("~/Desktop/GenePCG/R/weak_signal/ScaleFree.1_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=20)
      save(res, file=sprintf("~/Desktop/GenePCG/R/weak_signal/Results/ScaleFree.1_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}


#### Random

nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.02,0.03)){
      load(file=sprintf("~/Desktop/GenePCG/R/weak_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta))
      res=eval_models(X, omega, rep=20)
      save(res,sprintf("~/Desktop/GenePCG/R/weak_signal/Results/Random_n%d_p%d_eta%g.RData", n, p, eta))
    }
  }
}

#### BlockDiag
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in c(4,8,10)){
      load(file=sprintf("~/Desktop/GenePCG/R/weak_signal/BlockDiag_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=20)
      save(res, file=sprintf("~/Desktop/GenePCG/R/weak_signal/Results/BlockDiag_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}


#### Plots ####
dim(res$ROC$ROC.out)
colnames(res$ROC$ROC.out)
res$ROC$ROC.out %>% 
  tidyr::gather("methods", "TPR", ENF:FGGM) %>%
  filter(!methods %in% c("ENF")) %>% 
  mutate(approach=factor(c(rep("MLE",1000), rep("CLEVEL",1000), 
                    rep("PCGII_all",1000*3),
                    rep("PCGII_70perct",1000*3),
                    rep("PCGII_30perct",1000*3),
                    rep("FGGM",1000)), levels = c("PCGII_all","PCGII_70perct", "PCGII_30perct", "CLEVEL", "FGGM", "MLE"))) %>% 
  mutate(Accuracy=factor(c(rep("1",1000), rep("1",1000), 
                    rep("1",1000),rep(".7",1000),rep(".3",1000),
                    rep("1",1000),rep(".7",1000),rep(".3",1000),
                    rep("1",1000),rep(".7",1000),rep(".3",1000),
                    rep("1",1000)),levels = c("1",".7",".3"))) %>% 
  ggplot(aes(x=FPR, y=TPR,col=approach)) + 
  geom_line(aes(linetype=Accuracy)) + 
  scale_linetype_manual(values=c("solid","dashed","dotted"))+
  xlim(c(0,0.35))

res$fdr$fdr.table %>% 
  filter(!methods %in% c("ENF.p","MLE.p","ENF.q")) %>%
  mutate(approach=factor(c(rep("MLE",58), rep("CLEVEL",58), 
                           rep("PCGII_all",58*3),
                           rep("PCGII_70perct",58*3),
                           rep("PCGII_30perct",58*3),
                           rep("FGGM",58)), levels = c("PCGII_all","PCGII_70perct", "PCGII_30perct", "CLEVEL", "FGGM", "MLE"))) %>% 
  mutate(Accuracy=factor(c(rep("1",58), rep("1",58), 
                           rep("1",58),rep(".7",58),rep(".3",58),
                           rep("1",58),rep(".7",58),rep(".3",58),
                           rep("1",58),rep(".7",58),rep(".3",58),
                           rep("1",58)),levels = c("1",".7",".3"))) %>% 
  ggplot(aes(x=FDR, y=fdr, col=approach))+
  geom_line(aes(linetype=Accuracy)) + 
  scale_linetype_manual(values=c("solid","dashed","dotted"))+
  geom_abline(slope = 1, intercept = 0)


my_power=function(RESULT, alpha=0.05){ # TP/P=TP/(TP+FN)
  ref=seq(0.001, 0.2005, 0.0035)
  ind=which(ref==alpha)
  recall=cbind(RESULT$metrics$tp.enf.p[ind,]/(RESULT$metrics$tp.enf.p[ind,]+RESULT$metrics$fn.enf.p[ind,]), # TPR
               RESULT$metrics$tp.enf.q[ind,]/(RESULT$metrics$tp.enf.q[ind,]+RESULT$metrics$fn.enf.q[ind,]),
               RESULT$metrics$tp.mle.p[ind,]/(RESULT$metrics$tp.mle.p[ind,]+RESULT$metrics$fn.mle.p[ind,]),
               RESULT$metrics$tp.mle.q[ind,]/(RESULT$metrics$tp.mle.q[ind,]+RESULT$metrics$fn.mle.q[ind,]),
               
               RESULT$metrics$tp.cLevel[ind,]/(RESULT$metrics$tp.cLevel[ind,]+RESULT$metrics$fn.cLevel[ind,]),
               
               RESULT$metrics$tp.PCGII_all[ind,]/(RESULT$metrics$tp.PCGII_all[ind,]+RESULT$metrics$fn.PCGII_all[ind,]),
               RESULT$metrics$tp.PCGII_all.30[ind,]/(RESULT$metrics$tp.PCGII_all.30[ind,]+RESULT$metrics$fn.PCGII_all.30[ind,]),
               RESULT$metrics$tp.PCGII_all.70[ind,]/(RESULT$metrics$tp.PCGII_all.70[ind,]+RESULT$metrics$fn.PCGII_all.70[ind,]),
               
               RESULT$metrics$tp.PCGII_30perct[ind,]/(RESULT$metrics$tp.PCGII_30perct[ind,]+RESULT$metrics$fn.PCGII_30perct[ind,]),
               RESULT$metrics$tp.PCGII_30perct.30[ind,]/(RESULT$metrics$tp.PCGII_30perct.30[ind,]+RESULT$metrics$fn.PCGII_30perct.30[ind,]),
               RESULT$metrics$tp.PCGII_30perct.70[ind,]/(RESULT$metrics$tp.PCGII_30perct.70[ind,]+RESULT$metrics$fn.PCGII_30perct.70[ind,]),
               
               RESULT$metrics$tp.PCGII_70perct[ind,]/(RESULT$metrics$tp.PCGII_70perct[ind,]+RESULT$metrics$fn.PCGII_70perct[ind,]),
               RESULT$metrics$tp.PCGII_70perct.30[ind,]/(RESULT$metrics$tp.PCGII_70perct.30[ind,]+RESULT$metrics$fn.PCGII_70perct.30[ind,]),
               RESULT$metrics$tp.PCGII_70perct.70[ind,]/(RESULT$metrics$tp.PCGII_70perct.70[ind,]+RESULT$metrics$fn.PCGII_70perct.70[ind,]),
               
               RESULT$metrics$tp.FGGM[ind,]/(RESULT$metrics$tp.FGGM[ind,]+RESULT$metrics$fn.FGGM[ind,])
  )
  
  
  power=as.data.frame(recall)
  colnames(power)=c("ENF.p","ENF.q","MLE.p","MLE.q",
                    "cLevel",
                    "PCGII_all","PCGII_all.30","PCGII_all.70",
                    "PCGII_30perct","PCGII_30perct.30","PCGII_30perct.70",
                    "PCGII_70perct","PCGII_70perct.30","PCGII_70perct.70",
                    "FGGM")
  power
}

my_power(res) %>% gather(key="approach",value = "power") %>% 
  filter(approach %in% c("cLevel", "MLE.q","FGGM","PCGII_70perct.70")) %>%
  mutate(approach = factor(plyr::revalue(approach, c("cLevel" = "CLEVEL", "MLE.q"="E-TEST","FGGM"="F-GGM","PCGII_70perct.70"="PCGII")), levels=c("PCGII","CLEVEL","F-GGM", "E-TEST"))) %>% ggplot(aes(x=approach, y=power))+geom_boxplot()

my_power(res) %>% gather(key="approach",value = "power") %>% 
  filter(!approach %in% c("ENF.p", "ENF.q","MLE.p","MLE.q","FGGM")) %>%
  ggplot(aes(x=approach, y=power))+geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

my_fdr=function(RESULT, alpha=0.05){ #FP/Disc=FP/(FP+TP)
  ref=seq(0.001, 0.2005, 0.0035)
  ind=which(ref==alpha)
  fdr=cbind.data.frame(RESULT$metrics$fp.enf.p[ind,]/(RESULT$metrics$tp.enf.p[ind,]+RESULT$metrics$fp.enf.p[ind,]), 
                       RESULT$metrics$fp.enf.q[ind,]/(RESULT$metrics$tp.enf.q[ind,]+RESULT$metrics$fp.enf.q[ind,]),
                       RESULT$metrics$fp.mle.p[ind,]/(RESULT$metrics$tp.mle.p[ind,]+RESULT$metrics$fp.mle.p[ind,]),
                       RESULT$metrics$fp.mle.q[ind,]/(RESULT$metrics$tp.mle.q[ind,]+RESULT$metrics$fp.mle.q[ind,]),
                       
                       RESULT$metrics$fp.cLevel[ind,]/(RESULT$metrics$tp.cLevel[ind,]+RESULT$metrics$fp.cLevel[ind,]),
                       
                       RESULT$metrics$fp.PCGII_all[ind,]/(RESULT$metrics$tp.PCGII_all[ind,]+RESULT$metrics$fp.PCGII_all[ind,]),
                       RESULT$metrics$fp.PCGII_all.30[ind,]/(RESULT$metrics$tp.PCGII_all.30[ind,]+RESULT$metrics$fp.PCGII_all.30[ind,]),
                       RESULT$metrics$fp.PCGII_all.70[ind,]/(RESULT$metrics$tp.PCGII_all.70[ind,]+RESULT$metrics$fp.PCGII_all.70[ind,]),

                       RESULT$metrics$fp.PCGII_30perct[ind,]/(RESULT$metrics$tp.PCGII_30perct[ind,]+RESULT$metrics$fp.PCGII_30perct[ind,]),
                       RESULT$metrics$fp.PCGII_30perct.30[ind,]/(RESULT$metrics$tp.PCGII_30perct.30[ind,]+RESULT$metrics$fp.PCGII_30perct.30[ind,]),
                       RESULT$metrics$fp.PCGII_30perct.70[ind,]/(RESULT$metrics$tp.PCGII_30perct.70[ind,]+RESULT$metrics$fp.PCGII_30perct.70[ind,]),
                       RESULT$metrics$fp.PCGII_70perct[ind,]/(RESULT$metrics$tp.PCGII_70perct[ind,]+RESULT$metrics$fp.PCGII_70perct[ind,]),
                       RESULT$metrics$fp.PCGII_70perct.30[ind,]/(RESULT$metrics$tp.PCGII_70perct.30[ind,]+RESULT$metrics$fp.PCGII_70perct.30[ind,]),
                       RESULT$metrics$fp.PCGII_70perct.70[ind,]/(RESULT$metrics$tp.PCGII_70perct.70[ind,]+RESULT$metrics$fp.PCGII_70perct.70[ind,]),
                       
                       RESULT$metrics$fp.FGGM[ind,]/(RESULT$metrics$tp.FGGM[ind,]+RESULT$metrics$fp.FGGM[ind,])
                       
  )
  
  colnames(fdr)=c("ENF.p","ENF.q","MLE.p","MLE.q",
                  "cLevel",
                  "PCGII_all","PCGII_all.30","PCGII_all.70",
                  "PCGII_30perct","PCGII_30perct.30","PCGII_30perct.70",
                  "PCGII_70perct","PCGII_70perct.30","PCGII_70perct.70",
                  "FGGM")
  fdr[is.na(fdr)] <- 0
  fdr
}

my_fdr(res) %>% gather(key="approach",value = "fdr") %>% 
  filter(approach %in% c("cLevel", "MLE.q","FGGM","PCGII_70perct.70")) %>%
  mutate(approach = factor(plyr::revalue(approach, c("cLevel" = "CLEVEL", "MLE.q"="E-TEST","FGGM"="F-GGM","PCGII_70perct.70"="PCGII")), levels=c("PCGII","CLEVEL","F-GGM", "E-TEST"))) %>% 
  ggplot(aes(x=approach, y=fdr))+
  geom_boxplot() + 
  geom_abline(slope=0, intercept = 0.05)

my_fdr(res) %>% gather(key="approach",value = "fdr") %>% 
  filter(!approach %in% c("ENF.p", "ENF.q","MLE.p","MLE.q","FGGM")) %>%
  ggplot(aes(x=approach, y=fdr))+
  geom_boxplot()+ 
  geom_abline(slope=0, intercept = 0.05) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# computing time ####

library(microbenchmark)
microbenchmark(
  "PCGII" = {
    temp1=inference(PCGII(df=sim.data, prior=double_prior(pri), lambda=lam))
    },
  "FGGM" = {
    temp2=FastGGM(sim.data, lambda=lam) 
    },
  times = 100, check = NULL, setup=set.seed(21)
  )


# n=80, p=100, scalefree m=2
# Unit: milliseconds
# expr       min        lq      mean    median        uq      max neval cld
# PCGII 336.36597 347.26473 357.84420 352.84758 358.11552 463.2188   100   b
# FGGM  46.78391  50.02574  53.92761  51.88242  55.62025 161.6234   100  a 
# 
# n=80, p=200, scalefree m=2
# Unit: milliseconds
# expr       min       lq      mean    median        uq       max neval cld
# PCGII 1169.8620 1219.012 1261.7604 1246.3490 1309.4233 1451.2327   100   b
# FGGM  234.6765  246.209  261.5894  250.9319  259.7453  376.1432   100  a 
# 
# elapsed, user.self, sys.self, user.child, and sys.child are columns containing values reported by system.time; see Sec. 7.1 Operating system access in The R language definition, or see system.time.
# The ‘user time’ is the CPU time charged for the execution of user instructions of the calling process. The ‘system time’ is the CPU time charged for execution by the system on behalf of the calling process.

