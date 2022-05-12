############# Last updated on 20220508 #############;

# This simulation study generates data from less ill-conditioned precision matrix, i.e. eigenvalue-ratio is not too large, although the partial correlation strength would be weaker;
# This simulation study applies the new-developed PCGII algorithm, which estimates the covariance matrix of residuals first via no-info node-wise penalized regressions and then corrects the bias with coefficients estimated via information incorporated regressions; a two-step approach
# This simulation study also checks the effects of tunning parameters empirically.

library(igraph)
library(tidyverse)
library(GeneNet)
library(FastGGM)
library(corpcor)
library(glmnet)
library(mvtnorm)

source("~/Desktop/GenePCG/R/shrinkagefunctions.R") # https://github.com/V-Bernal/GGM-Shrinkage


############# Simulate Data, 2022/05/07 #############;

makeBlockDiag=function(blocksize=4, p=100, min.beta=0.1, max.beta=1){ # blocksize has to be a factor of p
  reps=p/blocksize
  S=list()
  for (i in 1:reps) {
    bd=matrix(runif(1,min.beta,max.beta),blocksize,blocksize)
    diag(bd)=runif(1,1,1.25)
    S[[i]]=bd
  }
  as.matrix(Matrix::bdiag(S))
}

#### Strong Signal ####

## Scale Free 1
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220420)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      g <- sample_pa(p, power=1, m=e, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>50) {
        diag(omega)=4+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/ScaleFree1_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Scale Free.5
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220420)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      g <- sample_pa(p, power=.5, m=e, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>50) {
        diag(omega)=4+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/ScaleFree.5_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Scale Free.1
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220421)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      g <- sample_pa(p, power=.1, m=e, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>50) {
        diag(omega)=4+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/ScaleFree.1_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Random
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220407)

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.02,0.03)){
      g <- sample_gnp(p, eta, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>50) {
        diag(omega)=4+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta)) 
    }
  }
}


## Block Diagonal
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220410)

for(n in nl){    
  for(p in pl){
    for(e in c(4,8,10)){
      Sigma=makeBlockDiag(blocksize=e, p=p, min.beta=.3, max.beta=.9)
      omega=solve(Sigma)
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      X = list()
      for(i in 1:10){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/BlockDiag_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}



#### Weak Signal ####

## Scale Free 1
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220420)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      g <- sample_pa(p, power=1, m=e, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>5) {
        diag(omega)=4+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/weak_signal/ScaleFree1_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Scale Free.5
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220420)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      g <- sample_pa(p, power=.5, m=e, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>5) {
        diag(omega)=4+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/weak_signal/ScaleFree.5_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Scale Free.1
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220421)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      g <- sample_pa(p, power=.1, m=e, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>5) {
        diag(omega)=4+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/weak_signal/ScaleFree.1_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Random
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220407)

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.02,0.03)){
      g <- sample_gnp(p, eta, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>5) {
        diag(omega)=4+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/weak_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta)) 
    }
  }
}


## Block Diagonal
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220410)

for(n in nl){    
  for(p in pl){
    for(e in c(4,8,10)){
      Sigma=makeBlockDiag(blocksize=e, p=p, min.beta=0.1, max.beta=0.6)
      omega=solve(Sigma)
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/weak_signal/BlockDiag_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}



#### Mixed Signal ####

## Scale Free 1
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220420)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      g <- sample_pa(p, power=1, m=e, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>10) {
        diag(omega)[1:(p/2)]=4+increment
        diag(omega)[(p/2 + 1):(p)]=7+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/mixed_signal/ScaleFree1_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Scale Free.5
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220420)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      g <- sample_pa(p, power=.5, m=e, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>10) {
        diag(omega)[1:(p/2)]=4+increment
        diag(omega)[(p/2 + 1):(p)]=7+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/mixed_signal/ScaleFree.5_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Scale Free.1
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220421)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      g <- sample_pa(p, power=.1, m=e, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>10) {
        diag(omega)[1:(p/2)]=4+increment
        diag(omega)[(p/2 + 1):(p)]=7+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/mixed_signal/ScaleFree.1_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

## Random
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220407)

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.02,0.03)){
      g <- sample_gnp(p, eta, directed = FALSE)
      omega=as_adjacency_matrix(g) %>% as.matrix()
      for(h1 in 1:(p-1)){
        for(h2 in (h1+1):p){
          if(omega[h1,h2]!=0){
            temp=runif(1, 0.5, 1.5)*sample(c(-1,1),size=1)
            omega[h1,h2]=temp
            omega[h2,h1]=temp
          }
        }
      }
      
      diag(omega)=4
      increment=.05
      temp=eigen(omega)$values
      while (min(temp)<0 | max(temp)/min(temp)>10) {
        diag(omega)[1:(p/2)]=4+increment
        diag(omega)[(p/2 + 1):(p)]=7+increment
        increment=increment+0.05
        temp=eigen(omega)$values
      }
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      Sigma=solve(omega)
      #is.positive.definite(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/mixed_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta)) 
    }
  }
}


## Block Diagonal
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20220410)

for(n in nl){    
  for(p in pl){
    for(e in c(4,8,10)){
      Sigma_weak=makeBlockDiag(blocksize=e, p=p/2, min.beta=0.1, max.beta=0.6)
      Sigma_strong=makeBlockDiag(blocksize=e, p=p/2, min.beta=0.3, max.beta=0.9)
      Sigma=as.matrix(Matrix::bdiag(list(Sigma_strong, Sigma_weak)))
      omega=solve(Sigma)
      ppc=-sm2vec(cov2cor(omega))
      ppc=ppc[which(ppc!=0)]
      sparsity=round(length(ppc)/(p*(p-1)/2),3)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(omega, sparsity, Sigma, X, file=sprintf("~/Desktop/GenePCG/R/mixed_signal/BlockDiag_simu_n%d_p%d_e%d.RData", n, p, e)) 
    }
  }
}

#### Functions ####

sigs2vec=function(sigs, P){
  require(corpcor)
  m=matrix(0,P,P)
  for (h in 1: dim(sigs)[1]){
    m[sigs[h,1],sigs[h,2]]=1
  }
  sm2vec(m)
}

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

PCGII2=function(df, prior, lambda){
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
  
  if(missing(lambda)){
    shat=sqrt(n/(log(p)^3))
    lambda=sqrt(2*(2+0.01)*log(p/shat)/n)    
  }
  
  default_penalty=rep(1,p-1)
  for (i in 1 : p){
    penalty_fac=default_penalty
    temp.node=prior[with(prior,row==i),'col']
    
    for(nds in temp.node){
      if (nds < i) {penalty_fac[nds]=0} else {penalty_fac[nds-1]=0.3}
    }
    
    out = glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda, penalty.factor=penalty_fac)
    Coef = out$beta
    CoefMatrix[i, ] = as.numeric(Coef) / apply(X[, -i], 2, sd)
    
    out = glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda)
    Predict = predict(out, XS[, -i], type = "link")
    Eresidual[, i] = X[, i] - Predict
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

clevel=function(df, lambda){
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
  
  if(missing(lambda)){
    shat=sqrt(n/(log(p)^3))
    lambda=sqrt(2*(2+0.01)*log(p/shat)/n)    
  }
  
  for (i in 1 : p){
    out = glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda)
    
    Coef = out$beta
    Predict = predict(out, XS[, -i], type = "link")
    CoefMatrix[i, ] = as.numeric(Coef) / apply(X[, -i], 2, sd)
    Eresidual[, i] = X[, i] - Predict
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

inference=function(list, alpha=0.05, c0=0.25){
  Est=list$Est
  tscore=list$tscore
  kappa=list$kappa
  EstThresh=list$EstThresh
  n=list$n; p=list$p; 
  t0=2; tau = seq(0, 3.5, 0.01); smax = n / 2; lentau = length(tau) 
  
  
  resprop = list() # selected edges with different tau's, a list of 351 elements
  rejectprop = c()
  for (i in 1 : lentau){ # tau vary from 0 to 3.50 by 0.01, length=351
    Threshold = tau[i] * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
    
    # c=0
    SRec = 1 * (abs(Est) > Threshold) # selected edge (matrix with 0 & 1) at tau[i]
    NoNSRec = 1 * (SRec == 0)
    resprop[[i]] = which(SRec == 1, arr.ind = TRUE) # select those significant edges at tau[i], off-diagonal elements, first columns, then second columns
    rejectprop = c(rejectprop, max(1, (sum(SRec) - p))) 
  }
  
  # c=0
  FDPprop = 2 * (p * (p - 1)) * ( 1 - pnorm( tau * sqrt(log(p)) ) ) / rejectprop # FDP corresponding to each tau (page 10)
  
  FDPresprop = c()
  
  # determine thresholding parameter tau by controling FDP 
  if (sum(FDPprop <= alpha) > 0) tauprop = min(c(2, tau[FDPprop <= alpha])) 
  if (sum(FDPprop <= alpha) == 0) tauprop = 2
  Threshold = tauprop * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
  SRec = 1 * (abs(Est) > Threshold); NoNSRec = 1 * (SRec == 0) # SRec is a matrix (0-1 matrix)
  FDPresprop = which(SRec == 1, arr.ind = TRUE) # selected edge location
  
  sigs=as.data.frame(FDPresprop[which(FDPresprop[,1]!=FDPresprop[,2]),])
  #colnames(sigs)=c("node1","node2")
  
  return(list(sigs=sigs))  
}

eval_models=function(X, omega, rep=3){
  
  p=dim(omega)[1] 
  n=dim(X[[1]])[1]
  truth<-sm2vec(omega, diag = F) 
  TP<- which(truth!=0) 
  CN=p*(p-1)/2-length(TP) 
  CP=length(TP)
  
  shat=sqrt(n/(log(p)^3))
  lam1=sqrt(2*(2+0.01)*log(p/shat)/n)   
  lam2=sqrt(2*log(p/sqrt(n))/n) ## FGGM
  lam3=sqrt(log(p)/n) ## smallest
    
  true.par=-cov2cor(omega)
  diag(true.par)=0
  
  prior.temp=cbind.data.frame(which(true.par!=0, arr.ind = TRUE), signal=round(true.par[which(true.par!=0)],6)) %>% 
    transform(row = pmin(row, col), col = pmax(row, col)) %>% 
    arrange(row, col) %>% 
    unique() %>%
    arrange(desc(abs(signal))) ## all
  
  if(dim(prior.temp)[1]!=CP) {stop("Prior information error!")}
  
  prior.all=prior.temp
  prior.top5=prior.temp[1:5,] 
  prior.top15=prior.temp[1:15,] 
  prior.top25=prior.temp[1:25,] 
  prior.random10=prior.temp[sample(1:dim(prior.temp)[1],size = 10),]
  
  prior.all=rbind.data.frame(prior.all, prior.all %>% transform(row = pmax(row, col), col = pmin(row, col))) %>% 
    arrange(row, col)
  prior.top5=rbind.data.frame(prior.top5, prior.top5 %>% transform(row = pmax(row, col), col = pmin(row, col))) %>% 
    arrange(row, col)
  prior.top15=rbind.data.frame(prior.top15, prior.top15 %>% transform(row = pmax(row, col), col = pmin(row, col))) %>% 
    arrange(row, col)
  prior.top25=rbind.data.frame(prior.top25, prior.top25 %>% transform(row = pmax(row, col), col = pmin(row, col))) %>% 
    arrange(row, col)
  prior.random10=rbind.data.frame(prior.random10, prior.random10 %>% transform(row = pmax(row, col), col = pmin(row, col))) %>% 
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
    
    cLevel_lambda1.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_lambda1.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_lambda1.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_lambda1.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    cLevel_lambda2.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_lambda2.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_lambda2.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_lambda2.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    cLevel_lambda3.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_lambda3.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_lambda3.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    cLevel_lambda3.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_all_lambda1.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all_lambda1.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all_lambda1.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all_lambda1.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_all_lambda2.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all_lambda2.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all_lambda2.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all_lambda2.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_all_lambda3.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all_lambda3.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all_lambda3.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_all_lambda3.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_top5_lambda1.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top5_lambda1.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top5_lambda1.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top5_lambda1.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_top5_lambda2.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top5_lambda2.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top5_lambda2.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top5_lambda2.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_top5_lambda3.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top5_lambda3.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top5_lambda3.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top5_lambda3.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_top15_lambda1.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top15_lambda1.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top15_lambda1.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top15_lambda1.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_top15_lambda2.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top15_lambda2.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top15_lambda2.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top15_lambda2.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_top15_lambda3.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top15_lambda3.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top15_lambda3.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top15_lambda3.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_top25_lambda1.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top25_lambda1.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top25_lambda1.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top25_lambda1.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_top25_lambda2.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top25_lambda2.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top25_lambda2.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top25_lambda2.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_top25_lambda3.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top25_lambda3.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top25_lambda3.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_top25_lambda3.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_random10_lambda1.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_random10_lambda1.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_random10_lambda1.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_random10_lambda1.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)  
    
    PCGII_random10_lambda2.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_random10_lambda2.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_random10_lambda2.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_random10_lambda2.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    PCGII_random10_lambda3.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_random10_lambda3.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_random10_lambda3.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    PCGII_random10_lambda3.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    
    FGGM_lambda1.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM_lambda1.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM_lambda1.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM_lambda1.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep) 
    
    FGGM_lambda2.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM_lambda2.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM_lambda2.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM_lambda2.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep) 
    
    FGGM_lambda3.ordered.TPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM_lambda3.ordered.FPR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM_lambda3.ordered.FDR=matrix(Inf, nrow=p*(p-1)/2, ncol=rep)
    FGGM_lambda3.ordered.PPV=matrix(Inf, nrow=p*(p-1)/2, ncol=rep) 
    
    tempstore.ENF=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.MLE=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points 
    tempstore.cLevel_lambda1=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.cLevel_lambda2=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.cLevel_lambda3=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_all_lambda1=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_all_lambda2=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_all_lambda3=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_top5_lambda1=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_top5_lambda2=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_top5_lambda3=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_top15_lambda1=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_top15_lambda2=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_top15_lambda3=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_top25_lambda1=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_top25_lambda2=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_top25_lambda3=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_random10_lambda1=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_random10_lambda2=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.PCGII_random10_lambda3=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.FGGM_lambda1=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.FGGM_lambda2=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    tempstore.FGGM_lambda3=matrix(Inf, nrow=1000, ncol=rep) # for ROC curve, 1000 cutting points
    
    ### Inference
    
    # Number of total correctly selected edges
    all.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    all.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.mle.p=matrix(Inf, nrow=length(al), ncol=rep) 
    all.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    all.cLevel_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    all.cLevel_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    all.cLevel_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_all_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_all_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_all_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_top5_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_top5_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_top5_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_top15_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_top15_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_top15_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_top25_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_top25_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_top25_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_random10_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_random10_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    all.PCGII_random10_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    all.FGGM_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    all.FGGM_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    all.FGGM_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    
    # true positives
    tp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tp.cLevel_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tp.cLevel_lambda2=matrix(Inf, nrow=length(al), ncol=rep)    
    tp.cLevel_lambda3=matrix(Inf, nrow=length(al), ncol=rep)    
    tp.PCGII_all_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_all_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_all_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_top5_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_top5_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_top5_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_top15_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_top15_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_top15_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_top25_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_top25_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_top25_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_random10_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_random10_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tp.PCGII_random10_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tp.FGGM_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tp.FGGM_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tp.FGGM_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false positives
    fp.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fp.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fp.cLevel_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fp.cLevel_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fp.cLevel_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_all_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_all_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_all_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_top5_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_top5_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_top5_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_top15_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_top15_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_top15_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_top25_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_top25_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_top25_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_random10_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_random10_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fp.PCGII_random10_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fp.FGGM_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fp.FGGM_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fp.FGGM_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    
    # true negatives
    tn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    tn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    tn.cLevel_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tn.cLevel_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tn.cLevel_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_all_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_all_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_all_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_top5_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_top5_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_top5_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_top15_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_top15_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_top15_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_top25_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_top25_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_top25_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_random10_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_random10_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tn.PCGII_random10_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    tn.FGGM_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    tn.FGGM_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    tn.FGGM_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    
    # false negatives
    fn.enf.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.enf.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.p=matrix(Inf, nrow=length(al), ncol=rep)
    fn.mle.q=matrix(Inf, nrow=length(al), ncol=rep)
    fn.cLevel_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fn.cLevel_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fn.cLevel_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_all_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_all_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_all_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_top5_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_top5_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_top5_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_top15_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_top15_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_top15_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_top25_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_top25_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_top25_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_random10_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_random10_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fn.PCGII_random10_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
    fn.FGGM_lambda1=matrix(Inf, nrow=length(al), ncol=rep)
    fn.FGGM_lambda2=matrix(Inf, nrow=length(al), ncol=rep)
    fn.FGGM_lambda3=matrix(Inf, nrow=length(al), ncol=rep)
  }
  
  
  for (k in 1:rep){
    print(paste0(k,"th rep starts!!!!!"))
    sim.data=X[[k]]
    
    ### Fast GGM
    out=FastGGM(sim.data, lambda=lam1)
    FGGM_lambda1.test=unMat(X_est=out$partialCor, X_p=out$p_partialCor) # Est, pvals
    FGGM_lambda1.test=cbind.data.frame(FGGM_lambda1.test, truth=truth, qval=p.adjust(FGGM_lambda1.test[,4], method="BH"))
    
    out=FastGGM(sim.data, lambda=lam2)
    FGGM_lambda2.test=unMat(X_est=out$partialCor, X_p=out$p_partialCor) # Est, pvals
    FGGM_lambda2.test=cbind.data.frame(FGGM_lambda2.test, truth=truth, qval=p.adjust(FGGM_lambda2.test[,4], method="BH"))
    
    out=FastGGM(sim.data, lambda=lam3)
    FGGM_lambda3.test=unMat(X_est=out$partialCor, X_p=out$p_partialCor) # Est, pvals
    FGGM_lambda3.test=cbind.data.frame(FGGM_lambda3.test, truth=truth, qval=p.adjust(FGGM_lambda3.test[,4], method="BH"))
    
    print("Fast GGM done")
    
    ### cLevel
    cLevel_lambda1=clevel(sim.data, lambda = lam1)
    cLevel_lambda2=clevel(sim.data, lambda = lam2)
    cLevel_lambda3=clevel(sim.data, lambda = lam3)
    # estimates by cLevel
    Est_cLevel_lambda1=sm2vec(cLevel_lambda1$Est, diag = F)
    Est_cLevel_lambda2=sm2vec(cLevel_lambda2$Est, diag = F)
    Est_cLevel_lambda3=sm2vec(cLevel_lambda3$Est, diag = F)
    # test statistics of cLevel
    tscore_cLevel_lambda1=sm2vec(cLevel_lambda1$tscore, diag = F)
    tscore_cLevel_lambda2=sm2vec(cLevel_lambda2$tscore, diag = F)
    tscore_cLevel_lambda3=sm2vec(cLevel_lambda3$tscore, diag = F)
    print("cLevel done")
    
    ### PCGII-prior.all
    PCGII_all_lambda1=PCGII2(df=sim.data, prior = prior.all, lambda = lam1)
    PCGII_all_lambda2=PCGII2(df=sim.data, prior = prior.all, lambda = lam2)
    PCGII_all_lambda3=PCGII2(df=sim.data, prior = prior.all, lambda = lam3)
    # estimates by PCGII
    Est_PCGII_all_lambda1=sm2vec(PCGII_all_lambda1$Est, diag = F)
    Est_PCGII_all_lambda2=sm2vec(PCGII_all_lambda2$Est, diag = F)
    Est_PCGII_all_lambda3=sm2vec(PCGII_all_lambda3$Est, diag = F)
    # test statistics of PCGII
    tscore_PCGII_all_lambda1=sm2vec(PCGII_all_lambda1$tscore, diag = F)
    tscore_PCGII_all_lambda2=sm2vec(PCGII_all_lambda2$tscore, diag = F)
    tscore_PCGII_all_lambda3=sm2vec(PCGII_all_lambda3$tscore, diag = F)
    
    ### PCGII-prior.top5
    PCGII_top5_lambda1=PCGII2(df=sim.data, prior = prior.top5, lambda = lam1)
    PCGII_top5_lambda2=PCGII2(df=sim.data, prior = prior.top5, lambda = lam2)
    PCGII_top5_lambda3=PCGII2(df=sim.data, prior = prior.top5, lambda = lam3)
    # estimates by PCGII
    Est_PCGII_top5_lambda1=sm2vec(PCGII_top5_lambda1$Est, diag = F)
    Est_PCGII_top5_lambda2=sm2vec(PCGII_top5_lambda2$Est, diag = F)
    Est_PCGII_top5_lambda3=sm2vec(PCGII_top5_lambda3$Est, diag = F)
    # test statistics of PCGII
    tscore_PCGII_top5_lambda1=sm2vec(PCGII_top5_lambda1$tscore, diag = F)
    tscore_PCGII_top5_lambda2=sm2vec(PCGII_top5_lambda2$tscore, diag = F)
    tscore_PCGII_top5_lambda3=sm2vec(PCGII_top5_lambda3$tscore, diag = F)
    
    ### PCGII-prior.top15
    PCGII_top15_lambda1=PCGII2(df=sim.data, prior = prior.top15, lambda = lam1)
    PCGII_top15_lambda2=PCGII2(df=sim.data, prior = prior.top15, lambda = lam2)
    PCGII_top15_lambda3=PCGII2(df=sim.data, prior = prior.top15, lambda = lam3)
    # estimates by PCGII
    Est_PCGII_top15_lambda1=sm2vec(PCGII_top15_lambda1$Est, diag = F)
    Est_PCGII_top15_lambda2=sm2vec(PCGII_top15_lambda2$Est, diag = F)
    Est_PCGII_top15_lambda3=sm2vec(PCGII_top15_lambda3$Est, diag = F)
    # test statistics of PCGII
    tscore_PCGII_top15_lambda1=sm2vec(PCGII_top15_lambda1$tscore, diag = F)
    tscore_PCGII_top15_lambda2=sm2vec(PCGII_top15_lambda2$tscore, diag = F)
    tscore_PCGII_top15_lambda3=sm2vec(PCGII_top15_lambda3$tscore, diag = F)
    
    ### PCGII-prior.top25
    PCGII_top25_lambda1=PCGII2(df=sim.data, prior = prior.top25, lambda = lam1)
    PCGII_top25_lambda2=PCGII2(df=sim.data, prior = prior.top25, lambda = lam2)
    PCGII_top25_lambda3=PCGII2(df=sim.data, prior = prior.top25, lambda = lam3)
    # estimates by PCGII
    Est_PCGII_top25_lambda1=sm2vec(PCGII_top25_lambda1$Est, diag = F)
    Est_PCGII_top25_lambda2=sm2vec(PCGII_top25_lambda2$Est, diag = F)
    Est_PCGII_top25_lambda3=sm2vec(PCGII_top25_lambda3$Est, diag = F)
    # test statistics of PCGII
    tscore_PCGII_top25_lambda1=sm2vec(PCGII_top25_lambda1$tscore, diag = F)
    tscore_PCGII_top25_lambda2=sm2vec(PCGII_top25_lambda2$tscore, diag = F)
    tscore_PCGII_top25_lambda3=sm2vec(PCGII_top25_lambda3$tscore, diag = F)
    
    ### PCGII-prior.random10
    PCGII_random10_lambda1=PCGII2(df=sim.data, prior = prior.random10, lambda = lam1)
    PCGII_random10_lambda2=PCGII2(df=sim.data, prior = prior.random10, lambda = lam2)
    PCGII_random10_lambda3=PCGII2(df=sim.data, prior = prior.random10, lambda = lam3)
    # estimates by PCGII
    Est_PCGII_random10_lambda1=sm2vec(PCGII_random10_lambda1$Est, diag = F)
    Est_PCGII_random10_lambda2=sm2vec(PCGII_random10_lambda2$Est, diag = F)
    Est_PCGII_random10_lambda3=sm2vec(PCGII_random10_lambda3$Est, diag = F)
    # test statistics of PCGII
    tscore_PCGII_random10_lambda1=sm2vec(PCGII_random10_lambda1$tscore, diag = F)
    tscore_PCGII_random10_lambda2=sm2vec(PCGII_random10_lambda2$tscore, diag = F)
    tscore_PCGII_random10_lambda3=sm2vec(PCGII_random10_lambda3$tscore, diag = F)
    
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
                             Est_cLevel_lambda1, Est_cLevel_lambda2, Est_cLevel_lambda3, 
                             Est_PCGII_all_lambda1, Est_PCGII_all_lambda2, Est_PCGII_all_lambda3, 
                             Est_PCGII_top5_lambda1, Est_PCGII_top5_lambda2, Est_PCGII_top5_lambda3, 
                             Est_PCGII_top15_lambda1, Est_PCGII_top15_lambda2, Est_PCGII_top15_lambda3, 
                             Est_PCGII_top25_lambda1, Est_PCGII_top25_lambda2, Est_PCGII_top25_lambda3, 
                             Est_PCGII_random10_lambda1, Est_PCGII_random10_lambda2, Est_PCGII_random10_lambda3, 
                             
                             Est_FGGM_lambda1=FGGM_lambda1.test$pcor, 
                             Est_FGGM_lambda2=FGGM_lambda2.test$pcor, 
                             Est_FGGM_lambda3=FGGM_lambda3.test$pcor,
                             shrunk_p, # Shrunk_pcor
                             
                             tscore_cLevel_lambda1, tscore_cLevel_lambda2, tscore_cLevel_lambda3,
                             tscore_PCGII_all_lambda1, tscore_PCGII_all_lambda2, tscore_PCGII_all_lambda3,
                             tscore_PCGII_top5_lambda1, tscore_PCGII_top5_lambda2, tscore_PCGII_top5_lambda3,
                             tscore_PCGII_top15_lambda1, tscore_PCGII_top15_lambda2, tscore_PCGII_top15_lambda3,
                             tscore_PCGII_top25_lambda1, tscore_PCGII_top25_lambda2, tscore_PCGII_top25_lambda3,
                             tscore_PCGII_random10_lambda1, tscore_PCGII_random10_lambda2, tscore_PCGII_random10_lambda3,
                             
                             FGGM_lambda1_pval=FGGM_lambda1.test$pval,
                             FGGM_lambda1_qval=FGGM_lambda1.test$qval,
                             FGGM_lambda2_pval=FGGM_lambda2.test$pval,
                             FGGM_lambda2_qval=FGGM_lambda2.test$qval,
                             FGGM_lambda3_pval=FGGM_lambda3.test$pval,
                             FGGM_lambda3_qval=FGGM_lambda3.test$qval,
                             
                             ENF_p=ENF.test$pval, ENF_q=p.adjust(ENF.test$pval, method="BH"),
                             ShrunkMLE_p=pval.shrunk, ShrunkMLE_q=p.adjust(pval.shrunk, method="BH"))
    
    
    ## Ranking
    FGGM_lambda1.ordered<-results[order(results$FGGM_lambda1_qval, results$FGGM_lambda1_pval, decreasing = F),c("truth","FGGM_lambda1_qval","FGGM_lambda1_pval")]
    FGGM_lambda2.ordered<-results[order(results$FGGM_lambda2_qval, results$FGGM_lambda2_pval, decreasing = F),c("truth","FGGM_lambda2_qval","FGGM_lambda2_pval")]
    FGGM_lambda3.ordered<-results[order(results$FGGM_lambda3_qval, results$FGGM_lambda3_pval, decreasing = F),c("truth","FGGM_lambda3_qval","FGGM_lambda3_pval")]
    
    
    cLevel_lambda1.ordered<-results[order(abs(results$tscore_cLevel_lambda1), decreasing = T),c("truth","Est_cLevel_lambda1","tscore_cLevel_lambda1")]
    cLevel_lambda2.ordered<-results[order(abs(results$tscore_cLevel_lambda2), decreasing = T),c("truth","Est_cLevel_lambda2","tscore_cLevel_lambda2")]
    cLevel_lambda3.ordered<-results[order(abs(results$tscore_cLevel_lambda3), decreasing = T),c("truth","Est_cLevel_lambda3","tscore_cLevel_lambda3")]
    
    PCGII_all_lambda1.ordered<-results[order(abs(results$tscore_PCGII_all_lambda1), decreasing = T),c("truth","Est_PCGII_all_lambda1","tscore_PCGII_all_lambda1")]
    PCGII_all_lambda2.ordered<-results[order(abs(results$tscore_PCGII_all_lambda2), decreasing = T),c("truth","Est_PCGII_all_lambda2","tscore_PCGII_all_lambda2")]
    PCGII_all_lambda3.ordered<-results[order(abs(results$tscore_PCGII_all_lambda3), decreasing = T),c("truth","Est_PCGII_all_lambda3","tscore_PCGII_all_lambda3")]
    
    PCGII_top5_lambda1.ordered<-results[order(abs(results$tscore_PCGII_top5_lambda1), decreasing = T),c("truth","Est_PCGII_top5_lambda1","tscore_PCGII_top5_lambda1")]
    PCGII_top5_lambda2.ordered<-results[order(abs(results$tscore_PCGII_top5_lambda2), decreasing = T),c("truth","Est_PCGII_top5_lambda2","tscore_PCGII_top5_lambda2")]
    PCGII_top5_lambda3.ordered<-results[order(abs(results$tscore_PCGII_top5_lambda3), decreasing = T),c("truth","Est_PCGII_top5_lambda3","tscore_PCGII_top5_lambda3")]
    
    PCGII_top15_lambda1.ordered<-results[order(abs(results$tscore_PCGII_top15_lambda1), decreasing = T),c("truth","Est_PCGII_top15_lambda1","tscore_PCGII_top15_lambda1")]
    PCGII_top15_lambda2.ordered<-results[order(abs(results$tscore_PCGII_top15_lambda2), decreasing = T),c("truth","Est_PCGII_top15_lambda2","tscore_PCGII_top15_lambda2")]
    PCGII_top15_lambda3.ordered<-results[order(abs(results$tscore_PCGII_top15_lambda3), decreasing = T),c("truth","Est_PCGII_top15_lambda3","tscore_PCGII_top15_lambda3")]
    
    PCGII_top25_lambda1.ordered<-results[order(abs(results$tscore_PCGII_top25_lambda1), decreasing = T),c("truth","Est_PCGII_top25_lambda1","tscore_PCGII_top25_lambda1")]
    PCGII_top25_lambda2.ordered<-results[order(abs(results$tscore_PCGII_top25_lambda2), decreasing = T),c("truth","Est_PCGII_top25_lambda2","tscore_PCGII_top25_lambda2")]
    PCGII_top25_lambda3.ordered<-results[order(abs(results$tscore_PCGII_top25_lambda3), decreasing = T),c("truth","Est_PCGII_top25_lambda3","tscore_PCGII_top25_lambda3")]
    
    PCGII_random10_lambda1.ordered<-results[order(abs(results$tscore_PCGII_random10_lambda1), decreasing = T),c("truth","Est_PCGII_random10_lambda1","tscore_PCGII_random10_lambda1")]
    PCGII_random10_lambda2.ordered<-results[order(abs(results$tscore_PCGII_random10_lambda2), decreasing = T),c("truth","Est_PCGII_random10_lambda2","tscore_PCGII_random10_lambda2")]
    PCGII_random10_lambda3.ordered<-results[order(abs(results$tscore_PCGII_random10_lambda3), decreasing = T),c("truth","Est_PCGII_random10_lambda3","tscore_PCGII_random10_lambda3")]
    
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
      
      cLevel_lambda1.ordered.TPR[loop, k]=sum(cLevel_lambda1.ordered[1:loop,]$truth!=0)/CP 
      cLevel_lambda1.ordered.FPR[loop, k]=sum(cLevel_lambda1.ordered[1:loop,]$truth==0)/CN 
      cLevel_lambda1.ordered.PPV[loop, k]=sum(cLevel_lambda1.ordered[1:loop,]$truth!=0)/loop
      cLevel_lambda1.ordered.FDR[loop, k]=sum(cLevel_lambda1.ordered[1:loop,]$truth==0)/loop
      
      cLevel_lambda2.ordered.TPR[loop, k]=sum(cLevel_lambda2.ordered[1:loop,]$truth!=0)/CP 
      cLevel_lambda2.ordered.FPR[loop, k]=sum(cLevel_lambda2.ordered[1:loop,]$truth==0)/CN 
      cLevel_lambda2.ordered.PPV[loop, k]=sum(cLevel_lambda2.ordered[1:loop,]$truth!=0)/loop
      cLevel_lambda2.ordered.FDR[loop, k]=sum(cLevel_lambda2.ordered[1:loop,]$truth==0)/loop
      
      cLevel_lambda3.ordered.TPR[loop, k]=sum(cLevel_lambda3.ordered[1:loop,]$truth!=0)/CP 
      cLevel_lambda3.ordered.FPR[loop, k]=sum(cLevel_lambda3.ordered[1:loop,]$truth==0)/CN 
      cLevel_lambda3.ordered.PPV[loop, k]=sum(cLevel_lambda3.ordered[1:loop,]$truth!=0)/loop
      cLevel_lambda3.ordered.FDR[loop, k]=sum(cLevel_lambda3.ordered[1:loop,]$truth==0)/loop
      
      PCGII_all_lambda1.ordered.TPR[loop, k]=sum(PCGII_all_lambda1.ordered[1:loop,]$truth!=0)/CP 
      PCGII_all_lambda1.ordered.FPR[loop, k]=sum(PCGII_all_lambda1.ordered[1:loop,]$truth==0)/CN 
      PCGII_all_lambda1.ordered.PPV[loop, k]=sum(PCGII_all_lambda1.ordered[1:loop,]$truth!=0)/loop
      PCGII_all_lambda1.ordered.FDR[loop, k]=sum(PCGII_all_lambda1.ordered[1:loop,]$truth==0)/loop
      
      PCGII_all_lambda2.ordered.TPR[loop, k]=sum(PCGII_all_lambda2.ordered[1:loop,]$truth!=0)/CP 
      PCGII_all_lambda2.ordered.FPR[loop, k]=sum(PCGII_all_lambda2.ordered[1:loop,]$truth==0)/CN 
      PCGII_all_lambda2.ordered.PPV[loop, k]=sum(PCGII_all_lambda2.ordered[1:loop,]$truth!=0)/loop
      PCGII_all_lambda2.ordered.FDR[loop, k]=sum(PCGII_all_lambda2.ordered[1:loop,]$truth==0)/loop
      
      PCGII_all_lambda3.ordered.TPR[loop, k]=sum(PCGII_all_lambda3.ordered[1:loop,]$truth!=0)/CP 
      PCGII_all_lambda3.ordered.FPR[loop, k]=sum(PCGII_all_lambda3.ordered[1:loop,]$truth==0)/CN 
      PCGII_all_lambda3.ordered.PPV[loop, k]=sum(PCGII_all_lambda3.ordered[1:loop,]$truth!=0)/loop
      PCGII_all_lambda3.ordered.FDR[loop, k]=sum(PCGII_all_lambda3.ordered[1:loop,]$truth==0)/loop
      
      PCGII_top5_lambda1.ordered.TPR[loop, k]=sum(PCGII_top5_lambda1.ordered[1:loop,]$truth!=0)/CP 
      PCGII_top5_lambda1.ordered.FPR[loop, k]=sum(PCGII_top5_lambda1.ordered[1:loop,]$truth==0)/CN 
      PCGII_top5_lambda1.ordered.PPV[loop, k]=sum(PCGII_top5_lambda1.ordered[1:loop,]$truth!=0)/loop
      PCGII_top5_lambda1.ordered.FDR[loop, k]=sum(PCGII_top5_lambda1.ordered[1:loop,]$truth==0)/loop
      
      PCGII_top5_lambda2.ordered.TPR[loop, k]=sum(PCGII_top5_lambda2.ordered[1:loop,]$truth!=0)/CP 
      PCGII_top5_lambda2.ordered.FPR[loop, k]=sum(PCGII_top5_lambda2.ordered[1:loop,]$truth==0)/CN 
      PCGII_top5_lambda2.ordered.PPV[loop, k]=sum(PCGII_top5_lambda2.ordered[1:loop,]$truth!=0)/loop
      PCGII_top5_lambda2.ordered.FDR[loop, k]=sum(PCGII_top5_lambda2.ordered[1:loop,]$truth==0)/loop
      
      PCGII_top5_lambda3.ordered.TPR[loop, k]=sum(PCGII_top5_lambda3.ordered[1:loop,]$truth!=0)/CP 
      PCGII_top5_lambda3.ordered.FPR[loop, k]=sum(PCGII_top5_lambda3.ordered[1:loop,]$truth==0)/CN 
      PCGII_top5_lambda3.ordered.PPV[loop, k]=sum(PCGII_top5_lambda3.ordered[1:loop,]$truth!=0)/loop
      PCGII_top5_lambda3.ordered.FDR[loop, k]=sum(PCGII_top5_lambda3.ordered[1:loop,]$truth==0)/loop
      
      PCGII_top15_lambda1.ordered.TPR[loop, k]=sum(PCGII_top15_lambda1.ordered[1:loop,]$truth!=0)/CP 
      PCGII_top15_lambda1.ordered.FPR[loop, k]=sum(PCGII_top15_lambda1.ordered[1:loop,]$truth==0)/CN 
      PCGII_top15_lambda1.ordered.PPV[loop, k]=sum(PCGII_top15_lambda1.ordered[1:loop,]$truth!=0)/loop
      PCGII_top15_lambda1.ordered.FDR[loop, k]=sum(PCGII_top15_lambda1.ordered[1:loop,]$truth==0)/loop
      
      PCGII_top15_lambda2.ordered.TPR[loop, k]=sum(PCGII_top15_lambda2.ordered[1:loop,]$truth!=0)/CP 
      PCGII_top15_lambda2.ordered.FPR[loop, k]=sum(PCGII_top15_lambda2.ordered[1:loop,]$truth==0)/CN 
      PCGII_top15_lambda2.ordered.PPV[loop, k]=sum(PCGII_top15_lambda2.ordered[1:loop,]$truth!=0)/loop
      PCGII_top15_lambda2.ordered.FDR[loop, k]=sum(PCGII_top15_lambda2.ordered[1:loop,]$truth==0)/loop
      
      PCGII_top15_lambda3.ordered.TPR[loop, k]=sum(PCGII_top15_lambda3.ordered[1:loop,]$truth!=0)/CP 
      PCGII_top15_lambda3.ordered.FPR[loop, k]=sum(PCGII_top15_lambda3.ordered[1:loop,]$truth==0)/CN 
      PCGII_top15_lambda3.ordered.PPV[loop, k]=sum(PCGII_top15_lambda3.ordered[1:loop,]$truth!=0)/loop
      PCGII_top15_lambda3.ordered.FDR[loop, k]=sum(PCGII_top15_lambda3.ordered[1:loop,]$truth==0)/loop
      
      PCGII_top25_lambda1.ordered.TPR[loop, k]=sum(PCGII_top25_lambda1.ordered[1:loop,]$truth!=0)/CP 
      PCGII_top25_lambda1.ordered.FPR[loop, k]=sum(PCGII_top25_lambda1.ordered[1:loop,]$truth==0)/CN 
      PCGII_top25_lambda1.ordered.PPV[loop, k]=sum(PCGII_top25_lambda1.ordered[1:loop,]$truth!=0)/loop
      PCGII_top25_lambda1.ordered.FDR[loop, k]=sum(PCGII_top25_lambda1.ordered[1:loop,]$truth==0)/loop
      
      PCGII_top25_lambda2.ordered.TPR[loop, k]=sum(PCGII_top25_lambda2.ordered[1:loop,]$truth!=0)/CP 
      PCGII_top25_lambda2.ordered.FPR[loop, k]=sum(PCGII_top25_lambda2.ordered[1:loop,]$truth==0)/CN 
      PCGII_top25_lambda2.ordered.PPV[loop, k]=sum(PCGII_top25_lambda2.ordered[1:loop,]$truth!=0)/loop
      PCGII_top25_lambda2.ordered.FDR[loop, k]=sum(PCGII_top25_lambda2.ordered[1:loop,]$truth==0)/loop
      
      PCGII_top25_lambda3.ordered.TPR[loop, k]=sum(PCGII_top25_lambda3.ordered[1:loop,]$truth!=0)/CP 
      PCGII_top25_lambda3.ordered.FPR[loop, k]=sum(PCGII_top25_lambda3.ordered[1:loop,]$truth==0)/CN 
      PCGII_top25_lambda3.ordered.PPV[loop, k]=sum(PCGII_top25_lambda3.ordered[1:loop,]$truth!=0)/loop
      PCGII_top25_lambda3.ordered.FDR[loop, k]=sum(PCGII_top25_lambda3.ordered[1:loop,]$truth==0)/loop
      
      PCGII_random10_lambda1.ordered.TPR[loop, k]=sum(PCGII_random10_lambda1.ordered[1:loop,]$truth!=0)/CP 
      PCGII_random10_lambda1.ordered.FPR[loop, k]=sum(PCGII_random10_lambda1.ordered[1:loop,]$truth==0)/CN 
      PCGII_random10_lambda1.ordered.PPV[loop, k]=sum(PCGII_random10_lambda1.ordered[1:loop,]$truth!=0)/loop
      PCGII_random10_lambda1.ordered.FDR[loop, k]=sum(PCGII_random10_lambda1.ordered[1:loop,]$truth==0)/loop
      
      PCGII_random10_lambda2.ordered.TPR[loop, k]=sum(PCGII_random10_lambda2.ordered[1:loop,]$truth!=0)/CP 
      PCGII_random10_lambda2.ordered.FPR[loop, k]=sum(PCGII_random10_lambda2.ordered[1:loop,]$truth==0)/CN 
      PCGII_random10_lambda2.ordered.PPV[loop, k]=sum(PCGII_random10_lambda2.ordered[1:loop,]$truth!=0)/loop
      PCGII_random10_lambda2.ordered.FDR[loop, k]=sum(PCGII_random10_lambda2.ordered[1:loop,]$truth==0)/loop
      
      PCGII_random10_lambda3.ordered.TPR[loop, k]=sum(PCGII_random10_lambda3.ordered[1:loop,]$truth!=0)/CP 
      PCGII_random10_lambda3.ordered.FPR[loop, k]=sum(PCGII_random10_lambda3.ordered[1:loop,]$truth==0)/CN 
      PCGII_random10_lambda3.ordered.PPV[loop, k]=sum(PCGII_random10_lambda3.ordered[1:loop,]$truth!=0)/loop
      PCGII_random10_lambda3.ordered.FDR[loop, k]=sum(PCGII_random10_lambda3.ordered[1:loop,]$truth==0)/loop
  
      
      FGGM_lambda1.ordered.TPR[loop, k]=sum(FGGM_lambda1.ordered[1:loop,]$truth!=0)/CP 
      FGGM_lambda1.ordered.FPR[loop, k]=sum(FGGM_lambda1.ordered[1:loop,]$truth==0)/CN 
      FGGM_lambda1.ordered.PPV[loop, k]=sum(FGGM_lambda1.ordered[1:loop,]$truth!=0)/loop
      FGGM_lambda1.ordered.FDR[loop, k]=sum(FGGM_lambda1.ordered[1:loop,]$truth==0)/loop
      
      FGGM_lambda2.ordered.TPR[loop, k]=sum(FGGM_lambda2.ordered[1:loop,]$truth!=0)/CP 
      FGGM_lambda2.ordered.FPR[loop, k]=sum(FGGM_lambda2.ordered[1:loop,]$truth==0)/CN 
      FGGM_lambda2.ordered.PPV[loop, k]=sum(FGGM_lambda2.ordered[1:loop,]$truth!=0)/loop
      FGGM_lambda2.ordered.FDR[loop, k]=sum(FGGM_lambda2.ordered[1:loop,]$truth==0)/loop
      
      FGGM_lambda3.ordered.TPR[loop, k]=sum(FGGM_lambda3.ordered[1:loop,]$truth!=0)/CP 
      FGGM_lambda3.ordered.FPR[loop, k]=sum(FGGM_lambda3.ordered[1:loop,]$truth==0)/CN 
      FGGM_lambda3.ordered.PPV[loop, k]=sum(FGGM_lambda3.ordered[1:loop,]$truth!=0)/loop
      FGGM_lambda3.ordered.FDR[loop, k]=sum(FGGM_lambda3.ordered[1:loop,]$truth==0)/loop
    }
    
    for (c in 1:1000){
      tempstore.ENF[c,k]=max(ENF.ordered.TPR[ENF.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.MLE[c,k]=max(MLE.ordered.TPR[MLE.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.cLevel_lambda1[c,k]=max(cLevel_lambda1.ordered.TPR[cLevel_lambda1.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.cLevel_lambda2[c,k]=max(cLevel_lambda2.ordered.TPR[cLevel_lambda2.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.cLevel_lambda3[c,k]=max(cLevel_lambda3.ordered.TPR[cLevel_lambda3.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      
      tempstore.PCGII_all_lambda1[c,k]=max(PCGII_all_lambda1.ordered.TPR[PCGII_all_lambda1.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_all_lambda2[c,k]=max(PCGII_all_lambda2.ordered.TPR[PCGII_all_lambda2.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_all_lambda3[c,k]=max(PCGII_all_lambda3.ordered.TPR[PCGII_all_lambda3.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_top5_lambda1[c,k]=max(PCGII_top5_lambda1.ordered.TPR[PCGII_top5_lambda1.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_top5_lambda2[c,k]=max(PCGII_top5_lambda2.ordered.TPR[PCGII_top5_lambda2.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_top5_lambda3[c,k]=max(PCGII_top5_lambda3.ordered.TPR[PCGII_top5_lambda3.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_top15_lambda1[c,k]=max(PCGII_top15_lambda1.ordered.TPR[PCGII_top15_lambda1.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_top15_lambda2[c,k]=max(PCGII_top15_lambda2.ordered.TPR[PCGII_top15_lambda2.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_top15_lambda3[c,k]=max(PCGII_top15_lambda3.ordered.TPR[PCGII_top15_lambda3.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_top25_lambda1[c,k]=max(PCGII_top25_lambda1.ordered.TPR[PCGII_top25_lambda1.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_top25_lambda2[c,k]=max(PCGII_top25_lambda2.ordered.TPR[PCGII_top25_lambda2.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_top25_lambda3[c,k]=max(PCGII_top25_lambda3.ordered.TPR[PCGII_top25_lambda3.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_random10_lambda1[c,k]=max(PCGII_random10_lambda1.ordered.TPR[PCGII_random10_lambda1.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_random10_lambda2[c,k]=max(PCGII_random10_lambda2.ordered.TPR[PCGII_random10_lambda2.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.PCGII_random10_lambda3[c,k]=max(PCGII_random10_lambda3.ordered.TPR[PCGII_random10_lambda3.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      
      
      tempstore.FGGM_lambda1[c,k]=max(FGGM_lambda1.ordered.TPR[FGGM_lambda1.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.FGGM_lambda2[c,k]=max(FGGM_lambda2.ordered.TPR[FGGM_lambda2.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
      tempstore.FGGM_lambda3[c,k]=max(FGGM_lambda3.ordered.TPR[FGGM_lambda3.ordered.FPR[,k]<=cutt[c],k]) # empirical TPR
    }
    print("ROC done")
    
    #### FDR 
    print("FDR starts")
    for (a in 1:length(al)){
      temp=inference(cLevel_lambda1, alpha = al[a])$sigs # c=0
      sigs_cLevel_lambda1=sigs2vec(temp, p) # significant edges
      
      temp=inference(cLevel_lambda2, alpha = al[a])$sigs # c=0
      sigs_cLevel_lambda2=sigs2vec(temp, p) # significant edges
      
      temp=inference(cLevel_lambda3, alpha = al[a])$sigs # c=0
      sigs_cLevel_lambda3=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_all_lambda1, alpha = al[a])$sigs # c=0
      sigs_PCGII_all_lambda1=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_all_lambda2, alpha = al[a])$sigs # c=0
      sigs_PCGII_all_lambda2=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_all_lambda3, alpha = al[a])$sigs # c=0
      sigs_PCGII_all_lambda3=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_top5_lambda1, alpha = al[a])$sigs # c=0
      sigs_PCGII_top5_lambda1=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_top5_lambda2, alpha = al[a])$sigs # c=0
      sigs_PCGII_top5_lambda2=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_top5_lambda3, alpha = al[a])$sigs # c=0
      sigs_PCGII_top5_lambda3=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_top15_lambda1, alpha = al[a])$sigs # c=0
      sigs_PCGII_top15_lambda1=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_top15_lambda2, alpha = al[a])$sigs # c=0
      sigs_PCGII_top15_lambda2=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_top15_lambda3, alpha = al[a])$sigs # c=0
      sigs_PCGII_top15_lambda3=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_top25_lambda1, alpha = al[a])$sigs # c=0
      sigs_PCGII_top25_lambda1=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_top25_lambda2, alpha = al[a])$sigs # c=0
      sigs_PCGII_top25_lambda2=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_top25_lambda3, alpha = al[a])$sigs # c=0
      sigs_PCGII_top25_lambda3=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_random10_lambda1, alpha = al[a])$sigs # c=0
      sigs_PCGII_random10_lambda1=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_random10_lambda2, alpha = al[a])$sigs # c=0
      sigs_PCGII_random10_lambda2=sigs2vec(temp, p) # significant edges
      
      temp=inference(PCGII_random10_lambda3, alpha = al[a])$sigs # c=0
      sigs_PCGII_random10_lambda3=sigs2vec(temp, p) # significant edges
      
      # Number of total selected edges
      all.enf.p[a,k]=sum(results$ENF_p<=al[a])
      all.enf.q[a,k]=sum(results$ENF_q<=al[a]) 
      all.mle.p[a,k]=sum(results$ShrunkMLE_p<=al[a]) 
      all.mle.q[a,k]=sum(results$ShrunkMLE_q<=al[a])
      all.cLevel_lambda1[a,k]=sum(sigs_cLevel_lambda1==1)
      all.cLevel_lambda2[a,k]=sum(sigs_cLevel_lambda2==1)
      all.cLevel_lambda3[a,k]=sum(sigs_cLevel_lambda3==1)
      
      all.PCGII_all_lambda1[a,k]=sum(sigs_PCGII_all_lambda1==1)
      all.PCGII_all_lambda2[a,k]=sum(sigs_PCGII_all_lambda2==1)
      all.PCGII_all_lambda3[a,k]=sum(sigs_PCGII_all_lambda3==1)
      
      all.PCGII_top5_lambda1[a,k]=sum(sigs_PCGII_top5_lambda1==1)
      all.PCGII_top5_lambda2[a,k]=sum(sigs_PCGII_top5_lambda2==1)
      all.PCGII_top5_lambda3[a,k]=sum(sigs_PCGII_top5_lambda3==1)
      
      all.PCGII_top15_lambda1[a,k]=sum(sigs_PCGII_top15_lambda1==1)
      all.PCGII_top15_lambda2[a,k]=sum(sigs_PCGII_top15_lambda2==1)
      all.PCGII_top15_lambda3[a,k]=sum(sigs_PCGII_top15_lambda3==1)
      
      all.PCGII_top25_lambda1[a,k]=sum(sigs_PCGII_top25_lambda1==1)
      all.PCGII_top25_lambda2[a,k]=sum(sigs_PCGII_top25_lambda2==1)
      all.PCGII_top25_lambda3[a,k]=sum(sigs_PCGII_top25_lambda3==1)
      
      all.PCGII_random10_lambda1[a,k]=sum(sigs_PCGII_random10_lambda1==1)
      all.PCGII_random10_lambda2[a,k]=sum(sigs_PCGII_random10_lambda2==1)
      all.PCGII_random10_lambda3[a,k]=sum(sigs_PCGII_random10_lambda3==1)
      
      
      all.FGGM_lambda1[a,k]=sum(results$FGGM_lambda1_qval<=al[a])
      all.FGGM_lambda2[a,k]=sum(results$FGGM_lambda2_qval<=al[a])
      all.FGGM_lambda3[a,k]=sum(results$FGGM_lambda3_qval<=al[a])
      
      # true positives
      tp.enf.p[a,k]=sum(which(results$ENF_p<=al[a]) %in% TP)
      tp.enf.q[a,k]=sum(which(results$ENF_q<=al[a]) %in% TP)
      tp.mle.p[a,k]=sum(which(results$ShrunkMLE_p<=al[a]) %in% TP)
      tp.mle.q[a,k]=sum(which(results$ShrunkMLE_q<=al[a]) %in% TP)
      tp.cLevel_lambda1[a,k]=sum(which(sigs_cLevel_lambda1==1) %in% TP)
      tp.cLevel_lambda2[a,k]=sum(which(sigs_cLevel_lambda2==1) %in% TP)
      tp.cLevel_lambda3[a,k]=sum(which(sigs_cLevel_lambda3==1) %in% TP)
      
      tp.PCGII_all_lambda1[a,k]=sum(which(sigs_PCGII_all_lambda1==1) %in% TP)
      tp.PCGII_all_lambda2[a,k]=sum(which(sigs_PCGII_all_lambda2==1) %in% TP)
      tp.PCGII_all_lambda3[a,k]=sum(which(sigs_PCGII_all_lambda3==1) %in% TP)
      tp.PCGII_top5_lambda1[a,k]=sum(which(sigs_PCGII_top5_lambda1==1) %in% TP)
      tp.PCGII_top5_lambda2[a,k]=sum(which(sigs_PCGII_top5_lambda2==1) %in% TP)
      tp.PCGII_top5_lambda3[a,k]=sum(which(sigs_PCGII_top5_lambda3==1) %in% TP)
      tp.PCGII_top15_lambda1[a,k]=sum(which(sigs_PCGII_top15_lambda1==1) %in% TP)
      tp.PCGII_top15_lambda2[a,k]=sum(which(sigs_PCGII_top15_lambda2==1) %in% TP)
      tp.PCGII_top15_lambda3[a,k]=sum(which(sigs_PCGII_top15_lambda3==1) %in% TP)
      tp.PCGII_top25_lambda1[a,k]=sum(which(sigs_PCGII_top25_lambda1==1) %in% TP)
      tp.PCGII_top25_lambda2[a,k]=sum(which(sigs_PCGII_top25_lambda2==1) %in% TP)
      tp.PCGII_top25_lambda3[a,k]=sum(which(sigs_PCGII_top25_lambda3==1) %in% TP)
      tp.PCGII_random10_lambda1[a,k]=sum(which(sigs_PCGII_random10_lambda1==1) %in% TP)
      tp.PCGII_random10_lambda2[a,k]=sum(which(sigs_PCGII_random10_lambda2==1) %in% TP)
      tp.PCGII_random10_lambda3[a,k]=sum(which(sigs_PCGII_random10_lambda3==1) %in% TP)
      
      tp.FGGM_lambda1[a,k]=sum(which(results$FGGM_lambda1_qval<=al[a]) %in% TP)
      tp.FGGM_lambda2[a,k]=sum(which(results$FGGM_lambda2_qval<=al[a]) %in% TP)
      tp.FGGM_lambda3[a,k]=sum(which(results$FGGM_lambda3_qval<=al[a]) %in% TP)
      
      # false positives
      fp.enf.p[a,k]=sum(!which(results$ENF_p<=al[a]) %in% TP)
      fp.enf.q[a,k]=sum(!which(results$ENF_q<=al[a]) %in% TP)
      fp.mle.p[a,k]=sum(!which(results$ShrunkMLE_p<=al[a]) %in% TP)
      fp.mle.q[a,k]=sum(!which(results$ShrunkMLE_q<=al[a]) %in% TP)
      fp.cLevel_lambda1[a,k]=sum(!which(sigs_cLevel_lambda1==1) %in% TP)
      fp.cLevel_lambda2[a,k]=sum(!which(sigs_cLevel_lambda2==1) %in% TP)
      fp.cLevel_lambda3[a,k]=sum(!which(sigs_cLevel_lambda3==1) %in% TP)

      fp.PCGII_all_lambda1[a,k]=sum(!which(sigs_PCGII_all_lambda1==1) %in% TP)
      fp.PCGII_all_lambda2[a,k]=sum(!which(sigs_PCGII_all_lambda2==1) %in% TP)
      fp.PCGII_all_lambda3[a,k]=sum(!which(sigs_PCGII_all_lambda3==1) %in% TP)
      fp.PCGII_top5_lambda1[a,k]=sum(!which(sigs_PCGII_top5_lambda1==1) %in% TP)
      fp.PCGII_top5_lambda2[a,k]=sum(!which(sigs_PCGII_top5_lambda2==1) %in% TP)
      fp.PCGII_top5_lambda3[a,k]=sum(!which(sigs_PCGII_top5_lambda3==1) %in% TP)
      fp.PCGII_top15_lambda1[a,k]=sum(!which(sigs_PCGII_top15_lambda1==1) %in% TP)
      fp.PCGII_top15_lambda2[a,k]=sum(!which(sigs_PCGII_top15_lambda2==1) %in% TP)
      fp.PCGII_top15_lambda3[a,k]=sum(!which(sigs_PCGII_top15_lambda3==1) %in% TP)
      fp.PCGII_top25_lambda1[a,k]=sum(!which(sigs_PCGII_top25_lambda1==1) %in% TP)
      fp.PCGII_top25_lambda2[a,k]=sum(!which(sigs_PCGII_top25_lambda2==1) %in% TP)
      fp.PCGII_top25_lambda3[a,k]=sum(!which(sigs_PCGII_top25_lambda3==1) %in% TP)
      fp.PCGII_random10_lambda1[a,k]=sum(!which(sigs_PCGII_random10_lambda1==1) %in% TP)
      fp.PCGII_random10_lambda2[a,k]=sum(!which(sigs_PCGII_random10_lambda2==1) %in% TP)
      fp.PCGII_random10_lambda3[a,k]=sum(!which(sigs_PCGII_random10_lambda3==1) %in% TP)

      fp.FGGM_lambda1[a,k]=sum(!which(results$FGGM_lambda1_qval<=al[a]) %in% TP)
      fp.FGGM_lambda2[a,k]=sum(!which(results$FGGM_lambda2_qval<=al[a]) %in% TP)
      fp.FGGM_lambda3[a,k]=sum(!which(results$FGGM_lambda3_qval<=al[a]) %in% TP)
      
      # true negatives
      tn.enf.p[a,k]=sum(!which(results$ENF_p>al[a]) %in% TP)
      tn.enf.q[a,k]=sum(!which(results$ENF.qval>al[a]) %in% TP)
      tn.mle.p[a,k]=sum(!which(results$ShrunkMLE_p>al[a]) %in% TP)
      tn.mle.q[a,k]=sum(!which(results$ShrunkMLE_q>al[a]) %in% TP)
      tn.cLevel_lambda1[a,k]=sum(!which(sigs_cLevel_lambda1!=1) %in% TP)
      tn.cLevel_lambda2[a,k]=sum(!which(sigs_cLevel_lambda2!=1) %in% TP)
      tn.cLevel_lambda3[a,k]=sum(!which(sigs_cLevel_lambda3!=1) %in% TP)
      
      tn.PCGII_all_lambda1[a,k]=sum(!which(sigs_PCGII_all_lambda1!=1) %in% TP)
      tn.PCGII_all_lambda2[a,k]=sum(!which(sigs_PCGII_all_lambda2!=1) %in% TP)
      tn.PCGII_all_lambda3[a,k]=sum(!which(sigs_PCGII_all_lambda3!=1) %in% TP)
      tn.PCGII_top5_lambda1[a,k]=sum(!which(sigs_PCGII_top5_lambda1!=1) %in% TP)
      tn.PCGII_top5_lambda2[a,k]=sum(!which(sigs_PCGII_top5_lambda2!=1) %in% TP)
      tn.PCGII_top5_lambda3[a,k]=sum(!which(sigs_PCGII_top5_lambda3!=1) %in% TP)
      tn.PCGII_top15_lambda1[a,k]=sum(!which(sigs_PCGII_top15_lambda1!=1) %in% TP)
      tn.PCGII_top15_lambda2[a,k]=sum(!which(sigs_PCGII_top15_lambda2!=1) %in% TP)
      tn.PCGII_top15_lambda3[a,k]=sum(!which(sigs_PCGII_top15_lambda3!=1) %in% TP)
      tn.PCGII_top25_lambda1[a,k]=sum(!which(sigs_PCGII_top25_lambda1!=1) %in% TP)
      tn.PCGII_top25_lambda2[a,k]=sum(!which(sigs_PCGII_top25_lambda2!=1) %in% TP)
      tn.PCGII_top25_lambda3[a,k]=sum(!which(sigs_PCGII_top25_lambda3!=1) %in% TP)
      tn.PCGII_random10_lambda1[a,k]=sum(!which(sigs_PCGII_random10_lambda1!=1) %in% TP)
      tn.PCGII_random10_lambda2[a,k]=sum(!which(sigs_PCGII_random10_lambda2!=1) %in% TP)
      tn.PCGII_random10_lambda3[a,k]=sum(!which(sigs_PCGII_random10_lambda3!=1) %in% TP)
      
      tn.FGGM_lambda1[a,k]=sum(!which(results$FGGM_lambda1_qval>al[a]) %in% TP)
      tn.FGGM_lambda2[a,k]=sum(!which(results$FGGM_lambda2_qval>al[a]) %in% TP)
      tn.FGGM_lambda3[a,k]=sum(!which(results$FGGM_lambda3_qval>al[a]) %in% TP)
      
      # false negatives
      fn.enf.p[a,k]=sum(which(results$ENF_p>al[a]) %in% TP)
      fn.enf.q[a,k]=sum(which(results$ENF_q>al[a]) %in% TP)
      fn.mle.p[a,k]=sum(which(results$ShrunkMLE_p>al[a]) %in% TP)
      fn.mle.q[a,k]=sum(which(results$ShrunkMLE_q>al[a]) %in% TP)
      fn.cLevel_lambda1[a,k]=sum(which(sigs_cLevel_lambda1!=1) %in% TP)
      fn.cLevel_lambda2[a,k]=sum(which(sigs_cLevel_lambda2!=1) %in% TP)
      fn.cLevel_lambda3[a,k]=sum(which(sigs_cLevel_lambda3!=1) %in% TP)
      
      fn.PCGII_all_lambda1[a,k]=sum(which(sigs_PCGII_all_lambda1!=1) %in% TP)
      fn.PCGII_all_lambda2[a,k]=sum(which(sigs_PCGII_all_lambda2!=1) %in% TP)
      fn.PCGII_all_lambda3[a,k]=sum(which(sigs_PCGII_all_lambda3!=1) %in% TP)
      fn.PCGII_top5_lambda1[a,k]=sum(which(sigs_PCGII_top5_lambda1!=1) %in% TP)
      fn.PCGII_top5_lambda2[a,k]=sum(which(sigs_PCGII_top5_lambda2!=1) %in% TP)
      fn.PCGII_top5_lambda3[a,k]=sum(which(sigs_PCGII_top5_lambda3!=1) %in% TP)
      fn.PCGII_top15_lambda1[a,k]=sum(which(sigs_PCGII_top15_lambda1!=1) %in% TP)
      fn.PCGII_top15_lambda2[a,k]=sum(which(sigs_PCGII_top15_lambda2!=1) %in% TP)
      fn.PCGII_top15_lambda3[a,k]=sum(which(sigs_PCGII_top15_lambda3!=1) %in% TP)
      fn.PCGII_top25_lambda1[a,k]=sum(which(sigs_PCGII_top25_lambda1!=1) %in% TP)
      fn.PCGII_top25_lambda2[a,k]=sum(which(sigs_PCGII_top25_lambda2!=1) %in% TP)
      fn.PCGII_top25_lambda3[a,k]=sum(which(sigs_PCGII_top25_lambda3!=1) %in% TP)
      fn.PCGII_random10_lambda1[a,k]=sum(which(sigs_PCGII_random10_lambda1!=1) %in% TP)
      fn.PCGII_random10_lambda2[a,k]=sum(which(sigs_PCGII_random10_lambda2!=1) %in% TP)
      fn.PCGII_random10_lambda3[a,k]=sum(which(sigs_PCGII_random10_lambda3!=1) %in% TP)
      
      fn.FGGM_lambda1[a,k]=sum(which(results$FGGM_lambda1_qval>al[a]) %in% TP)
      fn.FGGM_lambda2[a,k]=sum(which(results$FGGM_lambda2_qval>al[a]) %in% TP)
      fn.FGGM_lambda3[a,k]=sum(which(results$FGGM_lambda3_qval>al[a]) %in% TP)
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
    rank.cLevel_lambda1=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                      TPR=rowMeans(cLevel_lambda1.ordered.TPR),
                                      FPR=rowMeans(cLevel_lambda1.ordered.FPR),
                                      PPV=rowMeans(cLevel_lambda1.ordered.PPV),
                                      FDR=rowMeans(cLevel_lambda1.ordered.FDR))
    rank.cLevel_lambda2=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                    TPR=rowMeans(cLevel_lambda2.ordered.TPR),
                                    FPR=rowMeans(cLevel_lambda2.ordered.FPR),
                                    PPV=rowMeans(cLevel_lambda2.ordered.PPV),
                                    FDR=rowMeans(cLevel_lambda2.ordered.FDR))
    rank.cLevel_lambda3=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                         TPR=rowMeans(cLevel_lambda3.ordered.TPR),
                                         FPR=rowMeans(cLevel_lambda3.ordered.FPR),
                                         PPV=rowMeans(cLevel_lambda3.ordered.PPV),
                                         FDR=rowMeans(cLevel_lambda3.ordered.FDR))
    rank.PCGII_all_lambda1=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                     TPR=rowMeans(PCGII_all_lambda1.ordered.TPR),
                                     FPR=rowMeans(PCGII_all_lambda1.ordered.FPR),
                                     PPV=rowMeans(PCGII_all_lambda1.ordered.PPV),
                                     FDR=rowMeans(PCGII_all_lambda1.ordered.FDR))
    rank.PCGII_all_lambda2=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                   TPR=rowMeans(PCGII_all_lambda2.ordered.TPR),
                                   FPR=rowMeans(PCGII_all_lambda2.ordered.FPR),
                                   PPV=rowMeans(PCGII_all_lambda2.ordered.PPV),
                                   FDR=rowMeans(PCGII_all_lambda2.ordered.FDR))
    rank.PCGII_all_lambda3=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_all_lambda3.ordered.TPR),
                                            FPR=rowMeans(PCGII_all_lambda3.ordered.FPR),
                                            PPV=rowMeans(PCGII_all_lambda3.ordered.PPV),
                                            FDR=rowMeans(PCGII_all_lambda3.ordered.FDR))
    rank.PCGII_top5_lambda1=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_top5_lambda1.ordered.TPR),
                                            FPR=rowMeans(PCGII_top5_lambda1.ordered.FPR),
                                            PPV=rowMeans(PCGII_top5_lambda1.ordered.PPV),
                                            FDR=rowMeans(PCGII_top5_lambda1.ordered.FDR))
    rank.PCGII_top5_lambda2=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_top5_lambda2.ordered.TPR),
                                            FPR=rowMeans(PCGII_top5_lambda2.ordered.FPR),
                                            PPV=rowMeans(PCGII_top5_lambda2.ordered.PPV),
                                            FDR=rowMeans(PCGII_top5_lambda2.ordered.FDR))
    rank.PCGII_top5_lambda3=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_top5_lambda3.ordered.TPR),
                                            FPR=rowMeans(PCGII_top5_lambda3.ordered.FPR),
                                            PPV=rowMeans(PCGII_top5_lambda3.ordered.PPV),
                                            FDR=rowMeans(PCGII_top5_lambda3.ordered.FDR))
    rank.PCGII_top15_lambda1=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_top15_lambda1.ordered.TPR),
                                            FPR=rowMeans(PCGII_top15_lambda1.ordered.FPR),
                                            PPV=rowMeans(PCGII_top15_lambda1.ordered.PPV),
                                            FDR=rowMeans(PCGII_top15_lambda1.ordered.FDR))
    rank.PCGII_top15_lambda2=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_top15_lambda2.ordered.TPR),
                                            FPR=rowMeans(PCGII_top15_lambda2.ordered.FPR),
                                            PPV=rowMeans(PCGII_top15_lambda2.ordered.PPV),
                                            FDR=rowMeans(PCGII_top15_lambda2.ordered.FDR))
    rank.PCGII_top15_lambda3=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_top15_lambda3.ordered.TPR),
                                            FPR=rowMeans(PCGII_top15_lambda3.ordered.FPR),
                                            PPV=rowMeans(PCGII_top15_lambda3.ordered.PPV),
                                            FDR=rowMeans(PCGII_top15_lambda3.ordered.FDR))
    rank.PCGII_top25_lambda1=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_top25_lambda1.ordered.TPR),
                                            FPR=rowMeans(PCGII_top25_lambda1.ordered.FPR),
                                            PPV=rowMeans(PCGII_top25_lambda1.ordered.PPV),
                                            FDR=rowMeans(PCGII_top25_lambda1.ordered.FDR))
    rank.PCGII_top25_lambda2=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_top25_lambda2.ordered.TPR),
                                            FPR=rowMeans(PCGII_top25_lambda2.ordered.FPR),
                                            PPV=rowMeans(PCGII_top25_lambda2.ordered.PPV),
                                            FDR=rowMeans(PCGII_top25_lambda2.ordered.FDR))
    rank.PCGII_top25_lambda3=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_top25_lambda3.ordered.TPR),
                                            FPR=rowMeans(PCGII_top25_lambda3.ordered.FPR),
                                            PPV=rowMeans(PCGII_top25_lambda3.ordered.PPV),
                                            FDR=rowMeans(PCGII_top25_lambda3.ordered.FDR))
    rank.PCGII_random10_lambda1=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_random10_lambda1.ordered.TPR),
                                            FPR=rowMeans(PCGII_random10_lambda1.ordered.FPR),
                                            PPV=rowMeans(PCGII_random10_lambda1.ordered.PPV),
                                            FDR=rowMeans(PCGII_random10_lambda1.ordered.FDR))
    rank.PCGII_random10_lambda2=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_random10_lambda2.ordered.TPR),
                                            FPR=rowMeans(PCGII_random10_lambda2.ordered.FPR),
                                            PPV=rowMeans(PCGII_random10_lambda2.ordered.PPV),
                                            FDR=rowMeans(PCGII_random10_lambda2.ordered.FDR))
    rank.PCGII_random10_lambda3=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                            TPR=rowMeans(PCGII_random10_lambda3.ordered.TPR),
                                            FPR=rowMeans(PCGII_random10_lambda3.ordered.FPR),
                                            PPV=rowMeans(PCGII_random10_lambda3.ordered.PPV),
                                            FDR=rowMeans(PCGII_random10_lambda3.ordered.FDR))
    rank.FGGM_lambda1=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                TPR=rowMeans(FGGM_lambda1.ordered.TPR),
                                FPR=rowMeans(FGGM_lambda1.ordered.FPR),
                                PPV=rowMeans(FGGM_lambda1.ordered.PPV),
                                FDR=rowMeans(FGGM_lambda1.ordered.FDR))
    rank.FGGM_lambda2=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                       TPR=rowMeans(FGGM_lambda2.ordered.TPR),
                                       FPR=rowMeans(FGGM_lambda2.ordered.FPR),
                                       PPV=rowMeans(FGGM_lambda2.ordered.PPV),
                                       FDR=rowMeans(FGGM_lambda2.ordered.FDR))
    rank.FGGM_lambda3=cbind.data.frame(cut=seq(1,p*(p-1)/2,1),
                                       TPR=rowMeans(FGGM_lambda3.ordered.TPR),
                                       FPR=rowMeans(FGGM_lambda3.ordered.FPR),
                                       PPV=rowMeans(FGGM_lambda3.ordered.PPV),
                                       FDR=rowMeans(FGGM_lambda3.ordered.FDR))
    
    rank.table=cbind.data.frame(methods=c(rep("ENF",p*(p-1)/2),
                                          rep("MLE",p*(p-1)/2),
                                          rep("cLevel_lambda1",p*(p-1)/2),
                                          rep("cLevel_lambda2",p*(p-1)/2),
                                          rep("cLevel_lambda3",p*(p-1)/2),
                                          rep("PCGII_all_lambda1",p*(p-1)/2),
                                          rep("PCGII_all_lambda2",p*(p-1)/2),
                                          rep("PCGII_all_lambda3",p*(p-1)/2),
                                          rep("PCGII_top5_lambda1",p*(p-1)/2),
                                          rep("PCGII_top5_lambda2",p*(p-1)/2),
                                          rep("PCGII_top5_lambda3",p*(p-1)/2),
                                          rep("PCGII_top15_lambda1",p*(p-1)/2),
                                          rep("PCGII_top15_lambda2",p*(p-1)/2),
                                          rep("PCGII_top15_lambda3",p*(p-1)/2),
                                          rep("PCGII_top25_lambda1",p*(p-1)/2),
                                          rep("PCGII_top25_lambda2",p*(p-1)/2),
                                          rep("PCGII_top25_lambda3",p*(p-1)/2),
                                          rep("PCGII_random10_lambda1",p*(p-1)/2),
                                          rep("PCGII_random10_lambda2",p*(p-1)/2),
                                          rep("PCGII_random10_lambda3",p*(p-1)/2),
                                          rep("FGGM_lambda1",p*(p-1)/2),
                                          rep("FGGM_lambda2",p*(p-1)/2),
                                          rep("FGGM_lambda3",p*(p-1)/2)),
                                rbind.data.frame(rank.ENF, rank.MLE, 
                                                 rank.cLevel_lambda1, rank.cLevel_lambda2, rank.cLevel_lambda3, 
                                                 rank.PCGII_all_lambda1, rank.PCGII_all_lambda2, rank.PCGII_all_lambda3,
                                                 rank.PCGII_top5_lambda1, rank.PCGII_top5_lambda2, rank.PCGII_top5_lambda3,
                                                 rank.PCGII_top15_lambda1, rank.PCGII_top15_lambda2, rank.PCGII_top15_lambda3,
                                                 rank.PCGII_top25_lambda1, rank.PCGII_top25_lambda2, rank.PCGII_top25_lambda3,
                                                 rank.PCGII_random10_lambda1, rank.PCGII_random10_lambda2, rank.PCGII_random10_lambda3,
                                                 rank.FGGM_lambda1, rank.FGGM_lambda2, rank.FGGM_lambda3))
    
    
    ROC.out=cbind.data.frame(FPR=seq(0.001,1,0.001), # Sample FPR cut points
                             ENF=rowMeans(tempstore.ENF), # averaged max TPR at cutted FPR level
                             MLE=rowMeans(tempstore.MLE),
                             cLevel_lambda1=rowMeans(tempstore.cLevel_lambda1),
                             cLevel_lambda2=rowMeans(tempstore.cLevel_lambda2),
                             cLevel_lambda3=rowMeans(tempstore.cLevel_lambda3),
                             PCGII_all_lambda1=rowMeans(tempstore.PCGII_all_lambda1),
                             PCGII_all_lambda2=rowMeans(tempstore.PCGII_all_lambda2),
                             PCGII_all_lambda3=rowMeans(tempstore.PCGII_all_lambda3),
                             PCGII_top5_lambda1=rowMeans(tempstore.PCGII_top5_lambda1),
                             PCGII_top5_lambda2=rowMeans(tempstore.PCGII_top5_lambda2),
                             PCGII_top5_lambda3=rowMeans(tempstore.PCGII_top5_lambda3),
                             PCGII_top15_lambda1=rowMeans(tempstore.PCGII_top15_lambda1),
                             PCGII_top15_lambda2=rowMeans(tempstore.PCGII_top15_lambda2),
                             PCGII_top15_lambda3=rowMeans(tempstore.PCGII_top15_lambda3),
                             PCGII_top25_lambda1=rowMeans(tempstore.PCGII_top25_lambda1),
                             PCGII_top25_lambda2=rowMeans(tempstore.PCGII_top25_lambda2),
                             PCGII_top25_lambda3=rowMeans(tempstore.PCGII_top25_lambda3),
                             PCGII_random10_lambda1=rowMeans(tempstore.PCGII_random10_lambda1),
                             PCGII_random10_lambda2=rowMeans(tempstore.PCGII_random10_lambda2),
                             PCGII_random10_lambda3=rowMeans(tempstore.PCGII_random10_lambda3),
                             FGGM_lambda1=rowMeans(tempstore.FGGM_lambda1),
                             FGGM_lambda2=rowMeans(tempstore.FGGM_lambda2),
                             FGGM_lambda3=rowMeans(tempstore.FGGM_lambda3))
    #ROC.table=tidyr::gather(ROC.out, methods, TPR, ENF:cPCG_df)
    
    
    # Number of total correctly selected edges
    All.enf.p=rowMeans(all.enf.p)
    All.enf.q=rowMeans(all.enf.q)
    All.mle.p=rowMeans(all.mle.p)
    All.mle.q=rowMeans(all.mle.q)
    All.cLevel_lambda1=rowMeans(all.cLevel_lambda1)
    All.cLevel_lambda2=rowMeans(all.cLevel_lambda2)
    All.cLevel_lambda3=rowMeans(all.cLevel_lambda3)
    
    All.PCGII_all_lambda1=rowMeans(all.PCGII_all_lambda1)
    All.PCGII_all_lambda2=rowMeans(all.PCGII_all_lambda2)
    All.PCGII_all_lambda3=rowMeans(all.PCGII_all_lambda3)
    All.PCGII_top5_lambda1=rowMeans(all.PCGII_top5_lambda1)
    All.PCGII_top5_lambda2=rowMeans(all.PCGII_top5_lambda2)
    All.PCGII_top5_lambda3=rowMeans(all.PCGII_top5_lambda3)
    All.PCGII_top15_lambda1=rowMeans(all.PCGII_top15_lambda1)
    All.PCGII_top15_lambda2=rowMeans(all.PCGII_top15_lambda2)
    All.PCGII_top15_lambda3=rowMeans(all.PCGII_top15_lambda3)
    All.PCGII_top25_lambda1=rowMeans(all.PCGII_top25_lambda1)
    All.PCGII_top25_lambda2=rowMeans(all.PCGII_top25_lambda2)
    All.PCGII_top25_lambda3=rowMeans(all.PCGII_top25_lambda3)
    All.PCGII_random10_lambda1=rowMeans(all.PCGII_random10_lambda1)
    All.PCGII_random10_lambda2=rowMeans(all.PCGII_random10_lambda2)
    All.PCGII_random10_lambda3=rowMeans(all.PCGII_random10_lambda3)
    
    All.FGGM_lambda1=rowMeans(all.FGGM_lambda1)
    All.FGGM_lambda2=rowMeans(all.FGGM_lambda2)
    All.FGGM_lambda3=rowMeans(all.FGGM_lambda3)
    
    # true positives
    TP.enf.p=rowMeans(tp.enf.p)
    TP.enf.q=rowMeans(tp.enf.q)
    TP.mle.p=rowMeans(tp.mle.p)
    TP.mle.q=rowMeans(tp.mle.q)
    TP.cLevel_lambda1=rowMeans(tp.cLevel_lambda1)
    TP.cLevel_lambda2=rowMeans(tp.cLevel_lambda2)
    TP.cLevel_lambda3=rowMeans(tp.cLevel_lambda3)
    
    TP.PCGII_all_lambda1=rowMeans(tp.PCGII_all_lambda1)
    TP.PCGII_all_lambda2=rowMeans(tp.PCGII_all_lambda2)
    TP.PCGII_all_lambda3=rowMeans(tp.PCGII_all_lambda3)
    TP.PCGII_top5_lambda1=rowMeans(tp.PCGII_top5_lambda1)
    TP.PCGII_top5_lambda2=rowMeans(tp.PCGII_top5_lambda2)
    TP.PCGII_top5_lambda3=rowMeans(tp.PCGII_top5_lambda3)
    TP.PCGII_top15_lambda1=rowMeans(tp.PCGII_top15_lambda1)
    TP.PCGII_top15_lambda2=rowMeans(tp.PCGII_top15_lambda2)
    TP.PCGII_top15_lambda3=rowMeans(tp.PCGII_top15_lambda3)
    TP.PCGII_top25_lambda1=rowMeans(tp.PCGII_top25_lambda1)
    TP.PCGII_top25_lambda2=rowMeans(tp.PCGII_top25_lambda2)
    TP.PCGII_top25_lambda3=rowMeans(tp.PCGII_top25_lambda3)
    TP.PCGII_random10_lambda1=rowMeans(tp.PCGII_random10_lambda1)
    TP.PCGII_random10_lambda2=rowMeans(tp.PCGII_random10_lambda2)
    TP.PCGII_random10_lambda3=rowMeans(tp.PCGII_random10_lambda3)
    
    TP.FGGM_lambda1=rowMeans(tp.FGGM_lambda1)
    TP.FGGM_lambda2=rowMeans(tp.FGGM_lambda2)
    TP.FGGM_lambda3=rowMeans(tp.FGGM_lambda3)
    
    # false positives
    FP.enf.p=rowMeans(fp.enf.p)
    FP.enf.q=rowMeans(fp.enf.q)
    FP.mle.p=rowMeans(fp.mle.p)
    FP.mle.q=rowMeans(fp.mle.q)
    FP.cLevel_lambda1=rowMeans(fp.cLevel_lambda1)
    FP.cLevel_lambda2=rowMeans(fp.cLevel_lambda2)
    FP.cLevel_lambda3=rowMeans(fp.cLevel_lambda3)
    
    FP.PCGII_all_lambda1=rowMeans(fp.PCGII_all_lambda1)
    FP.PCGII_all_lambda2=rowMeans(fp.PCGII_all_lambda2)
    FP.PCGII_all_lambda3=rowMeans(fp.PCGII_all_lambda3)
    FP.PCGII_top5_lambda1=rowMeans(fp.PCGII_top5_lambda1)
    FP.PCGII_top5_lambda2=rowMeans(fp.PCGII_top5_lambda2)
    FP.PCGII_top5_lambda3=rowMeans(fp.PCGII_top5_lambda3)
    FP.PCGII_top15_lambda1=rowMeans(fp.PCGII_top15_lambda1)
    FP.PCGII_top15_lambda2=rowMeans(fp.PCGII_top15_lambda2)
    FP.PCGII_top15_lambda3=rowMeans(fp.PCGII_top15_lambda3)
    FP.PCGII_top25_lambda1=rowMeans(fp.PCGII_top25_lambda1)
    FP.PCGII_top25_lambda2=rowMeans(fp.PCGII_top25_lambda2)
    FP.PCGII_top25_lambda3=rowMeans(fp.PCGII_top25_lambda3)
    FP.PCGII_random10_lambda1=rowMeans(fp.PCGII_random10_lambda1)
    FP.PCGII_random10_lambda2=rowMeans(fp.PCGII_random10_lambda2)
    FP.PCGII_random10_lambda3=rowMeans(fp.PCGII_random10_lambda3)
    
    FP.FGGM_lambda1=rowMeans(fp.FGGM_lambda1)
    FP.FGGM_lambda2=rowMeans(fp.FGGM_lambda2)
    FP.FGGM_lambda3=rowMeans(fp.FGGM_lambda3)
    
    # true negatives
    TN.enf.p=rowMeans(tn.enf.p)
    TN.enf.q=rowMeans(tn.enf.q)
    TN.mle.p=rowMeans(tn.mle.p)
    TN.mle.q=rowMeans(tn.mle.q)
    TN.cLevel_lambda1=rowMeans(tn.cLevel_lambda1)
    TN.cLevel_lambda2=rowMeans(tn.cLevel_lambda2)
    TN.cLevel_lambda3=rowMeans(tn.cLevel_lambda3)
    
    TN.PCGII_all_lambda1=rowMeans(tn.PCGII_all_lambda1)
    TN.PCGII_all_lambda2=rowMeans(tn.PCGII_all_lambda2)
    TN.PCGII_all_lambda3=rowMeans(tn.PCGII_all_lambda3)
    TN.PCGII_top5_lambda1=rowMeans(tn.PCGII_top5_lambda1)
    TN.PCGII_top5_lambda2=rowMeans(tn.PCGII_top5_lambda2)
    TN.PCGII_top5_lambda3=rowMeans(tn.PCGII_top5_lambda3)
    TN.PCGII_top15_lambda1=rowMeans(tn.PCGII_top15_lambda1)
    TN.PCGII_top15_lambda2=rowMeans(tn.PCGII_top15_lambda2)
    TN.PCGII_top15_lambda3=rowMeans(tn.PCGII_top15_lambda3)
    TN.PCGII_top25_lambda1=rowMeans(tn.PCGII_top25_lambda1)
    TN.PCGII_top25_lambda2=rowMeans(tn.PCGII_top25_lambda2)
    TN.PCGII_top25_lambda3=rowMeans(tn.PCGII_top25_lambda3)
    TN.PCGII_random10_lambda1=rowMeans(tn.PCGII_random10_lambda1)
    TN.PCGII_random10_lambda2=rowMeans(tn.PCGII_random10_lambda2)
    TN.PCGII_random10_lambda3=rowMeans(tn.PCGII_random10_lambda3)
    
    TN.FGGM_lambda1=rowMeans(tn.FGGM_lambda1)
    TN.FGGM_lambda2=rowMeans(tn.FGGM_lambda2)
    TN.FGGM_lambda3=rowMeans(tn.FGGM_lambda3)
    
    # false negatives
    FN.enf.p=rowMeans(fn.enf.p)
    FN.enf.q=rowMeans(fn.enf.q)
    FN.mle.p=rowMeans(fn.mle.p)
    FN.mle.q=rowMeans(fn.mle.q)
    FN.cLevel_lambda1=rowMeans(fn.cLevel_lambda1)
    FN.cLevel_lambda2=rowMeans(fn.cLevel_lambda2)
    FN.cLevel_lambda3=rowMeans(fn.cLevel_lambda3)
    
    FN.PCGII_all_lambda1=rowMeans(fn.PCGII_all_lambda1)
    FN.PCGII_all_lambda2=rowMeans(fn.PCGII_all_lambda2)
    FN.PCGII_all_lambda3=rowMeans(fn.PCGII_all_lambda3)
    FN.PCGII_top5_lambda1=rowMeans(fn.PCGII_top5_lambda1)
    FN.PCGII_top5_lambda2=rowMeans(fn.PCGII_top5_lambda2)
    FN.PCGII_top5_lambda3=rowMeans(fn.PCGII_top5_lambda3)
    FN.PCGII_top15_lambda1=rowMeans(fn.PCGII_top15_lambda1)
    FN.PCGII_top15_lambda2=rowMeans(fn.PCGII_top15_lambda2)
    FN.PCGII_top15_lambda3=rowMeans(fn.PCGII_top15_lambda3)
    FN.PCGII_top25_lambda1=rowMeans(fn.PCGII_top25_lambda1)
    FN.PCGII_top25_lambda2=rowMeans(fn.PCGII_top25_lambda2)
    FN.PCGII_top25_lambda3=rowMeans(fn.PCGII_top25_lambda3)
    FN.PCGII_random10_lambda1=rowMeans(fn.PCGII_random10_lambda1)
    FN.PCGII_random10_lambda2=rowMeans(fn.PCGII_random10_lambda2)
    FN.PCGII_random10_lambda3=rowMeans(fn.PCGII_random10_lambda3)
    
    FN.FGGM_lambda1=rowMeans(fn.FGGM_lambda1)
    FN.FGGM_lambda2=rowMeans(fn.FGGM_lambda2)
    FN.FGGM_lambda3=rowMeans(fn.FGGM_lambda3)
    
    # empirical fdr
    fdr.enf.p=FP.enf.p/All.enf.p
    fdr.enf.p[is.na(fdr.enf.p)]=0
    fdr.enf.q=FP.enf.q/All.enf.q
    fdr.enf.q[is.na(fdr.enf.q)]=0
    fdr.mle.p=FP.mle.p/All.mle.p
    fdr.mle.p[is.na(fdr.mle.p)]=0
    fdr.mle.q=FP.mle.q/All.mle.q
    fdr.mle.q[is.na(fdr.mle.q)]=0
    fdr.cLevel_lambda1=FP.cLevel_lambda1/All.cLevel_lambda1
    fdr.cLevel_lambda1[is.na(fdr.cLevel_lambda1)]=0
    fdr.cLevel_lambda2=FP.cLevel_lambda2/All.cLevel_lambda2
    fdr.cLevel_lambda2[is.na(fdr.cLevel_lambda2)]=0
    fdr.cLevel_lambda3=FP.cLevel_lambda3/All.cLevel_lambda3
    fdr.cLevel_lambda3[is.na(fdr.cLevel_lambda3)]=0
    
    fdr.PCGII_all_lambda1=FP.PCGII_all_lambda1/All.PCGII_all_lambda1
    fdr.PCGII_all_lambda1[is.na(fdr.PCGII_all_lambda1)]=0
    fdr.PCGII_all_lambda2=FP.PCGII_all_lambda2/All.PCGII_all_lambda2
    fdr.PCGII_all_lambda2[is.na(fdr.PCGII_all_lambda2)]=0
    fdr.PCGII_all_lambda3=FP.PCGII_all_lambda3/All.PCGII_all_lambda3
    fdr.PCGII_all_lambda3[is.na(fdr.PCGII_all_lambda3)]=0
    
    fdr.PCGII_top5_lambda1=FP.PCGII_top5_lambda1/All.PCGII_top5_lambda1
    fdr.PCGII_top5_lambda1[is.na(fdr.PCGII_top5_lambda1)]=0
    fdr.PCGII_top5_lambda2=FP.PCGII_top5_lambda2/All.PCGII_top5_lambda2
    fdr.PCGII_top5_lambda2[is.na(fdr.PCGII_top5_lambda2)]=0
    fdr.PCGII_top5_lambda3=FP.PCGII_top5_lambda3/All.PCGII_top5_lambda3
    fdr.PCGII_top5_lambda3[is.na(fdr.PCGII_top5_lambda3)]=0
    
    fdr.PCGII_top15_lambda1=FP.PCGII_top15_lambda1/All.PCGII_top15_lambda1
    fdr.PCGII_top15_lambda1[is.na(fdr.PCGII_top15_lambda1)]=0
    fdr.PCGII_top15_lambda2=FP.PCGII_top15_lambda2/All.PCGII_top15_lambda2
    fdr.PCGII_top15_lambda2[is.na(fdr.PCGII_top15_lambda2)]=0
    fdr.PCGII_top15_lambda3=FP.PCGII_top15_lambda3/All.PCGII_top15_lambda3
    fdr.PCGII_top15_lambda3[is.na(fdr.PCGII_top15_lambda3)]=0
    
    fdr.PCGII_top25_lambda1=FP.PCGII_top25_lambda1/All.PCGII_top25_lambda1
    fdr.PCGII_top25_lambda1[is.na(fdr.PCGII_top25_lambda1)]=0
    fdr.PCGII_top25_lambda2=FP.PCGII_top25_lambda2/All.PCGII_top25_lambda2
    fdr.PCGII_top25_lambda2[is.na(fdr.PCGII_top25_lambda2)]=0
    fdr.PCGII_top25_lambda3=FP.PCGII_top25_lambda3/All.PCGII_top25_lambda3
    fdr.PCGII_top25_lambda3[is.na(fdr.PCGII_top25_lambda3)]=0
    
    fdr.PCGII_random10_lambda1=FP.PCGII_random10_lambda1/All.PCGII_random10_lambda1
    fdr.PCGII_random10_lambda1[is.na(fdr.PCGII_random10_lambda1)]=0
    fdr.PCGII_random10_lambda2=FP.PCGII_random10_lambda2/All.PCGII_random10_lambda2
    fdr.PCGII_random10_lambda2[is.na(fdr.PCGII_random10_lambda2)]=0
    fdr.PCGII_random10_lambda3=FP.PCGII_random10_lambda3/All.PCGII_random10_lambda3
    fdr.PCGII_random10_lambda3[is.na(fdr.PCGII_random10_lambda3)]=0
    
    fdr.FGGM_lambda1=FP.FGGM_lambda1/All.FGGM_lambda1
    fdr.FGGM_lambda1[is.na(fdr.FGGM_lambda1)]=0
    fdr.FGGM_lambda2=FP.FGGM_lambda2/All.FGGM_lambda2
    fdr.FGGM_lambda2[is.na(fdr.FGGM_lambda2)]=0
    fdr.FGGM_lambda3=FP.FGGM_lambda3/All.FGGM_lambda3
    fdr.FGGM_lambda3[is.na(fdr.FGGM_lambda3)]=0
    
    
    fdr.table=cbind.data.frame(
      FDR=rep(al,25), # nominal FDR
      methods=c(rep("ENF.p",length(al)),
                rep("ENF.q",length(al)),
                rep("MLE.p",length(al)),
                rep("MLE.q",length(al)),
                rep("cLevel_lambda1",length(al)),
                rep("cLevel_lambda2",length(al)),
                rep("cLevel_lambda3",length(al)),
                
                rep("PCGII_all_lambda1",length(al)),
                rep("PCGII_all_lambda2",length(al)),
                rep("PCGII_all_lambda3",length(al)),
                rep("PCGII_top5_lambda1",length(al)),
                rep("PCGII_top5_lambda2",length(al)),
                rep("PCGII_top5_lambda3",length(al)),
                rep("PCGII_top15_lambda1",length(al)),
                rep("PCGII_top15_lambda2",length(al)),
                rep("PCGII_top15_lambda3",length(al)),
                rep("PCGII_top25_lambda1",length(al)),
                rep("PCGII_top25_lambda2",length(al)),
                rep("PCGII_top25_lambda3",length(al)),
                rep("PCGII_random10_lambda1",length(al)),
                rep("PCGII_random10_lambda2",length(al)),
                rep("PCGII_random10_lambda3",length(al)),
                
                rep("FGGM_lambda1",length(al)),
                rep("FGGM_lambda2",length(al)),
                rep("FGGM_lambda3",length(al))),
      fdr=c(fdr.enf.p, fdr.enf.q, fdr.mle.p, fdr.mle.q, 
            fdr.cLevel_lambda1, fdr.cLevel_lambda2, fdr.cLevel_lambda3, 
            
            fdr.PCGII_all_lambda1, fdr.PCGII_all_lambda2,  fdr.PCGII_all_lambda3, 
            fdr.PCGII_top5_lambda1, fdr.PCGII_top5_lambda2,  fdr.PCGII_top5_lambda3, 
            fdr.PCGII_top15_lambda1, fdr.PCGII_top15_lambda2,  fdr.PCGII_top15_lambda3, 
            fdr.PCGII_top25_lambda1, fdr.PCGII_top25_lambda2,  fdr.PCGII_top25_lambda3, 
            fdr.PCGII_random10_lambda1, fdr.PCGII_random10_lambda2,  fdr.PCGII_random10_lambda3, 
            
            fdr.FGGM_lambda1,fdr.FGGM_lambda2,fdr.FGGM_lambda3))
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
          
          cLevel_lambda1.ordered.TPR=cLevel_lambda1.ordered.TPR, cLevel_lambda1.ordered.FPR=cLevel_lambda1.ordered.FPR, 
          cLevel_lambda1.ordered.PPV=cLevel_lambda1.ordered.PPV, cLevel_lambda1.ordered.FDR=cLevel_lambda1.ordered.FDR,
          cLevel_lambda2.ordered.TPR=cLevel_lambda2.ordered.TPR, cLevel_lambda2.ordered.FPR=cLevel_lambda2.ordered.FPR, 
          cLevel_lambda2.ordered.PPV=cLevel_lambda2.ordered.PPV, cLevel_lambda2.ordered.FDR=cLevel_lambda2.ordered.FDR,
          cLevel_lambda3.ordered.TPR=cLevel_lambda3.ordered.TPR, cLevel_lambda3.ordered.FPR=cLevel_lambda3.ordered.FPR, 
          cLevel_lambda3.ordered.PPV=cLevel_lambda3.ordered.PPV, cLevel_lambda3.ordered.FDR=cLevel_lambda3.ordered.FDR,
          
          PCGII_all_lambda1.ordered.TPR=PCGII_all_lambda1.ordered.TPR, PCGII_all_lambda1.ordered.FPR=PCGII_all_lambda1.ordered.FPR, 
          PCGII_all_lambda1.ordered.PPV=PCGII_all_lambda1.ordered.PPV, PCGII_all_lambda1.ordered.FDR=PCGII_all_lambda1.ordered.FDR,
          PCGII_all_lambda2.ordered.TPR=PCGII_all_lambda2.ordered.TPR, PCGII_all_lambda2.ordered.FPR=PCGII_all_lambda2.ordered.FPR, 
          PCGII_all_lambda2.ordered.PPV=PCGII_all_lambda2.ordered.PPV, PCGII_all_lambda2.ordered.FDR=PCGII_all_lambda2.ordered.FDR,
          PCGII_all_lambda3.ordered.TPR=PCGII_all_lambda3.ordered.TPR, PCGII_all_lambda3.ordered.FPR=PCGII_all_lambda3.ordered.FPR, 
          PCGII_all_lambda3.ordered.PPV=PCGII_all_lambda3.ordered.PPV, PCGII_all_lambda3.ordered.FDR=PCGII_all_lambda3.ordered.FDR,

          PCGII_top5_lambda1.ordered.TPR=PCGII_top5_lambda1.ordered.TPR, PCGII_top5_lambda1.ordered.FPR=PCGII_top5_lambda1.ordered.FPR, 
          PCGII_top5_lambda1.ordered.PPV=PCGII_top5_lambda1.ordered.PPV, PCGII_top5_lambda1.ordered.FDR=PCGII_top5_lambda1.ordered.FDR,
          PCGII_top5_lambda2.ordered.TPR=PCGII_top5_lambda2.ordered.TPR, PCGII_top5_lambda2.ordered.FPR=PCGII_top5_lambda2.ordered.FPR, 
          PCGII_top5_lambda2.ordered.PPV=PCGII_top5_lambda2.ordered.PPV, PCGII_top5_lambda2.ordered.FDR=PCGII_top5_lambda2.ordered.FDR,
          PCGII_top5_lambda3.ordered.TPR=PCGII_top5_lambda3.ordered.TPR, PCGII_top5_lambda3.ordered.FPR=PCGII_top5_lambda3.ordered.FPR, 
          PCGII_top5_lambda3.ordered.PPV=PCGII_top5_lambda3.ordered.PPV, PCGII_top5_lambda3.ordered.FDR=PCGII_top5_lambda3.ordered.FDR,
          PCGII_top15_lambda1.ordered.TPR=PCGII_top15_lambda1.ordered.TPR, PCGII_top15_lambda1.ordered.FPR=PCGII_top15_lambda1.ordered.FPR, 
          PCGII_top15_lambda1.ordered.PPV=PCGII_top15_lambda1.ordered.PPV, PCGII_top15_lambda1.ordered.FDR=PCGII_top15_lambda1.ordered.FDR,
          PCGII_top15_lambda2.ordered.TPR=PCGII_top15_lambda2.ordered.TPR, PCGII_top15_lambda2.ordered.FPR=PCGII_top15_lambda2.ordered.FPR, 
          PCGII_top15_lambda2.ordered.PPV=PCGII_top15_lambda2.ordered.PPV, PCGII_top15_lambda2.ordered.FDR=PCGII_top15_lambda2.ordered.FDR,
          PCGII_top15_lambda3.ordered.TPR=PCGII_top15_lambda3.ordered.TPR, PCGII_top15_lambda3.ordered.FPR=PCGII_top15_lambda3.ordered.FPR, 
          PCGII_top15_lambda3.ordered.PPV=PCGII_top15_lambda3.ordered.PPV, PCGII_top15_lambda3.ordered.FDR=PCGII_top15_lambda3.ordered.FDR,
          PCGII_top25_lambda1.ordered.TPR=PCGII_top25_lambda1.ordered.TPR, PCGII_top25_lambda1.ordered.FPR=PCGII_top25_lambda1.ordered.FPR, 
          PCGII_top25_lambda1.ordered.PPV=PCGII_top25_lambda1.ordered.PPV, PCGII_top25_lambda1.ordered.FDR=PCGII_top25_lambda1.ordered.FDR,
          PCGII_top25_lambda2.ordered.TPR=PCGII_top25_lambda2.ordered.TPR, PCGII_top25_lambda2.ordered.FPR=PCGII_top25_lambda2.ordered.FPR, 
          PCGII_top25_lambda2.ordered.PPV=PCGII_top25_lambda2.ordered.PPV, PCGII_top25_lambda2.ordered.FDR=PCGII_top25_lambda2.ordered.FDR,
          PCGII_top25_lambda3.ordered.TPR=PCGII_top25_lambda3.ordered.TPR, PCGII_top25_lambda3.ordered.FPR=PCGII_top25_lambda3.ordered.FPR, 
          PCGII_top25_lambda3.ordered.PPV=PCGII_top25_lambda3.ordered.PPV, PCGII_top25_lambda3.ordered.FDR=PCGII_top25_lambda3.ordered.FDR,
          PCGII_random10_lambda1.ordered.TPR=PCGII_random10_lambda1.ordered.TPR, PCGII_random10_lambda1.ordered.FPR=PCGII_random10_lambda1.ordered.FPR, 
          PCGII_random10_lambda1.ordered.PPV=PCGII_random10_lambda1.ordered.PPV, PCGII_random10_lambda1.ordered.FDR=PCGII_random10_lambda1.ordered.FDR,
          PCGII_random10_lambda2.ordered.TPR=PCGII_random10_lambda2.ordered.TPR, PCGII_random10_lambda2.ordered.FPR=PCGII_random10_lambda2.ordered.FPR, 
          PCGII_random10_lambda2.ordered.PPV=PCGII_random10_lambda2.ordered.PPV, PCGII_random10_lambda2.ordered.FDR=PCGII_random10_lambda2.ordered.FDR,
          PCGII_random10_lambda3.ordered.TPR=PCGII_random10_lambda3.ordered.TPR, PCGII_random10_lambda3.ordered.FPR=PCGII_random10_lambda3.ordered.FPR, 
          PCGII_random10_lambda3.ordered.PPV=PCGII_random10_lambda3.ordered.PPV, PCGII_random10_lambda3.ordered.FDR=PCGII_random10_lambda3.ordered.FDR,
          
          FGGM_lambda1.ordered.TPR=FGGM_lambda1.ordered.TPR, FGGM_lambda1.ordered.FPR=FGGM_lambda1.ordered.FPR, 
          FGGM_lambda1.ordered.PPV=FGGM_lambda1.ordered.PPV, FGGM_lambda1.ordered.FDR=FGGM_lambda1.ordered.FDR,
          FGGM_lambda2.ordered.TPR=FGGM_lambda2.ordered.TPR, FGGM_lambda2.ordered.FPR=FGGM_lambda2.ordered.FPR, 
          FGGM_lambda2.ordered.PPV=FGGM_lambda2.ordered.PPV, FGGM_lambda2.ordered.FDR=FGGM_lambda2.ordered.FDR,
          FGGM_lambda3.ordered.TPR=FGGM_lambda3.ordered.TPR, FGGM_lambda3.ordered.FPR=FGGM_lambda3.ordered.FPR, 
          FGGM_lambda3.ordered.PPV=FGGM_lambda3.ordered.PPV, FGGM_lambda3.ordered.FDR=FGGM_lambda3.ordered.FDR,
          rank.table=rank.table,
          
          # empirical TPR
          tempstore.ENF=tempstore.ENF, tempstore.MLE=tempstore.MLE, 
          tempstore.cLevel_lambda1=tempstore.cLevel_lambda1, 
          tempstore.cLevel_lambda2=tempstore.cLevel_lambda2, 
          tempstore.cLevel_lambda3=tempstore.cLevel_lambda3, 
          
          tempstore.PCGII_all_lambda1=tempstore.PCGII_all_lambda1, 
          tempstore.PCGII_all_lambda2=tempstore.PCGII_all_lambda2, 
          tempstore.PCGII_all_lambda3=tempstore.PCGII_all_lambda3, 
          tempstore.PCGII_top5_lambda1=tempstore.PCGII_top5_lambda1, 
          tempstore.PCGII_top5_lambda2=tempstore.PCGII_top5_lambda2, 
          tempstore.PCGII_top5_lambda3=tempstore.PCGII_top5_lambda3, 
          tempstore.PCGII_top15_lambda1=tempstore.PCGII_top15_lambda1, 
          tempstore.PCGII_top15_lambda2=tempstore.PCGII_top15_lambda2, 
          tempstore.PCGII_top15_lambda3=tempstore.PCGII_top15_lambda3, 
          tempstore.PCGII_top25_lambda1=tempstore.PCGII_top25_lambda1, 
          tempstore.PCGII_top25_lambda2=tempstore.PCGII_top25_lambda2, 
          tempstore.PCGII_top25_lambda3=tempstore.PCGII_top25_lambda3, 
          tempstore.PCGII_random10_lambda1=tempstore.PCGII_random10_lambda1, 
          tempstore.PCGII_random10_lambda2=tempstore.PCGII_random10_lambda2, 
          tempstore.PCGII_random10_lambda3=tempstore.PCGII_random10_lambda3, 
          
          tempstore.FGGM_lambda1=tempstore.FGGM_lambda1,
          tempstore.FGGM_lambda2=tempstore.FGGM_lambda2,
          tempstore.FGGM_lambda3=tempstore.FGGM_lambda3,
          
          ROC.out=ROC.out),
        
        # FDR
        metrics=list( 
          # Number of total correctly selected edges
          all.enf.p=all.enf.p, all.enf.q=all.enf.q,  
          all.mle.p=all.mle.p,  all.mle.q=all.mle.q,
          all.cLevel_lambda1=all.cLevel_lambda1, 
          all.cLevel_lambda2=all.cLevel_lambda2, 
          all.cLevel_lambda3=all.cLevel_lambda3, 
          
          all.PCGII_all_lambda1=all.PCGII_all_lambda1, 
          all.PCGII_all_lambda2=all.PCGII_all_lambda2, 
          all.PCGII_all_lambda3=all.PCGII_all_lambda3,
          all.PCGII_top5_lambda1=all.PCGII_top5_lambda1, 
          all.PCGII_top5_lambda2=all.PCGII_top5_lambda2, 
          all.PCGII_top5_lambda3=all.PCGII_top5_lambda3,
          all.PCGII_top15_lambda1=all.PCGII_top15_lambda1, 
          all.PCGII_top15_lambda2=all.PCGII_top15_lambda2, 
          all.PCGII_top15_lambda3=all.PCGII_top15_lambda3,
          all.PCGII_top25_lambda1=all.PCGII_top25_lambda1, 
          all.PCGII_top25_lambda2=all.PCGII_top25_lambda2, 
          all.PCGII_top25_lambda3=all.PCGII_top25_lambda3,
          all.PCGII_random10_lambda1=all.PCGII_random10_lambda1, 
          all.PCGII_random10_lambda2=all.PCGII_random10_lambda2, 
          all.PCGII_random10_lambda3=all.PCGII_random10_lambda3,
          
          
          all.FGGM_lambda1=all.FGGM_lambda1,
          all.FGGM_lambda2=all.FGGM_lambda2,
          all.FGGM_lambda3=all.FGGM_lambda3,
          
          # true positives
          tp.enf.p=tp.enf.p, tp.enf.q=tp.enf.q, 
          tp.mle.p=tp.mle.p, tp.mle.q=tp.mle.q,
          tp.cLevel_lambda1=tp.cLevel_lambda1, 
          tp.cLevel_lambda2=tp.cLevel_lambda2, 
          tp.cLevel_lambda3=tp.cLevel_lambda3,
          
          tp.PCGII_all_lambda1=tp.PCGII_all_lambda1, 
          tp.PCGII_all_lambda2=tp.PCGII_all_lambda2, 
          tp.PCGII_all_lambda3=tp.PCGII_all_lambda3,
          tp.PCGII_top5_lambda1=tp.PCGII_top5_lambda1, 
          tp.PCGII_top5_lambda2=tp.PCGII_top5_lambda2, 
          tp.PCGII_top5_lambda3=tp.PCGII_top5_lambda3,
          tp.PCGII_top15_lambda1=tp.PCGII_top15_lambda1, 
          tp.PCGII_top15_lambda2=tp.PCGII_top15_lambda2, 
          tp.PCGII_top15_lambda3=tp.PCGII_top15_lambda3,
          tp.PCGII_top25_lambda1=tp.PCGII_top25_lambda1, 
          tp.PCGII_top25_lambda2=tp.PCGII_top25_lambda2, 
          tp.PCGII_top25_lambda3=tp.PCGII_top25_lambda3,
          tp.PCGII_random10_lambda1=tp.PCGII_random10_lambda1, 
          tp.PCGII_random10_lambda2=tp.PCGII_random10_lambda2, 
          tp.PCGII_random10_lambda3=tp.PCGII_random10_lambda3,
          
          tp.FGGM_lambda1=tp.FGGM_lambda1,
          tp.FGGM_lambda2=tp.FGGM_lambda2,
          tp.FGGM_lambda3=tp.FGGM_lambda3,
          
          # false positives
          fp.enf.p=fp.enf.p, fp.enf.q=fp.enf.q, 
          fp.mle.p=fp.mle.p, fp.mle.q=fp.mle.q,
          fp.cLevel_lambda1=fp.cLevel_lambda1, 
          fp.cLevel_lambda2=fp.cLevel_lambda2, 
          fp.cLevel_lambda3=fp.cLevel_lambda3,
          
          fp.PCGII_all_lambda1=fp.PCGII_all_lambda1, 
          fp.PCGII_all_lambda2=fp.PCGII_all_lambda2, 
          fp.PCGII_all_lambda3=fp.PCGII_all_lambda3,
          fp.PCGII_top5_lambda1=fp.PCGII_top5_lambda1, 
          fp.PCGII_top5_lambda2=fp.PCGII_top5_lambda2, 
          fp.PCGII_top5_lambda3=fp.PCGII_top5_lambda3,
          fp.PCGII_top15_lambda1=fp.PCGII_top15_lambda1, 
          fp.PCGII_top15_lambda2=fp.PCGII_top15_lambda2, 
          fp.PCGII_top15_lambda3=fp.PCGII_top15_lambda3,
          fp.PCGII_top25_lambda1=fp.PCGII_top25_lambda1, 
          fp.PCGII_top25_lambda2=fp.PCGII_top25_lambda2, 
          fp.PCGII_top25_lambda3=fp.PCGII_top25_lambda3,
          fp.PCGII_random10_lambda1=fp.PCGII_random10_lambda1, 
          fp.PCGII_random10_lambda2=fp.PCGII_random10_lambda2, 
          fp.PCGII_random10_lambda3=fp.PCGII_random10_lambda3,
          
          fp.FGGM_lambda1=fp.FGGM_lambda1,
          fp.FGGM_lambda2=fp.FGGM_lambda2,
          fp.FGGM_lambda3=fp.FGGM_lambda3,
          
          # false negatives
          tn.enf.p=tn.enf.p, tn.enf.q=tn.enf.q, 
          tn.mle.p=tn.mle.p, tn.mle.q=tn.mle.q,
          tn.cLevel_lambda1=tn.cLevel_lambda1, 
          tn.cLevel_lambda2=tn.cLevel_lambda2, 
          tn.cLevel_lambda3=tn.cLevel_lambda3, 
          
          tn.PCGII_all_lambda1=tn.PCGII_all_lambda1, 
          tn.PCGII_all_lambda2=tn.PCGII_all_lambda2, 
          tn.PCGII_all_lambda3=tn.PCGII_all_lambda3,
          tn.PCGII_top5_lambda1=tn.PCGII_top5_lambda1, 
          tn.PCGII_top5_lambda2=tn.PCGII_top5_lambda2, 
          tn.PCGII_top5_lambda3=tn.PCGII_top5_lambda3,
          tn.PCGII_top15_lambda1=tn.PCGII_top15_lambda1, 
          tn.PCGII_top15_lambda2=tn.PCGII_top15_lambda2, 
          tn.PCGII_top15_lambda3=tn.PCGII_top15_lambda3,
          tn.PCGII_top25_lambda1=tn.PCGII_top25_lambda1, 
          tn.PCGII_top25_lambda2=tn.PCGII_top25_lambda2, 
          tn.PCGII_top25_lambda3=tn.PCGII_top25_lambda3,
          tn.PCGII_random10_lambda1=tn.PCGII_random10_lambda1, 
          tn.PCGII_random10_lambda2=tn.PCGII_random10_lambda2, 
          tn.PCGII_random10_lambda3=tn.PCGII_random10_lambda3,
          
          tn.FGGM_lambda1=tn.FGGM_lambda1,
          tn.FGGM_lambda2=tn.FGGM_lambda2,
          tn.FGGM_lambda3=tn.FGGM_lambda3,
          
          # false negatives
          fn.enf.p=fn.enf.p, fn.enf.q=fn.enf.q,
          fn.mle.p=fn.mle.p, fn.mle.q=fn.mle.q,
          fn.cLevel_lambda1=fn.cLevel_lambda1, 
          fn.cLevel_lambda2=fn.cLevel_lambda2, 
          fn.cLevel_lambda3=fn.cLevel_lambda3, 
          
          fn.PCGII_all_lambda1=fn.PCGII_all_lambda1, 
          fn.PCGII_all_lambda2=fn.PCGII_all_lambda2, 
          fn.PCGII_all_lambda3=fn.PCGII_all_lambda3,
          fn.PCGII_top5_lambda1=fn.PCGII_top5_lambda1, 
          fn.PCGII_top5_lambda2=fn.PCGII_top5_lambda2, 
          fn.PCGII_top5_lambda3=fn.PCGII_top5_lambda3,
          fn.PCGII_top15_lambda1=fn.PCGII_top15_lambda1, 
          fn.PCGII_top15_lambda2=fn.PCGII_top15_lambda2, 
          fn.PCGII_top15_lambda3=fn.PCGII_top15_lambda3,
          fn.PCGII_top25_lambda1=fn.PCGII_top25_lambda1, 
          fn.PCGII_top25_lambda2=fn.PCGII_top25_lambda2, 
          fn.PCGII_top25_lambda3=fn.PCGII_top25_lambda3,
          fn.PCGII_random10_lambda1=fn.PCGII_random10_lambda1, 
          fn.PCGII_random10_lambda2=fn.PCGII_random10_lambda2, 
          fn.PCGII_random10_lambda3=fn.PCGII_random10_lambda3,
          
          fn.FGGM_lambda1=fn.FGGM_lambda1,
          fn.FGGM_lambda2=fn.FGGM_lambda2,
          fn.FGGM_lambda3=fn.FGGM_lambda3,
          
          All.enf.p=All.enf.p, All.enf.q=All.enf.q, 
          All.mle.p=All.mle.p, All.mle.q=All.mle.q,
          All.cLevel_lambda1=All.cLevel_lambda1, 
          All.cLevel_lambda2=All.cLevel_lambda2, 
          All.cLevel_lambda3=All.cLevel_lambda3, 
          
          All.PCGII_all_lambda1=All.PCGII_all_lambda1, 
          All.PCGII_all_lambda2=All.PCGII_all_lambda2, 
          All.PCGII_all_lambda3=All.PCGII_all_lambda3,
          All.PCGII_top5_lambda1=All.PCGII_top5_lambda1, 
          All.PCGII_top5_lambda2=All.PCGII_top5_lambda2, 
          All.PCGII_top5_lambda3=All.PCGII_top5_lambda3,
          All.PCGII_top15_lambda1=All.PCGII_top15_lambda1, 
          All.PCGII_top15_lambda2=All.PCGII_top15_lambda2, 
          All.PCGII_top15_lambda3=All.PCGII_top15_lambda3,
          All.PCGII_top25_lambda1=All.PCGII_top25_lambda1, 
          All.PCGII_top25_lambda2=All.PCGII_top25_lambda2, 
          All.PCGII_top25_lambda3=All.PCGII_top25_lambda3,
          All.PCGII_random10_lambda1=All.PCGII_random10_lambda1, 
          All.PCGII_random10_lambda2=All.PCGII_random10_lambda2, 
          All.PCGII_random10_lambda3=All.PCGII_random10_lambda3,
          
          All.FGGM_lambda1=All.FGGM_lambda1,
          All.FGGM_lambda2=All.FGGM_lambda2,
          All.FGGM_lambda3=All.FGGM_lambda3,
          
          # true positives
          TP.enf.p=TP.enf.p,TP.enf.q=TP.enf.q, 
          TP.mle.p=TP.mle.p,TP.mle.q=TP.mle.q, 
          TP.cLevel_lambda1=TP.cLevel_lambda1,
          TP.cLevel_lambda2=TP.cLevel_lambda2,
          TP.cLevel_lambda3=TP.cLevel_lambda3,
          
          TP.PCGII_all_lambda1=TP.PCGII_all_lambda1, 
          TP.PCGII_all_lambda2=TP.PCGII_all_lambda2, 
          TP.PCGII_all_lambda3=TP.PCGII_all_lambda3,
          TP.PCGII_top5_lambda1=TP.PCGII_top5_lambda1, 
          TP.PCGII_top5_lambda2=TP.PCGII_top5_lambda2, 
          TP.PCGII_top5_lambda3=TP.PCGII_top5_lambda3,
          TP.PCGII_top15_lambda1=TP.PCGII_top15_lambda1, 
          TP.PCGII_top15_lambda2=TP.PCGII_top15_lambda2, 
          TP.PCGII_top15_lambda3=TP.PCGII_top15_lambda3,
          TP.PCGII_top25_lambda1=TP.PCGII_top25_lambda1, 
          TP.PCGII_top25_lambda2=TP.PCGII_top25_lambda2, 
          TP.PCGII_top25_lambda3=TP.PCGII_top25_lambda3,
          TP.PCGII_random10_lambda1=TP.PCGII_random10_lambda1, 
          TP.PCGII_random10_lambda2=TP.PCGII_random10_lambda2, 
          TP.PCGII_random10_lambda3=TP.PCGII_random10_lambda3,
          
          TP.FGGM_lambda1=TP.FGGM_lambda1,
          TP.FGGM_lambda2=TP.FGGM_lambda2,
          TP.FGGM_lambda3=TP.FGGM_lambda3,
          
          # false positives
          FP.enf.p=FP.enf.p,FP.enf.q=FP.enf.q, 
          FP.mle.p=FP.mle.p, FP.mle.q=FP.mle.q,
          FP.cLevel_lambda1=FP.cLevel_lambda1, 
          FP.cLevel_lambda2=FP.cLevel_lambda2,
          FP.cLevel_lambda3=FP.cLevel_lambda3,
          
          FP.PCGII_all_lambda1=FP.PCGII_all_lambda1, 
          FP.PCGII_all_lambda2=FP.PCGII_all_lambda2, 
          FP.PCGII_all_lambda3=FP.PCGII_all_lambda3,
          FP.PCGII_top5_lambda1=FP.PCGII_top5_lambda1, 
          FP.PCGII_top5_lambda2=FP.PCGII_top5_lambda2, 
          FP.PCGII_top5_lambda3=FP.PCGII_top5_lambda3,
          FP.PCGII_top15_lambda1=FP.PCGII_top15_lambda1, 
          FP.PCGII_top15_lambda2=FP.PCGII_top15_lambda2, 
          FP.PCGII_top15_lambda3=FP.PCGII_top15_lambda3,
          FP.PCGII_top25_lambda1=FP.PCGII_top25_lambda1, 
          FP.PCGII_top25_lambda2=FP.PCGII_top25_lambda2, 
          FP.PCGII_top25_lambda3=FP.PCGII_top25_lambda3,
          FP.PCGII_random10_lambda1=FP.PCGII_random10_lambda1, 
          FP.PCGII_random10_lambda2=FP.PCGII_random10_lambda2, 
          FP.PCGII_random10_lambda3=FP.PCGII_random10_lambda3,
          
          FP.FGGM_lambda1=FP.FGGM_lambda1,
          FP.FGGM_lambda2=FP.FGGM_lambda2,
          FP.FGGM_lambda3=FP.FGGM_lambda3,
          
          # true negatives
          TN.enf.p=TN.enf.p, TN.enf.q=TN.enf.q, 
          TN.mle.p=TN.mle.p, TN.mle.q=TN.mle.q,
          TN.cLevel_lambda1=TN.cLevel_lambda1, 
          TN.cLevel_lambda2=TN.cLevel_lambda2, 
          TN.cLevel_lambda3=TN.cLevel_lambda3, 
          
          TN.PCGII_all_lambda1=TN.PCGII_all_lambda1, 
          TN.PCGII_all_lambda2=TN.PCGII_all_lambda2, 
          TN.PCGII_all_lambda3=TN.PCGII_all_lambda3,
          TN.PCGII_top5_lambda1=TN.PCGII_top5_lambda1, 
          TN.PCGII_top5_lambda2=TN.PCGII_top5_lambda2, 
          TN.PCGII_top5_lambda3=TN.PCGII_top5_lambda3,
          TN.PCGII_top15_lambda1=TN.PCGII_top15_lambda1, 
          TN.PCGII_top15_lambda2=TN.PCGII_top15_lambda2, 
          TN.PCGII_top15_lambda3=TN.PCGII_top15_lambda3,
          TN.PCGII_top25_lambda1=TN.PCGII_top25_lambda1, 
          TN.PCGII_top25_lambda2=TN.PCGII_top25_lambda2, 
          TN.PCGII_top25_lambda3=TN.PCGII_top25_lambda3,
          TN.PCGII_random10_lambda1=TN.PCGII_random10_lambda1, 
          TN.PCGII_random10_lambda2=TN.PCGII_random10_lambda2, 
          TN.PCGII_random10_lambda3=TN.PCGII_random10_lambda3,
          
          TN.FGGM_lambda1=TN.FGGM_lambda1,
          TN.FGGM_lambda2=TN.FGGM_lambda2,
          TN.FGGM_lambda3=TN.FGGM_lambda3,
          
          # false negatives
          FN.enf.p=FN.enf.p, FN.enf.q=FN.enf.q, 
          FN.mle.p=FN.mle.p,FN.mle.q=FN.mle.q,
          FN.cLevel_lambda1=FN.cLevel_lambda1,
          FN.cLevel_lambda2=FN.cLevel_lambda2, 
          fn.cLevel_lambda3=FN.cLevel_lambda3, 
          
          FN.PCGII_all_lambda1=FN.PCGII_all_lambda1, 
          FN.PCGII_all_lambda2=FN.PCGII_all_lambda2, 
          FN.PCGII_all_lambda3=FN.PCGII_all_lambda3,
          FN.PCGII_top5_lambda1=FN.PCGII_top5_lambda1, 
          FN.PCGII_top5_lambda2=FN.PCGII_top5_lambda2, 
          FN.PCGII_top5_lambda3=FN.PCGII_top5_lambda3,
          FN.PCGII_top15_lambda1=FN.PCGII_top15_lambda1, 
          FN.PCGII_top15_lambda2=FN.PCGII_top15_lambda2, 
          FN.PCGII_top15_lambda3=FN.PCGII_top15_lambda3,
          FN.PCGII_top25_lambda1=FN.PCGII_top25_lambda1, 
          FN.PCGII_top25_lambda2=FN.PCGII_top25_lambda2, 
          FN.PCGII_top25_lambda3=FN.PCGII_top25_lambda3,
          FN.PCGII_random10_lambda1=FN.PCGII_random10_lambda1, 
          FN.PCGII_random10_lambda2=FN.PCGII_random10_lambda2, 
          FN.PCGII_random10_lambda3=FN.PCGII_random10_lambda3,
          
          FN.FGGM_lambda1=FN.FGGM_lambda1,
          FN.FGGM_lambda2=FN.FGGM_lambda2,
          FN.FGGM_lambda3=FN.FGGM_lambda3
        ),
        # empirical fdr
        fdr=list(
          fdr.enf.p=fdr.enf.p, fdr.enf.q=fdr.enf.q, 
          fdr.mle.p=fdr.mle.p, fdr.mle.q=fdr.mle.q,
          fdr.cLevel_lambda1=fdr.cLevel_lambda1, 
          fdr.cLevel_lambda2=fdr.cLevel_lambda2, 
          fdr.cLevel_lambda3=fdr.cLevel_lambda3, 
          
          fdr.PCGII_all_lambda1=fdr.PCGII_all_lambda1, 
          fdr.PCGII_all_lambda2=fdr.PCGII_all_lambda2, 
          fdr.PCGII_all_lambda3=fdr.PCGII_all_lambda3,
          fdr.PCGII_top5_lambda1=fdr.PCGII_top5_lambda1, 
          fdr.PCGII_top5_lambda2=fdr.PCGII_top5_lambda2, 
          fdr.PCGII_top5_lambda3=fdr.PCGII_top5_lambda3,
          fdr.PCGII_top15_lambda1=fdr.PCGII_top15_lambda1, 
          fdr.PCGII_top15_lambda2=fdr.PCGII_top15_lambda2, 
          fdr.PCGII_top15_lambda3=fdr.PCGII_top15_lambda3,
          fdr.PCGII_top25_lambda1=fdr.PCGII_top25_lambda1, 
          fdr.PCGII_top25_lambda2=fdr.PCGII_top25_lambda2, 
          fdr.PCGII_top25_lambda3=fdr.PCGII_top25_lambda3,
          fdr.PCGII_random10_lambda1=fdr.PCGII_random10_lambda1, 
          fdr.PCGII_random10_lambda2=fdr.PCGII_random10_lambda2, 
          fdr.PCGII_random10_lambda3=fdr.PCGII_random10_lambda3,
          
          fdr.FGGM_lambda1=fdr.FGGM_lambda1,
          fdr.FGGM_lambda2=fdr.FGGM_lambda2,
          fdr.FGGM_lambda3=fdr.FGGM_lambda3,
          
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
      load(file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/ScaleFree1_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=25)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/Results/ScaleFree1_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}


#### Scale Free .5
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      load(file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/ScaleFree.5_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=25)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/Results/ScaleFree.5_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}

#### Scale Free .1
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in 1:3){
      load(file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/ScaleFree.1_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=25)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/Results/ScaleFree.1_n%d_p%d_e%d.RData", n, p, e))
    }
  }
}


#### Random

nl = c(60) # Sample Size
pl = c(100, 200)  # Number of Genes

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.02,0.03)){
      load(file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta))
      res=eval_models(X, omega, rep=25)
      save(res,file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/Results/Random_n%d_p%d_eta%g.RData", n, p, eta))
    }
  }
}

nl = c(80) # Sample Size
pl = c(100, 200)  # Number of Genes

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.02,0.03)){
      load(file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/Random_simu_n%d_p%d_eta%g.RData", n, p, eta))
      res=eval_models(X, omega, rep=25)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/Results/Random_n%d_p%d_eta%g.RData", n, p, eta))
    }
  }
}

#### BlockDiag
nl=c(60,80)
pl=c(100,200)

for(n in nl){    
  for(p in pl){
    for(e in c(4,8,10)){
      load(file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/BlockDiag_simu_n%d_p%d_e%d.RData", n, p, e))
      res=eval_models(X, omega, rep=25)
      save(res, file=sprintf("~/Desktop/GenePCG/R/20220507/strong_signal/Results/BlockDiag_n%d_p%d_e%d.RData", n, p, e))
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
res$ROC$ROC.out %>% 
  tidyr::gather("methods", "TPR", ENF:FGGM_lambda3) %>%
  filter(!methods %in% c("ENF")) %>% 
  mutate(group=c(rep("MLE",1000),rep(c(rep("lam1",1000),rep("lam2",1000),rep("lam3",1000)),7))) %>%
  mutate(approach=c(rep("MLE",1000), rep("CLEVEL",3000), rep("PCGII_all",3000),rep("PCGII_top5",3000),rep("PCGII_top15",3000),rep("PCGII_top25",3000),rep("PCGII_random10",3000), rep("FGGM",3000))) %>% 
  ggplot(aes(x=FPR, y=TPR,col=approach)) + 
  geom_line(aes(linetype=group))

res$fdr$fdr.table %>% 
  filter(!methods %in% c("ENF.p","MLE.p","ENF.q")) %>%
  mutate(group=c(rep("MLE",58),rep(c(rep("lam1",58),rep("lam2",58),rep("lam3",58)),7))) %>%
  mutate(approach=c(rep("MLE",58), rep("CLEVEL",58*3), rep("PCGII_all",58*3),rep("PCGII_top5",58*3),rep("PCGII_top15",58*3),rep("PCGII_top25",58*3),rep("PCGII_random10",58*3), rep("FGGM",58*3))) %>% 
  ggplot(aes(x=FDR, y=fdr, col=approach))+
  geom_line(aes(linetype=group)) + 
  geom_abline(slope = 1, intercept = 0)


my_power=function(RESULT, alpha=0.05){ # TP/P=TP/(TP+FN)
  ref=seq(0.001, 0.2005, 0.0035)
  ind=which(ref==alpha)
  recall=cbind(RESULT$metrics$tp.enf.p[ind,]/(RESULT$metrics$tp.enf.p[ind,]+RESULT$metrics$fn.enf.p[ind,]), # TPR
               RESULT$metrics$tp.enf.q[ind,]/(RESULT$metrics$tp.enf.q[ind,]+RESULT$metrics$fn.enf.q[ind,]),
               RESULT$metrics$tp.mle.p[ind,]/(RESULT$metrics$tp.mle.p[ind,]+RESULT$metrics$fn.mle.p[ind,]),
               RESULT$metrics$tp.mle.q[ind,]/(RESULT$metrics$tp.mle.q[ind,]+RESULT$metrics$fn.mle.q[ind,]),
               
               RESULT$metrics$tp.cLevel_lambda1[ind,]/(RESULT$metrics$tp.cLevel_lambda1[ind,]+RESULT$metrics$fn.cLevel_lambda1[ind,]),
               RESULT$metrics$tp.cLevel_lambda2[ind,]/(RESULT$metrics$tp.cLevel_lambda2[ind,]+RESULT$metrics$fn.cLevel_lambda2[ind,]),
               RESULT$metrics$tp.cLevel_lambda3[ind,]/(RESULT$metrics$tp.cLevel_lambda3[ind,]+RESULT$metrics$fn.cLevel_lambda3[ind,]),
               
               RESULT$metrics$tp.PCGII_all_lambda1[ind,]/(RESULT$metrics$tp.PCGII_all_lambda1[ind,]+RESULT$metrics$fn.PCGII_all_lambda1[ind,]),
               RESULT$metrics$tp.PCGII_all_lambda2[ind,]/(RESULT$metrics$tp.PCGII_all_lambda2[ind,]+RESULT$metrics$fn.PCGII_all_lambda2[ind,]),
               RESULT$metrics$tp.PCGII_all_lambda3[ind,]/(RESULT$metrics$tp.PCGII_all_lambda3[ind,]+RESULT$metrics$fn.PCGII_all_lambda3[ind,]),
               
               RESULT$metrics$tp.PCGII_top5_lambda1[ind,]/(RESULT$metrics$tp.PCGII_top5_lambda1[ind,]+RESULT$metrics$fn.PCGII_top5_lambda1[ind,]),
               RESULT$metrics$tp.PCGII_top5_lambda2[ind,]/(RESULT$metrics$tp.PCGII_top5_lambda2[ind,]+RESULT$metrics$fn.PCGII_top5_lambda2[ind,]),
               RESULT$metrics$tp.PCGII_top5_lambda3[ind,]/(RESULT$metrics$tp.PCGII_top5_lambda3[ind,]+RESULT$metrics$fn.PCGII_top5_lambda3[ind,]),
               
               RESULT$metrics$tp.PCGII_top15_lambda1[ind,]/(RESULT$metrics$tp.PCGII_top15_lambda1[ind,]+RESULT$metrics$fn.PCGII_top15_lambda1[ind,]),
               RESULT$metrics$tp.PCGII_top15_lambda2[ind,]/(RESULT$metrics$tp.PCGII_top15_lambda2[ind,]+RESULT$metrics$fn.PCGII_top15_lambda2[ind,]),
               RESULT$metrics$tp.PCGII_top15_lambda3[ind,]/(RESULT$metrics$tp.PCGII_top15_lambda3[ind,]+RESULT$metrics$fn.PCGII_top15_lambda3[ind,]),
               
               RESULT$metrics$tp.PCGII_top25_lambda1[ind,]/(RESULT$metrics$tp.PCGII_top25_lambda1[ind,]+RESULT$metrics$fn.PCGII_top25_lambda1[ind,]),
               RESULT$metrics$tp.PCGII_top25_lambda2[ind,]/(RESULT$metrics$tp.PCGII_top25_lambda2[ind,]+RESULT$metrics$fn.PCGII_top25_lambda2[ind,]),
               RESULT$metrics$tp.PCGII_top25_lambda3[ind,]/(RESULT$metrics$tp.PCGII_top25_lambda3[ind,]+RESULT$metrics$fn.PCGII_top25_lambda3[ind,]),
               
               RESULT$metrics$tp.PCGII_random10_lambda1[ind,]/(RESULT$metrics$tp.PCGII_random10_lambda1[ind,]+RESULT$metrics$fn.PCGII_random10_lambda1[ind,]),
               RESULT$metrics$tp.PCGII_random10_lambda2[ind,]/(RESULT$metrics$tp.PCGII_random10_lambda2[ind,]+RESULT$metrics$fn.PCGII_random10_lambda2[ind,]),
               RESULT$metrics$tp.PCGII_random10_lambda3[ind,]/(RESULT$metrics$tp.PCGII_random10_lambda3[ind,]+RESULT$metrics$fn.PCGII_random10_lambda3[ind,]),
               
               RESULT$metrics$tp.FGGM_lambda1[ind,]/(RESULT$metrics$tp.FGGM_lambda1[ind,]+RESULT$metrics$fn.FGGM_lambda1[ind,]),
               RESULT$metrics$tp.FGGM_lambda2[ind,]/(RESULT$metrics$tp.FGGM_lambda2[ind,]+RESULT$metrics$fn.FGGM_lambda2[ind,]),
               RESULT$metrics$tp.FGGM_lambda3[ind,]/(RESULT$metrics$tp.FGGM_lambda3[ind,]+RESULT$metrics$fn.FGGM_lambda3[ind,])
  )
  
  
  power=as.data.frame(recall)
  colnames(power)=c("ENF.p","ENF.q","MLE.p","MLE.q",
                    "cLevel_lambda1","cLevel_lambda2","cLevel_lambda3",
                    "PCGII_all_lambda1","PCGII_all_lambda2","PCGII_all_lambda3",
                    "PCGII_top5_lambda1","PCGII_top5_lambda2","PCGII_top5_lambda3",
                    "PCGII_top15_lambda1","PCGII_top15_lambda2","PCGII_top15_lambda3",
                    "PCGII_top25_lambda1","PCGII_top25_lambda2","PCGII_top25_lambda3",
                    "PCGII_random10_lambda1","PCGII_random10_lambda2","PCGII_random10_lambda3",
                    "FGGM_lambda1","FGGM_lambda2","FGGM_lambda3")
  power
}

my_fdr=function(RESULT, alpha=0.05){ #FP/Disc=FP/(FP+TP)
  ref=seq(0.001, 0.2005, 0.0035)
  ind=which(ref==alpha)
  fdr=cbind.data.frame(RESULT$metrics$fp.enf.p[ind,]/(RESULT$metrics$tp.enf.p[ind,]+RESULT$metrics$fp.enf.p[ind,]), 
                       RESULT$metrics$fp.enf.q[ind,]/(RESULT$metrics$tp.enf.q[ind,]+RESULT$metrics$fp.enf.q[ind,]),
                       RESULT$metrics$fp.mle.p[ind,]/(RESULT$metrics$tp.mle.p[ind,]+RESULT$metrics$fp.mle.p[ind,]),
                       RESULT$metrics$fp.mle.q[ind,]/(RESULT$metrics$tp.mle.q[ind,]+RESULT$metrics$fp.mle.q[ind,]),
                       
                       RESULT$metrics$fp.cLevel_lambda1[ind,]/(RESULT$metrics$tp.cLevel_lambda1[ind,]+RESULT$metrics$fp.cLevel_lambda1[ind,]),
                       RESULT$metrics$fp.cLevel_lambda2[ind,]/(RESULT$metrics$tp.cLevel_lambda2[ind,]+RESULT$metrics$fp.cLevel_lambda2[ind,]),
                       RESULT$metrics$fp.cLevel_lambda3[ind,]/(RESULT$metrics$tp.cLevel_lambda3[ind,]+RESULT$metrics$fp.cLevel_lambda3[ind,]),
                       
                       RESULT$metrics$fp.PCGII_all_lambda1[ind,]/(RESULT$metrics$tp.PCGII_all_lambda1[ind,]+RESULT$metrics$fp.PCGII_all_lambda1[ind,]),
                       RESULT$metrics$fp.PCGII_all_lambda2[ind,]/(RESULT$metrics$tp.PCGII_all_lambda2[ind,]+RESULT$metrics$fp.PCGII_all_lambda2[ind,]),
                       RESULT$metrics$fp.PCGII_all_lambda3[ind,]/(RESULT$metrics$tp.PCGII_all_lambda3[ind,]+RESULT$metrics$fp.PCGII_all_lambda3[ind,]),
                       RESULT$metrics$fp.PCGII_top5_lambda1[ind,]/(RESULT$metrics$tp.PCGII_top5_lambda1[ind,]+RESULT$metrics$fp.PCGII_top5_lambda1[ind,]),
                       RESULT$metrics$fp.PCGII_top5_lambda2[ind,]/(RESULT$metrics$tp.PCGII_top5_lambda2[ind,]+RESULT$metrics$fp.PCGII_top5_lambda2[ind,]),
                       RESULT$metrics$fp.PCGII_top5_lambda3[ind,]/(RESULT$metrics$tp.PCGII_top5_lambda3[ind,]+RESULT$metrics$fp.PCGII_top5_lambda3[ind,]),
                       RESULT$metrics$fp.PCGII_top15_lambda1[ind,]/(RESULT$metrics$tp.PCGII_top15_lambda1[ind,]+RESULT$metrics$fp.PCGII_top15_lambda1[ind,]),
                       RESULT$metrics$fp.PCGII_top15_lambda2[ind,]/(RESULT$metrics$tp.PCGII_top15_lambda2[ind,]+RESULT$metrics$fp.PCGII_top15_lambda2[ind,]),
                       RESULT$metrics$fp.PCGII_top15_lambda3[ind,]/(RESULT$metrics$tp.PCGII_top15_lambda3[ind,]+RESULT$metrics$fp.PCGII_top15_lambda3[ind,]),
                       RESULT$metrics$fp.PCGII_top25_lambda1[ind,]/(RESULT$metrics$tp.PCGII_top25_lambda1[ind,]+RESULT$metrics$fp.PCGII_top25_lambda1[ind,]),
                       RESULT$metrics$fp.PCGII_top25_lambda2[ind,]/(RESULT$metrics$tp.PCGII_top25_lambda2[ind,]+RESULT$metrics$fp.PCGII_top25_lambda2[ind,]),
                       RESULT$metrics$fp.PCGII_top25_lambda3[ind,]/(RESULT$metrics$tp.PCGII_top25_lambda3[ind,]+RESULT$metrics$fp.PCGII_top25_lambda3[ind,]),
                       RESULT$metrics$fp.PCGII_random10_lambda1[ind,]/(RESULT$metrics$tp.PCGII_random10_lambda1[ind,]+RESULT$metrics$fp.PCGII_random10_lambda1[ind,]),
                       RESULT$metrics$fp.PCGII_random10_lambda2[ind,]/(RESULT$metrics$tp.PCGII_random10_lambda2[ind,]+RESULT$metrics$fp.PCGII_random10_lambda2[ind,]),
                       RESULT$metrics$fp.PCGII_random10_lambda3[ind,]/(RESULT$metrics$tp.PCGII_random10_lambda3[ind,]+RESULT$metrics$fp.PCGII_random10_lambda3[ind,]),
                       
                       RESULT$metrics$fp.FGGM_lambda1[ind,]/(RESULT$metrics$tp.FGGM_lambda1[ind,]+RESULT$metrics$fp.FGGM_lambda1[ind,]),
                       RESULT$metrics$fp.FGGM_lambda2[ind,]/(RESULT$metrics$tp.FGGM_lambda2[ind,]+RESULT$metrics$fp.FGGM_lambda2[ind,]),
                       RESULT$metrics$fp.FGGM_lambda3[ind,]/(RESULT$metrics$tp.FGGM_lambda3[ind,]+RESULT$metrics$fp.FGGM_lambda3[ind,])
                       
  )
  
  colnames(fdr)=c("ENF.p","ENF.q","MLE.p","MLE.q",
                    "cLevel_lambda1","cLevel_lambda2","cLevel_lambda3",
                  "PCGII_all_lambda1","PCGII_all_lambda2","PCGII_all_lambda3",
                  "PCGII_top5_lambda1","PCGII_top5_lambda2","PCGII_top5_lambda3",
                  "PCGII_top15_lambda1","PCGII_top15_lambda2","PCGII_top15_lambda3",
                  "PCGII_top25_lambda1","PCGII_top25_lambda2","PCGII_top25_lambda3",
                  "PCGII_random10_lambda1","PCGII_random10_lambda2","PCGII_random10_lambda3",
                  "FGGM_lambda1","FGGM_lambda2","FGGM_lambda3")
  fdr[is.na(fdr)] <- 0
  fdr
}
