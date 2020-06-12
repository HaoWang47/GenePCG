# Last updated: 2020/06/08

# reference: https://rpubs.com/lgadar/generate-graphs

# Gene PCG -----

library(MASS)
library(igraph)
library(corpcor)
library(mvtnorm)

# Functions -----

makeSymm=function(m){
  m[upper.tri(m)]=t(m)[upper.tri(m)]
  return(m)
}

makeBlockDiag=function(blocksize=4, p=100, min.beta=0.1, max.beta=1){ # blocksize has to be a factor of p
  reps=p/blocksize
  Max=matrix(0, p, p)
  while (!is.positive.definite(Max)) {
    S=list()
    for (i in 1:reps) {
      bd=matrix(runif(1,min.beta,max.beta)*sample(c(1,-1),1),blocksize,blocksize)
      diag(bd)=runif(1,1.25,3)
      S[[i]]=bd
    }
    Max=as.matrix(Matrix::bdiag(S))
  }
  Max
}

# Following code generate specific network structures and generate multivariate normal random data from them. 

# Scale Free ------------------------------------------------------------------------
# simulate data from BA model (Barabási–Albert model) (Scale Free Network, m=1 or 2)
# BA model: https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model
#
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta=0.3
set.seed(20200520)

for(n in nl){    
  for(p in pl){
    for(e in 1:2){
      
      omega = matrix(0, p, p) 
      while (! is.positive.definite(omega)){
        omega=as(igraph.to.graphNEL(barabasi.game(n=p,m=e)), "matrix")
        w=which(omega==1 & lower.tri(omega))
        omega[w]=runif(n = length(w), min = min.beta, max = 1)*sample(x=c(-1,1), size=length(w), replace = T)
        omega=makeSymm(omega)
        diag(omega)=4.5
      }
      
      Sigma=solve(omega)
      pcor=diag(p)
      for(i in 2:p){
        for (j in (i-1):1){
          pcor[i,j]=-omega[i,j]/sqrt(omega[i,i]*omega[j,j])
          pcor[j,i]=pcor[i,j]
        }
      }
      
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(pcor, Sigma, X, file=sprintf("simu_data/ScaleFree_simu_n%d_p%d_e%d_min_beta%g.RData", n, p, e, min.beta)) 
      }
  }
}


# Block Diagonal ---------------------------------------------------------
# simulate data from Block Diagonal model (Block size = 4, 10)
# 
#
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta=0.3
max.beta=1
set.seed(20200520)

for(n in nl){    
  for(p in pl){
    for(e in c(4,10)){
    Sigma=makeBlockDiag(blocksize = e, p=p, min.beta = min.beta, max.beta = max.beta)
      pcor=cor2pcor(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(pcor, Sigma, X, file=sprintf("simu_data/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta)) 
    }
  }
}



# Random Graph ---------------------------------------------------------------
# simulate data from Random Graph model (sparsity eta=0.01, 0.02, 0.03)
# 
# 
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta=0.3
set.seed(20200520)

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.03,0.05)){
      omega = matrix(0, p, p) 
      while (! is.positive.definite(omega)){
        omega=as(igraph.to.graphNEL(erdos.renyi.game(n=p,p=eta,type="gnp")), "matrix")
        w=which(omega==1 & lower.tri(omega))
        omega[w]=runif(n = length(w), min = min.beta, max = 1)*sample(x=c(-1,1), size=length(w), replace = T)
        omega=makeSymm(omega)
        diag(omega)=4.5
      }
      Sigma=solve(omega)
      pcor=diag(p)
      for(i in 2:p){
        for (j in (i-1):1){
          pcor[i,j]=-omega[i,j]/sqrt(omega[i,i]*omega[j,j])
          pcor[j,i]=pcor[i,j]
        }
      }
      
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(pcor, Sigma, X, file=sprintf("simu_data/Random_simu_n%d_p%d_eta%g_min_beta%g.RData", n, p, eta, min.beta)) 
    }
  }
}

