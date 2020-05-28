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
  S=list()
  for (i in 1:reps) {
    bd=matrix(runif(1,min.beta,max.beta)*sample(1,c(1,-1)),blocksize,blocksize)
    diag(bd)=runif(1,1.25,3)
    S[[i]]=bd
  }
  as.matrix(Matrix::bdiag(S))
}

# Following code generate specific network structures and generate multivariate normal random data from them. 

# Scale Free ------------------------------------------------------------------------
# simulate data from BA model (Barabási–Albert model) (Scale Free Network, m=1 or 2)
# BA model: https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model
#
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta=0.1
set.seed(20200520)

for(n in nl){    
  for(p in pl){
    for(e in 1:2){
      PCOR = list()        
      X = list()
      for(i in 1:100){            
        pcor=as(igraph.to.graphNEL(barabasi.game(n=p,m=e)), "matrix")
        w=length(which(pcor==1 & lower.tri(pcor)))
        pcor[(pcor==1 & lower.tri(pcor))]=runif(n = w, min = min.beta, max = 1)*sample(x=c(-1,1), size=w, replace = T)
        pcor=makeSymm(pcor)
        diag(pcor)=1
        D=diag(1.25, p, p)
        Sigma=solve(D %*% pcor %*% D)
        #is.positive.definite(Sigma)
        #pppp=round(cor2pcor(Sigma),5)
        id.shuffle = sample(1:p, p, replace = FALSE)                        
        Xi = rmvnorm(n = n, sigma = Sigma)            
        Xi = Xi[,id.shuffle]            
        pcor = pcor[id.shuffle, id.shuffle]                        
        PCOR[[i]] = pcor           
        X[[i]] = Xi
      }
      save(PCOR, X, file=sprintf("simu_data/ScaleFree_simu_n%d_p%d_e%d_min_beta%g.RData", n, p, e, min.beta)) 
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
      PCOR = list()        
      X = list()
      for(i in 1:100){            
        Sigma=makeBlockDiag(blocksize = e, p=p, min.beta = min.beta, max.beta = max.beta)
        pcor=cor2pcor(Sigma)
        #is.positive.definite(Sigma)
        #pppp=round(cor2pcor(Sigma),5)
        id.shuffle = sample(1:p, p, replace = FALSE)                        
        Xi = rmvnorm(n = n, sigma = Sigma)            
        Xi = Xi[,id.shuffle]            
        pcor = pcor[id.shuffle, id.shuffle]                        
        PCOR[[i]] = pcor           
        X[[i]] = Xi
      }
      save(PCOR, X, file=sprintf("simu_data/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta)) 
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
    for(eta in c(0.01,0.02,0.03)){
      PCOR = list()        
      X = list()
      for(i in 1:100){            
        pcor = matrix(0, p, p)     
        w  = which(lower.tri(pcor))      
        ruf = runif(n = length(w), min = min.beta, max = 1)      
        pcor[w] = rbinom(n = length(w), size = 1, prob = eta) * ruf
        pcor=makeSymm(pcor)
        diag(pcor)=1
        D=diag(1.25, p, p)
        Sigma=solve(D %*% pcor %*% D)
        #is.positive.definite(Sigma)
        #pppp=round(cor2pcor(Sigma),5)
        id.shuffle = sample(1:p, p, replace = FALSE)                        
        Xi = rmvnorm(n = n, sigma = Sigma)            
        Xi = Xi[,id.shuffle]            
        pcor = pcor[id.shuffle, id.shuffle]                        
        PCOR[[i]] = pcor           
        X[[i]] = Xi
      }
      save(PCOR, X, file=sprintf("simu_data/Random_simu_n%d_p%d_eta%g_min_beta%g.RData", n, p, eta, min.beta)) 
    }
  }
}

