# Last updated: 2020/06/11

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

sigma2pcor=function(sigma){
  if (is.positive.definite(sigma)) {
    omega=solve(sigma)
    p=dim(omega)[1]
    pcor=diag(p)
    for(i in 2:p){
      for (j in (i-1):1){
        pcor[i,j]=-omega[i,j]/sqrt(omega[i,i]*omega[j,j])
        pcor[j,i]=pcor[i,j]
      }
    }
    return(pcor)
  } else stop('sigma not positive definite')
}


# Following code generate specific network structures and generate multivariate normal random data from them. 

# Scale Free ------------------------------------------------------------------------
# simulate data from BA model (Barabási–Albert model) (Scale Free Network, m=1 or 2)
# BA model: https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model
#
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta=0.3
set.seed(20200611)

for(n in nl){    
  for(p in pl){
    for(e in 1:2){
      A=as(igraph.to.graphNEL(barabasi.game(n=p,m=e)), "matrix")
      w=which(A!=0)
      A[w]=runif(n = length(w), min = min.beta, max = 1)*sample(x=c(-1,1), size=length(w), replace = T)      
      Sigh = solve(diag(p) - A)            
      Sigma = Sigh %*% t(Sigh)
      pcor=round(sigma2pcor(Sigma),5)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(A, pcor, Sigma, X, file=sprintf("simu_data_v2/ScaleFree_simu_n%d_p%d_e%d_min_beta%g.RData", n, p, e, min.beta)) 
      }
  }
}


# Block Diagonal ---------------------------------------------------------
# simulate data from Block Diagonal model (Block size = 4, 10)
# 
#
makeBlockDiag=function(blocksize=4, p=100, min.beta=0.1, max.beta=1){ # blocksize has to be a factor of p
  reps=p/blocksize
  S=list()
  for (i in 1:reps) {
    bd=matrix(runif(1,min.beta,max.beta)*sample(c(1,-1),1),blocksize,blocksize)
    S[[i]]=bd
  }
  mat=as.matrix(Matrix::bdiag(S))
  diag(mat)=1
  mat[upper.tri(mat)]=0
  mat
}

nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta=0.05
max.beta=.5
set.seed(20200611)

for(n in nl){    
  for(p in pl){
    for(e in c(4,8,10)){
      A=makeBlockDiag(blocksize=e, p=p, min.beta=min.beta, max.beta=max.beta)
      Sigh = solve(A)            
      Sigma = Sigh %*% t(Sigh)
      pcor=round(sigma2pcor(Sigma),5)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(A, pcor, Sigma, X, file=sprintf("simu_data_v2/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta)) 
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
      A = matrix(0, p, p)     
      w  = which(lower.tri(A))      
      ruf = runif(n = length(w), min = .1, max = 1)      
      A[w] = rbinom(n = length(w), size = 1, prob = 4/p) * ruf
      Sigh = solve(diag(p) - A)          
      Sigma = Sigh %*% t(Sigh)
      pcor=round(sigma2pcor(Sigma),5)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(pcor, Sigma, X, file=sprintf("simu_data_v2/Random_simu_n%d_p%d_eta%g_min_beta%g.RData", n, p, eta, min.beta)) 
    }
  }
}

which(A!=0)

which(round(omega,3)!=0 & lower.tri(omega))

which(round(pcor,3)!=0 & lower.tri(pcor))
