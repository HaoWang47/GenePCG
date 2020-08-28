# Last updated: 2020/07/09

# These code used for final data simulation *******************************

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
#set.seed(20200611)
set.seed(20200627)

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
      print(paste("Setting: n=",n, ", p=",p, ", e=", e,"; True Sparsity:", length(which(sm2vec(pcor)!=0))/length(which(lower.tri(pcor)))) )
    }
  }
}




# Block Diagonal (all positives) --------------------------------------------
# simulate data from Block Diagonal model (Block size = 4, 10)
# 
#
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

nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta=0.3
max.beta=.9
#set.seed(20200611)
set.seed(20200627)

for(n in nl){    
  for(p in pl){
    for(e in c(4,8,10)){
      Sigma=makeBlockDiag(blocksize=e, p=p, min.beta=min.beta, max.beta=max.beta)
      pcor=sigma2pcor(Sigma)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(pcor, Sigma, X, file=sprintf("simu_data_v4/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta)) 
      print(paste("Setting: n=",n, ", p=",p, ", e=", e,"; True Sparsity:", length(which(sm2vec(pcor)!=0))/length(which(lower.tri(pcor))), "; df=", ceiling(p*length(which(sm2vec(pcor)!=0))/length(which(lower.tri(pcor)))*2)) )
    }
  }
}


# Scale Free ------------------------------------------------------------------------
# simulate data from BA model (Barabási–Albert model) (Scale Free Network, m=1 or 2)
# BA model: https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model
#
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
min.beta=0.3
set.seed(20200627)

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
      for(i in 1:3){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(A, pcor, Sigma, X, file=sprintf("simu_data_v4/ScaleFree_simu_n%d_p%d_e%d_min_beta%g.RData", n, p, e, min.beta)) 
      print(paste("Setting: n=",n, ", p=",p, ", e=", e,"; True Sparsity:", length(which(sm2vec(pcor)!=0))/length(which(lower.tri(pcor))), "; df=", ceiling(p*length(which(sm2vec(pcor)!=0))/length(which(lower.tri(pcor)))*2)) )
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
#set.seed(20200624)
set.seed(20200627)

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.02,0.03)){
      A = matrix(0, p, p)     
      w  = which(lower.tri(A))      
      ruf = runif(n = length(w), min = .1, max = 1)      
      A[w] = rbinom(n = length(w), size = 1, prob = eta) * ruf
      Sigh = solve(diag(p) - A)          
      Sigma = Sigh %*% t(Sigh)
      pcor=round(sigma2pcor(Sigma),5)
      X = list()
      for(i in 1:100){            
        Xi = rmvnorm(n = n, sigma = Sigma)            
        X[[i]] = Xi
      }
      save(pcor, Sigma, X, file=sprintf("simu_data_v4/Random_simu_n%d_p%d_eta%g_min_beta%g.RData", n, p, eta, min.beta)) 
      print(paste("Setting: n=",n, ", p=",p, ", eta=", eta,"; True Sparsity:", length(which(sm2vec(pcor)!=0))/length(which(lower.tri(pcor))), "; df=", ceiling(p*length(which(sm2vec(pcor)!=0))/length(which(lower.tri(pcor)))*2)) )
    }
  }
}


#----         
Blockdiag
[1] "Setting: n= 60 , p= 100 , e= 4 ; True Sparsity: 0.0303030303030303 ; df= 7"
[1] "Setting: n= 60 , p= 100 , e= 8 ; True Sparsity: 0.0736842105263158 ; df= 15"
[1] "Setting: n= 60 , p= 100 , e= 10 ; True Sparsity: 0.0909090909090909 ; df= 19"
[1] "Setting: n= 60 , p= 200 , e= 4 ; True Sparsity: 0.0150753768844221 ; df= 7"
[1] "Setting: n= 60 , p= 200 , e= 8 ; True Sparsity: 0.0351758793969849 ; df= 15"
[1] "Setting: n= 60 , p= 200 , e= 10 ; True Sparsity: 0.0452261306532663 ; df= 19"
[1] "Setting: n= 80 , p= 100 , e= 4 ; True Sparsity: 0.0303030303030303 ; df= 7"
[1] "Setting: n= 80 , p= 100 , e= 8 ; True Sparsity: 0.0736842105263158 ; df= 15"
[1] "Setting: n= 80 , p= 100 , e= 10 ; True Sparsity: 0.0909090909090909 ; df= 19"
[1] "Setting: n= 80 , p= 200 , e= 4 ; True Sparsity: 0.0150753768844221 ; df= 7"
[1] "Setting: n= 80 , p= 200 , e= 8 ; True Sparsity: 0.0351758793969849 ; df= 15"
[1] "Setting: n= 80 , p= 200 , e= 10 ; True Sparsity: 0.0452261306532663 ; df= 19" 

Scale free
[1] "Setting: n= 60 , p= 100 , e= 1 ; True Sparsity: 0.02 ; df= 4"
[1] "Setting: n= 60 , p= 100 , e= 2 ; True Sparsity: 0.0482828282828283 ; df= 10"
[1] "Setting: n= 60 , p= 200 , e= 1 ; True Sparsity: 0.01 ; df= 4"
[1] "Setting: n= 60 , p= 200 , e= 2 ; True Sparsity: 0.0261809045226131 ; df= 11"
[1] "Setting: n= 80 , p= 100 , e= 1 ; True Sparsity: 0.02 ; df= 4"
[1] "Setting: n= 80 , p= 100 , e= 2 ; True Sparsity: 0.0513131313131313 ; df= 11"
[1] "Setting: n= 80 , p= 200 , e= 1 ; True Sparsity: 0.01 ; df= 4"
[1] "Setting: n= 80 , p= 200 , e= 2 ; True Sparsity: 0.0257788944723618 ; df= 11"

Random
[1] "Setting: n= 60 , p= 100 , eta= 0.01 ; True Sparsity: 0.00787878787878788 ; df= 2"
[1] "Setting: n= 60 , p= 100 , eta= 0.02 ; True Sparsity: 0.0327272727272727 ; df= 7"
[1] "Setting: n= 60 , p= 100 , eta= 0.03 ; True Sparsity: 0.0533333333333333 ; df= 11"
[1] "Setting: n= 60 , p= 200 , eta= 0.01 ; True Sparsity: 0.0145226130653266 ; df= 6"
[1] "Setting: n= 60 , p= 200 , eta= 0.02 ; True Sparsity: 0.0383417085427136 ; df= 16"
[1] "Setting: n= 60 , p= 200 , eta= 0.03 ; True Sparsity: 0.0911557788944724 ; df= 37"
[1] "Setting: n= 80 , p= 100 , eta= 0.01 ; True Sparsity: 0.0133333333333333 ; df= 3"
[1] "Setting: n= 80 , p= 100 , eta= 0.02 ; True Sparsity: 0.0367676767676768 ; df= 8"
[1] "Setting: n= 80 , p= 100 , eta= 0.03 ; True Sparsity: 0.0511111111111111 ; df= 11"
[1] "Setting: n= 80 , p= 200 , eta= 0.01 ; True Sparsity: 0.0142211055276382 ; df= 6"
[1] "Setting: n= 80 , p= 200 , eta= 0.02 ; True Sparsity: 0.0437688442211055 ; df= 18"
[1] "Setting: n= 80 , p= 200 , eta= 0.03 ; True Sparsity: 0.0823618090452261 ; df= 33"

which(A!=0)

which(round(omega,3)!=0 & lower.tri(omega))

which(round(pcor,3)!=0 & lower.tri(pcor))
