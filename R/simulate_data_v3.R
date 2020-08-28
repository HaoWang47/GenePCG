# Last updated: 2020/06/11

# reference: https://rpubs.com/lgadar/generate-graphs

# Gene PCG -----

library(MASS)
library(igraph)
library(corpcor)
library(mvtnorm)
library(sgnesR)

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
set.seed(20200611)

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
      save(pcor, Sigma, X, file=sprintf("simu_data_v3/BlockDiag_simu_n%d_p%d_e%d_min_beta%g_max_beta%g.RData", n, p, e, min.beta, max.beta)) 
    }
  }
}




# Scale Free ------------------------------------------------------------------------
# simulate data from BA model (Barabási–Albert model) (Scale Free Network, m=1 or 2)
# BA model: https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model
#
nl = c(60) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20200611)

for(n in nl){    
  for(p in pl){
    for(e in 1:2){
      g <- barabasi.game(n=p, m=1.5, directed=F)
      V(g)$Ppop <- (sample(p,vcount(g), rep=T))
      V(g)$Rpop <- (sample(p, vcount(g), rep=T))
      rp<-new("rsgns.param",time=0,stop_time=1000,readout_interval=500)
      rc <- c(0.002, 0.005, 0.005, 0.005, 0.01, 0.02)
      rsg <- new('rsgns.data',network=g, rconst=rc)
      pcor=as(igraph.to.graphNEL(g), "matrix")
      X = list()
      for(i in 1:10){            
        Xi<- rsgns.rn(rsg, rp, timeseries=F, sample=n) 
        X[[i]] = t(Xi$expression)
        print(paste0("n=", n, ", p=", p, ", e=", e, ", i=", i, ", done!!!"))
      }
      save(pcor, X, file=sprintf("simu_data_v3/ScaleFree_simu_n%d_p%d_e%d.RData", n, p, e)) 
      }
  }
}

# Random Graph ---------------------------------------------------------------
# simulate data from Random Graph model (sparsity eta=0.01, 0.02, 0.03)
# 
# 
nl = c(60) # Sample Size
pl = c(100, 200)  # Number of Genes
set.seed(20200520)

for(n in nl){    
  for(p in pl){
    for(eta in c(0.01,0.03,0.05)){
      g <- erdos.renyi.game(n=p, p=eta, type="gnp", directed=F)
      V(g)$Ppop <- (sample(p,vcount(g), rep=T))
      V(g)$Rpop <- (sample(p, vcount(g), rep=T))
      rp<-new("rsgns.param",time=0,stop_time=1000,readout_interval=500)
      rc <- c(0.002, 0.005, 0.005, 0.005, 0.01, 0.02)
      rsg <- new('rsgns.data',network=g, rconst=rc)
      pcor=as(igraph.to.graphNEL(g), "matrix")
      X = list()
      for(i in 1:10){            
        Xi<- rsgns.rn(rsg, rp, timeseries=F, sample=n) 
        X[[i]] = t(Xi$expression)
        print(paste0("n=", n, ", p=", p, ", eta=", eta, ", i=", i, ", done!!!"))
      }
      save(pcor, X, file=sprintf("simu_data_v3/Random_simu_n%d_p%d_eta%g.RData", n, p, eta))
      
    }
  }
}


