# Utility functions


# This function takes PCGII inference results as input, generates a 
# adjacecy matrix corresponding to the significant partial correlations
# and returns a vector of all off-diagonal elements. The function does 
# require package corpcor. 
sigs2vec=function(sigs, P){
  require(corpcor)
  m=matrix(0,P,P)
  for (h in 1: dim(sigs)[1]){
    m[sigs[h,1],sigs[h,2]]=1
  }
  sm2vec(m)
}


# This function takes PCGII inference results as input, generates a 
# adjacecy matrix corresponding to the significant partial correlations
# and returns the matrix. The function does require package corpcor. 
sigs2mat=function(sigs, P){
  require(corpcor)
  m=matrix(0,P,P)
  for (h in 1: dim(sigs)[1]){
    m[sigs[h,1],sigs[h,2]]=1
    m[sigs[h,2],sigs[h,1]]=1
  }
  m
}


# This function takes the estimated p by p partial correlations matrix and 
# corresponding element-wise pvalue matrix as input, and return a matrix where
# each row corresponds to a pair of nodes with the estimated partial correlation
# and p-value. 
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


## An utility function to pre-process the input prior set. This function will ensure the input prior set corresponds to an undirected prior network. If the prior network is believed to be directed, no pre-processing of the prior set is needed. 
## Input:
# prior: the prior set, a k by 2 matrix, in which each row corresponds to a pair of nodes (any omics features) that are connected under prior belief. 
## Output:
# This function returns a pre-processed prior set, in which the connection between any pair of nodes is undirected.
# Remark: this function is not necessary. Prior set should be considered carefully before running the network analysis. If the prior network connections are believed to be undirected while the prior set only includes one way connections for simplicity, this function will duplicate the connections and swap the direction automactically.  
double_prior=function(prior){
  rbind.data.frame(prior, prior %>% transform(row = pmax(row, col), col = pmin(row, col))) %>% 
    arrange(row, col) %>% unique()
}