# Example

# load required libraries and functions
library(igraph)
library(tidyverse)
library(corpcor)
library(glmnet)
library(mvtnorm)

source("./R/PCGII.R")
source("./R/Utility.R")

# Simulating data
set.seed(1234567)
n=50 # sample size
p=30 # number of nodes

g=sample_pa(p, power=1, m=1, directed = FALSE) # undirected scale-free network with the power of the preferential attachment set as 1, the number of edges to add in each time step set as 2.
plot(g, vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5) # visulize simulated network structure
omega=as_adjacency_matrix(g) %>% as.matrix() # precision matrix structure corresponding to the simulated scale-free networ
for(h1 in 1:(p-1)){
  for(h2 in (h1+1):p){
    if(omega[h1,h2]!=0){
      temp=runif(1, 0.2, 0.5)*sample(c(-1,1),size=1) # randomly assign connection strength, i.e. partial correlations
      omega[h1,h2]=temp
      omega[h2,h1]=temp
    }
  }
}

diag(omega)=rowSums(abs(omega)) # make sure precision matrix is positive definite
diag(omega)=diag(omega)+0.10

Sigma=solve(omega) # population covariance matrix, which is used to generate data
X = rmvnorm(n = n, sigma = Sigma) # simulate expression data  

lam=2*sqrt(log(p)/n) ## fixed lambda

# directed prior network
prior_set=matrix(data=c(9,15, 3,4, 5,24, 16,20, 25,22, 28,8, 11,4), nrow=7, ncol=2, byrow = TRUE)
colnames(prior_set)=c("row", "col")
PCGII_out1=PCGII(df=X, prior=as.data.frame(prior_set), lambda = lam)
inference_out1=inference(list=PCGII_out1)
inference_out1$sigs # a data frame of pairs of nodes with significant partial correlations  
inference_out1$sigs %>% sigs2mat(P=p)  %>% 
  graph_from_adjacency_matrix(mode = "undirected") %>%
  plot(vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)

# stronger assumptions, undirected prior network
PCGII_out2=PCGII(df=X, prior=double_prior(prior_set), lambda = lam)
inference_out2=inference(list=PCGII_out2)
inference_out2$sigs # a data frame of pairs of nodes with significant partial correlations  
inference_out2$sigs %>% sigs2mat(P=p)  %>% 
  graph_from_adjacency_matrix(mode = "undirected") %>%
  plot(vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)


# More examples

## Centering within groups
load("./Data/simulated_data_twogroups.RData")
X %>% group_by(treatment) %>% summarise_all(mean)
n=nrow(X)
p=X %>% select(-treatment) %>% ncol()
temp=X %>% group_by(treatment) %>% summarise_all(mean) %>% select(-treatment) %>% as.matrix()
temp=temp[rep(1:nrow(temp), each=n/2),]
X_centered_by_group=X[,1:p]-temp
## network analysis of centered data
prior_set=matrix(data=c(9,15, 3,4, 5,24, 16,20, 25,22, 28,8, 11,4), nrow=7, ncol=2, byrow = TRUE)
colnames(prior_set)=c("row", "col")
lam=2*sqrt(log(p)/n) ## fixed lambda
PCGII_out2=PCGII(df=X_centered_by_group, prior=double_prior(prior_set), lambda = lam)
inference_out2=inference(list=PCGII_out2)
inference_out2$sigs # a data frame of pairs of nodes with significant partial correlations  
inference_out2$sigs %>% sigs2mat(P=p)  %>% 
  graph_from_adjacency_matrix(mode = "undirected") %>%
  plot(vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)


## Visualization
load("./Data/simulated_data_twogroups.RData")
n=nrow(X)
p=X %>% select(-treatment) %>% ncol()
temp=X %>% group_by(treatment) %>% summarise_all(mean) %>% select(-treatment) %>% as.matrix()
temp=temp[rep(1:nrow(temp), each=n/2),]
X_centered_by_group=X[,1:p]-temp
prior_set=matrix(data=c(9,15, 3,4, 5,24, 16,20, 25,22, 28,8, 11,4), nrow=7, ncol=2, byrow = TRUE) # prior set 
colnames(prior_set)=c("row", "col")
lam=2*sqrt(log(p)/n) ## fixed lambda
PCGII_out=PCGII(df=X_centered_by_group, prior=double_prior(prior_set), lambda = lam)
set.seed(333)
nodenames=colnames(X[,1:p])
inference_out=inference(PCGII_out)$sigs %>% mutate(weight=0) 
# assign weights to each connected pair of nodes, weight is equal to corresponding estimated partial correlatio 
for (k in 1:nrow(inference_out)) {
  inference_out$weight[k]=PCGII_out$Est[inference_out$row[k],inference_out$col[k]]
}
my_link=inference_out  %>%
  transform(row = pmin(row, col), col = pmax(row, col)) %>% 
  arrange(row, col) %>% 
  unique() %>%
  mutate(type=ifelse(weight<0,2,1))
colnames(my_link)[1:2]=c("from","to")
my_node=cbind.data.frame(id=1:p, gene=nodenames, types=c(rep(1,10),rep(2,10),rep(3,10)),type.label=c(rep("gene set 1",10),rep("gene set 2",10),rep("gene set 3",10)))
my_net <- graph_from_data_frame(d=my_link, vertices=my_node, directed=F) 

Ecolrs=c("skyblue","pink")
Vcolrs=c("gray80", "tomato", "gold")
V(my_net)$color <- Vcolrs[V(my_net)$types]
E(my_net)$width <- abs(E(my_net)$weight)*2
plot(my_net, edge.arrow.size=.2, edge.color=Ecolrs[E(my_net)$type],
     vertex.frame.color="#ffffff",
     vertex.label=V(my_net)$gene, vertex.label.color="black",
     layout=layout_in_circle(my_net)) 
