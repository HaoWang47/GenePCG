# GenePCGII

Gene network reconstruction by Partial Correlation Graph with Information Incorporation (PCGII)

### Authors:
> Hao Wang, Yumou Qiu and Peng Liu.

### Contact:
> [halewang@iastate.edu] (Hao Wang)

### Citation:
> Wang, H., Qiu, Y.\*, Guo, H., Yin, Y., Liu, P.\*, 2022+. Constructing Large Scale Gene Networks by Partial Correlation
Graphs with Information Incorporation. To be submitted.

> Qiu, Y., & Zhou, X. H. (2020). Estimating c-level partial correlation graphs with application to brain imaging. Biostatistics (Oxford, England), 21(4), 641–658. https://doi.org/10.1093/biostatistics/kxy076
---

This is a tutorial script for researchers who are interested in applying PCGII on omics data to learn the direct association structure of omics features. The main function `PGCII()` takes a biologically pre-processed expression data matrix as input, and returns a list of statistics (estimates and test statistics). The function `inference()` takes the list returned by `PGCII()` as input and conduct simultaneous test to identify significant partial correlations with False Discovery Rate (FDR) controlled at a pre-determined nominal level (0.05 by default). 

# Installation and Package loading
```r
# R version is required >= 4.1.2
# When the first time to use the package, please make sure all required packages are installed under your R environment
> library(igraph)
> library(tidyverse)
> library(corpcor)
> library(glmnet)
> source("./R/PCGII.R")
> source("./R/Utility.R")
```

# Data Preparation

To apply PCGII on omics dataset, some data pre-processing is necessary. For example, gene expression data is supposed to be normalized. Depending on the biological assumptions, treatment effects (group means) can be subtracted ahead if the covariance structures are believe to be the same across treatment groups. Mathematical standardization is not required.

Example:

```r
> load("./Data/simulated_data_twogroups.RData")
> X %>% group_by(treatment) %>% summarise_all(mean)
> n=nrow(X) # sample size
> p=X %>% select(-treatment) %>% ncol() # number of nodes
> temp=X %>% group_by(treatment) %>% summarise_all(mean) %>% select(-treatment) %>% as.matrix() # mean expression matrix by treatment groups
> temp=temp[rep(1:nrow(temp), each=n/2),] # duplicating the mean matrix
> X_centered_by_group=X[,1:p]-temp # expression data centered within treatment groups, ready for network analysis
```

# Visualization

To visualize the reconstructed biological network, we apply package `igraph`.

Example: 

```r
> load("./Data/simulated_data_twogroups.RData")
> n=nrow(X)
> p=X %>% select(-treatment) %>% ncol()
> temp=X %>% group_by(treatment) %>% summarise_all(mean) %>% select(-treatment) %>% as.matrix()
> temp=temp[rep(1:nrow(temp), each=n/2),]
> X_centered_by_group=X[,1:p]-temp
> prior_set=matrix(data=c(9,15, 3,4, 5,24, 16,20, 25,22, 28,8, 11,4), nrow=7, ncol=2, byrow = TRUE) # prior set 
> colnames(prior_set)=c("row", "col")
> lam=2*sqrt(log(p)/n) ## fixed lambda
> PCGII_out=PCGII(df=X_centered_by_group, prior=double_prior(prior_set), lambda = lam)
> set.seed(333)
> nodenames=colnames(X[,1:p])
> inference_out=inference(PCGII_out)$sigs %>% mutate(weight=0) 
> # assign weights to each connected pair of nodes, weight is equal to corresponding estimated partial correlatio 
> for (k in 1:nrow(inference_out)) {
+   inference_out$weight[k]=PCGII_out$Est[inference_out$row[k],inference_out$col[k]]
+ }
> my_link=inference_out  %>%
+   transform(row = pmin(row, col), col = pmax(row, col)) %>% 
+   arrange(row, col) %>% 
+   unique() %>%
+   mutate(type=ifelse(weight<0,2,1))
> colnames(my_link)[1:2]=c("from","to")
> my_node=cbind.data.frame(id=1:p, gene=nodenames, types=c(rep(1,10),rep(2,10),rep(3,10)),type.label=c(rep("gene set 1",10),rep("gene set 2",10),rep("gene set 3",10)))
> my_net <- graph_from_data_frame(d=my_link, vertices=my_node, directed=F) 
> 
> Ecolrs=c("skyblue","pink")
> Vcolrs=c("gray80", "tomato", "gold")
> V(my_net)$color <- Vcolrs[V(my_net)$types]
> E(my_net)$width <- abs(E(my_net)$weight)*2
> plot(my_net, edge.arrow.size=.2, edge.color=Ecolrs[E(my_net)$type],
+      vertex.frame.color="#ffffff",
+      vertex.label=V(my_net)$gene, vertex.label.color="black",
+      layout=layout_in_circle(my_net)) 
```

![Visualization Result](./figs/Example_Visualization.png)

# Example
```r
> # Simulating data
> set.seed(1234567)
> n=50 # sample size
> p=30 # number of nodes
> 
> g=sample_pa(p, power=1, m=1, directed = FALSE) # undirected scale-free network with the power of the preferential attachment set as 1, the number of edges to add in each time step set as 2.
> plot(g, vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5) # visulize simulated network structure
> omega=as_adjacency_matrix(g) %>% as.matrix() # precision matrix structure corresponding to the simulated scale-free networ
> for(h1 in 1:(p-1)){
+   for(h2 in (h1+1):p){
+     if(omega[h1,h2]!=0){
+       temp=runif(1, 0.2, 0.5)*sample(c(-1,1),size=1) # randomly assign connection strength, i.e. partial correlations
+       omega[h1,h2]=temp
+       omega[h2,h1]=temp
+     }
+   }
+ }
> 
> diag(omega)=rowSums(abs(omega)) # make sure precision matrix is positive definite
> diag(omega)=diag(omega)+0.10
> 
> Sigma=solve(omega) # population covariance matrix, which is used to generate data
> X = rmvnorm(n = n, sigma = Sigma) # simulate expression data  
> 
> lam=2*sqrt(log(p)/n) ## fixed lambda
> 
> # directed prior network
> prior_set=matrix(data=c(9,15, 3,4, 5,24, 16,20, 25,22, 28,8, 11,4), nrow=7, ncol=2, byrow = TRUE)
> colnames(prior_set)=c("row", "col")
> PCGII_out1=PCGII(df=X, prior=as.data.frame(prior_set), lambda = lam)
> inference_out1=inference(list=PCGII_out1)
> inference_out1$sigs # a data frame of pairs of nodes with significant partial correlations  
> inference_out1$sigs %>% sigs2mat(P=p)  %>% 
+   graph_from_adjacency_matrix(mode = "undirected") %>%
+   plot(vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)
> 
> # stronger assumptions, undirected prior network
> PCGII_out2=PCGII(df=X, prior=double_prior(prior_set), lambda = lam)
> inference_out2=inference(list=PCGII_out2)
> inference_out2$sigs # a data frame of pairs of nodes with significant partial correlations
> inference_out2$sigs %>% sigs2mat(P=p)  %>% 
+   graph_from_adjacency_matrix(mode = "undirected") %>%
+   plot(vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)
```

For more examples, please see R script `demo.R`.

