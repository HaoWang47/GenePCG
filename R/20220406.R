library(tidyverse)
library(GeneNet)
library(FastGGM)
library(corpcor)
library(glmnet)

source("R/Qiu.R")
source("R/PCGII.R")
source("R/shrinkagefunctions.R") # https://github.com/V-Bernal/GGM-Shrinkage

# Scale Free ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
el = c(1, 2)
for (n in nl) {
  for (p in pl) {
    for (e in el) {
      load( file=sprintf("simu_data/ScaleFree_simu_n%d_p%d_e%d.RData", n, p, e) )
      degf=n/5
      RESULT=model.eval(X, omega, rep = 3, degfree = degf, prior_prop=0.1)
      save(RESULT, file=sprintf("results/ScaleFree_n%d_p%d_e%d.RData", n, p, e) )
    }
  }
}


# Random ----
nl = c(60, 80) # Sample Size
pl = c(100, 200)  # Number of Genes
etal = c(.01, .02, .03)
for (n in nl) {
  for (p in pl) {
    for (eta in etal) {
      load( file=sprintf("simu_data/Random_simu_n%d_p%d_eta%g.RData", n, p, eta) )
      degf=20
      RESULT=model.eval(X, omega, rep = 3, degfree = degf, prior_prop=0.1)
      save(RESULT, file=sprintf("results/ScaleFree_n%d_p%d_e%d.RData", n, p, e) )
    }
  }
}




