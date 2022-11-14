# GenePCGII

Gene network reconstruction by Partial Correlation Graph with Information Incorporation (PCGII)

### Authors:
> Hao Wang, Yumou Qiu and Peng Liu.

### Contact:
> [halewang@iastate.edu] (Hao Wang)

### Citation:
> Wang, H., Qiu, Y.\*, Guo, H., Yin, Y., Liu, P.\*, 2022+. Constructing Large Scale Gene Networks by Partial Correlation
Graphs with Information Incorporation
> Qiu, Y., & Zhou, X. H. (2020). Estimating c-level partial correlation graphs with application to brain imaging. Biostatistics (Oxford, England), 21(4), 641–658. https://doi.org/10.1093/biostatistics/kxy076
---

This is a tutorial script for researchers who are interested in applying PCGII on omics data to learn the direct association structure of omics features. 

# Installation and Package loading
```r
# R version is required >= 4.1.2
# When the first time to use the package, please make sure all required packages are installed under your R environment, if not, please use commands below to install
> library(igraph)
> library(tidyverse)
> library(GeneNet)
> library(FastGGM)
> library(corpcor)
> library(glmnet)
> install.packages("leaps")
> install.packages("MASS")
# install "devtools" package in your R environment
> devtools::install_github("shanwai1234/GPWAS")
> library(GPWAS)
```
# Data Preparation
