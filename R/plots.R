ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>%
  ggplot(aes(x=FPR, y=TPR,col=methods)) + geom_line()+
  geom_abline(slope = 1, intercept = 0, size=0.3)+ggtitle("N=60, P=100, e=4")

fdr.table %>% ggplot(aes(x=FDR, y=fdr, col=methods))+geom_line()+ggtitle("N=60, P=100, e=4")


# plots 

library(ggplot2)
library(tidyverse)
library(patchwork)

load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n60_p100_e4_min_beta0.3_max_beta0.9.RData")

