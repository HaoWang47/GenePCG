ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>%
  ggplot(aes(x=FPR, y=TPR,col=methods)) + geom_line()+
  geom_abline(slope = 1, intercept = 0, size=0.3)+ggtitle("N=60, P=100, e=4")

fdr.table %>% ggplot(aes(x=FDR, y=fdr, col=methods))+geom_line()+ggtitle("N=60, P=100, e=4")


# plots 

# Library----

library(tidyverse)
library(patchwork)

# Functions----
area_roc=function(x){
  sum=0.001*x[1]/2
  for (i in 1:999){
    sum=(.001*(x[i+1]+x[i])/2)+sum
  }
  sum
}


AUC=function(RESULT){
  ENF=apply(RESULT$ROC$tempstore.ENF, 2, area_roc)
  MLE=apply(RESULT$ROC$tempstore.MLE, 2, area_roc)
  cLevel_theo=apply(RESULT$ROC$tempstore.cLevel_theo, 2, area_roc)
  F_GGM=apply(RESULT$ROC$tempstore.F_GGM, 2, area_roc)
  cLevel_cv=apply(RESULT$ROC$tempstore.cLevel_cv, 2, area_roc)
  result=rbind(
    cbind.data.frame(AUC=ENF, methods="ENF"),
    cbind.data.frame(AUC=MLE, methods="MLE"),
    cbind.data.frame(AUC=cLevel_theo, methods="cPCG.theo"),
    cbind.data.frame(AUC=F_GGM, methods="F_GGM"),
    cbind.data.frame(AUC=cLevel_df, methods="cPCG.df")
  )
  return(result)
}

f1score=function(RESULT){ 
  recall=cbind(RESULT$metrics$tp.enf.p[15,]/(RESULT$metrics$tp.enf.p[15,]+RESULT$metrics$fn.enf.p[15,]), # TPR
               RESULT$metrics$tp.enf.q[15,]/(RESULT$metrics$tp.enf.q[15,]+RESULT$metrics$fn.enf.q[15,]),
               RESULT$metrics$tp.mle.p[15,]/(RESULT$metrics$tp.mle.p[15,]+RESULT$metrics$fn.mle.p[15,]),
               RESULT$metrics$tp.mle.q[15,]/(RESULT$metrics$tp.mle.q[15,]+RESULT$metrics$fn.mle.q[15,]),
               RESULT$metrics$tp.clevel_theo[15,]/(RESULT$metrics$tp.clevel_theo[15,]+RESULT$metrics$fn.clevel_theo[15,]),
               RESULT$metrics$tp.clevel_cv[15,]/(RESULT$metrics$tp.clevel_cv[15,]+RESULT$metrics$fn.clevel_cv[15,]),
               RESULT$metrics$tp.clevel_df[15,]/(RESULT$metrics$tp.clevel_df[15,]+RESULT$metrics$fn.clevel_df[15,])
  )
  precision=cbind( # power
    RESULT$metrics$tp.enf.p[15,]/(RESULT$metrics$tp.enf.p[15,]+RESULT$metrics$fp.enf.p[15,]),
    RESULT$metrics$tp.enf.q[15,]/(RESULT$metrics$tp.enf.q[15,]+RESULT$metrics$fp.enf.q[15,]),
    RESULT$metrics$tp.mle.p[15,]/(RESULT$metrics$tp.mle.p[15,]+RESULT$metrics$fp.mle.p[15,]),
    RESULT$metrics$tp.mle.q[15,]/(RESULT$metrics$tp.mle.q[15,]+RESULT$metrics$fp.mle.q[15,]),
    RESULT$metrics$tp.clevel_theo[15,]/(RESULT$metrics$tp.clevel_theo[15,]+RESULT$metrics$fp.clevel_theo[15,]),
    RESULT$metrics$tp.clevel_cv[15,]/(RESULT$metrics$tp.clevel_cv[15,]+RESULT$metrics$fp.clevel_cv[15,]),
    RESULT$metrics$tp.clevel_df[15,]/(RESULT$metrics$tp.clevel_df[15,]+RESULT$metrics$fp.clevel_df[15,])
  )
  
  precision[which(is.na(precision))]=0
  
  Fscore=2*(precision*recall/(precision+recall)) 
  Fscore[which(is.na(Fscore))]=0
  colnames(Fscore)=c("ENF.p","ENF.q","MLE.p","MLE.q","cPCG.theo","cPCG.cv","cPCG.df")
  Fscore %>% as.data.frame() %>% gather("methods","F1score")
}

my_power=function(RESULT){
  recall=cbind(RESULT$metrics$tp.enf.p[15,]/(RESULT$metrics$tp.enf.p[15,]+RESULT$metrics$fn.enf.p[15,]), # TPR
               RESULT$metrics$tp.enf.q[15,]/(RESULT$metrics$tp.enf.q[15,]+RESULT$metrics$fn.enf.q[15,]),
               RESULT$metrics$tp.mle.p[15,]/(RESULT$metrics$tp.mle.p[15,]+RESULT$metrics$fn.mle.p[15,]),
               RESULT$metrics$tp.mle.q[15,]/(RESULT$metrics$tp.mle.q[15,]+RESULT$metrics$fn.mle.q[15,]),
               RESULT$metrics$tp.clevel_theo[15,]/(RESULT$metrics$tp.clevel_theo[15,]+RESULT$metrics$fn.clevel_theo[15,]),
               RESULT$metrics$tp.clevel_cv[15,]/(RESULT$metrics$tp.clevel_cv[15,]+RESULT$metrics$fn.clevel_cv[15,]),
               RESULT$metrics$tp.clevel_df[15,]/(RESULT$metrics$tp.clevel_df[15,]+RESULT$metrics$fn.clevel_df[15,])
  )
  
  
  power=as.data.frame(recall)
  colnames(power)=c("ENF.p","ENF.q","MLE.p","MLE.q","cPCG.theo","cPCG.cv","cPCG.df")
  power %>% gather("methods","power")
}

# Data----

# BD----

{
### BD 60,100 
load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n60_p100_e4_min_beta0.3_max_beta0.9.RData")

BD60_100_4_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:F_GGM) %>% mutate(Structure="BD4",P=100,N=60)

BD60_100_4_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD4",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_100_4_power=my_power(RESULT) %>% mutate(Structure="BD4",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_100_4_fscore=f1score(RESULT) %>% mutate(Structure="BD4",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_100_4_AUC=AUC(RESULT) %>% mutate(Structure="BD4",P=100,N=60)

load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n60_p100_e8_min_beta0.3_max_beta0.9.RData")

BD60_100_8_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD8",P=100,N=60)

BD60_100_8_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD8",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_100_8_power=my_power(RESULT) %>% mutate(Structure="BD8",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_100_8_fscore=f1score(RESULT) %>% mutate(Structure="BD8",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_100_8_AUC=AUC(RESULT) %>% mutate(Structure="BD8",P=100,N=60)

load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n60_p100_e10_min_beta0.3_max_beta0.9.RData")

BD60_100_10_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD10",P=100,N=60)

BD60_100_10_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD10",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_100_10_power=my_power(RESULT) %>% mutate(Structure="BD10",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_100_10_fscore=f1score(RESULT) %>% mutate(Structure="BD10",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_100_10_AUC=AUC(RESULT) %>% mutate(Structure="BD10",P=100,N=60)

### BD 60,200
load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n60_p200_e4_min_beta0.3_max_beta0.9.RData")

BD60_200_4_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD4",P=200,N=60)

BD60_200_4_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD4",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_200_4_power=my_power(RESULT) %>% mutate(Structure="BD4",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_200_4_fscore=f1score(RESULT) %>% mutate(Structure="BD4",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_200_4_AUC=AUC(RESULT) %>% mutate(Structure="BD4",P=200,N=60)

load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n60_p200_e8_min_beta0.3_max_beta0.9.RData")

BD60_200_8_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD8",P=200,N=60)

BD60_200_8_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD8",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_200_8_power=my_power(RESULT) %>% mutate(Structure="BD8",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_200_8_fscore=f1score(RESULT) %>% mutate(Structure="BD8",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_200_8_AUC=AUC(RESULT) %>% mutate(Structure="BD8",P=200,N=60)

load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n60_p200_e10_min_beta0.3_max_beta0.9.RData")

BD60_200_10_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD10",P=200,N=60)

BD60_200_10_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD10",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_200_10_power=my_power(RESULT) %>% mutate(Structure="BD10",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_200_10_fscore=f1score(RESULT) %>% mutate(Structure="BD10",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD60_200_10_AUC=AUC(RESULT) %>% mutate(Structure="BD10",P=200,N=60)

### BD 80,100
load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n80_p100_e4_min_beta0.3_max_beta0.9.RData")

BD80_100_4_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD4",P=100,N=80)

BD80_100_4_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD4",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_100_4_power=my_power(RESULT) %>% mutate(Structure="BD4",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_100_4_fscore=f1score(RESULT) %>% mutate(Structure="BD4",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_100_4_AUC=AUC(RESULT) %>% mutate(Structure="BD4",P=100,N=80)

load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n80_p100_e8_min_beta0.3_max_beta0.9.RData")

BD80_100_8_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD8",P=100,N=80)

BD80_100_8_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD8",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_100_8_power=my_power(RESULT) %>% mutate(Structure="BD8",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_100_8_fscore=f1score(RESULT) %>% mutate(Structure="BD8",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_100_8_AUC=AUC(RESULT) %>% mutate(Structure="BD8",P=100,N=80)

load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n80_p100_e10_min_beta0.3_max_beta0.9.RData")

BD80_100_10_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD10",P=100,N=80)

BD80_100_10_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD10",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_100_10_power=my_power(RESULT) %>% mutate(Structure="BD10",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_100_10_fscore=f1score(RESULT) %>% mutate(Structure="BD10",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_100_10_AUC=AUC(RESULT) %>% mutate(Structure="BD10",P=100,N=80)

### BD 80,200
load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n80_p200_e4_min_beta0.3_max_beta0.9.RData")

BD80_200_4_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD4",P=200,N=80)

BD80_200_4_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD4",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_200_4_power=my_power(RESULT) %>% mutate(Structure="BD4",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_200_4_fscore=f1score(RESULT) %>% mutate(Structure="BD4",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_200_4_AUC=AUC(RESULT) %>% mutate(Structure="BD4",P=200,N=80)

load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n80_p200_e8_min_beta0.3_max_beta0.9.RData")

BD80_200_8_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD8",P=200,N=80)

BD80_200_8_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD8",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_200_8_power=my_power(RESULT) %>% mutate(Structure="BD8",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_200_8_fscore=f1score(RESULT) %>% mutate(Structure="BD8",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_200_8_AUC=AUC(RESULT) %>% mutate(Structure="BD8",P=200,N=80)

load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/BlockDiag_simu_n80_p200_e10_min_beta0.3_max_beta0.9.RData")

BD80_200_10_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="BD10",P=200,N=80)

BD80_200_10_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="BD10",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_200_10_power=my_power(RESULT) %>% mutate(Structure="BD10",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_200_10_fscore=f1score(RESULT) %>% mutate(Structure="BD10",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))

BD80_200_10_AUC=AUC(RESULT) %>% mutate(Structure="BD10",P=200,N=80)
}

# ROC----
rbind(BD60_100_4_ROC,BD60_100_8_ROC,BD60_100_10_ROC,
      BD60_200_4_ROC,BD60_200_8_ROC,BD60_200_10_ROC,
      BD80_100_4_ROC,BD80_100_8_ROC,BD80_100_10_ROC,
      BD80_200_4_ROC,BD80_200_8_ROC,BD80_200_10_ROC) %>% 
  mutate(Structure=factor(Structure, levels=c("BD4","BD8","BD10"))) %>% 
  filter(methods %in% c("cPCG_df","cPCG_theo","MLE"),Structure %in% c("BD8","BD10")) %>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv","cPCG.theo","ENF/MLE"))) %>%
  ggplot(aes(x=FPR, y=TPR,col=methods))+
  geom_line()+
  geom_abline(slope = 1, intercept = 0, size=0.3)+
  ggtitle("Block Diagonal")+ 
  facet_grid(Structure~P*N)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/BD_ROC.jpg",width = 6, height = 5, dpi = 300)

# fdr----
rbind(BD60_100_4_fdr,BD60_100_8_fdr,BD60_100_10_fdr,
      BD60_200_4_fdr,BD60_200_8_fdr,BD60_200_10_fdr,
      BD80_100_4_fdr,BD80_100_8_fdr,BD80_100_10_fdr,
      BD80_200_4_fdr,BD80_200_8_fdr,BD80_200_10_fdr) %>% 
  mutate(Structure=factor(Structure, levels=c("BD4","BD8","BD10"))) %>% 
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"),Structure %in% c("BD8","BD10"))%>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv", "cPCG.theo", "ENF.q", "MLE.q"))) %>%
  ggplot(aes(x=FDR, y=fdr, col=methods))+
  geom_line()+   
  geom_abline(slope = 1, intercept = 0, size=0.3)+
  ggtitle("Block Diagonal")+ 
  facet_grid(Structure~P*N)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/BD_fdr.jpg",width = 6, height = 5, dpi = 300)

# fscore----
rbind(BD60_100_4_fscore,BD60_100_8_fscore,BD60_100_10_fscore,
      BD60_200_4_fscore,BD60_200_8_fscore,BD60_200_10_fscore,
      BD80_100_4_fscore,BD80_100_8_fscore,BD80_100_10_fscore,
      BD80_200_4_fscore,BD80_200_8_fscore,BD80_200_10_fscore) %>% 
  mutate(Structure=factor(Structure, levels=c("BD4","BD8","BD10"))) %>% 
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"),Structure %in% c("BD8","BD10"))%>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv", "cPCG.theo", "ENF.q", "MLE.q"))) %>%
  ggplot(aes(y=F1score, x=methods))+
  geom_boxplot() + 
  ggtitle("Block Diagonal")+ 
  facet_grid(Structure~P*N, scales="free_y")+
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  #theme( axis.text.x = element_blank())+
  theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/BD_fscore.jpg",width = 6, height = 6, dpi = 150)

# power----
rbind(BD60_100_4_power,BD60_100_8_power,BD60_100_10_power,
      BD60_200_4_power,BD60_200_8_power,BD60_200_10_power,
      BD80_100_4_power,BD80_100_8_power,BD80_100_10_power,
      BD80_200_4_power,BD80_200_8_power,BD80_200_10_power) %>% 
  mutate(Structure=factor(Structure, levels=c("BD4","BD8","BD10"))) %>% 
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"),Structure %in% c("BD8","BD10")) %>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv", "cPCG.theo", "ENF.q", "MLE.q"))) %>%
  ggplot(aes(y=power, x=methods))+
  geom_boxplot() + 
  ggtitle("Block Diagonal")+ 
  facet_grid(Structure~P*N, scales="free_y")+
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  #theme( axis.text.x = element_blank())+
  theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/BD_power.jpg",width = 6, height = 6, dpi = 300)

# AUC----
rbind(BD60_100_4_AUC,BD60_100_8_AUC,BD60_100_10_AUC,
      BD60_200_4_AUC,BD60_200_8_AUC,BD60_200_10_AUC,
      BD80_100_4_AUC,BD80_100_8_AUC,BD80_100_10_AUC,
      BD80_200_4_AUC,BD80_200_8_AUC,BD80_200_10_AUC) %>% 
  mutate(Structure=factor(Structure, levels=c("BD4","BD8","BD10")), methods=as.character(methods)) %>% 
  filter(methods %in% c("cPCG.df","cPCG.theo","MLE"),Structure %in% c("BD8","BD10")) %>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv","cPCG.theo","ENF/MLE"))) %>%
  ggplot(aes(y=AUC, x=methods))+
  geom_boxplot() + 
  ggtitle("Block Diagonal")+ 
  facet_grid(Structure~P*N, scales="free_y")+
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  #theme( axis.text.x = element_blank())+
  theme(legend.position="bottom", legend.box = "horizontal")
#theme(legend.position = "none", axis.text.x = element_blank(),strip.text.x = element_text(angle = 90, vjust = 0.5))

ggsave("figs/BD_AUC.jpg",width = 6, height = 6, dpi = 300)

# AUC----
rbind(BD60_100_4_AUC,BD60_100_8_AUC,BD60_100_10_AUC,
      BD60_200_4_AUC,BD60_200_8_AUC,BD60_200_10_AUC,
      BD80_100_4_AUC,BD80_100_8_AUC,BD80_100_10_AUC,
      BD80_200_4_AUC,BD80_200_8_AUC,BD80_200_10_AUC) %>% 
  mutate(Structure=factor(Structure, levels=c("BD4","BD8","BD10")), methods=as.character(methods)) %>% 
  filter(methods %in% c("cPCG.df","cPCG.theo","MLE"),Structure %in% c("BD8","BD10")) %>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv","cPCG.theo","ENF/MLE"))) %>%
  ggplot(aes(y=AUC,col=methods)) + 
  geom_boxplot() + 
  ggtitle("Block Diagonal")+ 
  coord_flip()+
  facet_grid(Structure*P*N~.)+
  theme( axis.text.y = element_blank(),strip.text.y = element_text(angle = 90, vjust = 0.5))

ggsave("figs/BD_AUC_new.jpg",width = 6, height = 14, dpi = 300)


# Random ----

{
  ### Ram 60,100 
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n60_p100_eta0.01_min_beta0.3.RData")
  
  Ram60_100_0.01_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.01",P=100,N=60)
  
  Ram60_100_0.01_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.01",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_100_0.01_power=my_power(RESULT) %>% mutate(Structure="Ram0.01",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_100_0.01_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.01",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_100_0.01_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.01",P=100,N=60)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n60_p100_eta0.02_min_beta0.3.RData")
  
  Ram60_100_0.02_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.02",P=100,N=60)
  
  Ram60_100_0.02_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.02",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_100_0.02_power=my_power(RESULT) %>% mutate(Structure="Ram0.02",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_100_0.02_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.02",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_100_0.02_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.02",P=100,N=60)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n60_p100_eta0.03_min_beta0.3.RData")
  
  Ram60_100_0.03_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.03",P=100,N=60)
  
  Ram60_100_0.03_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.03",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_100_0.03_power=my_power(RESULT) %>% mutate(Structure="Ram0.03",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_100_0.03_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.03",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_100_0.03_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.03",P=100,N=60)
  
  ### Ram 60,200
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n60_p200_eta0.01_min_beta0.3.RData")
  
  Ram60_200_0.01_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.01",P=200,N=60)
  
  Ram60_200_0.01_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.01",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_200_0.01_power=my_power(RESULT) %>% mutate(Structure="Ram0.01",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_200_0.01_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.01",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_200_0.01_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.01",P=200,N=60)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n60_p200_eta0.02_min_beta0.3.RData")
  
  Ram60_200_0.02_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.02",P=200,N=60)
  
  Ram60_200_0.02_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.02",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_200_0.02_power=my_power(RESULT) %>% mutate(Structure="Ram0.02",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_200_0.02_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.02",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_200_0.02_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.02",P=200,N=60)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n60_p200_eta0.03_min_beta0.3.RData")
  
  Ram60_200_0.03_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.03",P=200,N=60)
  
  Ram60_200_0.03_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.03",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_200_0.03_power=my_power(RESULT) %>% mutate(Structure="Ram0.03",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_200_0.03_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.03",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram60_200_0.03_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.03",P=200,N=60)
  
  ### Ram 80,100
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n80_p100_eta0.01_min_beta0.3.RData")
  
  Ram80_100_0.01_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.01",P=100,N=80)
  
  Ram80_100_0.01_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.01",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_100_0.01_power=my_power(RESULT) %>% mutate(Structure="Ram0.01",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_100_0.01_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.01",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_100_0.01_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.01",P=100,N=80)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n80_p100_eta0.02_min_beta0.3.RData")
  
  Ram80_100_0.02_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.02",P=100,N=80)
  
  Ram80_100_0.02_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.02",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_100_0.02_power=my_power(RESULT) %>% mutate(Structure="Ram0.02",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_100_0.02_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.02",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_100_0.02_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.02",P=100,N=80)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n80_p100_eta0.03_min_beta0.3.RData")
  
  Ram80_100_0.03_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.03",P=100,N=80)
  
  Ram80_100_0.03_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.03",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_100_0.03_power=my_power(RESULT) %>% mutate(Structure="Ram0.03",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_100_0.03_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.03",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_100_0.03_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.03",P=100,N=80)
  
  ### Ram 80,200
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n80_p200_eta0.01_min_beta0.3.RData")
  
  Ram80_200_0.01_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.01",P=200,N=80)
  
  Ram80_200_0.01_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.01",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_200_0.01_power=my_power(RESULT) %>% mutate(Structure="Ram0.01",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_200_0.01_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.01",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_200_0.01_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.01",P=200,N=80)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n80_p200_eta0.02_min_beta0.3.RData")
  
  Ram80_200_0.02_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.02",P=200,N=80)
  
  Ram80_200_0.02_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.02",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_200_0.02_power=my_power(RESULT) %>% mutate(Structure="Ram0.02",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_200_0.02_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.02",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_200_0.02_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.02",P=200,N=80)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/Random_simu_n80_p200_eta0.03_min_beta0.3.RData")
  
  Ram80_200_0.03_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="Ram0.03",P=200,N=80)
  
  Ram80_200_0.03_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="Ram0.03",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_200_0.03_power=my_power(RESULT) %>% mutate(Structure="Ram0.03",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_200_0.03_fscore=f1score(RESULT) %>% mutate(Structure="Ram0.03",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  Ram80_200_0.03_AUC=AUC(RESULT) %>% mutate(Structure="Ram0.03",P=200,N=80)
}

# ROC----
rbind(Ram60_100_0.01_ROC,Ram60_100_0.02_ROC,Ram60_100_0.03_ROC,
      Ram60_200_0.01_ROC,Ram60_200_0.02_ROC,Ram60_200_0.03_ROC,
      Ram80_100_0.01_ROC,Ram80_100_0.02_ROC,Ram80_100_0.03_ROC,
      Ram80_200_0.01_ROC,Ram80_200_0.02_ROC,Ram80_200_0.03_ROC) %>% 
  mutate(Structure=factor(Structure, levels=c("Ram0.01","Ram0.02","Ram0.03"))) %>% 
  filter(methods %in% c("cPCG_df","cPCG_theo","MLE"),Structure %in% c("Ram0.01","Ram0.02")) %>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv","cPCG.theo","ENF/MLE"))) %>%
  ggplot(aes(x=FPR, y=TPR,col=methods))+
  geom_line()+
  geom_abline(slope = 1, intercept = 0, size=0.3)+
  ggtitle("Random Graph")+ 
  facet_grid(Structure~P*N)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/Ram_ROC.jpg",width = 6, height = 5, dpi = 300)

# fdr----
rbind(Ram60_100_0.01_fdr,Ram60_100_0.02_fdr,Ram60_100_0.03_fdr,
      Ram60_200_0.01_fdr,Ram60_200_0.02_fdr,Ram60_200_0.03_fdr,
      Ram80_100_0.01_fdr,Ram80_100_0.02_fdr,Ram80_100_0.03_fdr,
      Ram80_200_0.01_fdr,Ram80_200_0.02_fdr,Ram80_200_0.03_fdr) %>% 
  mutate(Structure=factor(Structure, levels=c("Ram0.01","Ram0.02","Ram0.03"))) %>% 
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"),Structure %in% c("Ram0.01","Ram0.02"))%>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv", "cPCG.theo", "ENF.q", "MLE.q"))) %>%
  ggplot(aes(x=FDR, y=fdr, col=methods))+
  geom_line()+   
  geom_abline(slope = 1, intercept = 0, size=0.3)+
  ggtitle("Random Graph")+ 
  facet_grid(Structure~P*N)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/Ram_fdr.jpg",width = 6, height = 5, dpi = 300)

# fscore----
rbind(Ram60_100_0.01_fscore,Ram60_100_0.02_fscore,Ram60_100_0.03_fscore,
      Ram60_200_0.01_fscore,Ram60_200_0.02_fscore,Ram60_200_0.03_fscore,
      Ram80_100_0.01_fscore,Ram80_100_0.02_fscore,Ram80_100_0.03_fscore,
      Ram80_200_0.01_fscore,Ram80_200_0.02_fscore,Ram80_200_0.03_fscore) %>% 
  mutate(Structure=factor(Structure, levels=c("Ram0.01","Ram0.02","Ram0.03"))) %>% 
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"),Structure %in% c("Ram0.01","Ram0.02"))%>%
  mutate(methods=factor(methods, labels = c("cPCG.cv", "cPCG.theo", "ENF.q", "MLE.q"))) %>%
  ggplot(aes(y=F1score, x=methods))+
  geom_boxplot() + 
  ggtitle("Random Graph")+ 
  facet_grid(Structure~P*N, scales="free_y")+
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  #theme( axis.text.x = element_blank())+
  theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/Ram_fscore.jpg",width = 6, height = 6, dpi = 150)

# power----
rbind(Ram60_100_0.01_power,Ram60_100_0.02_power,Ram60_100_0.03_power,
      Ram60_200_0.01_power,Ram60_200_0.02_power,Ram60_200_0.03_power,
      Ram80_100_0.01_power,Ram80_100_0.02_power,Ram80_100_0.03_power,
      Ram80_200_0.01_power,Ram80_200_0.02_power,Ram80_200_0.03_power) %>% 
  mutate(Structure=factor(Structure, levels=c("Ram0.01","Ram0.02","Ram.03"))) %>%
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"),Structure %in% c("Ram0.01","Ram0.02")) %>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv", "cPCG.theo", "ENF.q", "MLE.q"))) %>%
  ggplot(aes(y=power, x=methods))+
  geom_boxplot() + 
  ggtitle("Random Graph")+ 
  facet_grid(Structure~P*N, scales="free_y")+
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  #theme( axis.text.x = element_blank())+
  theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/Ram_power.jpg",width = 6, height = 6, dpi = 300)

# AUC----
rbind(Ram60_100_0.01_AUC,Ram60_100_0.02_AUC,Ram60_100_0.03_AUC,
      Ram60_200_0.01_AUC,Ram60_200_0.02_AUC,Ram60_200_0.03_AUC,
      Ram80_100_0.01_AUC,Ram80_100_0.02_AUC,Ram80_100_0.03_AUC,
      Ram80_200_0.01_AUC,Ram80_200_0.02_AUC,Ram80_200_0.03_AUC) %>% 
  mutate(Structure=factor(Structure, levels=c("Ram0.01","Ram0.02","Ram0.03"))) %>% 
  mutate(methods=factor(methods, labels = c("ENF.q" ,"MLE.q", "cPCG.theo", "cPCG.cv", "cPCG.df"))) %>%
  filter(!methods %in% c("cPCG.cv"),Structure %in% c("Ram0.01","Ram0.02")) %>% 
  mutate(methods=factor(methods, levels=c("cPCG.df","cPCG.theo","ENF.q","MLE.q"))) %>%
  ggplot(aes(y=AUC, x=methods))+
  geom_boxplot() + 
  ggtitle("Random Graph")+ 
  facet_grid(Structure~P*N, scales="free_y")+
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  #theme( axis.text.x = element_blank())+
  theme(legend.position="bottom", legend.box = "horizontal")
  #theme(legend.position = "none", axis.text.x = element_blank(),strip.text.x = element_text(angle = 90, vjust = 0.5))

ggsave("figs/Ram_AUC.jpg",width = 6, height = 6, dpi = 300)

# AUC----
rbind(Ram60_100_0.01_AUC,Ram60_100_0.02_AUC,Ram60_100_0.03_AUC,
      Ram60_200_0.01_AUC,Ram60_200_0.02_AUC,Ram60_200_0.03_AUC,
      Ram80_100_0.01_AUC,Ram80_100_0.02_AUC,Ram80_100_0.03_AUC,
      Ram80_200_0.01_AUC,Ram80_200_0.02_AUC,Ram80_200_0.03_AUC) %>% 
  mutate(Structure=factor(Structure, levels=c("Ram0.01","Ram0.02","Ram0.03"))) %>% 
  mutate(methods=factor(methods, labels = c("ENF.q" ,"MLE.q", "cPCG.theo", "cPCG.cv", "cPCG.df"))) %>%
  filter(!methods %in% c("cPCG.cv")) %>% 
  mutate(methods=factor(methods, levels=c("cPCG.df","cPCG.theo","ENF.q","MLE.q"))) %>%
  ggplot(aes(y=AUC,col=methods)) + 
  geom_boxplot() + 
  ggtitle("Random")+ 
  coord_flip()+
  facet_grid(Structure*P*N~.)+
  theme( axis.text.y = element_blank(),strip.text.y = element_text(angle = 90, vjust = 0.5))

ggsave("figs/Ram_AUC_new.jpg",width = 6, height = 14, dpi = 300)



# Scale ----

{
  ### scalef 60,100 
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/ScaleFree_simu_n60_p100_e1_min_beta0.3.RData")
  
  scalef60_100_1_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="scalef1",P=100,N=60)
  
  scalef60_100_1_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="scalef1",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_100_1_power=my_power(RESULT) %>% mutate(Structure="scalef1",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_100_1_fscore=f1score(RESULT) %>% mutate(Structure="scalef1",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_100_1_AUC=AUC(RESULT) %>% mutate(Structure="scalef1",P=100,N=60)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/ScaleFree_simu_n60_p100_e2_min_beta0.3.RData")
  
  scalef60_100_2_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="scalef2",P=100,N=60)
  
  scalef60_100_2_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="scalef2",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_100_2_power=my_power(RESULT) %>% mutate(Structure="scalef2",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_100_2_fscore=f1score(RESULT) %>% mutate(Structure="scalef2",P=100,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_100_2_AUC=AUC(RESULT) %>% mutate(Structure="scalef2",P=100,N=60)
 
  
  ### scalef 60,200
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/ScaleFree_simu_n60_p200_e1_min_beta0.3.RData")
  
  scalef60_200_1_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="scalef1",P=200,N=60)
  
  scalef60_200_1_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="scalef1",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_200_1_power=my_power(RESULT) %>% mutate(Structure="scalef1",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_200_1_fscore=f1score(RESULT) %>% mutate(Structure="scalef1",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_200_1_AUC=AUC(RESULT) %>% mutate(Structure="scalef1",P=200,N=60)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/ScaleFree_simu_n60_p200_e2_min_beta0.3.RData")
  
  scalef60_200_2_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="scalef2",P=200,N=60)
  
  scalef60_200_2_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="scalef2",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_200_2_power=my_power(RESULT) %>% mutate(Structure="scalef2",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_200_2_fscore=f1score(RESULT) %>% mutate(Structure="scalef2",P=200,N=60) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef60_200_2_AUC=AUC(RESULT) %>% mutate(Structure="scalef2",P=200,N=60)

  ### scalef 80,100
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/ScaleFree_simu_n80_p100_e1_min_beta0.3.RData")
  
  scalef80_100_1_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="scalef1",P=100,N=80)
  
  scalef80_100_1_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="scalef1",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_100_1_power=my_power(RESULT) %>% mutate(Structure="scalef1",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_100_1_fscore=f1score(RESULT) %>% mutate(Structure="scalef1",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_100_1_AUC=AUC(RESULT) %>% mutate(Structure="scalef1",P=100,N=80)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/ScaleFree_simu_n80_p100_e2_min_beta0.3.RData")
  
  scalef80_100_2_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="scalef2",P=100,N=80)
  
  scalef80_100_2_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="scalef2",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_100_2_power=my_power(RESULT) %>% mutate(Structure="scalef2",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_100_2_fscore=f1score(RESULT) %>% mutate(Structure="scalef2",P=100,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_100_2_AUC=AUC(RESULT) %>% mutate(Structure="scalef2",P=100,N=80)
  

  
  ### scalef 80,200
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/ScaleFree_simu_n80_p200_e1_min_beta0.3.RData")
  
  scalef80_200_1_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="scalef1",P=200,N=80)
  
  scalef80_200_1_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="scalef1",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_200_1_power=my_power(RESULT) %>% mutate(Structure="scalef1",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_200_1_fscore=f1score(RESULT) %>% mutate(Structure="scalef1",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_200_1_AUC=AUC(RESULT) %>% mutate(Structure="scalef1",P=200,N=80)
  
  load("/Users/Wanghao/GitHub/GenePCG/temp/6_27/ScaleFree_simu_n80_p200_e2_min_beta0.3.RData")
  
  scalef80_200_2_ROC=RESULT$ROC$ROC.out %>% tidyr::gather("methods", "TPR", ENF:cPCG_df) %>% mutate(Structure="scalef2",P=200,N=80)
  
  scalef80_200_2_fdr=RESULT$fdr$fdr.table %>% mutate(Structure="scalef2",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_200_2_power=my_power(RESULT) %>% mutate(Structure="scalef2",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_200_2_fscore=f1score(RESULT) %>% mutate(Structure="scalef2",P=200,N=80) #%>% filter(!methods %in% c("ENF.p","MLE.p"))
  
  scalef80_200_2_AUC=AUC(RESULT) %>% mutate(Structure="scalef2",P=200,N=80)
  
}

# ROC----
rbind(scalef60_100_1_ROC,scalef60_100_2_ROC,
      scalef60_200_1_ROC,scalef60_200_2_ROC,
      scalef80_100_1_ROC,scalef80_100_2_ROC,
      scalef80_200_1_ROC,scalef80_200_2_ROC) %>% 
  mutate(Structure=factor(Structure, levels=c("scalef1","scalef2"))) %>% 
  filter(methods %in% c("cPCG_df","cPCG_theo","MLE")) %>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv","cPCG.theo","ENF.MLE"))) %>%
  ggplot(aes(x=FPR, y=TPR,col=methods))+
  geom_line()+
  geom_abline(slope = 1, intercept = 0, size=0.3)+
  ggtitle("Scale Free")+ 
  facet_grid(Structure~P*N)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/scalef_ROC.jpg",width = 6, height = 5, dpi = 300)

# fdr----
rbind(scalef60_100_1_fdr,scalef60_100_2_fdr,
      scalef60_200_1_fdr,scalef60_200_2_fdr,
      scalef80_100_1_fdr,scalef80_100_2_fdr,
      scalef80_200_1_fdr,scalef80_200_2_fdr) %>% 
  mutate(Structure=factor(Structure, levels=c("scalef1","scalef2"))) %>% 
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"))%>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv", "cPCG.theo", "ENF.q", "MLE.q"))) %>%
  ggplot(aes(x=FDR, y=fdr, col=methods))+
  geom_line()+   
  geom_abline(slope = 1, intercept = 0, size=0.3)+
  ggtitle("Scale Free")+ 
  facet_grid(Structure~P*N)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/scalef_fdr.jpg",width = 6, height = 5, dpi = 300)

# fscore----
rbind(scalef60_100_1_fscore,scalef60_100_2_fscore,
      scalef60_200_1_fscore,scalef60_200_2_fscore,
      scalef80_100_1_fscore,scalef80_100_2_fscore,
      scalef80_200_1_fscore,scalef80_200_2_fscore) %>% 
  mutate(Structure=factor(Structure, levels=c("scalef1","scalef2"))) %>% 
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"))%>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv", "cPCG.theo", "ENF.q", "MLE.q"))) %>%
  ggplot(aes(y=F1score, x=methods))+
  geom_boxplot() + 
  ggtitle("Scale Free")+ 
  facet_grid(Structure~P*N, scales="free_y")+
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  #theme( axis.text.x = element_blank())+
  theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/scalef_fscore.jpg",width = 6, height = 6, dpi = 150)

# power----
rbind(scalef60_100_1_power,scalef60_100_2_power,
      scalef60_200_1_power,scalef60_200_2_power,
      scalef80_100_1_power,scalef80_100_2_power,
      scalef80_200_1_power,scalef80_200_2_power) %>% 
  mutate(Structure=factor(Structure, levels=c("scalef1","scalef2"))) %>%
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv")) %>% 
  mutate(methods=factor(methods, labels = c("cPCG.cv", "cPCG.theo", "ENF.q", "MLE.q"))) %>%
  ggplot(aes(y=power, x=methods))+
  geom_boxplot() + 
  ggtitle("Scale Free")+ 
  facet_grid(Structure~P*N, scales="free_y")+
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  #theme( axis.text.x = element_blank())+
  theme(legend.position="bottom", legend.box = "horizontal")

ggsave("figs/scalef_power.jpg",width = 6, height = 6, dpi = 300)

# AUC----
rbind(scalef60_100_1_AUC,scalef60_100_2_AUC,
      scalef60_200_1_AUC,scalef60_200_2_AUC,
      scalef80_100_1_AUC,scalef80_100_2_AUC,
      scalef80_200_1_AUC,scalef80_200_2_AUC) %>% 
  mutate(Structure=factor(Structure, levels=c("scalef1","scalef2"))) %>% 
  mutate(methods=factor(methods, labels = c("ENF.q" ,"MLE.q", "cPCG.theo", "cPCG.cv", "cPCG.df"))) %>%
  filter(!methods %in% c("cPCG.cv")) %>% 
  mutate(methods=factor(methods, levels=c("cPCG.df","cPCG.theo","ENF.q","MLE.q"))) %>%
  ggplot(aes(y=AUC,x=methods)) + 
  geom_boxplot() + 
  ggtitle("Scale Free")+ 
  facet_grid(Structure~P*N, scales="free_y")+
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  #theme( axis.text.x = element_blank())+
  theme(legend.position="bottom", legend.box = "horizontal")
#theme(legend.position = "none", axis.text.x = element_blank(),strip.text.x = element_text(angle = 90, vjust = 0.5))

ggsave("figs/scalef_AUC.jpg",width = 6, height = 6, dpi = 300)

# AUC----
rbind(scalef60_100_1_AUC,scalef60_100_2_AUC,
      scalef60_200_1_AUC,scalef60_200_2_AUC,
      scalef80_100_1_AUC,scalef80_100_2_AUC,
      scalef80_200_1_AUC,scalef80_200_2_AUC) %>% 
  mutate(Structure=factor(Structure, levels=c("scalef1","scalef2"))) %>% 
  mutate(methods=factor(methods, labels = c("ENF.q" ,"MLE.q", "cPCG.theo", "cPCG.cv", "cPCG.df"))) %>%
  filter(!methods %in% c("cPCG.cv")) %>% 
  mutate(methods=factor(methods, levels=c("cPCG.df","cPCG.theo","ENF.q","MLE.q"))) %>%
  ggplot(aes(y=AUC,col=methods)) + 
  geom_boxplot() + 
  ggtitle("Scale Free")+ 
  coord_flip()+
  facet_grid(Structure*P*N~.)+
  theme(axis.text.y = element_blank(),strip.text.y = element_text(angle = 90, vjust = 0.5))

ggsave("figs/scalef_AUC_new.jpg",width = 6, height = 12, dpi = 300)

# power table----
## BD
t1=rbind(BD60_100_4_power,BD60_100_8_power,BD60_100_10_power,
      BD60_200_4_power,BD60_200_8_power,BD60_200_10_power,
      BD80_100_4_power,BD80_100_8_power,BD80_100_10_power,
      BD80_200_4_power,BD80_200_8_power,BD80_200_10_power) %>% 
  mutate(Structure=factor(Structure, levels=c("BD4","BD8","BD10"))) %>% 
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"),Structure %in% c("BD4","BD8","BD10")) %>%
  group_by(methods, Structure, P, N) %>%
  summarize(Mean=mean(power), SD=sd(power))

## Random
t2=rbind(Ram60_100_0.01_power,Ram60_100_0.02_power,Ram60_100_0.03_power,
      Ram60_200_0.01_power,Ram60_200_0.02_power,Ram60_200_0.03_power,
      Ram80_100_0.01_power,Ram80_100_0.02_power,Ram80_100_0.03_power,
      Ram80_200_0.01_power,Ram80_200_0.02_power,Ram80_200_0.03_power) %>% 
  mutate(Structure=factor(Structure, levels=c("Ram0.01","Ram0.02","Ram0.03"))) %>%
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv"),Structure %in% c("Ram0.01","Ram0.02","Ram0.03")) %>%
  group_by(methods, Structure, P, N) %>%
  summarize(Mean=mean(power), SD=sd(power))

## Scale free
t3=rbind(scalef60_100_1_power,scalef60_100_2_power,
      scalef60_200_1_power,scalef60_200_2_power,
      scalef80_100_1_power,scalef80_100_2_power,
      scalef80_200_1_power,scalef80_200_2_power) %>% 
  mutate(Structure=factor(Structure, levels=c("scalef1","scalef2"))) %>%
  filter(!methods %in% c("ENF.p","MLE.p","cPCG.cv")) %>% 
  group_by(methods, Structure, P, N) %>%
  summarize(Mean=mean(power), SD=sd(power))
