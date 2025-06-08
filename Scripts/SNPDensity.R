setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/ACR_Divergence/")
library(tidyverse)
library(ggplot2)

Summit_CF<-read.delim("Data/PopGen/CF_peaks_WinCount_mod.bed",header = F)
colnames(Summit_CF)<-c("Chr","Start","End","PeakID","BinIndex","SNPCount")
quantile(Summit_CF$SNPCount,probs = c(0,0.9,0.95,0.98,0.99,1))

as_tibble(Summit_CF) %>% 
    select(BinIndex,SNPCount) %>% 
    group_by(BinIndex) %>% 
    summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
    select(BinIndex, mean) %>% 
    mutate(
      SNPCountNorm= ((mean -min(mean))/(
        max(mean)-min(mean)))) %>% 
    ggplot(aes(x=BinIndex))+
    geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
    geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
    scale_x_continuous(breaks = c(1,10.5,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(1,20))+
    scale_y_continuous("Normalized SNP Density",breaks = c(0,0.5,1))+
    geom_vline(xintercept = 10.5,linetype="dashed",color="red",size=0.6)+
    labs(title="Coastal",x="")+ #"#0072B2","#D55E00"
    theme_classic(base_size = 16)+
    theme(
      axis.ticks.x = element_blank(),
      axis.text = element_text(size=16,colour = "black"),
      axis.line.y.right = element_line(color = "gray60"),
      axis.title.y = element_text(size=16,colour = "black"),
      plot.title = element_text(hjust = 0.5, colour = "#0072B2",face = "bold"))
  
 

Summit_CH<-read.delim("Data/PopGen/CH_peaks_WinCount_mod.bed",header = F)
colnames(Summit_CH)<-c("Chr","Start","End","PeakID","BinIndex","SNPCount")
quantile(Summit_CH$SNPCount,probs = c(0,0.9,0.95,0.98,0.99,1))

as_tibble(Summit_CH) %>% 
  select(BinIndex,SNPCount) %>% 
  group_by(BinIndex) %>% 
  summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
  select(BinIndex, mean) %>% 
  mutate(
    SNPCountNorm= ((mean -min(mean))/(
      max(mean)-min(mean)))) %>% 
  ggplot(aes(x=BinIndex))+
  geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
  geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
  scale_x_continuous(breaks = c(1,10.5,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(1,20))+
  scale_y_continuous("Normalized SNP Density",breaks = c(0,0.5,1))+
  geom_vline(xintercept = 10.5,linetype="dashed",color="red",size=0.6)+
  labs(title="Inland",x="")+ #"#0072B2","#D55E00"
  theme_classic(base_size = 16)+
  theme(
    axis.ticks.x = element_blank(),
    axis.text = element_text(size=16,colour = "black"),
    axis.line.y.right = element_line(color = "gray60"),
    axis.title.y = element_text(size=16,colour = "black"),
    plot.title = element_text(hjust = 0.5, colour = "#D55E00",face = "bold"))

 

