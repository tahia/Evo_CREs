setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/ACR_Divergence/")

library(ggplot2)
library(xlsx)
library(tidyverse)
options(java.parameters = "-Xmx6g" )

########## SNP Diversity
########## SNP And TAG Count
## Coastal
Summit_CF<-read.delim("Data/FIL/Catalog/CF_peaks_TAGCount.bed",header = F)
colnames(Summit_CF)<-c("Chr","Start","End","PeakID","BinIndex","TAGCount")
SNP_CF<-read.delim("Data/FIL/Catalog/CF_peaks_WinCount.bed",header = F)
Summit_CF$SNPCount<-SNP_CF$V6

#colnames(Summit_CF)<-c("PeakID","BinIndex","TAGCount","SNPCount")
Summit_CF$BinIndex<-gsub("Bin_","",Summit_CF$BinIndex)
Summit_CF$BinIndex<-as.numeric(as.character(Summit_CF$BinIndex))
quantile(Summit_CF$TAGCount,probs = c(0,0.9,0.95,0.98,0.99,1))
quantile(Summit_CF$SNPCount,probs = c(0,0.9,0.95,0.98,0.99,1))

    
(SumCoastal<-as_tibble(Summit_CF) %>% 
    select(BinIndex,TAGCount,SNPCount) %>% 
  group_by(BinIndex) %>% 
  summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
  select(BinIndex,TAGCount_mean, SNPCount_mean) %>% 
  mutate(TAGCountNorm= ((TAGCount_mean -min(TAGCount_mean))/(
    max(TAGCount_mean)-min(TAGCount_mean))),
    SNPCountNorm= ((SNPCount_mean -min(SNPCount_mean))/(
      max(SNPCount_mean)-min(SNPCount_mean)))) %>% 
  ggplot(aes(x=BinIndex))+
  geom_point(aes(y=TAGCountNorm),size=0.75,color="seagreen")+
  geom_line(aes(y=TAGCountNorm),size=0.85,color="seagreen")+
  geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
  geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
  scale_x_continuous(breaks = c(0,10.5,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(0,20))+
  scale_y_continuous("Normalized Tn5 Integration",breaks = c(0,0.5,1),
                     sec.axis = sec_axis(~.*1,name="Normalized SNP Density",
                                         breaks = c(0,0.5,1)))+
  geom_vline(xintercept = 10.5,linetype="dashed",color="red",size=0.6)+
  labs(title="Coastal",x="")+ #"#0072B2","#D55E00"
  theme_classic(base_size = 16)+
  theme(
    axis.ticks.x = element_blank(),
    axis.text = element_text(size=16,colour = "black"),
    axis.line.y.right = element_line(color = "gray60"),
    axis.line.y.left = element_line(color = "seagreen"),
    axis.title.y = element_text(size=16,colour = "black"),
    plot.title = element_text(hjust = 0.5, colour = "#0072B2",face = "bold"))

)  

### Inland  
Summit_CH<-read.delim("Data/Redo_2024_V2/HAL/Catalog/CH_peaks_TAGCount.bed",header = F)
colnames(Summit_CH)<-c("Chr","Start","End","PeakID","BinIndex","TAGCount")
SNP_CH<-read.delim("Data/Redo_2024_V2/HAL/Catalog/CH_peaks_WinCount.bed",header = F)
Summit_CH$SNPCount<-SNP_CH$V6

#colnames(Summit_CF)<-c("PeakID","BinIndex","TAGCount","SNPCount")
Summit_CH$BinIndex<-gsub("Bin_","",Summit_CH$BinIndex)
Summit_CH$BinIndex<-as.numeric(as.character(Summit_CH$BinIndex))
quantile(Summit_CH$TAGCount,probs = c(0,0.9,0.95,0.98,0.99,1))
quantile(Summit_CH$SNPCount,probs = c(0,0.9,0.95,0.98,0.99,1))

(SumInland<-  as_tibble(Summit_CH) %>% 
  select(BinIndex,TAGCount,SNPCount) %>% 
  group_by(BinIndex) %>% 
  summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
  select(BinIndex, TAGCount_mean, SNPCount_mean) %>% 
  mutate(TAGCountNorm= ((TAGCount_mean -min(TAGCount_mean))/(
    max(TAGCount_mean)-min(TAGCount_mean))),
    SNPCountNorm= ((SNPCount_mean -min(SNPCount_mean))/(
      max(SNPCount_mean)-min(SNPCount_mean)))) %>% 
  ggplot(aes(x=BinIndex))+
  geom_point(aes(y=TAGCountNorm),size=0.75,color="seagreen")+
  geom_line(aes(y=TAGCountNorm),size=0.85,color="seagreen")+
  geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
  geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
  scale_x_continuous(breaks = c(0,10.5,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(0,20))+
  scale_y_continuous("Normalized Tn5 Integration",breaks = c(0,0.5,1),
                     sec.axis = sec_axis(~.*1,name="Normalized SNP Density",
                                         breaks = c(0,0.5,1)))+
  geom_vline(xintercept = 10.5,linetype="dashed",color="red",size=0.6)+
  labs(title="Inland",x="")+ #"#0072B2","#D55E00"
  theme_classic(base_size = 16)+
  theme(
    axis.ticks.x = element_blank(),
    axis.text = element_text(size=16,colour = "black"),
    axis.line.y.right = element_line(color = "gray60"),
    axis.line.y.left = element_line(color = "seagreen"),
    axis.title.y = element_text(size=16,colour = "black"),
    plot.title = element_text(hjust = 0.5, colour = "#D55E00",face = "bold"))

)


library(cowplot)
Fig3 <-
  ggdraw() +
  draw_plot(SumInland,x = 0, y = 0,width = 0.5,height = 1) +
  draw_plot(SumCoastal,x = 0.5, y = 0,width = 0.5,height = 1) +
  draw_label(x=0.02,y=0.98,label = "A)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.52,y=0.98,label = "B)", color = "black", size = 16, fontface = "bold")




tiff("Plots/Fig3.tiff",width=12,height=6,units="in",res=300)
Fig3
dev.off()


######## Supplementary Figure
# Do some wranglling to merge annotation and Tn5 integration
Coastal_Narrow<-read.delim("Data/FIL/Catalog/CF_peaks.clean2.catalog.narrowPeak",
                           header = F) %>% 
  as_tibble() %>% 
  `colnames<-`(c("Chr","Start","End","PeakID", "Score","Strand","FE","pval","qval","Summit") ) %>% 
  select(Chr,Start,End,PeakID,Summit) %>% 
  mutate(Peak_ID=paste(Chr,Start,End,Summit,sep = "_")) %>% 
  select(PeakID , Peak_ID)

CF_Annot<-read.delim("Data/FIL/Catalog_Annot/CF_peaks.clean2.catalog.ACR.bed",
                          header = F) %>% 
  as_tibble() %>% 
  `colnames<-`(c("Chr","Start","End","Peak_ID", "FE","qval", "Summit", "GeneID","Category")) %>% 
  select(Peak_ID,GeneID,Category) %>% 
  left_join(.,Coastal_Narrow, 
            by=c("Peak_ID"="Peak_ID"))

Inland_Narrow<-read.delim("Data/Redo_2024_V2/HAL/Catalog/CH_peaks.clean2.catalog.narrowPeak",
                           header = F) %>% 
  as_tibble() %>% 
  `colnames<-`(c("Chr","Start","End","PeakID", "Score","Strand","FE","pval","qval","Summit") ) %>% 
  select(Chr,Start,End,PeakID,Summit) %>% 
  mutate(Peak_ID=paste(Chr,Start,End,Summit,sep = "_")) %>% 
  select(PeakID , Peak_ID)

CH_Annot<-read.delim("Data/Redo_2024_V2/HAL/Catalog_Annot/CH_peaks.clean2.catalog.ACR.bed",
                     header = F) %>% 
  as_tibble() %>% 
  `colnames<-`(c("Chr","Start","End","Peak_ID", "FE","qval", "Summit", "GeneID","Category")) %>% 
  select(Peak_ID,GeneID,Category) %>% 
  left_join(.,Inland_Narrow, 
            by=c("Peak_ID"="Peak_ID"))


(SumgCoastal<-as_tibble(Summit_CF) %>% 
    left_join(.,y=CF_Annot,
              #by=c("PeakID"="Peak_ID")) %>% 
              by=c("PeakID"="PeakID")) %>% 
    select(BinIndex,TAGCount,SNPCount,Category) %>% 
    filter(Category=="gACR") %>% 
    group_by(BinIndex) %>% 
    summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
    select(BinIndex,TAGCount_mean, SNPCount_mean) %>% 
    mutate(TAGCountNorm= ((TAGCount_mean -min(TAGCount_mean))/(
      max(TAGCount_mean)-min(TAGCount_mean))),
      SNPCountNorm= ((SNPCount_mean -min(SNPCount_mean))/(
        max(SNPCount_mean)-min(SNPCount_mean)))) %>% 
    ggplot(aes(x=BinIndex))+
    geom_point(aes(y=TAGCountNorm),size=0.75,color="seagreen")+
    geom_line(aes(y=TAGCountNorm),size=0.85,color="seagreen")+
    geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
    geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
    scale_x_continuous(breaks = c(0,10,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(0,20))+
    scale_y_continuous("Normalized Tn5 Integration",breaks = c(0,0.5,1),
                       sec.axis = sec_axis(~.*1,name="Normalized SNP Density",
                                           breaks = c(0,0.5,1)))+
    geom_vline(xintercept = 10,linetype="dashed",color="red",size=0.6)+
    labs(title="gACR Coastal",x="")+ #"#0072B2","#D55E00"
    theme_classic(base_size = 16)+
    theme(
      axis.ticks.x = element_blank(),
      axis.text = element_text(size=16,colour = "black"),
      axis.line.y.right = element_line(color = "gray60"),
      axis.line.y.left = element_line(color = "seagreen"),
      axis.title.y = element_text(size=16,colour = "black"),
      plot.title = element_text(hjust = 0.5, colour = "#0072B2",face = "bold"))
  
)  

(SumpCoastal<-as_tibble(Summit_CF) %>% 
    left_join(.,y=CF_Annot,
              by=c("PeakID"="PeakID")) %>% 
    select(BinIndex,TAGCount,SNPCount,Category) %>% 
    filter(Category=="pACR") %>% 
    group_by(BinIndex) %>% 
    summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
    select(BinIndex,TAGCount_mean, SNPCount_mean) %>% 
    mutate(TAGCountNorm= ((TAGCount_mean -min(TAGCount_mean))/(
      max(TAGCount_mean)-min(TAGCount_mean))),
      SNPCountNorm= ((SNPCount_mean -min(SNPCount_mean))/(
        max(SNPCount_mean)-min(SNPCount_mean)))) %>% 
    ggplot(aes(x=BinIndex))+
    geom_point(aes(y=TAGCountNorm),size=0.75,color="seagreen")+
    geom_line(aes(y=TAGCountNorm),size=0.85,color="seagreen")+
    geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
    geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
    scale_x_continuous(breaks = c(0,10,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(0,20))+
    scale_y_continuous("Normalized Tn5 Integration",breaks = c(0,0.5,1),
                       sec.axis = sec_axis(~.*1,name="Normalized SNP Density",
                                           breaks = c(0,0.5,1)))+
    geom_vline(xintercept = 10,linetype="dashed",color="red",size=0.6)+
    labs(title="pACR Coastal",x="")+ #"#0072B2","#D55E00"
    theme_classic(base_size = 16)+
    theme(
      axis.ticks.x = element_blank(),
      axis.text = element_text(size=16,colour = "black"),
      axis.line.y.right = element_line(color = "gray60"),
      axis.line.y.left = element_line(color = "seagreen"),
      axis.title.y = element_text(size=16,colour = "black"),
      plot.title = element_text(hjust = 0.5, colour = "#0072B2",face = "bold"))
  
)  

(SumdCoastal<-as_tibble(Summit_CF) %>% 
    left_join(.,y=CF_Annot,
              by=c("PeakID"="PeakID")) %>% 
    select(BinIndex,TAGCount,SNPCount,Category) %>% 
    filter(Category=="dACR") %>% 
    group_by(BinIndex) %>% 
    summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
    select(BinIndex,TAGCount_mean, SNPCount_mean) %>% 
    mutate(TAGCountNorm= ((TAGCount_mean -min(TAGCount_mean))/(
      max(TAGCount_mean)-min(TAGCount_mean))),
      SNPCountNorm= ((SNPCount_mean -min(SNPCount_mean))/(
        max(SNPCount_mean)-min(SNPCount_mean)))) %>% 
    ggplot(aes(x=BinIndex))+
    geom_point(aes(y=TAGCountNorm),size=0.75,color="seagreen")+
    geom_line(aes(y=TAGCountNorm),size=0.85,color="seagreen")+
    geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
    geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
    scale_x_continuous(breaks = c(0,10,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(0,20))+
    scale_y_continuous("Normalized Tn5 Integration",breaks = c(0,0.5,1),
                       sec.axis = sec_axis(~.*1,name="Normalized SNP Density",
                                           breaks = c(0,0.5,1)))+
    geom_vline(xintercept = 10,linetype="dashed",color="red",size=0.6)+
    labs(title="dACR Coastal",x="")+ #"#0072B2","#D55E00"
    theme_classic(base_size = 16)+
    theme(
      axis.ticks.x = element_blank(),
      axis.text = element_text(size=16,colour = "black"),
      axis.line.y.right = element_line(color = "gray60"),
      axis.line.y.left = element_line(color = "seagreen"),
      axis.title.y = element_text(size=16,colour = "black"),
      plot.title = element_text(hjust = 0.5, colour = "#0072B2",face = "bold"))
  
)  

(SumgInland<-as_tibble(Summit_CH) %>% 
    left_join(.,y=CH_Annot,
              by=c("PeakID"="PeakID")) %>% 
    select(BinIndex,TAGCount,SNPCount,Category) %>% 
    filter(Category=="gACR") %>%  
    group_by(BinIndex) %>% 
    summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
    select(BinIndex, TAGCount_mean, SNPCount_mean) %>% 
    mutate(TAGCountNorm= ((TAGCount_mean -min(TAGCount_mean))/(
      max(TAGCount_mean)-min(TAGCount_mean))),
      SNPCountNorm= ((SNPCount_mean -min(SNPCount_mean))/(
        max(SNPCount_mean)-min(SNPCount_mean)))) %>% 
    ggplot(aes(x=BinIndex))+
    geom_point(aes(y=TAGCountNorm),size=0.75,color="seagreen")+
    geom_line(aes(y=TAGCountNorm),size=0.85,color="seagreen")+
    geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
    geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
    scale_x_continuous(breaks = c(0,10,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(0,20))+
    scale_y_continuous("Normalized Tn5 Integration",breaks = c(0,0.5,1),
                       sec.axis = sec_axis(~.*1,name="Normalized SNP Density",
                                           breaks = c(0,0.5,1)))+
    geom_vline(xintercept = 10,linetype="dashed",color="red",size=0.6)+
    labs(title="gACR Inland",x="")+ #"#0072B2","#D55E00"
    theme_classic(base_size = 16)+
    theme(
      axis.ticks.x = element_blank(),
      axis.text = element_text(size=16,colour = "black"),
      axis.line.y.right = element_line(color = "gray60"),
      axis.line.y.left = element_line(color = "seagreen"),
      axis.title.y = element_text(size=16,colour = "black"),
      plot.title = element_text(hjust = 0.5, colour = "#D55E00",face = "bold"))
  
)

(SumpInland<-as_tibble(Summit_CH) %>% 
    left_join(.,y=CH_Annot,
              by=c("PeakID"="PeakID")) %>% 
    select(BinIndex,TAGCount,SNPCount,Category) %>% 
    filter(Category=="pACR") %>%  
    group_by(BinIndex) %>% 
    summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
    select(BinIndex, TAGCount_mean, SNPCount_mean) %>% 
    mutate(TAGCountNorm= ((TAGCount_mean -min(TAGCount_mean))/(
      max(TAGCount_mean)-min(TAGCount_mean))),
      SNPCountNorm= ((SNPCount_mean -min(SNPCount_mean))/(
        max(SNPCount_mean)-min(SNPCount_mean)))) %>% 
    ggplot(aes(x=BinIndex))+
    geom_point(aes(y=TAGCountNorm),size=0.75,color="seagreen")+
    geom_line(aes(y=TAGCountNorm),size=0.85,color="seagreen")+
    geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
    geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
    scale_x_continuous(breaks = c(0,10,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(0,20))+
    scale_y_continuous("Normalized Tn5 Integration",breaks = c(0,0.5,1),
                       sec.axis = sec_axis(~.*1,name="Normalized SNP Density",
                                           breaks = c(0,0.5,1)))+
    geom_vline(xintercept = 10,linetype="dashed",color="red",size=0.6)+
    labs(title="pACR Inland",x="")+ #"#0072B2","#D55E00"
    theme_classic(base_size = 16)+
    theme(
      axis.ticks.x = element_blank(),
      axis.text = element_text(size=16,colour = "black"),
      axis.line.y.right = element_line(color = "gray60"),
      axis.line.y.left = element_line(color = "seagreen"),
      axis.title.y = element_text(size=16,colour = "black"),
      plot.title = element_text(hjust = 0.5, colour = "#D55E00",face = "bold"))
  
)

(SumdInland<-as_tibble(Summit_CH) %>% 
    left_join(.,y=CH_Annot,
              by=c("PeakID"="PeakID")) %>% 
    select(BinIndex,TAGCount,SNPCount,Category) %>% 
    filter(Category=="dACR") %>%  
    group_by(BinIndex) %>% 
    summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
    select(BinIndex, TAGCount_mean, SNPCount_mean) %>% 
    mutate(TAGCountNorm= ((TAGCount_mean -min(TAGCount_mean))/(
      max(TAGCount_mean)-min(TAGCount_mean))),
      SNPCountNorm= ((SNPCount_mean -min(SNPCount_mean))/(
        max(SNPCount_mean)-min(SNPCount_mean)))) %>% 
    ggplot(aes(x=BinIndex))+
    geom_point(aes(y=TAGCountNorm),size=0.75,color="seagreen")+
    geom_line(aes(y=TAGCountNorm),size=0.85,color="seagreen")+
    geom_point(aes(y=SNPCountNorm),size=0.75,color="gray60")+
    geom_line(aes(y=SNPCountNorm),size=0.85,color="gray60")+
    scale_x_continuous(breaks = c(0,10,20),labels = c("-2Kb","Summit","+2Kb"),limits = c(0,20))+
    scale_y_continuous("Normalized Tn5 Integration",breaks = c(0,0.5,1),
                       sec.axis = sec_axis(~.*1,name="Normalized SNP Density",
                                           breaks = c(0,0.5,1)))+
    geom_vline(xintercept = 10,linetype="dashed",color="red",size=0.6)+
    labs(title="dACR Inland",x="")+ #"#0072B2","#D55E00"
    theme_classic(base_size = 16)+
    theme(
      axis.ticks.x = element_blank(),
      axis.text = element_text(size=16,colour = "black"),
      axis.line.y.right = element_line(color = "gray60"),
      axis.line.y.left = element_line(color = "seagreen"),
      axis.title.y = element_text(size=16,colour = "black"),
      plot.title = element_text(hjust = 0.5, colour = "#D55E00",face = "bold"))
  
)


library(cowplot)
SF1 <-
  ggdraw() +
  draw_plot(SumgInland,x = 0, y = 0.5,width = 0.33,height = 0.5) +
  draw_plot(SumpInland,x = 0.33, y = 0.5,width = 0.33,height = 0.5) +
  draw_plot(SumdInland,x = 0.66, y = 0.5,width = 0.33,height = 0.5) +
  draw_plot(SumgCoastal,x = 0, y = 0,width = 0.33,height = 0.5) +
  draw_plot(SumpCoastal,x = 0.33, y = 0,width = 0.33,height = 0.5) +
  draw_plot(SumdCoastal,x = 0.66, y = 0,width = 0.33,height = 0.5) +
  draw_label(x=0.02,y=0.98,label = "A)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.35,y=0.98,label = "B)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.67,y=0.98,label = "C)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.02,y=0.48,label = "D)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.35,y=0.48,label = "E)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.67,y=0.48,label = "F)", color = "black", size = 16, fontface = "bold")



tiff("Plots/Supplementary_S1.tiff",width=10,height=9,units="in",res=300)
SF1
dev.off()

