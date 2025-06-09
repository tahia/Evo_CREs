setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/ACR_Divergence/")
library(ggplot2)
library(tidyverse)


##### Read per bin count, normalize and write the data frame
#Run only once
#datMatCov<-read.delim("Data/FIL/Catalog/CF_Summit.tab",header = F,sep = "\t",skip = 3,na.strings = c(NA,"nan"))
datMatCov<-read.delim("Data/HAL/Catalog/CH_Summit.tab",header = F,sep = "\t",skip = 3,na.strings = c(NA,"nan"))

NormdatMatCov<-matrix(data = NA,nrow = nrow(datMatCov),ncol = ncol(datMatCov))

# offset<-as.numeric(quantile(melt(datMatCov)[,2],probs = 0.0))
# scaleF<-as.numeric(quantile(melt(datMatCov)[,2],probs = 0.98))

for (i in 1:nrow(datMatCov)) {
  vect<-as.vector(unlist(datMatCov[i,]))
  offset<-as.numeric(quantile(vect[!is.na(vect)],probs = 0.02))
  scaleF<-as.numeric(quantile(vect[!is.na(vect)],probs = 0.98))
  vect<- (vect - offset)/scaleF
  NormdatMatCov[i,]<-vect
}

NormdatMatCov<-as.data.frame(NormdatMatCov)
plot(colMeans(NormdatMatCov,na.rm = T))

CF<-as.data.frame(cbind(Count=colMeans(NormdatMatCov,na.rm = T),
                          Pos=seq(-4999,4999,10),
                         Ecotype=rep("Coastal",1000)
))

write.csv(CF,file="Data/FIL/Catalog/CF_Summit_Normalized.csv",row.names = F)

CH<-as.data.frame(cbind(Count=colMeans(NormdatMatCov,na.rm = T),
                        Pos=seq(-4999,4999,10),
                        Ecotype=rep("Inland",1000)
))

write.csv(CH,file="Data/HAL/Catalog/CH_Summit_Normalized.csv",row.names = F)

#################
CF<-read.csv("Data/FIL/Catalog/CF_Summit_Normalized.csv")
CH<-read.csv("Data/HAL/Catalog/CH_Summit_Normalized.csv")

datMeta<-rbind(CF,CH) %>% 
  mutate(Ecotype=as.factor(Ecotype))

(MetaPlot<-ggplot(data=datMeta,aes(y=Count,x=Pos,color=Ecotype))+
  geom_line(size=1.1,linetype=6,alpha=0.95,show.legend = F)+
  geom_vline(xintercept = 0,linetype=2,size=0.3)+
  scale_x_continuous(limits = c(-2000,2000),breaks = c(-2000,0,2000),
                       labels = c("-2K","Summit","2K"))+
  scale_y_continuous(limits = c(0.1,0.9),
                     breaks = c(0.25,0.5,0.75),
                     labels = c(0.25,0.5,0.75))+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  labs(title = "", x="",y="Normalized Count")+
  #facet_wrap(~Ecotype)+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1,color="black",face="bold") )
)

###### ACR length distribution
FIL_ACR<-read.delim("Data/FIL/Catalog/CF_peaks.clean2.catalog.narrowPeak",header = F) %>% 
  mutate(Start=V2,End=V3) %>% 
  mutate(Length=abs(End-Start+1)) %>% 
  mutate(Ecotype="Coastal") %>% 
  select(Length,Ecotype)
  
HAL_ACR<-read.delim("Data/HAL/Catalog/CH_peaks.clean2.catalog.narrowPeak",header = F) %>% 
  mutate(Start=V2,End=V3) %>% 
  mutate(Length=abs(End-Start+1)) %>% 
  mutate(Ecotype="Inlnad") %>% 
  select(Length,Ecotype)

ACR_Length<-rbind(FIL_ACR,HAL_ACR) %>% 
  mutate(Length = ifelse(Length >2000, 2000, Length))

(LenDenPlot<-ggplot(data=ACR_Length,aes(x=Length,color=Ecotype))+
  geom_density(size=1.1,linetype=6,alpha=0.7,show.legend = F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  labs(title = "", x="Length (bp)",y="Density")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1,color="black",face="bold") )
)

Fetdat<-as.data.frame(cbind(Count=c(12618, 6989, 8692, 6311, 4087, 8348),
                            Percent = c("44%","25%","31%","34%","22%","44%"),
                            Pos = c(6000,16000,23000,3000,8000,14000),
                            Category=c("gACR","pACR","dACR","gACR","pACR","dACR"),
                            Ecotype=c(rep("Coastal",3),rep("Inland",3)))) %>% 
  mutate(Category=factor(Category,levels=c("dACR","pACR","gACR"))) %>% 
  mutate(Count=as.numeric(Count)) %>% 
  mutate(Pos=as.numeric(Pos))

(FetPlot<-ggplot(Fetdat, aes(fill=Category, y=Count, x=Ecotype)) + 
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(x=Ecotype,y=Pos,label=Percent),size=6)+
  scale_fill_manual(values = c("#999999", "#E69F00", "#009E73"))+
  scale_y_continuous(breaks = c(0,5000,10000,15000, 20000, 25000),
                       labels = c("0","5K","10K","15K", "20K","25K"))+
  labs(y="Number of ACRs")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y=element_text(size=12,vjust=1,color="black",face="bold"),
        axis.text.x=element_text(size=12,vjust=1,color=c("#0072B2","#D55E00"),face="bold"),
        legend.position=c(0.75, 0.88),
        legend.title  = element_text(size=14,vjust=1,color="black",face="bold"),
        legend.text = element_text(size=12,vjust=1,color="black",face="bold"))
)

CF_Near<-read.delim("Data/FIL/Catalog/CF_dist.bed",
                    header = F) %>% select(V17) %>% rename(Distance=V17) %>% 
                    mutate(Ecotype="Coastal")
  

CH_Near<-read.delim("Data/HAL/Catalog/CH_dist.bed",
                    header = F) %>% select(V17) %>% rename(Distance=V17) %>% 
                    mutate(Ecotype="Inland")

Closestdat<-rbind(CF_Near,CH_Near) %>% mutate(Distance=ifelse(Distance>10000,10000,Distance))

(ClosePlot<-ggplot(data=Closestdat,aes(x=Distance,color=Ecotype))+
  geom_density(size=1.1,linetype=6,alpha=0.7)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
    scale_x_continuous(breaks = c(0,2500,5000,7500,10000),labels = c(0,2.5,5,7.5,10))+
  labs(title = "", x="Distance from the nearest gene (Kb)",y="Density")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1,color="black",face="bold"),
        legend.position=c(0.75, 0.85),
        legend.title  = element_text(size=14,vjust=1,color="black",face="bold"),
        legend.text = element_text(size=12,vjust=1,color="black",face="bold")
        )
)

GC_CF<-read.delim("Data/FIL/Catalog/CF_Summit_5K_W100_S50_Avg_PerGC_mod.tab",
                  header = F) %>% select(V3,V5) %>% 
                  rename(Index=V3,GC=V5) %>% 
                  mutate(Position= (Index-201)*25 ) %>% 
                  mutate(Ecotype="Coastal")

GC_CH<-read.delim("Data/HAL/Catalog/CH_Summit_5K_W100_S50_Avg_PerGC_mod.tab",
                  header = F) %>% select(V3,V5) %>% 
  rename(Index=V3,GC=V5) %>% 
  mutate(Position= (Index-201)*25 ) %>% 
  mutate(Ecotype="Inland")

GCdat<-rbind(GC_CF,GC_CH)

(GCPlot<-ggplot(data=GCdat,aes(x=Position,y=GC*100,color=Ecotype))+
    geom_line(size=1.1,linetype=6,alpha=0.95,show.legend = F)+
    geom_vline(xintercept = -47,linetype=2,size=0.3)+
    scale_x_continuous(limits = c(-2000,2000),breaks = c(-2000,0,2000),
                       labels = c("-2K","Summit","2K"))+
    scale_y_continuous(limits = c(10,90),
                       breaks = c(25,50,75),
                       labels = c(25,50,75))+
    scale_color_manual(values = c("#0072B2","#D55E00"))+
    labs(title = "", x="",y="GC Content (%)")+
    theme_bw()+
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
          axis.title = element_text(size=14,face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=12,vjust=1,color="black",face="bold") )

)

CF_TSSEn<-read.delim("Data/FIL/Catalog/CF_TSS_Avg_NACR.tab",
                     header = F) %>% select(V3,V5) %>% 
                      rename(Index=V3,ACR=V5) %>% 
                      mutate(Position= ((Index-100)*50)-150 )%>% 
                      mutate(Ecotype="Coastal") %>% 
                      mutate(RACR=ACR/max(ACR)  )

CH_TSSEn<-read.delim("Data/HAL/Catalog/CH_TSS_Avg_NACR.tab",
                     header = F) %>% select(V3,V5) %>% 
  rename(Index=V3,ACR=V5) %>% 
  mutate(Position= ((Index-100)*50)-150 )%>% 
  mutate(Ecotype="Inland") %>% 
  mutate(RACR=ACR/max(ACR)  )
  
TSSEndat<-rbind(CF_TSSEn,CH_TSSEn)

(TSSEnPlot<-ggplot(data=TSSEndat,aes(x=Position,y=RACR,color=Ecotype))+
    geom_line(size=1.1,linetype=6,alpha=0.95,show.legend = F)+
    geom_vline(xintercept = 0,linetype=2,size=0.3)+
    scale_x_continuous(limits = c(-2000,2000),breaks = c(-2000,0,2000),
                       labels = c("-2K","TSS","2K"))+
    scale_y_continuous(limits = c(.10,1),
                       breaks = c(.25,.50,.75),
                       labels = c(.25,.50,.75))+
    scale_color_manual(values = c("#0072B2","#D55E00"))+
    labs(title = "", x="",y="Relaive ACR frequency")+
    theme_bw()+
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
          axis.title = element_text(size=14,face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=12,vjust=1,color="black",face="bold") )
  
)

library(cowplot)
#legends <- get_legend(p3legend)

poly.plot <-
  ggdraw() +
  draw_plot(ClosePlot,x = 0, y = 0.6,width = 0.3,height = 0.4) +
  draw_plot(LenDenPlot,x = 0.3, y = 0.6,width = 0.3,height = 0.4) +
  draw_plot(FetPlot,x = 0.59, y = 0.6,width = 0.4,height = 0.4) +
  draw_label(x=0.02,y=0.98,label = "A)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.32,y=0.98,label = "B)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.62,y=0.98,label = "C)", color = "black", size = 16, fontface = "bold")+
  draw_plot(MetaPlot,x = 0.0, y = 0.3,width = 0.3,height = 0.3) +
  draw_plot(GCPlot,x = 0.3, y = 0.3,width = 0.3,height = 0.3) +
  draw_plot(TSSEnPlot,x = 0.6, y = 0.3,width = 0.3,height = 0.3) +
  draw_label(x=0.02,y=0.58,label = "D)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.32,y=0.58,label = "E)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.62,y=0.58,label = "F)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.02,y=0.25,label = "I)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.52,y=0.25,label = "G)", color = "black", size = 16, fontface = "bold")
  
  
tiff("Plots/Fig2.tiff",width=10,height=11,units="in",res=300)
poly.plot
dev.off()

####### Supplement Fig S2
datMatCov<-read.delim("Data/SELECTED_MACSPool_FE3_Clean_Intersect/Replicate_TSS_TAG/R1_CH.dedup_b10.bw_Summit.tab",header = F,sep = "\t",skip = 3,na.strings = c(NA,"nan"))


CH<-as.data.frame(cbind(Count=colSums(datMatCov,na.rm = T),
                        Pos=seq(-4999,4999,10),
                        
                        Sample=rep("R1_CH",1000)
))

CH$Dens<-as.numeric(as.character(CH$Count))/sum(as.numeric(as.character(CH$Count)))
write.csv(CH,file="Data/SELECTED_MACSPool_FE3_Clean_Intersect/Replicate_TSS_TAG/R1_CH_TSSTAG_Normalized_v3.csv",row.names = F)

### Coastal
CRep1<-read.csv("Data/SELECTED_MACSPool_FE3_Clean_Intersect/Replicate_TSS_TAG/R1_CF_TSSTAG_Normalized_v3.csv")
CRep2<-read.csv("Data/SELECTED_MACSPool_FE3_Clean_Intersect/Replicate_TSS_TAG/R2_CF_TSSTAG_Normalized_v3.csv")
CRep3<-read.csv("Data/SELECTED_MACSPool_FE3_Clean_Intersect/Replicate_TSS_TAG/R3_CF_TSSTAG_Normalized_v3.csv") 

CRep1$REP<-rep("REP1",nrow(CRep1))
CRep2$REP<-rep("REP2",nrow(CRep2))
CRep3$REP<-rep("REP3",nrow(CRep3))

CREP<-rbind(CRep1,CRep2,CRep3)
CREP$Ecotype<-rep("Coastal",nrow(CREP))

#####Inland
IRep1<-read.csv("Data/SELECTED_MACSPool_FE3_Clean_Intersect/Replicate_TSS_TAG/R1_CH_TSSTAG_Normalized_v3.csv")
IRep2<-read.csv("Data/SELECTED_MACSPool_FE3_Clean_Intersect/Replicate_TSS_TAG/R2_CH_TSSTAG_Normalized_v3.csv")
IRep3<-read.csv("Data/SELECTED_MACSPool_FE3_Clean_Intersect/Replicate_TSS_TAG/R3_CH_TSSTAG_Normalized_v3.csv") 

IRep1$REP<-rep("REP1",nrow(IRep1))
IRep2$REP<-rep("REP2",nrow(IRep2))
IRep3$REP<-rep("REP3",nrow(IRep3))

IREP<-rbind(IRep1,IRep2,IRep3)
IREP$Ecotype<-rep("Inland",nrow(IREP))

datREPTSS<-rbind(CREP,IREP)
datREPTSS$Count[which(is.infinite(datREPTSS$Count))]<-NA

(SupMetaPlot<-as_tibble(datREPTSS) %>% 
    drop_na() %>% 
    #ggplot(aes(y=Count,x=Pos,color=Ecotype))+
    ggplot(aes(y=Dens,x=Pos,color=Ecotype))+
    geom_line(size=1.1,linetype=1,alpha=0.95,show.legend = F)+
    geom_vline(xintercept = 0,linetype=2,size=0.3)+
    scale_x_continuous(limits = c(-4991,4999),breaks = c(-4991,0,4999),
                       labels = c("-5K","TSS","5K"))+
    # scale_y_continuous(limits = c(0.075,0.33),
    #                    breaks = c(0.2,0.3),
    #                    labels = c(0.2,0.3))+
    scale_color_manual(values = c("#0072B2","#D55E00"))+
    #labs(title = "", x="",y="Normalized Count")+
    labs(title = "", x="",y="Density")+
    facet_wrap(Ecotype~REP)+
    #           scales = "free")+
    theme_bw()+
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.title = element_text(hjust = 0.5,size=12, face = "bold"),
          axis.title = element_text(size=14,face="bold"),
          strip.text.x = element_text(size = 14,face="bold"),
          strip.background  = element_rect(colour = "transparent"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=14,vjust=1,color="black",face="bold") )
)

tiff("Plots/Supplementary_FigS2.tiff",width=10,height=10,units="in",res=100)
SupMetaPlot
dev.off()
