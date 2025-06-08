setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/ACR_Divergence/")
library(ggplot2)
library(tidyverse)
library(limma)
library(DESeq2)
library(readxl)
library(xlsx)


###### Test for DEG between ecotype and if these are enriched with Nonr-ACR
Design<-read.csv("Data/TAGSeq_2021_Design.csv")
Design<-Design[which(Design$Treatment=="Control" ),]
rownames(Design)<-Design$SampleID

countH<-read.csv("Data/TAGSeq_Count_REFHAL.csv")
colnames(countH)[1]<-"HAL"
countH<-countH[,c(1,which(colnames(countH) %in% Design$SampleID))]


countF<-read.csv("Data/TAGSeq_Count_REFFIL.csv")
colnames(countF)[1]<-"FIL"
countF<-countF[,c(1,which(colnames(countF) %in% Design$SampleID))]

Ortho<-read.csv("Data/FIL_HAL_1-1_OrthologPairs_OrthoFinder_m50.csv")
countOrtho<-merge(merge(Ortho,countH,by="HAL",
                        all.x=F,all.y=F), countF,by="FIL",
                  all.x=F,all.y=F )
rownames(countOrtho)<-countOrtho$HAL
countOrtho<-countOrtho[,-c(1:2)]
countTable<-countOrtho[,rownames(Design)]
countTable<-countTable[-which(as.vector(rowMeans(countTable,na.rm = F)) <5),]

dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = Design,
                              design = ~ 1)
dds

dds = estimateSizeFactors(dds)
sizeFactors(dds)

####### DEG Analysis

Design$Genotype<-factor(Design$Genotype,levels = c("FIL","HAL"))
length(levels(Design$Genotype))

se<-SummarizedExperiment(assays = data.matrix(countTable),
                         colData = DataFrame(Design))

##### G
ddsG<- DESeqDataSet(se = se, design = as.formula("~Genotype"))
cat(paste(c("Start Running \"G Gene\" Analysis at",date(), "........\n")))
desG<-DESeq(ddsG,test="LRT", reduced= as.formula("~1"), parallel = T)

resultsNames(desG)

GL_G<-results(desG,test = "LRT",
              parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)
GL_G<-GL_G[which(GL_G$padj < 0.05),]

GL_G_M<-merge(GL_G,Ortho,by.x="row",by.y="HAL",
              all.x=F,all.y=F)

write.csv(GL_G_M,"Results/DEG_Analysis_Control.csv",row.names = F)

#### Get genewise expression
vst<-vst(dds,blind = F)
dat<-as.data.frame(assay(vst))
tdat<-as.data.frame(t(dat))
tdat$LIBID<-rownames(tdat)
Design$LIBID<-rownames(Design)
VSTDes<-merge(Design,tdat,by="LIBID")

gene<-"PhHAL.5G264600"
gene<-"PhHAL.2G476000" #WRKY
gene<-"PhHAL.1G090900" #F-box like
gene<-"PhHAL.4G259400" #MATE efflux family protein

p<-as_tibble(VSTDes) %>% dplyr::select(Genotype, gene) %>% group_by(Genotype) %>%
  summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>%
  mutate(Group=paste(Tissue,Ecotype,sep = ":"))%>% ggplot()+
  geom_point(aes(x=Treatment, y=mean,color=Ecotype,shape=Tissue),alpha=1,size=3,show.legend = T)+
  geom_line(aes(x=Treatment, y=mean,group=Group,color=Ecotype,linetype=Tissue),alpha=0.75,size=1,show.legend = T)+
  geom_linerange(aes(x=Treatment,ymin=mean-se,ymax=mean+se,color=Ecotype,linetype=Tissue),alpha=0.75,size=0.5,show.legend = F)+
  scale_color_manual(values = c( "#3949AB","#CB4335" ))+
  scale_x_discrete(expand = c(0.1,0.1))+
  labs(x="Treatment",y= "Normalized Count",title = Descript,color="Ecotype")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12,face="bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1),
        strip.background = element_rect(fill = "transparent", color = NA),
        strip.text = element_text(size=12,face = "bold")
        #axis.ticks.x = element_blank()
  )
print(p)

p<-plotCounts(dds,gene =gene , intgroup = "Genotype",returnData = T) %>% group_by(Genotype) %>%
  summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>%
  mutate(Genotype=if_else(Genotype=="FIL","Coastal","Inland")) %>% 
  ggplot()+
  geom_col(aes(x=Genotype, y=mean,color=Genotype,fill=Genotype),alpha=1,show.legend = T)+
  geom_errorbar(aes(x=Genotype,ymin=mean-se,ymax=mean+se),
                 color="black",alpha=0.75,size=0.5,show.legend = F,width=0.1)+
  scale_color_manual(values = c( "#3949AB","#CB4335" ))+
  scale_fill_manual(values = c( "#3949AB","#CB4335" ))+
  #scale_x_discrete(expand = c(0.1,0.1))+
  labs(x="Genotype",y= "Expression Count",color="Genotype")+
  theme_bw(base_size = 18)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12,face="bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=18,vjust=1),
        strip.background = element_rect(fill = "transparent", color = NA),
        strip.text = element_text(size=12,face = "bold"),
        legend.position = c(0.8,0.8)
        #axis.ticks.x = element_blank()
  )
print(p)


