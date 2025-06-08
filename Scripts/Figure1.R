setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/ACR_Divergence/")
library(gridExtra)
library(ggplot2)
library(xlsx)
library(tidyverse)

###### Get CNS relationship
library("ggtree")
library(treeio)
Tree<-read.tree("Data/PH_CNS/CNS_Tree")
Tree$tip.label<-c("PHC","PHI","PVK","SV")

(plotTreeCore<-ggplot(Tree, aes(x, y)) + 
    geom_tree(linetype=2) + 
    theme_tree()+
    #geom_tiplab(parse=T,size=4, offset=-0.9,aes(angle=0,vjust=-0.25),color="blue")+
    #geom_text(aes(x=branch,label=Tree$nodelab),
    #          show.legend = F, color="red",size=4,fontface="bold",angle=45, vjust=2)+
    #geom_hilight(mapping=aes(subset = node %in% c(1, 2,3), fill="steelblue", alpha=.6))+
    annotate("rect", xmin = -0.5, xmax = 17 , ymin = 0.5, ymax = 4.5,
             alpha = .2,fill = NA,color="#009E73",size=1.5)+
    annotate("rect", xmin = 1, xmax = 17, ymin = 1.8, ymax = 4.5,
             alpha = .4,color = "#D55E00", fill=NA,size=1.5)+
    annotate("text", x = 13.5, y = 0.75, label = "core-Panicoid",
             fontface =1 ,color="#009E73",
               parse = TRUE,size=4)+
    annotate("text", x = 13.5, y = 4.25, label = "core-Panicum",
             fontface =1 ,color="#D55E00",
             parse = TRUE,size=4)+
    geom_text( x = 14, y = 1, label = paste0("italic('S. viridis')", "") ,
             parse = TRUE,size=3)+
    geom_text( x = 14, y = 2, label = "italic('P. virgatum K-subgenome')" ,
               parse = TRUE,size=3)+
    geom_text( x = 14, y = 3, label = "italic('P. hallii var. filipes')",
               parse = TRUE,size=3)+
    geom_text( x = 14, y = 4, label = "italic('P. hallii var. hallii')",
               parse = TRUE,size=3)+
    # geom_text( x = 14.5, y = 5, label = "italic('P. virgatum K-subgenome')",
    #            parse = TRUE,size=3)+
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          axis.text=element_blank(),
          #legend.position=c(1.05, 0.4),
          legend.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_blank())
)

(plotTreeLInland<-ggplot(Tree, aes(x, y)) + 
    geom_tree(linetype=2) + 
    theme_tree()+
    #geom_tiplab(parse=T,size=4, offset=-0.9,aes(angle=0,vjust=-0.25),color="blue")+
    #geom_text(aes(x=branch,label=Tree$nodelab),
    #          show.legend = F, color="red",size=4,fontface="bold",angle=45, vjust=2)+
    #geom_hilight(mapping=aes(subset = node %in% c(1, 2,3), fill="steelblue", alpha=.6))+
    annotate("rect", xmin = -0.5, xmax = 17, ymin = 0.5, ymax = 3.6,
             alpha = .2,fill = NA,color="#3949AB",size=1.5)+
    annotate("text", x = 13.5, y = 0.75, label = "pan-lost-Inland",
             fontface =1 ,color="#3949AB",
             parse = TRUE,size=4)+
    geom_text( x = 14, y = 1, label = paste0("italic('S. viridis')", "") ,
               parse = TRUE,size=3)+
    geom_text( x = 14, y = 2, label = "italic('P. virgatum K-subgenome')" ,
               parse = TRUE,size=3)+
    geom_text( x = 14, y = 3, label = "italic('P. hallii var. filipes')",
               parse = TRUE,size=3)+
    geom_text( x = 14, y = 4, label = "italic('P. hallii var. hallii')",
               parse = TRUE,size=3)+
    # geom_text( x = 14.5, y = 5, label = "italic('P. virgatum K-subgenome')",
    #            parse = TRUE,size=3)+
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          axis.text=element_blank(),
          #legend.position=c(1.05, 0.4),
          legend.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_blank())
)

(plotTreeLCoastal<-ggplot(Tree, aes(x, y)) + 
    geom_tree(linetype=2) + 
    theme_tree()+
    #geom_tiplab(parse=T,size=4, offset=-0.9,aes(angle=0,vjust=-0.25),color="blue")+
    #geom_text(aes(x=branch,label=Tree$nodelab),
    #          show.legend = F, color="red",size=4,fontface="bold",angle=45, vjust=2)+
    #geom_hilight(mapping=aes(subset = node %in% c(1, 2,3), fill="steelblue", alpha=.6))+
    annotate("rect", xmin = -0.5, xmax = 17, ymin = 0.5, ymax = 3.6,
             alpha = .2,fill = NA,color="#CB4335",size=1.5)+
    annotate("text", x = 13.5, y = 0.75, label = "pan-lost-Coastal",
             fontface =1 ,color="#CB4335",
             parse = TRUE,size=4)+
    geom_text( x = 14, y = 1, label = paste0("italic('S. viridis')", "") ,
               parse = TRUE,size=3)+
    geom_text( x = 14, y = 2, label = "italic('P. virgatum K-subgenome')" ,
               parse = TRUE,size=3)+
    geom_text( x = 14, y = 4, label = "italic('P. hallii var. filipes')",
               parse = TRUE,size=3)+
    geom_text( x = 14, y = 3, label = "italic('P. hallii var. hallii')",
               parse = TRUE,size=3)+
    # geom_text( x = 14.5, y = 5, label = "italic('P. virgatum K-subgenome')",
    #            parse = TRUE,size=3)+
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          axis.text=element_blank(),
          #legend.position=c(1.05, 0.4),
          legend.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_blank())
)

########## CNS stats
SumTab<-as.data.frame(cbind(Class=c("core-Panicoid", "core-Panicum", 
                                    "pan-lost-Coastal","pan-lost-Inland"),
                            Total_Number=c(56512,35115,2538,2565),
                            Ang_Length=c(195,247,161,149)),
                      row.names = F)


(PlotNum<-  SumTab  %>% mutate (Total_Number=as.numeric(Total_Number))  %>% 
    ggplot()+
    geom_col(aes(x=Class,y=Total_Number,fill=Class),alpha=0.9,show.legend = F)+
    scale_fill_manual(values = c("#009E73","#D55E00","#3949AB","#CB4335"))+
    labs(title = "", y="Number of CNS",
         x="",
         fill="Class")+
    theme_bw()+
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.title = element_text(hjust = 0.5,size=16, face = "bold"),
          axis.title.x = element_text(size=10,face="bold"),
          axis.title.y = element_text(size=12,face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(size=10,vjust=0.3,angle = 45,color="black",face="bold"),
          axis.text.y=element_text(size=12,vjust=1,color="black",face="bold"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0.5,size=16, face = "bold"),
          legend.position=c(1.05, 0.4),
          legend.background = element_rect(fill = "transparent"), 
          legend.title = element_text(size=14,vjust=1,color="black",face="bold"),
          legend.text = element_text(size=12,vjust=1,color="black",face="bold")
    )
)

Distance<-read.delim("Data/PH_CNS/Dist/ALL_CNS_GENEID_Dist.clean.bed",header = F,
                     col.names = c("Chr","Start","End","GeneID","Distance","Category") ) %>% 
  select(Distance, Category) %>% 
  mutate(Distance=abs(Distance)) %>% 
  mutate(Distance=ifelse(Distance > 2000, 2000, Distance)) %>% 
  mutate(Category=ifelse(Category=="pan-lost-Inland","pan-lost-Inland",Category)) %>% 
  mutate(Category=ifelse(Category=="pan-lost-Coastal","pan-lost-Coastal",Category))

(PlotDist<-ggplot(Distance) +
  geom_density(aes(x=Distance,color=Category),size=1.2)+
  scale_x_continuous(limits = c(0,2000))+
  scale_color_manual(values = c("#009E73","#D55E00","#3949AB","#CB4335"))+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5,size=16, face = "bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1,color="black",face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0.5,size=16, face = "bold"),
        legend.position=c(0.67, 0.8),
        legend.background = element_rect(fill = "transparent"), 
        legend.title = element_text(size=14,vjust=1,color="black",face="bold"),
        legend.text = element_text(size=14,vjust=1,color="black",face="bold")
  )
)

library(cowplot)
#legends <- get_legend(p3legend)

Fig1 <-
  ggdraw() +
  draw_plot(plotTreeCore,x = 0, y = 0.65,width = 0.45,height = 0.33) +
  draw_plot(plotTreeLInland,x = 0.,y=0.33,width = 0.45,height = 0.33) +
  draw_plot(plotTreeLCoastal,x = 0,y=0,width = 0.45,height = 0.33) +
  draw_plot(PlotNum,x = 0.5,y=0.5,width = 0.5,height = 0.5) +
  draw_plot(PlotDist,x = 0.5,y=0,width = 0.5,height = 0.5) +
  draw_label(x=0.02,y=0.98,label = "A)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.52,y=0.98,label = "B)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.52,y=0.5,label = "C)", color = "black", size = 16, fontface = "bold")



tiff("Plots/Fig1.tiff",width=10,height=8,units="in",res=300)
Fig1
dev.off()

