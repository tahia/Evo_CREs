setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/ACR_Divergence/")
library(gridExtra)
library(ggplot2)
library(xlsx)
library(tidyverse)
library(ggplot2)
library(cowplot)

################ Coastal Loss
# PhHAL.1G090900 Pahal.1G093100 :F-box protein lost CNS in coastal
# CH_peak: Chr01:5977481-5977873
# pan-lost-coastal-CNS: Chr01:5977682-5977718
# PhHAL.1G090900: Chr01:5974189-5977828
# Pahal.1G093100: Chr01:6625456-6628905

FILLostBED<-read.csv("Data/BEDGraph/Pahal.1G093100.count.bed",
                 header = F,sep = "\t",col.names =c("Chr","Start","End","count"))
FILLostBED$Index<-seq(1:nrow(FILLostBED))

(FILLost<-ggplot(FILLostBED)+geom_col(aes(x=Index,y=count),width=1,color="#3949AB",fill="#3949AB")+
  scale_y_continuous(limits = c(-50,50),labels = NULL, breaks = NULL,name = "")+
  scale_x_continuous(labels = NULL, breaks = NULL, name = "")+
  geom_segment(xend = 9, yend = -42, x = 187, y = -42,
               arrow = arrow(length=unit(0.30,"cm")))+
  annotate("rect", xmin = 14, xmax = 44, ymin = -44, ymax = -40,
           fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
  annotate("rect", xmin = 58, xmax = 68, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
  annotate("rect", xmin = 113, xmax = 187, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
  annotate("text", x = 100, y = -48,label = "Pahal.1G093100",color="black",size=5)+
  annotate("text", x = 3, y = -48,label = "+1 Kb",color="black",size=5)+
  annotate("text", x = 217, y = -48,label = "-1 Kb",color="black",size=5)+
  #annotate("rect", xmin = 39, xmax = 176, ymin = -12, ymax = -4,
  #         alpha = .2,fill = NA,color="black",linewidth=1,linetype=1)+
  #annotate("text", x = 15, y = -30,color="black",label="ACR peak",size=5)+
  #annotate("text", x = 40, y = -43,color="black",label="pan-lost-Coastal-CNS",size=5)+
  annotate("text", x = -4, y = 16,label = "ACR",color="black",size=6,angle=90)+
  annotate("text", x = -4, y = -25,label = "CNS",color="black",size=6,angle=90)+
  geom_hline(yintercept = -17,color="grey40", lty=2)+
  geom_hline(yintercept = -32,color="grey40", lty=2)+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5,size=16, face = "bold"),
        axis.line = element_blank(), 
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=12,vjust=1,color="black",face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0.5,size=16, face = "bold"),
        legend.position=c(1.05, 0.4),
        legend.background = element_rect(fill = "transparent"), 
        legend.title = element_text(size=14,vjust=1,color="black",face="bold"),
        legend.text = element_text(size=12,vjust=1,color="black",face="bold")
  )

)

# CH_peak: Chr01:5977481-5977873
# PhHAL.1G090900: Chr01:5974189-5977828
# pan-lost-coastal-CNS: Chr01:5977682-5977718


HALGainBED<-read.csv("Data/BEDGraph/PhHAL.1G090900.count.bed",
                 header = F,sep = "\t",col.names =c("Chr","Start","End","count"))
HALGainBED$Index<-seq(1:nrow(HALGainBED))

(HALGain<-ggplot(HALGainBED)+geom_col(aes(x=Index,y=count),width=1,color="#CB4335",fill="#CB4335")+
  scale_y_continuous(limits = c(-50,50),labels = NULL, breaks = NULL,name = "")+
  scale_x_continuous(labels = NULL, breaks = NULL,name = "")+
  # annotate("rect", xmin = 39, xmax = 184, ymin = -12, ymax = -4,
  #          alpha = .2,fill = NA,color="black",linewidth=1,linetype=1)+
  geom_hline(yintercept = -17,color="grey40", lty=2)+
  geom_hline(yintercept = -32,color="grey40", lty=2)+
  geom_segment(xend = 9, yend = -42, x = 187, y = -42,
                 arrow = arrow(length=unit(0.30,"cm")))+
  annotate("rect", xmin = 14, xmax = 44, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
  annotate("rect", xmin = 56, xmax = 68, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
  annotate("rect", xmin = 113, xmax = 187, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
  annotate("rect", xmin = 172, xmax = 185, ymin = -12, ymax = -6,
           alpha = .9,fill = "#CB4335",color="#CB4335",linewidth=1,linetype=1)+
  #annotate("text", x = 15, y = -30,color="black",label="ACR peak",size=5)+
  annotate("rect", xmin = 178, xmax = 179, ymin = -28, ymax = -22,
           alpha = .9,fill = "orange2",color="orange2",linewidth=1,linetype=1)+
  annotate("text", x = 115, y = -48,label = "PhHAL.1G090900",color="black",size=5)+
  annotate("text", x = 3, y = -48,label = "+1 Kb",color="black",size=5)+
  annotate("text", x = 225, y = -48,label = "-1 Kb",color="black",size=5)+
  annotate("rect", xmin = 165, xmax = 190, ymin = -30, ymax = 48,alpha=0.1,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
  #annotate("text", x = 44, y = -43,color="black",label="pan-lost-Coastal-CNS",size=5)+
  annotate("text", x = -4, y = 16,label = "ACR",color="black",size=6,angle=90)+
  annotate("text", x = -4, y = -25,label = "CNS",color="black",size=6,angle=90)+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5,size=16, face = "bold"),
        axis.line = element_blank(), 
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=12,vjust=1,color="black",face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0.5,size=16, face = "bold"),
        legend.position=c(1.05, 0.4),
        legend.background = element_rect(fill = "transparent"), 
        legend.title = element_text(size=14,vjust=1,color="black",face="bold"),
        legend.text = element_text(size=12,vjust=1,color="black",face="bold")
  )

)

########## Inland Loss
#PhHAL.4G259400 Pahal.4G245800 MATE efflux family protein
# CF_peak: Chr04_44448238_44449042_400
# pan-lost-coastal-CNS: Chr04:44449002-44449020
# PhHAL.4G259400: Chr04:42764717-42771566
# Pahal.4G245800: Chr04:44439278-44457198


FILGainBED<-read.csv("Data/BEDGraph/Pahal.4G245800.count.bed",
                     header = F,sep = "\t",col.names =c("Chr","Start","End","count"))
FILGainBED$Index<-seq(1:nrow(FILGainBED))

#scaled to 1/3 to match with A
(FILGain<-ggplot(FILGainBED)+geom_col(aes(x=Index,y=count/3),width=1,color="#3949AB",fill="#3949AB")+
    scale_y_continuous(limits = c(-50,50),labels = NULL, breaks = NULL,name = "")+
    scale_x_continuous(labels = NULL, breaks = NULL,name = "")+
    # annotate("rect", xmin = 39, xmax = 678, ymin = -12, ymax = -4,
    #          alpha = .2,fill = NA,color="black",linewidth=1,linetype=1)+
    geom_hline(yintercept = -17,color="grey40", lty=2)+
    geom_hline(yintercept = -32,color="grey40", lty=2)+
    geom_segment(xend = 22, yend = -42, x = 679, y = -42,
                 arrow = arrow(length=unit(0.30,"cm")))+
    annotate("rect", xmin = 44, xmax = 75, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 95, xmax = 105, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 125, xmax = 135, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 155, xmax = 170, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 190, xmax = 200, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 220, xmax = 230, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 255, xmax = 285, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 635, xmax = 679, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 660, xmax = 689, ymin = -12, ymax = -6,
             alpha = .9,fill = "#3949AB",color="#3949AB",linewidth=1,linetype=1)+
    annotate("rect", xmin = 357, xmax = 389, ymin = -12, ymax = -6,
             alpha = .9,fill = "#3949AB",color="#3949AB",linewidth=1,linetype=1)+
    #annotate("text", x = 40, y = -30,color="black",label="ACR peak",size=5)+
    annotate("rect", xmin = 387, xmax = 388, ymin = -28, ymax = -22,
             alpha = .9,fill = "orange2",color="orange2",linewidth=1,linetype=1)+
    annotate("text", x = 340, y = -48,label = "Pahal.4G245800",color="black",size=5)+
    annotate("text", x = 3, y = -48,label = "+1 Kb",color="black",size=5)+
    annotate("text", x = 715, y = -48,label = "-1 Kb",color="black",size=5)+
    annotate("rect", xmin = 351, xmax = 395, ymin = -30, ymax = 45,alpha=0.1,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("text", x = -4, y = 16,label = "ACR",color="black",size=6,angle=90)+
    annotate("text", x = -4, y = -25,label = "CNS",color="black",size=6,angle=90)+
    #annotate("text", x = 120, y = -43,color="black",label="pan-lost-Inland-CNS",size=5)+
    theme_bw()+
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.title = element_text(hjust = 0.5,size=16, face = "bold"),
          axis.line = element_blank(), 
          axis.title.x = element_text(size=10,face="bold"),
          axis.title.y = element_text(size=12,face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.y=element_text(size=12,vjust=1,color="black",face="bold"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0.5,size=16, face = "bold"),
          legend.position=c(1.05, 0.4),
          legend.background = element_rect(fill = "transparent"), 
          legend.title = element_text(size=14,vjust=1,color="black",face="bold"),
          legend.text = element_text(size=12,vjust=1,color="black",face="bold")
    )
  
)

#PhHAL.4G259400: Chr04:42764717-42771566

HALLostBED<-read.csv("Data/BEDGraph/PhHAL.4G259400.bed.count.bed",
                     header = F,sep = "\t",col.names =c("Chr","Start","End","count"))
HALLostBED$Index<-seq(1:nrow(HALLostBED))

(HALLost<-ggplot(HALLostBED)+geom_col(aes(x=Index,y=count/3),width=1,color="#CB4335",fill="#CB4335")+
    scale_y_continuous(limits = c(-50,50),labels = NULL, breaks = NULL,name = "")+
    scale_x_continuous(labels = NULL, breaks = NULL,name = "")+
    geom_hline(yintercept = -17,color="grey40", lty=2)+
    geom_hline(yintercept = -32,color="grey40", lty=2)+
    geom_segment(xend = 22, yend = -42, x = 240, y = -42,
                 arrow = arrow(length=unit(0.30,"cm")))+
    annotate("rect", xmin = 32, xmax = 52, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 61, xmax = 66, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 76, xmax = 81, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 91, xmax = 100, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 110, xmax = 115, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 125, xmax = 130, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 140, xmax = 160, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    annotate("rect", xmin = 220, xmax = 240, ymin = -44, ymax = -40,
             fill = "gray70",color="gray70",linewidth=0.5,linetype=1)+
    # annotate("rect", xmin = 39, xmax = 236, ymin = -12, ymax = -4,
    #          alpha = .2,fill = NA,color="black",linewidth=1,linetype=1)+
    annotate("rect", xmin = 223, xmax = 240, ymin = -12, ymax = -6,
             alpha = .9,fill = "#CB4335",color="#CB4335",linewidth=1,linetype=1)+
    annotate("text", x = 135, y = -48,label = "PhHAL.4G259400",color="black",size=5)+
    annotate("text", x = 3, y = -48,label = "+1 Kb",color="black",size=5)+
    annotate("text", x = 270, y = -48,label = "-1 Kb",color="black",size=5)+
    # annotate("text", x = 20, y = -30,color="black",label="ACR peak",size=5)+
    # annotate("text", x = 50, y = -43,color="black",label="pan-lost-Inland-CNS",size=5)+
    annotate("text", x = -4, y = 16,label = "ACR",color="black",size=6,angle=90)+
    annotate("text", x = -4, y = -25,label = "CNS",color="black",size=6,angle=90)+
    theme_bw()+
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.title = element_text(hjust = 0.5,size=16, face = "bold"),
          axis.line = element_blank(), 
          axis.title.x = element_text(size=10,face="bold"),
          axis.title.y = element_text(size=12,face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.y=element_text(size=12,vjust=1,color="black",face="bold"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0.5,size=16, face = "bold"),
          legend.position=c(1.05, 0.4),
          legend.background = element_rect(fill = "transparent"), 
          legend.title = element_text(size=14,vjust=1,color="black",face="bold"),
          legend.text = element_text(size=12,vjust=1,color="black",face="bold")
    )
  
)

Fig4 <-
  ggdraw() +
  draw_plot(FILLost,x = 0, y = 0.5,width = 0.5,height = 0.5) +
  draw_plot(HALGain,x = 0.,y=0,width = 0.5,height = 0.5) +
  draw_plot(FILGain,x = 0.5, y = 0.5,width = 0.5,height = 0.5) +
  draw_plot(HALLost,x = 0.5,y=0,width = 0.5,height = 0.5) +
  draw_label(x=0.02,y=0.98,label = "A)", color = "black", size = 16, fontface = "bold")+
  draw_label(x=0.52,y=0.98,label = "B)", color = "black", size = 16, fontface = "bold")



tiff("Plots/Fig4.tiff",width=10,height=8,units="in",res=300)
Fig4
dev.off()



