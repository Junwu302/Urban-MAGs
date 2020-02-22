## plot the geographic distribution maps
library(plyr)
library(maptools)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(ggpubr )
library(cowplot)
load("./data/Sample_Info.RData")
Sample_Info = Sample_Info[!duplicated(Sample_Info$uuid),c(1,7,11,20,23:27,29,31)]
df =  ddply(Sample_Info, .variables = c("city","city_latitude","city_longitude"),nrow)
df = df[!is.na(df$city_latitude),]
colnames(df)[1:3] = c("city","latitude","longitude")
df$sampleNum = "< 10"
df$sampleNum[df$V1>=10 & df$V1<100] = "10 ~ 100"
df$sampleNum[df$V1>=100 & df$V1<500] = "100 ~ 500"
df$sampleNum[df$V1>=500] = "> 500"
df$sampleNum = factor(df$sampleNum, levels = c("> 500","100 ~ 500","10 ~ 100","< 10"))
df = df[order(df$V1),]
#pdf("./Figures/Fig_1A.pdf", width = 12, height = 9)
mp  = NULL #定义一个空的地图
mapworld = borders("world",colour = "gray40",fill="white") #绘制基本地图
g1_A = ggplot(data = df) + mapworld + coord_cartesian() +
  geom_point(aes(x=longitude,y=latitude,color=sampleNum),alpha = 0.7,shape=16,size = 5)+
  labs(color = "Sample Number")  + scale_color_brewer(palette = "Set1") + 
  scale_y_continuous(breaks = (-2:2) * 45) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() + 
  theme(axis.title = element_text(size=18, color = "black"),
        axis.text = element_text(size=16, color = "black"),
        plot.margin = unit(c(.5, 1, .5, .5), "cm"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=18),
        legend.position = c(0.15,0.25))
#dev.off()


load("./data/checkM_Res.RData")
checkM_Res$Quality[checkM_Res$Quality == ""] = "Low"
checkM_Res$Quality = factor(checkM_Res$Quality, levels = c("Low","Medium","High"))

#pdf("./Figures/Fig_1B.pdf", width = 12, height = 3)
g1_B = ggplot(data = checkM_Res, mapping = aes(x=QS)) +
  geom_histogram(aes(y=..density..),binwidth=3,colour="black", fill="white") +
  geom_density(alpha=.2, fill="gray60") + xlim(-50,120) + ylim(0,0.02) +
  geom_vline(mapping = aes(xintercept=50),color="red",size=1, linetype=2) +
  xlab("QS = Completeness - 5*Contamination") + ylab("Density") +
  theme_classic() +
  theme(axis.title = element_text(size=16, color = "black"),
        axis.text = element_text(size=14, color = "black"),
        plot.margin = unit(c(.5, 1, .5, 1), "cm"),)
#dev.off()


df1 = checkM_Res[,c("Completeness","Quality")]
df2 = checkM_Res[,c("Contamination","Quality")]
g1_C = ggplot(data = df1, mapping = aes(x=Quality, y=Completeness)) +
  geom_boxplot(mapping = aes(fill = Quality), outlier.shape = NA) +
  geom_jitter(colour = "gray30",size = 0.1, alpha=0.7) +
  scale_fill_npg() +
  theme_classic() +
  theme(axis.title = element_text(size=16, color = "black"),
        axis.text = element_text(size=14, color = "black"),
        plot.margin = unit(c(.5, 1, .5, 1), "cm"),
        legend.position = "none")

g1_D = ggplot(data = df2, mapping = aes(x=Quality, y=Contamination)) +
  geom_boxplot(mapping = aes(fill = Quality), outlier.shape = NA) +
  geom_jitter(colour = "gray30",size = 0.1, alpha=0.7) +
  scale_fill_npg() + theme_classic() +
  theme(axis.title = element_text(size=16, color = "black"),
        axis.text = element_text(size=14, color = "black"),
        plot.margin = unit(c(.5, 1, .5, 1), "cm"),
        legend.position = "none")

g1 = as_ggplot(arrangeGrob(g1_A, g1_B, g1_C, g1_D, ncol=4,nrow=3, 
                           layout_matrix = rbind(c(1,1,2,2),c(1,1,3,4),c(1,1,3,4)))) + 
  draw_plot_label(label = c("(A)", "(B)", "(C)"), size = 16,x = c(0, 0.5, 0.5), y = c(1, 1, 0.66)) # Add labels
pdf("Figures/Fig_1.pdf",width = 18, height = 9)
g1
dev.off()



