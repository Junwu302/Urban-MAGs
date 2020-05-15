## plot the geographic distribution maps
library(plyr)
library(maptools)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(ggpubr )
library(cowplot)
library(plotly)
source("theme_basic.R")
load("./data/Sample_Info.RData")
df =  ddply(Sample_Info, .variables = c("City_Name","city_latitude","city_longitude"),nrow)
df = df[!is.na(df$city_latitude),]
colnames(df)[1:3] = c("city","latitude","longitude")
df$sampleNum = "< 10"
df$sampleNum[df$V1>=10 & df$V1<100] = "10 ~ 100"
df$sampleNum[df$V1>=100 & df$V1<500] = "100 ~ 500"
df$sampleNum[df$V1>=500] = "> 500"
df$sampleNum = factor(df$sampleNum, levels = c("> 500","100 ~ 500","10 ~ 100","< 10"))
df = df[order(df$V1),]

mp  = NULL #定义一个空的地图
mapworld = borders("world",colour = "gray40",fill="white",size=.25) #绘制基本地图
g1_A = ggplot(data = df) + mapworld + coord_cartesian() +
  geom_point(aes(x=longitude,y=latitude,color=sampleNum),alpha = 0.7,shape=16,size = 5)+
  labs(color = "Sample Number")  + scale_color_brewer(palette = "Set1") + 
  scale_y_continuous(breaks = (-2:2) * 45) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  xlab("Longitude") + ylab("Latitude") + theme_0 +
  theme(legend.position = c(0.15,0.25))

load("./data/checkM_Res.RData")
checkM_Res$Quality[checkM_Res$Quality == ""] = "Low"
checkM_Res$Quality = factor(checkM_Res$Quality, levels = c("Low","Medium","High"))
df = checkM_Res[checkM_Res$Contamination <10,c("Completeness","Contamination","Strain_heterogeneity","Quality")]
df = df[order(df$Quality),]
fig = plot_ly(df, x = ~Completeness, z = ~Contamination, y = ~Strain_heterogeneity, 
              color = ~Quality, size = .5,
              colors = c('#BEBEBE','#0C4B8E','#BF382A'))
fig = fig %>% add_markers()
fig = fig %>% layout(scene = list(xaxis = list(title = 'Completeness'),
                                  yaxis = list(title = 'Strain heterogeneity'),
                                  zaxis = list(title = 'Contamination')))
fig = plot_ly(df, x = ~Completeness, z = ~Contamination, y = ~Strain_heterogeneity, 
              color = ~Quality, size = .5,
              colors = c('#BEBEBE','#0C4B8E','#BF382A'))
fig = fig %>% add_markers()
fig = fig %>% layout(scene = list(xaxis = list(title = 'Completeness'),
                                  yaxis = list(title = 'Strain heterogeneity'),
                                  zaxis = list(title = 'Contamination')))
fig

load("data/MAGs_NCBIRef_ANI_AQ.RData")
load("data/MAGs_IGG_ANI_AQ.RData")
all_MAGs = gsub(".fa$","",MAGs_IGG_ANI_AQ$Query)
NCBIRef_MAGs = gsub(".fa$","",MAGs_NCBIRef_ANI_AQ$Query[MAGs_NCBIRef_ANI_AQ$ANI >=95 &
                                                          MAGs_NCBIRef_ANI_AQ$Ref.AQ >= 60 &
                                                          MAGs_NCBIRef_ANI_AQ$Query.AQ >= 60])
IGG_MAGs = gsub(".fa$","",MAGs_IGG_ANI_AQ$Query[MAGs_IGG_ANI_AQ$ANI >=95 &
                                                  MAGs_IGG_ANI_AQ$Ref.AQ >= 60 &
                                                  MAGs_IGG_ANI_AQ$Query.AQ >= 60])
length(all_MAGs)
length(NCBIRef_MAGs)
length(IGG_MAGs)
sum(NCBIRef_MAGs %in% IGG_MAGs)

library(ggplot2)
library(ggimage)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(plyr)
library(reshape2)

## figures/Fig_2A.pdf is the venn plot generated usig Visio
# df = data.frame(x=1, y = 1)
# g1_E = ggplot(data = df, mapping = aes(x,y)) + geom_image(image = "figures/Fig_2A.pdf",size=0.8) +
#   xlab("") + ylab("") + theme_nothing()+
#   theme(axis.text = element_blank(),plot.margin = unit(c(.5, 1, .5, 1), "cm"))
# 

load("./data/Sample_Info.RData")
MAGs_df = data.frame(MAGs = all_MAGs,stringsAsFactors = F)
MAGs_df$uuid = gsub("_bin.*","",MAGs_df$MAGs)
MAGs_df$Matched = "None"
MAGs_df$Matched[MAGs_df$MAGs %in% NCBIRef_MAGs] = "NCBI Ref"
MAGs_df$Matched[MAGs_df$MAGs %in% IGG_MAGs] = "IGG"


MAGs_df = merge(MAGs_df, Sample_Info[,c("City_Name","uuid","Continent_Name")])
MAGs_df = ddply(MAGs_df, .variables = c("uuid","Continent_Name"), .fun = function(df){
  Total_Num = nrow(df)
  unmatched_num = sum(df$Matched == "None") 
  unmatched_frac = 100*sum(df$Matched == "None")/Total_Num
  return(c(Total_Num, unmatched_num, unmatched_frac))
})
colnames(MAGs_df)[-c(1,2)] = c("Total_Num","unmatched_num", "unmatched_frac")

MAGs_df = MAGs_df[MAGs_df$Total_Num >= 5,]
ind = ddply(MAGs_df, .variables = "Continent_Name", .fun = function(x){
  median(x$unmatched_frac)
})
ind = ind[order(ind$V1, decreasing = T),]
MAGs_df$Continent_Name = factor(MAGs_df$Continent_Name, levels = ind$Continent_Name)
g1_D = ggplot(data = MAGs_df, mapping = aes(x=Continent_Name, y=unmatched_frac)) +
  geom_boxplot(mapping = aes(fill = Continent_Name), outlier.shape = NA) +
  geom_jitter(colour = "gray30",size = 0.1, alpha=0.7) +
  ylab("Fraction of unmatched MAGs") + xlab("") +
  scale_fill_npg() + theme_1

fig_wd = 8
pdf("figures/Fig1_A.pdf", width = fig_wd, height = 0.75*fig_wd)
g1_A
dev.off()
pdf("figures/Fig1_B.pdf", width = fig_wd, height = 0.6*fig_wd)
g1_B
dev.off()

# pdf("figures/Fig1_E.pdf", width = 450, height = 450)
# g1_E
# dev.off()
pdf("figures/Fig1_D.pdf", width = fig_wd, height = 0.6*fig_wd)
g1_D
dev.off()


