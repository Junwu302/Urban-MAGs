library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(gridExtra)
library(ggpubr)
load("data/COG_city_adonis2.RData")
city_name = read.csv("data/city_name.csv",stringsAsFactors = F)
df = COG_city_adonis2[,c(2,1,3:7)]
colnames(df)[1:2] = c("Var1","Var2")
df = rbind(COG_city_adonis2, df)
df = merge(df, city_name[,c(1,2)],by.x = "Var1", by.y = "city")
df = merge(df, city_name[,c(1,2)],by.x = "Var2", by.y = "city")
df$Var1 = df$City_Name.x
df$Var2 = df$City_Name.y
all_cities = unique(c(df$Var1, df$Var2))
dismat = matrix(0, nrow = length(all_cities), ncol = length(all_cities))
pmat = matrix(1, nrow = length(all_cities), ncol = length(all_cities))
colnames(dismat) = rownames(dismat) = rownames(pmat) = colnames(pmat) = all_cities
for(city in rownames(dismat)){
  x = df[df$Var1==city,c("Var2","R2","pval")]
  dismat[city,x$Var2] = x$R2
  pmat[city,x$Var2] = x$pval
}
dismat[upper.tri(dismat)] = NA
pmat[upper.tri(pmat)] = NA

df1 = reshape2::melt(dismat, na.rm = TRUE)
df2 = reshape2::melt(pmat, na.rm = TRUE)
df = cbind(df1, pval=df2$value)
df$label = ""
df$label[df$pval <= 0.05 & df$pval >0.01] = "*"
df$label[df$pval <= 0.01 & df$pval >0.001] = "**"
df$label[df$pval <= 0.001] = "***"
df$Var1 = factor(df$Var1, levels = all_cities)
df$Var2 = factor(df$Var2, levels = all_cities)
colnames(df)[3] = "R2"
g5_A = ggplot(data = df, aes(Var1, Var2, fill = R2))+
  geom_tile(color = "white",na.rm = TRUE)+
  scale_fill_distiller(palette="Spectral") +
  geom_text(aes(Var1, Var2, label = label), size = 5) +
  ylab("") + xlab("") + theme_minimal()+ 
  theme(axis.title = element_text(size=14, color = "black"),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=14, color = "black")) +
  coord_fixed()


library(plyr)
library(umap)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(ggpubr)
load("data/KEGG_prevalance.RData")
load("data/KEGG_abundance.RData")
load("data/Sample_Info.RData")

KEGG_city = merge(Sample_Info[,c("City_Name","uuid","Continent_Name")],KEGG_prevalance, by="uuid")
ind = KEGG_city$City_Name %in% names(table(KEGG_city$City_Name))[table(KEGG_city$City_Name)>=20]

KEGG_city = ddply(KEGG_city[ind,], .variables = "City_Name", .fun = function(df){
  n = length(unique(df$uuid))
  m = apply(df[,-c(1:3)], 2, function(x,uuid){
    return(length(unique(uuid[x])))
  },df$uuid)
  return(m/n)
})

KEGG_sample = apply(KEGG_prevalance[,-1], 2, function(x, uuid){
  return(length(unique(uuid[x])))
},KEGG_prevalance$uuid)
KEGG_sample = KEGG_sample[order(KEGG_sample, decreasing = T)]

df2 = data.frame(t(KEGG_city[,-1]))
colnames(df2) = KEGG_city$City_Name

ind = apply(df2, 1, function(x){
  sum(x<0.5) >= length(x)
})
g5_B = as_ggplot(pheatmap(df2[!ind,], show_rownames = F, silent = T,fontsize = 14,
                          colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(100))[[4]])

## UMAP for cities
KEGG_df = merge(KEGG_abundance, Sample_Info[,c("uuid","City_Name","Continent_Name")],by="uuid")
city_col =  c(brewer.pal(n = 9, name = "Set1"),"#006745","#00f9ff","#00ff04","#e73838","#08465c")
city_num = data.frame(table(KEGG_df$City_Name), stringsAsFactors = F)
KEGG_df = KEGG_df[KEGG_df$City_Name %in% city_num$Var1[city_num$Freq >= 20],]
kegg_umap = umap(KEGG_df[,grepl("^K[0-9].*",colnames(KEGG_df))])
KEGG_df = cbind(KEGG_df[,c("City_Name","Continent_Name")], data.frame(kegg_umap$layout))

city_num = data.frame(table(KEGG_df$City_Name), stringsAsFactors = F)
KEGG_df$City_Name = factor(KEGG_df$City_Name, levels = city_num$Var1[order(city_num$Freq, decreasing = T)])

g5_C = ggplot(data = KEGG_df, mapping = aes(x=X1, y=X2, color=City_Name, shape=Continent_Name)) +
  geom_point(size=1.5) + xlab("") + ylab("")+
  scale_color_manual(values = city_col)+ theme_bw() +
  theme(axis.title = element_text(size=14, color = "black"),
        axis.text = element_text(size=12, color = "black"),
        legend.title = element_text(size = 14, color="black"),
        legend.text = element_text(size=12, color="black"))

load("data/KEGG_surface_adonis2.RData")
KEGG_surface_adonis2$Var1 = str_to_title(KEGG_surface_adonis2$Var1)
KEGG_surface_adonis2$Var2 = str_to_title(KEGG_surface_adonis2$Var2)
all_sur = unique(c(KEGG_surface_adonis2$Var1, KEGG_surface_adonis2$Var2))
dismat = matrix(0, nrow = length(all_sur), ncol = length(all_sur))
pmat = matrix(1, nrow = length(all_sur), ncol = length(all_sur))
colnames(dismat) = rownames(dismat) = rownames(pmat) = colnames(pmat) = all_sur
for(sur in rownames(dismat)){
  x = KEGG_surface_adonis2[KEGG_surface_adonis2$Var1==sur,c("Var2","R2","pval")]
  dismat[sur,x$Var2] = x$R2
  dismat[x$Var2,sur] = x$R2
  pmat[sur,x$Var2] = x$pval
  pmat[x$Var2,sur] = x$pval
}
dismat[upper.tri(dismat)] = NA
pmat[upper.tri(pmat)] = NA

df1 = reshape2::melt(dismat, na.rm = TRUE)
df2 = reshape2::melt(pmat, na.rm = TRUE)
df = cbind(df1, pval=df2$value)
df$label = ""
df$label[df$pval <= 0.05 & df$pval >0.01] = "*"
df$label[df$pval <= 0.01 & df$pval >0.001] = "**"
df$label[df$pval <= 0.001] = "***"
df$Var1 = factor(df$Var1, levels = all_sur)
df$Var2 = factor(df$Var2, levels = all_sur)
colnames(df)[3] = "R2"
g5_D = ggplot(data = df, aes(Var1, Var2, fill = R2))+
  geom_tile(color = "white",na.rm = TRUE)+
  scale_fill_distiller(palette="Spectral") +
  geom_text(aes(Var1, Var2, label = label), size = 5) +
  ylab("") + xlab("") + theme_minimal()+ 
  theme(axis.title = element_text(size=14, color = "black"),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, color = "black")) + coord_fixed()


library(plyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(Hmisc)
load("data/antisMash_Res.RData")
#load("city_info.RData")
#load("MIBIG_info.RData")
load("data/OTU_tax.RData")
load("data/Sample_Info.RData")
load("data/OTU_Abundance.RData")

antisMash_Res  = antisMash_Res[order(antisMash_Res$OTU, antisMash_Res$cluster, antisMash_Res$MIBIG_hits),]
antisMash_Res = antisMash_Res[!duplicated(antisMash_Res[,c("cluster","OTU")]),]


df = ddply(antisMash_Res, .variables = "OTU", .fun = function(x){
  Total_BGC = nrow(x)
  Known_BGC = sum(!is.na(x$MIBIG_hits))
  return(c(Total_BGC, Known_BGC))
})
colnames(df) = c("OTU","Total_BGC", "Known_BGC")
df = merge(df, OTU_tax, by="OTU",all.x = T)

novel_BGC = c(1-sum(df$Known_BGC[df$Species!= "s__"])/sum(df$Total_BGC[df$Species!= "s__"]),
              1-sum(df$Known_BGC[df$Species== "s__"])/sum(df$Total_BGC[df$Species== "s__"]))
novel_BGC = data.frame(x=c("Annotated OTU", "Novel OTU"), y=novel_BGC*100,
                       label=paste(round(novel_BGC*100,2),"%",sep=''),stringsAsFactors = F)
ggplot(novel_BGC, aes(x=x, y=y, fill=x))+
  geom_bar(stat="identity", width = 0.6)+
  geom_text(aes(y=y, label=label), hjust=-0.2, color="black", size=4) +
  theme_bw() + coord_flip() + ylab("% Novel BGC") + ylim(0,100)+ 
  scale_fill_npg() +
  theme(axis.title.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, color = "black"),
        legend.position = "none")

df = merge(antisMash_Res, OTU_tax, by = "OTU", all.x = T)
df$OTU_annotated = "Novel OTU"
df$OTU_annotated[df$Species != "s__"] = "Annotated OTU"
df = ddply(df, .variables = c("OTU_annotated","cluster_type"), .fun = function(x){
  return(nrow(x))
})
df = df[order(df$OTU_annotated, -df$V1),]
ind = ddply(df,.variables = "cluster_type", .fun = function(x){sum(x$V1)})


ind = ind[order(ind$V1, decreasing = T),]
df = df[df$cluster_type %in% ind$cluster_type[1:15],]
df$cluster_type = factor(df$cluster_type, levels = ind$cluster_type[1:15])
df$OTU_annotated = factor(df$OTU_annotated, levels = c("Annotated OTU", "Novel OTU"))
df$V1[df$OTU_annotated=="Novel OTU"] = -df$V1[df$OTU_annotated=="Novel OTU"]
g5_E = ggplot(data = df, aes(x =cluster_type , y = V1, fill = OTU_annotated)) +
  geom_bar(stat = "identity",position = "identity",color="black",size=0.25) +
  scale_y_continuous(labels = abs, limits = c(-700, 400), breaks = seq(-700, 400, 100)) +
  ylab("Number of BGCs") + coord_flip() + scale_fill_nejm() + theme_classic()+ 
  theme(plot.margin = unit(c(.5, 1, .5, 1), "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        legend.title = element_blank(),
        legend.text =  element_text(size = 14,color = "black"),
        legend.position = c(0.2,0.9))

antisMash_Res  = antisMash_Res[order(antisMash_Res$OTU, antisMash_Res$cluster, antisMash_Res$MIBIG_hits),]
antisMash_Res = antisMash_Res[!duplicated(antisMash_Res[,c("cluster","OTU")]),]

antisMash_Res = merge(antisMash_Res[,c("cluster_type","MIBIG_hits", "OTU")], 
                      OTU_Abundance[,c("uuid","relative_abundance","prevalence","OTU")])
antisMash_Res = merge(Sample_Info[!duplicated(Sample_Info$uuid),c("uuid","City_Name","Continent_Name")], antisMash_Res)

## uuid
antisMash_df = ddply(antisMash_Res, .variables = c("Continent_Name","City_Name"), .fun = function(df){
  N = length(unique(df$uuid))
  df = ddply(df, .variables = "cluster_type",.fun = function(x){
    length(unique(x$uuid))
  })
  df[,2] = df[,2]/N
  df$sampleNum = N
  return(df)
})
antisMash_df = antisMash_df[antisMash_df$sampleNum >= 20,]
library(ggplot2)
all_cities = unique(antisMash_df$City_Name)
all_type = unique(antisMash_df$cluster_type)
df3 = data.frame(matrix(0, nrow = length(all_type), ncol = length(all_cities)))
rownames(df3) = all_type
colnames(df3) = all_cities
for(i in 1:length(all_type)){
  bgc = all_type[i]
  x = antisMash_df[antisMash_df$cluster_type == bgc,]
  x = ddply(x, .variables = "City_Name", .fun = function(y){
    mean(y$V1)
  })
  df3[i,x$City_Name] = x$V1
}
ind = rowMeans(df3)
ind = ind[order(ind, decreasing = T)]
df3 = df3[names(ind)[1:15],]

g5_F = as_ggplot(pheatmap(df3, fontsize = 14, silent = T)[[4]])

fig_wd = 8
pdf("figures/Fig5_A.pdf", width = fig_wd, height = 0.6*fig_wd)
g5_A
dev.off()
pdf("figures/Fig5_B.pdf", width = fig_wd, height = 0.6*fig_wd)
g5_B
dev.off()
pdf("figures/Fig5_C.pdf", width = fig_wd, height = 0.6*fig_wd)
g5_C
dev.off()
pdf("figures/Fig5_D.pdf", width = fig_wd, height = 0.6*fig_wd)
g5_D
dev.off()
pdf("figures/Fig5_E.pdf", width = fig_wd, height = 0.6*fig_wd)
g5_E
dev.off()
pdf("figures/Fig5_F.pdf", width = fig_wd, height = 0.6*fig_wd)
g5_F
dev.off()


