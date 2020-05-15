load("data/COG_conti_adonis2.RData")
city_name = read.csv("data/city_name.csv",stringsAsFactors = F)
df = COG_conti_adonis2[,c(2,1,3:7)]
colnames(df)[1:2] = c("Var1","Var2")
df = rbind(COG_conti_adonis2, df)
df = merge(df, city_name[,c(3,4)],by.x = "Var1", by.y = "continent")
df = merge(df, city_name[,c(3,4)],by.x = "Var2", by.y = "continent")
df$Var1 = df$Continent_Name.x
df$Var2 = df$Continent_Name.y
all_conti = unique(c(df$Var1, df$Var2))
dismat = matrix(0, nrow = length(all_conti), ncol = length(all_conti))
pmat = matrix(1, nrow = length(all_conti), ncol = length(all_conti))
colnames(dismat) = rownames(dismat) = rownames(pmat) = colnames(pmat) = all_conti
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
df$Var1 = factor(df$Var1, levels = all_conti)
df$Var2 = factor(df$Var2, levels = all_conti)
colnames(df)[3] = "R2"

gS4_A = ggplot(data = df, aes(Var1, Var2, fill = R2))+
  geom_tile(color = "white",na.rm = TRUE)+
  scale_fill_distiller(palette="Spectral") +
  geom_text(aes(Var1, Var2, label = label), size = 5) +
  ylab("") + xlab("") + theme_minimal()+ 
  theme(axis.title = element_text(size=14, color = "black"),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black")) + coord_fixed()

library(ggplot2)
library(ggpubr)
load("data/KEGG_city_adonis2.RData")
load("data/KEGG_conti_adonis2.RData")
city_name = read.csv("data/city_name.csv",stringsAsFactors = F)
KEGG_city_adonis2 = merge(merge(KEGG_city_adonis2, city_name[,c(1,2)],by.x = "Var1", by.y="city"),
                          city_name[,c(1,2)], by.x="Var2", by.y="city")
KEGG_city_adonis2$Var1 = KEGG_city_adonis2$City_Name.x
KEGG_city_adonis2$Var2 = KEGG_city_adonis2$City_Name.y
all_cities = unique(c(KEGG_city_adonis2$Var1, KEGG_city_adonis2$Var2))
all_cities = all_cities[order(all_cities)]
dismat = matrix(0, nrow=length(all_cities), ncol=length(all_cities))
pmat = matrix(1, nrow=length(all_cities), ncol=length(all_cities))
colnames(dismat) = rownames(dismat) = rownames(pmat) = colnames(pmat)=all_cities
for(i in 1:nrow(KEGG_city_adonis2)){
  dismat[KEGG_city_adonis2$Var1[i], KEGG_city_adonis2$Var2[i]] = KEGG_city_adonis2$R2[i]
  dismat[KEGG_city_adonis2$Var2[i], KEGG_city_adonis2$Var1[i]] = KEGG_city_adonis2$R2[i]
  pmat[KEGG_city_adonis2$Var1[i], KEGG_city_adonis2$Var2[i]] = KEGG_city_adonis2$pval[i]
  pmat[KEGG_city_adonis2$Var2[i], KEGG_city_adonis2$Var1[i]] = KEGG_city_adonis2$pval[i]
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
gS4_B = ggplot(data = df, aes(Var1, Var2, fill = R2))+
  geom_tile(color = "white",na.rm = TRUE)+
  scale_fill_distiller(palette="Spectral") +
  geom_text(aes(Var1, Var2, label = label), size = 5) +
  ylab("") + xlab("") + theme_minimal()+ 
  theme(axis.title = element_text(size=14, color = "black"),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black")) + coord_fixed()


KEGG_conti_adonis2 = merge(merge(KEGG_conti_adonis2, city_name[!duplicated(city_name$continent),c(3,4)],by.x = "Var1", by.y="continent"),
                           city_name[!duplicated(city_name$continent),c(3,4)], by.x="Var2", by.y="continent")
KEGG_conti_adonis2$Var1 = KEGG_conti_adonis2$Continent_Name.x
KEGG_conti_adonis2$Var2 = KEGG_conti_adonis2$Continent_Name.y
all_conti = unique(c(KEGG_conti_adonis2$Var1, KEGG_conti_adonis2$Var2))
all_conti = all_conti[order(all_conti)]
dismat = matrix(0, nrow=length(all_conti), ncol=length(all_conti))
pmat = matrix(1, nrow=length(all_conti), ncol=length(all_conti))
colnames(dismat) = rownames(dismat) = rownames(pmat) = colnames(pmat)=all_conti
for(i in 1:nrow(KEGG_conti_adonis2)){
  dismat[KEGG_conti_adonis2$Var1[i], KEGG_conti_adonis2$Var2[i]] = KEGG_conti_adonis2$R2[i]
  dismat[KEGG_conti_adonis2$Var2[i], KEGG_conti_adonis2$Var1[i]] = KEGG_conti_adonis2$R2[i]
  pmat[KEGG_conti_adonis2$Var1[i], KEGG_conti_adonis2$Var2[i]] = KEGG_conti_adonis2$pval[i]
  pmat[KEGG_conti_adonis2$Var2[i], KEGG_conti_adonis2$Var1[i]] = KEGG_conti_adonis2$pval[i]
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
df$Var1 = factor(df$Var1, levels = all_conti)
df$Var2 = factor(df$Var2, levels = all_conti)
colnames(df)[3] = "R2"
gS4_C = ggplot(data = df, aes(Var1, Var2, fill = R2))+
  geom_tile(color = "white",na.rm = TRUE)+
  scale_fill_distiller(palette="Spectral") +
  geom_text(aes(Var1, Var2, label = label), size = 5) +
  ylab("") + xlab("") + theme_minimal()+ 
  theme(axis.title = element_text(size=14, color = "black"),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black")) + coord_fixed()


pdf("figures/FigS4_A.pdf",width = 6, height = 6)
gS4_A
dev.off()

pdf("figures/FigS4_B.pdf",width = 6, height = 6)
gS4_B
dev.off()

pdf("figures/FigS4_C.pdf",width = 6, height = 6)
gS4_C
dev.off()   

