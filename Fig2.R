###== diversity improvement
library(phyloseq)
library(ape)
library(reshape2)
library(ggplot2)
library(ggsci)


load("data/OTU_tax.RData")
tree = read.tree("data/gtdbtk.bac120.user_msa.tree")
OTU_annoted = OTU_tax[OTU_tax$Species != "s__",1:7]
OTU_novel = OTU_tax[OTU_tax$Species == "s__",1:7]
tax = rbind(OTU_annoted, OTU_novel)
OTU_tax = OTU_tax[OTU_tax$OTU %in% tree$tip.label,]
phyla = unique(OTU_tax$Phylum)
pd = data.frame(matrix(0, ncol=5, nrow=length(phyla)))
colnames(pd) = c("Phylum", "Annotated_OTU_PD", "All_OTUs_PD","Annoted_OTU_Num","Novel_OTU_Num")
for (i in 1:length(phyla)){
  phylum = as.character(phyla[i])
  annot.otu = OTU_annoted$OTU[OTU_annoted$Phylum == phylum]
  novel.otu = OTU_novel$OTU[OTU_novel$Phylum == phylum]
  all.out = c(annot.otu, novel.otu)
  
  if (length(annot.otu) > 0){
    subs.tre.annot = drop.tip(tree, tree$tip.label[-match(annot.otu, tree$tip.label)])
    pd[i,2] = sum(subs.tre.annot$edge.length)
  } 
  subs.tre.all = drop.tip(tree, tree$tip.label[-match(all.out, tree$tip.label)])
  pd[i,1] = phylum
  pd[i,3] = sum(subs.tre.all$edge.length)
  pd[i,4] = length(annot.otu)
  pd[i,5] = length(novel.otu)
}
pd$Proportion = (pd$All_OTUs_PD-pd$Annotated_OTU_PD)/pd$All_OTUs_PD *100
pd$Improvement = pd$All_OTUs_PD - pd$Annotated_OTU_PD

# reorder phylum
pd.fi = pd[order(pd$Improvement,decreasing=TRUE),]
pd.fi$Phylum = paste(gsub("^p__","",pd.fi$Phylum),
                     paste(paste("(",paste(pd.fi$Annoted_OTU_Num, 
                                           pd.fi$Novel_OTU_Num,sep="/"),sep = ''),")",sep=''))
pd.fi = pd.fi[1:10,]
pd.fi$Phylum = factor(pd.fi$Phylum, levels = pd.fi$Phylum)
df = melt(pd.fi[,c(1,6,7)])

g2_B = ggplot(data = df, mapping = aes(x=Phylum, y = value))+
  geom_bar(stat = "identity",mapping = aes(fill=Phylum)) +
  facet_wrap(variable~., scales = "free",nrow=2) +
  coord_flip()+ scale_fill_npg() + theme_bw()+
  theme(axis.title = element_blank(),
        strip.text = element_text(size=16, color = "black"),
        axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_blank(),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        legend.position = "right")

library(plyr)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(ggsci)
library(ggplot2)
library(phyloseq)
library(ape)
library(vegan)
library(picante)

load("data/OTU_tax.RData")
load("data/OTU_Abundance.RData")
load("data/Sample_Info.RData")

Novel_OTUs = merge(OTU_Abundance[,c("uuid","OTU")], OTU_tax[,c("OTU","Species")])
Novel_OTUs = merge(Novel_OTUs, Sample_Info)
df = ddply(Novel_OTUs, .variables = "City_Name",.fun = function(x){
  latitude = mean(x$city_latitude)
  longitude  = mean(x$city_longitude)
  TotalCount = length(unique(x$OTU))
  NovelCount = length(unique(x$OTU[x$Species == "s__"]))
  NovelRatio = NovelCount/TotalCount
  SampleNum = length(unique(x$uuid))
  Pop_density = mean(as.numeric(x$city_population_density))
  return(c(latitude, longitude, TotalCount, NovelCount, NovelRatio, SampleNum ,Pop_density))
})[-1,]
colnames(df) = c("city","latitude", "longitude", "TotalCount","NovelCount", "NovelRatio", "SampleNum", "Pop_density")
df = df[!is.nan(df$Pop_density),]
df = df[order(df$NovelCount, decreasing = T),]
df1 = df[1:15,]
df1$city = factor(df1$city, levels = df1$city)
df1 = data.frame(city = rep(df1$city,2), Count=c(df1$NovelCount, df1$TotalCount - df1$NovelCount),
                 labels = rep(c("Novel OTU","Annotated OTU"),each=nrow(df1)))
df1$labels = factor(df1$labels, levels = c("Annotated OTU","Novel OTU"))

g2_C = ggplot(data = df1, mapping = aes(x=city,y=Count, fill=labels)) +
  geom_bar(stat="identity",width = .9, position=position_dodge()) + ylim(0, 300) +
  ylab("Novel OTU Number") + xlab("") + scale_fill_npg() +theme_classic() + 
  theme(axis.title = element_text(size=14, color = "black"),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        legend.text = element_text(size=16, color = "black"),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.9))


g2_D = ggplot(data = df, mapping = aes(x=SampleNum, y = NovelCount)) +
  geom_point(size = 2) + geom_smooth(method = lm, fill="#ffc0cb") +
  stat_cor(method = "pearson", size=5, color='red') +
  scale_x_log10() + scale_y_log10()+ 
  xlab("Sample Number") + ylab("Novel OTU Number") + theme_classic() +
  theme(axis.title = element_text(size=16, color = "black"),
        axis.text = element_text(size=12, color = "black"))



load("data/OTU_Abundance.RData")
load("data/Sample_Info.RData")
load("data/OTU_tax.RData")

df = merge(Sample_Info[,c("uuid","City_Name","Continent_Name")], OTU_Abundance)
conti_otus = ddply(df[df$prevalence==1,], .variables = c("City_Name","Continent_Name"), .fun = function(x){
  return(length(unique(x$OTU)))
})
df = df[df$Continent_Name %in% conti_otus$Continent_Name[conti_otus$V1>=20],]
all_otus = unique(df$OTU)
top5_OTUs = data.frame(ddply(df, .variables = "Continent_Name", .fun = function(x, all_otus){
  res = rep("None",length(all_otus))
  names(res) = all_otus
  x = ddply(x,.variables = "OTU",.fun = function(a){
    sum(a$prevalence)
  })
  x = x[order(x$V1, decreasing = T),]
  res[x$OTU] = "Observed"
  res[x$OTU[1:5]] = "Top 5"
  return(res)
}, all_otus),stringsAsFactors = F)
top5_OTUs = top5_OTUs[,c(1, which(apply(top5_OTUs, 2, function(x){
  sum(x=="Top 5")>0
})))]
top5_OTUs = melt(top5_OTUs, id.vars = 1)
colnames(top5_OTUs) = c("continent","OTU","Prevalance")
top5_OTUs$Prevalance = factor(top5_OTUs$Prevalance,level=c("Top 5","Observed","None"))
df = merge(top5_OTUs, OTU_tax)
df$Species = gsub("^[a-z]__","",df$Species)
df$Species[df$Species==""] = "Novel"
df$OTU = paste(gsub("_","",df$OTU),paste(paste("(",df$Species,sep=""),")",sep=""))  
id = ddply(df, .variables = "OTU", .fun = function(x){
  sum(x$Prevalance=="Top 5")
})
id = id[order(id$V1),]
df$OTU = factor(df$OTU,level=id$OTU)

g2_E = ggplot(df, aes(x=continent, y = OTU, fill=Prevalance)) + 
  geom_tile(color="white", size=0.1) +
  scale_fill_manual(values=c('#e35604','#98c4ff','#e6e6e6'))+labs(x="",y="")+
  scale_y_discrete(position = "right") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.x = element_text(angle = 45,size=8, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=6, color = "black"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "top")



load("data/OTU_Abundance.RData")
load("data/Sample_Info.RData")
load("data/OTU_tax.RData")
tree = read.tree("data/gtdbtk.bac120.user_msa.tree")

OTU_Abundance = merge(Sample_Info, OTU_Abundance)
city_uuid_Num = ddply(OTU_Abundance, .variables = "City_Name", function(x){
  length(unique(x$uuid))
})

OTU_Abundance = OTU_Abundance[OTU_Abundance$City_Name %in% city_uuid_Num$City_Name[city_uuid_Num$V1>=20],]

all_otu = unique(OTU_Abundance$OTU)
all_city = unique(OTU_Abundance$City_Name)

otu_tbl = data.frame(matrix(0, length(all_city), length(all_otu)),stringsAsFactors = F)
rownames(otu_tbl) = all_city
colnames(otu_tbl) = all_otu
for(city in all_city){
  df = OTU_Abundance[OTU_Abundance$City_Name==city,] 
  df = ddply(df, .variables = "OTU",function(x){sum(x$prevalence)})
  otu_tbl[city,df$OTU] = df$V1
}

diversity_index = data.frame(city = all_city,shannon_index = diversity(otu_tbl, index = 'shannon'),
                             simpson_index = diversity(otu_tbl, index = 'simpson'), 
                             Richness = specnumber(otu_tbl),stringsAsFactors = F)
diversity_index = merge(diversity_index, Sample_Info[!duplicated(Sample_Info$City_Name),c("City_Name","Continent_Name")],
                        by.x = "city", by.y="City_Name")

diversity_index = melt(diversity_index, id.vars=c("city","Continent_Name"))
diversity_index = diversity_index[diversity_index$variable == "shannon_index",]
diversity_index = diversity_index[order(diversity_index$Continent_Name, diversity_index$value),]
diversity_index$city = factor(diversity_index$city, levels = diversity_index$city[!duplicated(diversity_index$city)])

g2_F = ggplot(data = diversity_index, mapping = aes(x=city,y=value,fill=Continent_Name)) +
  geom_bar(stat="identity",width = .9) + scale_fill_npg() + 
  theme_bw()+ labs(fill = "Continent") + xlab("") + ylab("Shannon Index") +
  theme(axis.title.y = element_text(size=18, color = "black"),
        strip.text = element_text(size=16, color = "black"),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=16, color = "black"),
        legend.text = element_text(size=16), 
        legend.title = element_text(size=18))

fig_wd = 8

pdf("figures/Fig2_B.pdf", width = fig_wd,height = 0.75*fig_wd)
g2_B
dev.off()
pdf("figures/Fig2_C.pdf", width = fig_wd,height = 0.5*fig_wd)
g2_C
dev.off()
pdf("figures/Fig2_D.pdf", width = fig_wd,height = 0.5*fig_wd)
g2_D
dev.off()
pdf("figures/Fig2_E.pdf", width = fig_wd,height = 0.75*fig_wd)
g2_E
dev.off()
pdf("figures/Fig2_F.pdf", width = fig_wd,height = 0.5*fig_wd)
g2_F
dev.off()



