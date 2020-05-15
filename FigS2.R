library(WGCNA)
library(plyr)
library(ape)
library(ggplot2)
library(reshape)
library(ggsci)
library(gridExtra)
library(cowplot)
load("data/OTU_tax.RData")
load("data/Protein_Cooccurrence.RData")
Protein_Cooccurrence = Protein_Cooccurrence[order(Protein_Cooccurrence$query, Protein_Cooccurrence$subject),]
OTU_jaccard = matrix(0, 1304,1304)
rownames(OTU_jaccard) = colnames(OTU_jaccard) = unique(c(Protein_Cooccurrence$query, Protein_Cooccurrence$subject))
OTU_jaccard[lower.tri(OTU_jaccard)] = Protein_Cooccurrence$jaccard_coef
OTU_jaccard = OTU_jaccard + t(OTU_jaccard) + diag(1,1304,1304)
OTU_TOM = TOMsimilarity(OTU_jaccard)
rownames(OTU_TOM) = colnames(OTU_TOM) = rownames(OTU_jaccard)
annotation_col = OTU_tax[,c(3:5)]
rownames(annotation_col) = OTU_tax$OTU
OTU_Tree = hclust(as.dist(1-OTU_TOM), method = "ward.D")
#write.tree(as.phylo.hclust(OT_Tree), tree.names == TRUE, file="OTU_PC_reAnnotation.newick")
cut_list = seq(0.5,0.99,0.001)
Cluster_score = data.frame()
reAnnot_df = data.frame()
unannot_genus = OTU_tax$OTU[OTU_tax$Genus == "g__"]
genus_reannote = data.frame()
for(i in 1:length(cut_list)){
  th = cut_list[i]
  dynamicMods = cutreeDynamicTree(OTU_Tree, maxTreeHeight = th, minModuleSize=2)
  dynamicMods = data.frame(OTU=OTU_Tree$labels, Mod = dynamicMods, stringsAsFactors = F)
  dynamicMods = merge(dynamicMods, OTU_tax)
  dynamicMods = dynamicMods[order(dynamicMods$OTU),]
  reAnnot_df = rbind(reAnnot_df, data.frame(t(dynamicMods$Mod)))
  Purity = apply(dynamicMods[dynamicMods$Mod>0,4:9], 2, function(x, y){
    ind = grep("__$",x, invert = T)
    x = x[ind]
    y = y[ind]
    return(sum(apply(table(x, y), 2, max))/length(y))
  },dynamicMods$Mod[dynamicMods$Mod>0])
  genus_reannote = rbind(genus_reannote, data.frame(Number = sum(dynamicMods$Mod[dynamicMods$OTU %in% unannot_genus]>0),
                                                    th = th, Accuracy = Purity["Genus"]))
  Module_Cov = sum(dynamicMods$Mod>0)/nrow(dynamicMods)
  Cluster_score = rbind(Cluster_score, data.frame(t(Purity),Module_Cov = Module_Cov))
}
Cluster_score$th =  cut_list
Cluster_score = Cluster_score[!is.nan(Cluster_score$Phylum),]
colnames(reAnnot_df) = OTU_tax$OTU
reAnnot_df$th = cut_list
genus_reannote = genus_reannote[genus_reannote$Number >0,]
genus_reannote=genus_reannote[order(genus_reannote$Accuracy, decreasing = T),]


df1 = melt(Cluster_score[,-c(6,7)], id.vars="th")
df2 = Cluster_score[,c(7,8)]
gS2_A = ggplot(data = df1, mapping = aes(x=th,y=value,color=variable)) +
  geom_line(size=1) + ylab("Module Purity ") + xlab("Cut Height") +
  geom_hline(yintercept=0.74,colour='black',linetype=2,size=1)+
  scale_colour_brewer(palette = "Dark2",name="Taxonomy Level") + theme_classic() +
  annotate("text", x=0.6, y=0.76, label= "Module Purity = 0.74",colour="red",size = 5) +
  theme(axis.title = element_text(size=16, color = "black"),
        axis.text = element_text(size=14, color = "black"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        legend.position = c(0.85,0.5))

gS2_B = ggplot(data = df2, mapping = aes(x=th,y=Module_Cov)) +
  geom_line(size=1) + ylab("Fraction of OUTs in modules") + xlab("Cut Height") +
  theme_classic() +
  theme(axis.title = element_text(size=16, color = "black"),
        axis.text = element_text(size=14, color = "black"))



pdf("figures/Fig_S2A.pdf", width = 8, height = 6)
gS2_A
dev.off()
pdf("figures/Fig_S2B.pdf", width = 8, height = 6)
gS2_B
dev.off()
