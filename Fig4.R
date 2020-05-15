## Heatmap TOM
library(WGCNA)
library(pheatmap)
library(ggplot2)
library(ggimage)
library(ggpubr)
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

g4_A = as_ggplot(pheatmap(OTU_TOM, annotation_col = annotation_col,annotation_legend=F,
                          annotation_names_col = T,show_rownames=F, show_colnames = F,
                          clustering_method="ward.D", fontsize = 14, silent = T)[[4]])+ 
  theme(legend.position = c(0.9,0.5),
        plot.margin = unit(c(.5, 1, .5, 1), "cm"))
pdf("figures/Fig4_A.pdf", height = 8, width = 8)
g4_A
dev.off()

## Annotation performance at Genues level
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
genus_reannote$Frac = genus_reannote$Number/length(unannot_genus)

g4_B = ggplot(data = genus_reannote, mapping = aes(x=Frac,y=Accuracy)) +
  geom_point(size = 1) + geom_smooth(method = "loess",span = 0.5) + 
  ylab("Estimated accuracy") + xlab("Faction of reannotaed") + 
  geom_vline(xintercept=0.39,colour='black',linetype=2,size=1)+
  geom_hline(yintercept=0.8,colour='black',linetype=2,size=1)+
  annotate("text", x=0.5, y=0.81, label= "Faction = 39%",colour="red",size = 5) +
  theme_classic() +
  theme(axis.title = element_text(size=16, color = "black"),
        axis.text = element_text(size=14, color = "black"))
pdf("figures/Fig4_B.pdf", height = 4, width = 8)
g4_B
dev.off()
# network visualization
library(igraph)
library(RColorBrewer)
cut_heigth = genus_reannote$th[1]
selected_mods = cutreeDynamicTree(OTU_Tree, maxTreeHeight = cut_heigth, minModuleSize=2)
selected_mods = data.frame(OTU=OTU_Tree$labels, Mod = selected_mods, stringsAsFactors = F)
ind = selected_mods$Mod[selected_mods$OTU %in% unannot_genus & selected_mods$Mod >0]
selected_mods = selected_mods[selected_mods$Mod %in% ind,]
df1 = merge(OTU_tax[,c("OTU","Genus")],selected_mods)
igraph1 = graph_from_adjacency_matrix(OTU_TOM[df1$OTU,df1$OTU],
                                      mode="undirected",weighted=TRUE,diag=FALSE)
V(igraph1)$label <- V(igraph1)$name
V(igraph1)$color = brewer.pal(length(unique(df1$Genus)), "Set1")[as.numeric(factor(df1$Genus))]
node.col1 = brewer.pal(length(unique(df1$Genus)), "Set1")
names(node.col1) = levels(factor(df1$Genus))
E(igraph1)$color = "darkgray"
names(node.col1) = c("Uannotated","Sphingomonas","UBA1936")
node.col1 = node.col1[c(2,3,1)]

cut_heigth = 0.851
selected_mods = cutreeDynamicTree(OTU_Tree, maxTreeHeight = cut_heigth, minModuleSize=2)
selected_mods = data.frame(OTU=OTU_Tree$labels, Mod = selected_mods, stringsAsFactors = F)
ind = selected_mods$Mod[selected_mods$OTU %in% unannot_genus & selected_mods$Mod >0]
selected_mods = selected_mods[selected_mods$Mod %in% ind,]
df2 = merge(OTU_tax[,c("OTU","Genus")],selected_mods)

selected_tom = OTU_TOM[df2$OTU,df2$OTU]
for(i in 1:nrow(selected_tom)){
  selected_tom[i,df2$Mod != df2$Mod[i]] = 0
}
igraph2 = graph_from_adjacency_matrix(selected_tom,mode="undirected",weighted=TRUE,diag=FALSE)
V(igraph2)$label <- V(igraph2)$name
node.col2 = rep(brewer.pal(3, "Set1")[2],nrow(df2))
node.col2[df2$Genus=="g__"] = brewer.pal(3, "Set1")[1]
E(igraph2)$color = "darkgray"
V(igraph2)$color = node.col2
node.col2 = brewer.pal(3, "Set1")[c(2,1)]
names(node.col2) = c("Annotated","Unannotated")


pdf("figures/Fig4_C.pdf", width = 4, height = 4)
plot(igraph2,vertex.label=NA,edge.width=1,vertex.size=5, edge.curved=TRUE)
legend("bottomleft", legend=names(node.col2), col = node.col2 , bty = "n", pch=20 ,
       pt.cex = 3, cex = 1, text.col=node.col2, horiz = FALSE)
dev.off()

pdf("figures/Fig4_D.pdf", width = 4, height = 4)
plot(igraph1,vertex.label=NA,edge.width=1,vertex.size=8, edge.curved=FALSE)
legend("bottomleft", legend=names(node.col1), col = node.col1 , bty = "n", pch=20 ,
       pt.cex = 3, cex = 1, text.col=node.col1, horiz = FALSE)
dev.off()

