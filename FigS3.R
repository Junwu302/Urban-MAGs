## Unannotated gene Numbers
load("data/MetaSUB_eggNOG_res.RData")
load("data/OTU_tax.RData")
OTU_geneNumbers = read.table("data/OTU_geneNumbers.tsv",header = F, stringsAsFactors = F)
OTU_geneNumbers = OTU_geneNumbers[,c(2,1)]
colnames(OTU_geneNumbers) = c("OTU","Total_Num")
OTU_geneNumbers$OTU = gsub(".gff$","",OTU_geneNumbers$OTU)
annoted_GeneNum = unlist(lapply(eggNOG_res, function(x){
  return(sum(x$COG_cat !="" & x$COG_cat != "S" & x$COG_cat != "R"))
}))
annoted_GeneNum = data.frame(OTU= names(annoted_GeneNum), Annoted_Num = annoted_GeneNum, stringsAsFactors = F)
OTU_geneNumbers = merge(OTU_geneNumbers, annoted_GeneNum)
OTU_geneNumbers = merge(OTU_geneNumbers, OTU_tax)
OTU_geneNumbers$Annot_ratio = OTU_geneNumbers$Annoted_Num/OTU_geneNumbers$Total_Num
OTU_geneNumbers$Annot_status = "Annotated"
OTU_geneNumbers$Annot_status[OTU_geneNumbers$Species == "s__"] = "Species_Unannotated"
OTU_geneNumbers$Annot_status[OTU_geneNumbers$Genus == "g__"] = "Genus_Unannotated"
OTU_geneNumbers$Annot_status  = factor(OTU_geneNumbers$Annot_status,
                                       levels = c("Genus_Unannotated","Species_Unannotated","Annotated"))
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(gridExtra)
#pdf("Gene_AnnotatedRatio.pdf",width = 10, height = 7)
df = data.frame(th = seq(0,1,0.01), frac = 1)
for(i in 1:nrow(df)){
  n = sum(OTU_geneNumbers$Annot_ratio <= df$th[i])/nrow(OTU_geneNumbers)
  df$frac[i] = n
}
gS3_A = ggplot(data = df, mapping = aes(x=th, y=frac)) +
  geom_line() + geom_vline(xintercept=0.56,colour='black',linetype=2,size=1)+
  geom_hline(yintercept=0.5,colour='black',linetype=2,size=1)+
  annotate("text", x=0.25, y=0.55, label= "y=0.5",colour="red",size = 5) +
  annotate("text", x=0.49, y=0.81, label= "x=0.56",colour="red",size = 5) +
  xlab("proportion of annotated genes") +
  ylab("Cumulative frequency (Proportion)") + theme_bw()+ 
  theme(axis.title = element_text(size=16, color = "black"),
        axis.text = element_text(size=14, color = "black"),
        plot.margin = unit(c(.5, 1, .5, 1), "cm"))

gS3_B = ggplot(data = OTU_geneNumbers, mapping = aes(x=Annot_status, y=Annot_ratio,fill = Annot_status)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_npg() + theme_bw() +
  xlab("Taxonomy annoation status") + ylab("Proportion of annotated gene") +
  theme(axis.title.y = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14, color = "black"),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2,0.9),
        plot.margin = unit(c(.5, 1, .5, 1), "cm"))

##== KEGG taxonomy
library(vegan)
library(umap)
load("data/KEGG_df.RData")
load("data/OTU_tax.RData")
KEGG_tax = merge(OTU_tax, KEGG_df)
# Phylum.anosim = anosim(KEGG_tax[,-c(1:8)], grouping = factor(KEGG_tax$Phylum))
# Class.anosim = anosim(KEGG_tax[,-c(1:8)], grouping = factor(KEGG_tax$Class))
# Order.anosim = anosim(KEGG_tax[,-c(1:8)], grouping = factor(KEGG_tax$Order))
# Family.anosim = anosim(KEGG_tax[KEGG_tax$Family!="f__",-c(1:8)], grouping = factor(KEGG_tax$Family[KEGG_tax$Family!="f__"]))
# Genus.anosim = anosim(KEGG_tax[KEGG_tax$Genus!="g__",-c(1:8)], grouping = factor(KEGG_tax$Genus[KEGG_tax$Genus!="g__"]))
# anosim_kegg_tax = list(Phylum=Phylum.anosim, Class = Class.anosim, 
#                        Order=Order.anosim, Family=Family.anosim,Genus=Genus.anosim)
# save(anosim_kegg_tax, file="anosim_kegg_tax.RData")
# load("data/anosim_kegg_tax.RData")
df = umap(KEGG_tax[,-c(1:8)])
library(ggplot2)
library(ggsci)
df = data.frame(df$layout)
df = cbind(KEGG_tax[,1:6], df)
ind = df$Phylum %in% names(table(KEGG_tax$Phylum))[table(KEGG_tax$Phylum) >= 10]
df = df[ind,c("Phylum","X1","X2")]
df$Phylum = gsub("^p__","",df$Phylum)
gS3_C = ggplot(data = df, mapping = aes(x=X1, y=X2, color=Phylum)) +
  geom_point(size = 1) + xlab("") + ylab("") +
  scale_colour_brewer(palette="Set1") + theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=14, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=16, color = "black"))




pdf("figures/FigS3_A.pdf",width = 6, height = 6)
gS3_A
dev.off()

pdf("figures/FigS3_B.pdf",width = 6, height = 6)
gS3_B
dev.off()

pdf("figures/FigS3_C.pdf",width = 6, height = 6)
gS3_C
dev.off()


