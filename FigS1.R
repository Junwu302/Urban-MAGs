library(gridExtra)
library(ggplot2)
library(ggsci)

load("data/OTU_tax.RData")
df = OTU_tax[OTU_tax$Kingdom == "d__Bacteria",-c(1,2)]
ranks = colnames(df)[1:6]
df = data.frame(do.call(rbind, lapply(ranks, function(x,df){
  x = df[,x]
  x = grep("__$",x, invert = T, value = T)
  x = as.data.frame(table(x),stringsAsFactors = F)
  x = x[order(x$Freq, decreasing = T),]
  x[6,1] = "Other Classified taxa"
  x[6,2] = sum(x$Freq[-c(1:5)])
  x$x = gsub("[a-z]__","",x$x)
  return(x[1:6,])
},df)),stringsAsFactors = F)
df$Rank = factor(rep(ranks, each=6),levels = ranks)
colnames(df) = c("Taxon","Counts","Rank")
df$Freq = 100*df$Counts/1302
df = df[order(df$Rank, -df$Freq),]
df$Taxon = factor(df$Taxon, levels = c(df$Taxon[df$Taxon != "Other Classified taxa"],
                                       "Other Classified taxa"))
cols = c(rep(rev(brewer.pal(5,"Set1")),6),"#A9A9A9")
pdf("figures/Fig_S1.pdf",width = 12, height = 7)
ggplot(data = df, mapping = aes(x=Rank,y=Freq)) +
  geom_bar(mapping = aes(fill=Taxon), stat="identity") +
  scale_fill_manual(values=cols) + 
  theme_classic() + ylab("Proportion (%)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16, color = "black"),
        axis.text = element_text(size=14, color = "black"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16))
dev.off()