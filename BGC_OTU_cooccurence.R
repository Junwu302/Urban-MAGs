library(plyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(reshape2)
library(Hmisc)
load("data/antisMash_Res.RData")
#load("city_info.RData")
#load("MIBIG_info.RData")
load("data/OTU_tax.RData")
load("data/Sample_Info.RData")
load("data/OTU_Abundance.RData")

antisMash_Res  = antisMash_Res[order(antisMash_Res$OTU, antisMash_Res$cluster, antisMash_Res$MIBIG_hits),]
antisMash_Res = antisMash_Res[!duplicated(antisMash_Res[,c("cluster","OTU")]),]

antisMash_Res = merge(antisMash_Res[,c("cluster_type", "OTU")], 
                      OTU_Abundance[,c("uuid","OTU")])
antisMash_Res = merge(Sample_Info[!duplicated(Sample_Info$uuid),c("uuid","City_Name","Continent_Name")], antisMash_Res)
antisMash_Res = merge(antisMash_Res, OTU_tax[,c("OTU","Genus")], by="OTU")
antisMash_Res$Genus[antisMash_Res$Genus == "g__"] = antisMash_Res$OTU[antisMash_Res$Genus == "g__"]
antisMash_Res$Genus = gsub("^g__","",antisMash_Res$Genus)
all_genus = unique(antisMash_Res$Genus)
all_bgc = unique(antisMash_Res$cluster_type[antisMash_Res$cluster_type != "Other"])
all_uuid = unique(antisMash_Res$uuid)
res = data.frame()

for(genus in all_genus){
  for(bgc in all_bgc){
    tmp = data.frame(bgc=rep(0, length(all_uuid)), genus=rep(0, length(all_uuid)))
    rownames(tmp) = all_uuid
    tmp[unique(antisMash_Res$uuid[antisMash_Res$Genus == genus]),"bgc"] = 1
    tmp[unique(antisMash_Res$uuid[antisMash_Res$cluster_type == bgc]),"genus"] = 1
    a = sum(tmp$genus == 1)
    b = sum(tmp$bgc == 1)
    n = nrow(tmp)
    x = sum(tmp$genus == 1 & tmp$bgc == 1)
    jaccard_coef = a/(n-a-b+x)
    pval = phyper(min(a,b), a,n-a,b) - phyper(x-1, a,n-a,b)
    res = rbind(res, data.frame(geuns=genus, bgc=bgc,jaccard_coef=jaccard_coef, 
                                pval = pval,stringsAsFactors = F))
  }
}

res = res[order(res$pval, -res$jaccard_coef),]
res$qval = p.adjust(res$pval, method = "BH")
df = res[res$qval < 0.01,]

net_df = data.frame(source=df$bgc, target=df$geuns, weigth=df$jaccard_coef)
attr_df = rbind(data.frame(node=unique(df$geuns),type="genus",stringsAsFactors = F),
                data.frame(node=unique(df$bgc),type="bgc",stringsAsFactors = F))
attr_df = merge(attr_df, rbind(data.frame(table(df$bgc)),data.frame(Var1=unique(df$geuns),Freq=10)),by.x="node",by.y="Var1")

OTU_tax = OTU_tax[,c("OTU","Genus","Phylum")]
OTU_tax$Genus[OTU_tax$Genus == "g__"] = OTU_tax$OTU[OTU_tax$Genus == "g__"]
OTU_tax$Genus = gsub("^g__","",OTU_tax$Genus)
OTU_tax$Phylum = gsub("p__","",OTU_tax$Phylum)
attr_df = merge(attr_df, OTU_tax[!duplicated(OTU_tax$Genus),c("Genus","Phylum")],
                by.x = "node",by.y = "Genus",all.x = T)
attr_df$Phylum[is.na(attr_df$Phylum)] = "BGC"
attr_df$label = attr_df$node
attr_df$label[attr_df$type == "genus"] = ""
write.csv(net_df, file="data/bgc_otu_net.csv",row.names = F, quote = F)
write.csv(attr_df, file="data/bgc_otu_attr.csv",row.names = F, quote = F)


top_bgc = attr_df[attr_df$type=="bgc",]
top_bgc = top_bgc[order(top_bgc$Freq, decreasing = T),]

for(bgc in top_bgc$node[1:10]){
  sub_net = net_df[net_df$source == bgc,]
  sub_attr = rbind(data.frame(node=bgc,type="bgc",size=max(sub_net$weigth),stringsAsFactors = F),
                             data.frame(node=sub_net$target,type="genus",size=sub_net$weigth,stringsAsFactors = F))
  net_file = paste0(c("data/",bgc,"_bgc_otu_net.csv"),collapse = "")
  attr_file = paste0(c("data/",bgc,"_bgc_otu_attr.csv"),collapse = "")
  write.csv(sub_net, file=net_file,row.names = F, quote = F)
  write.csv(sub_attr, file=attr_file,row.names = F, quote = F)
}




