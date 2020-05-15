library(plyr)
library(reshape2)
load("data/OTU_tax.RData")
load("data/OTU_Abundance.RData")
load("data/Sample_Info.RData")
df = merge(OTU_Abundance[,c("uuid","OTU","relative_abundance")], Sample_Info[,c("uuid","Continent_Name")], by="uuid")
df = ddply(df, .variables = c("OTU","Continent_Name"),.fun = function(x){
  x = ddply(x, .variables = "uuid", .fun = function(a){
    sum(a$relative_abundance)
  })
  mean(x$V1)
})
colnames(df)[3] = "abundance"
mat = data.frame(matrix(0, nrow = length(unique(df$OTU)), ncol=length(unique(df$Continent_Name))))
rownames(mat) = unique(df$OTU)
colnames(mat) = unique(df$Continent_Name)
for(conti in unique(df$Continent_Name)){
  mat[df$OTU[df$Continent_Name==conti],conti] = df$abundance[df$Continent_Name==conti]
}
mat$OTU = rownames(mat)
df = merge(OTU_tax, mat, by="OTU")
df = df[df$Kingdom == "d__Bacteria",-c(1,2)]
df$Phylum = gsub("p__","",df$Phylum)
df$Class = gsub("c__","",df$Class)
df$Order = gsub("o__","",df$Order)
df$Family = gsub("f__","",df$Family)
df$Family[df$Family == ""] = paste("Fnovel",1:sum(df$Family == ""),sep='')
df$Genus = gsub("g__","",df$Genus)
df$Genus[df$Genus == ""] = paste("Gnovel",1:sum(df$Genus == ""),sep='')
df$Species = gsub("s__","",df$Species)
df$Species[df$Species == ""] = paste("Snovel",1:sum(df$Species == ""),sep='')
### "OTU_326"
guide = apply(df[,1:6], 1, function(x){
  return(paste0(x, collapse = "."))
})
df$guide = guide
write.table(guide, file="guide.txt",row.names = F, col.names = F, 
            sep='\t',quote = F)

# node color
annot_1 = data.frame(tax=df$Family[grep("novel",df$Family)],
                     term = "clade_marker_color",col="r",stringsAsFactors = F)
annot_1 = rbind(annot_1, data.frame(tax=df$Genus[grep("novel",df$Genus)],
                     term = "clade_marker_color",col="r",stringsAsFactors = F))
annot_1 = rbind(annot_1, data.frame(tax=df$Species[grep("novel",df$Species)],
                                    term = "clade_marker_color",col="r",stringsAsFactors = F))
write.table(annot_1, file="annot_1.txt",row.names = F, col.names = F, sep='\t',quote = F)

# background
phylums = table(df$Phylum)
phylums = phylums[order(phylums, decreasing = T)]
annot_2 = data.frame(tax=names(phylums[1:4]), term = "annotation_background_color",
                     value=c("#03A8F4","#95FE71","#BF3A3A","#E58D30"),stringsAsFactors = F)
annot_2 = rbind(annot_2, data.frame(tax=names(phylums[1:4]), term="annotation",value=names(phylums[1:4]),
                                    stringsAsFactors = F))
write.table(annot_2, file="annot_2.txt",row.names = F, col.names = F, sep='\t',quote = F)

## Ring
## ring_height/ring_separator_color/ring_internal_separator_thickness
df = df[,c(1:6,10,11,13,9,8,12,7,14:ncol(df))]
annot_3 = data.frame(term="ring_edge_color",rank = 1:7,value='None',stringsAsFactors = F)
annot_3 = rbind(annot_3, data.frame(term="ring_internal_separator_thickness",rank = 1:7,value=0.5,stringsAsFactors = F))
annot_3 = rbind(annot_3, data.frame(term="ring_label",rank = 1:7,value=colnames(df)[7:13],stringsAsFactors = F))
annot_3 = rbind(annot_3, data.frame(term="ring_label_font_size",rank = 1:7,value=10,stringsAsFactors = F))
annot_3 = rbind(annot_3,data.frame(term="ring_height",rank = 1:7,value=1.5,stringsAsFactors = F))

annot_4 = data.frame()
cols = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628")
for(i in 1:7){
  x = df[df[,i+6]>0,c(i+6,14)]
  annot_4 = rbind(annot_4, data.frame(tax=x[,2],term="ring_alpha",rank=i,value=sqrt(sqrt(x[,1])),stringsAsFactors = F))
  annot_4 = rbind(annot_4, data.frame(tax=x[,2],term="ring_color",rank=i, value=toupper(cols[i]),stringsAsFactors = F))
}
write.table(annot_3, file="annot_3.txt",row.names = F, col.names = F, sep='\t',quote = F)
write.table(annot_4, file="annot_3.txt",append = T,row.names = F, col.names = F, sep='\t',quote = F)



