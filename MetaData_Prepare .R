library(plyr)
## merge the contigs, reads and metadata
contigs_files = read.table("all_contigs.tsv", sep='\t',stringsAsFactors = F)[,1]
#contigs_files = list.files("assemblies",pattern=".fasta$",recursive=T, full.names=T)
contigs_id = basename(gsub(".metaspades/contigs.fasta","",contigs_files))
contigs = data.frame(contigs_id = contigs_id, contigs_files = contigs_files, stringsAsFactors = F)

reads_files = read.table("all_reads.tsv", sep='\t',stringsAsFactors = F)[,1]

## 1. merge contigs and reads
Sample_Info = data.frame()
for(i in 1:length(contigs_files)){
  uuid = basename(gsub(".metaspades/contigs.fasta","",contigs_files[i]))
  reads = grep(uuid, reads_files, value = T)
  n = length(reads)
  if(n > 0){
    df = data.frame(uuid = rep(uuid,n), contigs = rep(contigs_files[i],n), 
                    reads = reads,stringsAsFactors = F)
    Sample_Info = rbind(Sample_Info, df)
  }
}
## 2. merged with metadata
meta = read.csv("./MetaSUB-metadata-master/complete_metadata.csv", header = T, stringsAsFactors = F)
meta$city_latitude[is.nan(meta$city_latitude)] = meta$latitude[is.nan(meta$city_latitude)]
meta$city_longitude[is.nan(meta$city_longitude)] = meta$longitude[is.nan(meta$city_longitude)]
meta$city_latitude[is.na(meta$city_latitude)] = meta$latitude[is.na(meta$city_latitude)]
meta$city_longitude[is.na(meta$city_longitude)] = meta$longitude[is.na(meta$city_longitude)]

meta$latitude = meta$city_latitude
meta$longitude = meta$city_longitude

Sample_Info = merge(Sample_Info, meta, by="uuid")

## Add the miss data manually
reads_unmapped = reads_files[!(reads_files %in% Sample_Info$reads)]
contigs_unmapped = contigs_files[!(contigs_files %in% Sample_Info$contigs)]
meta_unmapped = meta[!(meta$uuid %in% Sample_Info$uuid),]
#uuids3 = grep("5106-CEM-0507", contigs_unmapped, value = T)
#uuids4 = grep("5106-CEM-0507", contigs_unmapped, value = T)
contigs_unmapped = grep("Shanghai", contigs_unmapped, value = T)
contigs_unmapped = data.frame(uuid = paste(paste("pilot_Shanghai-China_MS-0",c("05","08","10"),sep=''),
                                            "-Shanghai",sep=''),contigs =contigs_unmapped,stringsAsFactors = F)

reads_unmapped = grep("Shanghai", reads_unmapped, value = T)
reads_unmapped = data.frame(uuid = paste(paste("pilot_Shanghai-China_MS-0",c("02","08","02","05","10"),sep=''),
                                                "-Shanghai",sep=''),reads=reads_unmapped, stringsAsFactors=F)
reads_unmapped = reads_unmapped[!duplicated(reads_unmapped$uuid),]

df = merge(contigs_unmapped, reads_unmapped)
df = merge(df, meta_unmapped, by="uuid")
Sample_Info = rbind(Sample_Info, df)

### City information correct
Sample_Info = Sample_Info[!grepl("(control|other)", Sample_Info$city),]
Sample_Info$city = gsub("\\s","_",Sample_Info$city)

city_name = read.csv("city_name.csv",stringsAsFactors = F)
Sample_Info = merge(Sample_Info, city_name[,c("city","City_Name","Continent_Name")], by.x= "city",by.y = "city")
save(Sample_Info, file="MetaSUB_Sample_Info.RData")



# city_meta = read.csv("./MetaSUB-metadata-master/city_metadata.csv", header = T, stringsAsFactors = F)
# city_involved = data.frame(city=unique(Sample_Info$city),stringsAsFactors = F)
# city_involved = merge(city_involved, city_meta, all.x=T)
# city_involved = city_involved[,c("city","latitude","longitude","coastal", "total_population",
#                                  "population_density_km2","continent")]
# write.csv(city_involved, file="city_involved.csv",row.names = F, quote = F)
# 

city_info = Sample_Info[!duplicated(Sample_Info$uuid),c("city","uuid","latitude","longitude",
                                                        "country","city_total_population","city_population_density","continent")]
save(city_info, file="city_info.RData")
## plot the geographic distr","ibution maps
library(plyr)
library(maptools)
library(ggplot2)
library(RColorBrewer)

df = Sample_Info[!duplicated(Sample_Info$uuid),c("uuid","city","latitude","longitude")]
df = ddply(df, .variables = c("city","latitude","longitude"),nrow)
#df$V1 = log10(df$V1)
pdf("SampleNumber_Geographic.pdf", width = 10, height = 7)
mp  = NULL #定义一个空的地图
mapworld = borders("world",colour = "gray50",fill="white") #绘制基本地图
ggplot(data = df) + mapworld + coord_cartesian() +
  geom_point(aes(x=longitude,y=latitude, colour = V1),alpha = 0.7,size = 3)+
  scale_y_continuous(breaks = (-2:2) * 45) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  scale_size(range=c(2,9)) +
  theme_bw() +
  theme(axis.text = element_text(size=14),axis.title = element_text(size = 16),
        legend.position="top",legend.direction = "horizontal",panel.grid=element_blank())
dev.off()





