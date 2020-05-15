library(plyr)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(stringr)
library(ggsci)
source("theme_basic.R")
load("data/OTU_tax.RData")
load("data/OTU_Abundance.RData")
load("data/Sample_Info.RData")
OTU_factors = merge(merge(OTU_Abundance[,c(1,2,3)], OTU_tax[,c("OTU","Species")], by="OTU"), Sample_Info, by="uuid")
OTU_factors$city_population_density = as.numeric(OTU_factors$city_population_density)

#Total Population vs. Novel OTU numer
df = OTU_factors[,c("uuid","OTU","Species","num_reads","city_total_population")]
df = ddply(df, .variables = c("uuid","city_total_population"), .fun = function(x){
  num = (sum(x$Species == "s__")+0.1)/x$num_reads[1] * 1e6
  return(num)
})
df$city_total_population = df$city_total_population/1e6


g3_A = ggplot(data = df, mapping = aes(x=city_total_population, y = V1)) +
  geom_point(size = 1) + geom_smooth(method = lm, fill="#ffc0cb",span =13) +
  stat_cor(method = "pearson", size=5, color='red') +
  ylab("# Novel OTUs per million reads") + xlab("City Population (Million)") + 
  theme_0

#Population density vs. Novel OTU numer
df = OTU_factors[,c("uuid","OTU","Species","num_reads","city_population_density")]
df = df[!is.na(df$city_population_density),]
df = ddply(df, .variables = c("uuid","city_population_density"), .fun = function(x){
  num = (sum(x$Species == "s__")+0.1)/x$num_reads[1] * 1e6
  return(num)
})
g3_B = ggplot(data = df, mapping = aes(x=city_population_density, y = V1)) +
  geom_point(size = 1) + geom_smooth(method = lm, fill="#ffc0cb",span =13) +
  stat_cor(method = "pearson", size=5, color='red') +
  ylab("# Novel OTUs per million reads") + xlab("City Population Density") + 
  theme_0

# Surface
df = OTU_factors[,c("uuid","OTU","Species","num_reads","surface_material")]
df = df[!(df$surface_material=="" | df$surface_material=="-"),]
df = ddply(df, .variables = c("uuid", "surface_material"), .fun = function(x){
  total_num = (length(unique(x$OTU))+0.1)/x$num_reads[1] * 1e6
  novel_num = (sum(x$Species == "s__")+0.1)/x$num_reads[1] * 1e6
  return(c(total_num, novel_num))
})
colnames(df)[3:4] = c("total_num","novel_num")
df = df[df$surface_material %in% names(table(df$surface_material))[table(df$surface_material)>=50],]
df$surface_material = str_to_title(df$surface_material)
df = melt(df[,-1], id.vars = "surface_material")
df$variable = as.character(df$variable)
df$variable[df$variable == "total_num"] = "All OTUs identified"
df$variable[df$variable == "novel_num"] = "Novel OTUs identified"
df$surface_material = factor(df$surface_material, levels = c("Skin","Glass","Wood","Metal","Plastic","Cement"))

g3_C = ggplot(data = df, mapping = aes(x=surface_material, y=value)) +
  geom_boxplot(aes(fill=variable),outlier.shape = NA) + 
  facet_wrap(~variable) + xlab("Surface Material") +
  ylab("Number per million reads") + scale_fill_npg() + theme_0 +
  theme(strip.text = element_text(size = 16, color = "black"))


# Coastal
df = OTU_factors[,c("uuid","OTU","Species","num_reads","coastal_city")]
df = df[!(df$coastal_city=="" | df$coastal_city == "nan"),]
df = ddply(df, .variables = c("uuid", "coastal_city"), .fun = function(x){
  total_num = (length(unique(x$OTU))+0.1)/x$num_reads[1] * 1e6
  novel_num = (sum(x$Species == "s__")+0.1)/x$num_reads[1] * 1e6
  return(c(total_num, novel_num))
})
df$coastal_city[df$coastal_city == "yes"] = "Coastal"
df$coastal_city[df$coastal_city == "no"] = "Inland"
colnames(df)[3:4] = c("total_num","novel_num")
df = melt(df[,-1], id.vars = "coastal_city")
df$variable = as.character(df$variable)
df$variable[df$variable == "total_num"] = "All OTUs identified"
df$variable[df$variable == "novel_num"] = "Novel OTUs identified"
df$coastal_city = factor(df$coastal_city, levels = c("Coastal","Inland"))

g3_D = ggplot(data = df, mapping = aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=coastal_city),outlier.shape = NA) + xlab("") +
  stat_compare_means(aes(group=coastal_city), label.y = c(2.0, 1.0),size=5) + ylim(0,2.5)+
  ylab("Number per million reads") + scale_fill_npg() + theme_classic() +
  theme_2 +theme(legend.title = element_blank(),
                 legend.position = c(0.85,0.9))
fig_wd = 8
pdf("figures/Fig3_A.pdf", width = fig_wd, height = 0.6*fig_wd)
g3_A
dev.off()
pdf("figures/Fig3_B.pdf", width = fig_wd, height = 0.6*fig_wd)
g3_B
dev.off()
pdf("figures/Fig3_C.pdf", width = fig_wd, height = 0.6*fig_wd)
g3_C
dev.off()
pdf("figures/Fig3_D.pdf", width = fig_wd, height = 0.6*fig_wd)
g3_D
dev.off()

