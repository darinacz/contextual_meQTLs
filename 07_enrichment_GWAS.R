library(data.table)
library(ggplot2)
library(ChIPseeker)
library(GenomicRanges)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("readxl")

#full list of results
meta<-read.table("meta.txt",sep="\t",header=T) #this is the full results list
all<-read.table("results_combined_all.txt",sep="\t",header=T) #this is the list of significant contmeQTLs

#load genome-wide significant results from GWAS, this is based on the IEU GWAS platform and filtered for genome-wide significant associations
load("genome_wide_hits_traits.Rdata") # this is called final

#all in meta
index<-(which(final$rsid %in% meta$rsID))
gwas_top<-final[index,]
gwas_top<-unique(gwas_top) 
length(unique(gwas_top$rsid)) #n=45,993 SNPs
length(unique(meta$rsID)) #n=133,121 => 34.5% in GWAS top hits
length(which(meta$rsID %in% gwas_top$rsid)) #n=61740 => 61740/167600 =36%   

length(unique(gwas_top$rsid)) 
length(unique(gwas_top$id)) 
write.table(unique(gwas_top$id),"overlap_ieu_gwas_catalogue.txt",sep="\t",quote=F,row.names=F)
write.table(gwas_top,"overlap_ieu_gwas_catalogue_5e-08_all.txt",sep="\t",quote=F,row.names=F)

#all in contextual meQTLs
index<-which(meta$SNP %in% all$SNP)
test<-meta[index,]
index<-(which(gwas_top$rsid %in% test$rsID))
gwas_initial<-gwas_top[index,]
gwas_initial<-unique(gwas_initial) 
length(unique(gwas_initial$rsid)) #n=2,307
length(unique(test$rsID)) #n=5,466 => 42.4 % in GWAS top hits
length(which(test$rsID %in% gwas_top$rsid)) #n=3778 => 3778/8435 =44.7%   

prop.test(x = c(61740,3778), n = c(167600,8435)) #p=4.11e-49=> overall GWAS hits more frequent in contextual meQTLs
length(unique(gwas_initial$rsid)) 
length(unique(gwas_initial$id))

#compare distributions across the categories
#categories where created using ChatGPT 3.5 and manually curated afterwards
cat<-read_excel("GWAS_cat_ChatGPT.xlsx")
gwas_top<-merge(gwas_top,cat,by.x="trait",by.y="Trait")
gwas_initial<-merge(gwas_initial,cat,by.x="trait",by.y="Trait")
table(gwas_top$Category)
table(gwas_inital$Category)

#comparisons for each domain via prop.test


