library(data.table)
library(ggplot2)
library(ChIPseeker)
library(GenomicRanges)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("readxl")

#full list of results
meta<-read.table("meta.txt",sep="\t",header=T) #this is the full results list
all<-read.table("results_combined_all.txt",sep="\t",header=T) #this is the list of significant contmeQTLs

#####overlap with EWAS results, based on EWAScatalgue###
##results from the EWAS catalogue were downloaded on Feb 19th 2024
##results with p<e-04 are available here, studies of at least 100 individuals, epi-genome wide
##https://www.ewascatalog.org/download/
##all results with p<1e-04 are available###
ewas<-fread("/Users/darina/ewascatalog-results_1902024.txt") #2,299,164 results
ewas_top<-ewas[ewas$P<9e-8,] #1,257,086 results epigenome-wide significant results

#all in meta
index<-(which(ewas_top$CpG %in% meta$CpG))
ewas_top<-ewas_top[index,]
ewas_top<-unique(ewas_top) 
length(unique(ewas_top$CpG)) #n=146,053
length(unique(meta$CpG)) #n=155,640 
length(which(meta$CpG %in% ewas_top$CpG)) #n=157,447 => 157447/167600=93.9%

#all in contextual meQTLs
index<-(which(ewas_top$CpG %in% all$CpG))
ewas_initial<-ewas_top[index,]
ewas_initial<-unique(ewas_initial) 
length(unique(ewas_initial$CpG)) #n=4,815
length(unique(all$CpG)) #n=5,115 => 94.1% in EWAS top hits
length(which(all$CpG %in% ewas_top$CpG)) #n=4,820 => 4820/5120=94%

prop.test(x = c(157447, 4820), n = c(167600, 5120)) #p=0.5775 => n.s. diff. in overall EWAS hits between all and meta 

#compare distributions across the categories
#categories where created using ChatGPT 3.5 and manually curated afterwards
cat<-read_excel("GWAS_cat_ChatGPT.xlsx")
gwas_top<-merge(gwas_top,cat,by.x="trait",by.y="Trait")
gwas_initial<-merge(gwas_initial,cat,by.x="trait",by.y="Trait")
table(gwas_top$Category)
table(gwas_inital$Category)

#compare distributions across the categories
#categories where created using ChatGPT 3.5 and manually curated afterwards
cat<-read_excel("EWAS_cat_ChatGPT.xlsx", header=F)

ewas_top<-merge(ewas_top,cat,by.x="StudyID",by.y="Trait")
ewas_initial<-merge(ewas_initial,cat,by.x="StudyID",by.y="Trait")

table(ewas_top$Category)
table(ewas_inital$Category)

#comparisons for each domain were performed via Fisher's enrichment tests

all<-as.data.frame(table(ewas_top$trait_domain))
names(all)[1]<-'Cat'
all$freq<-all$Freq/(dim(ewas_top)[1])
all[,1]<-as.character(all[,1])
all$Freq <- as.numeric(all$Freq)


meta<-as.data.frame(table(ewas_initial$trait_domain))
names(meta)[1]<-'Cat'
meta$freq<-meta$Freq/(dim(gwas_initial)[1])
meta<-as.data.frame(meta)
meta$Cat<-as.character(meta$Cat)

meta[33,]<-c("Anthropometrics","0","0")
meta[34,]<-c("Dental","0","0")
meta[35,]<-c("Metabolomics","0","0")
meta[36,]<-c("Microbiome","0","0")
meta[37,]<-c("Mortality","0","0")
meta[38,]<-c("Nutritional","0","0")
meta[39,]<-c("Social Support","0","0")
meta<-meta[order(meta$Cat),]
meta$Freq <- as.numeric(meta$Freq)

all.equal(all[,1],meta[,1]) #TRUE

#count events per group
results<-matrix(nrow=39, ncol=3)

events<-c(dim(ewas_initial)[1], dim(ewas_top)[1] - dim(ewas_initial)[1]) #this is always the same, number of overall obs per group

for (i in 1:39)
  
{exposure<-c(meta$Freq[i],all$Freq[i]-meta$Freq[i]) #make disjoint groups, so trait association in contmeQTLs and trait associations in non contmeQTLs
mod<-poisson.test(exposure, events, alternative="greater")
results[i,1] <- meta$Cat[i]
results[i,2] <- mod$p.value
results[i,3] <- mod$estimate}

results<-as.data.frame(results)
names(results)<-c("cat","p","rate_ratio")
for (i in 2:3)
{results[,i]<-as.numeric(results[,i])}

#with freq > 1%
index<-which(meta$freq<0.01)
remove_cat<-meta$Cat[index]
index<-which(results$cat %in% remove_cat)
results[-index,]

