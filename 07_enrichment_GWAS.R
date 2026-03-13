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

#comparisons for each domain were performed via Fisher's enrichment tests
all<-as.data.frame(table(gwas_top$domain))
names(all)[1]<-'Cat'
all$freq<-all$Freq/(dim(gwas_top)[1])
all[,1]<-as.character(all[,1])
all$Freq <- as.numeric(all$Freq)

meta<-as.data.frame(table(gwas_initial$domain))
names(meta)[1]<-'Cat'
meta$freq<-meta$Freq/(dim(gwas_initial)[1])
meta<-as.data.frame(meta)
meta$Cat<-as.character(meta$Cat)
#fill up the categories missing in the contmeQTLs
meta[34,]<-c("Geographical","0","0")
meta[35,]<-c("Trauma","0","0")
meta<-meta[order(meta$Cat),]
meta$Freq <- as.numeric(meta$Freq)

all.equal(all[,1],meta[,1]) #TRUE

#if relative frequency less than 1%, put these into the "other" category
#compare also to categories available in the MR-analysis as we want to put these results next to each other
which(meta$freq<0.01 & ! meta$Cat=="Puberty related" & ! meta$Cat=="Reproductive" & !meta$Cat=="Gastrointestinal") -> index
#leave puberty & reproductive & gastro in as this is available in the MR analysis
remove_cat<-meta$Cat[index]
meta_adapt<-meta
index<-which(meta_adapt$Cat %in% remove_cat)
meta_adapt[36,]<-c("other",0,0)
meta_adapt[36,2] <- sum(as.numeric(meta_adapt$Freq[index]))
meta_adapt[36,3] <- as.numeric(meta_adapt[36,2])/(dim(gwas_initial)[1])
index<-which(meta_adapt$Cat %in% remove_cat)
meta_adapt <- meta_adapt[-index,]
meta_adapt$Freq <- as.numeric(meta_adapt$Freq)

all_adapt<-all
index<-which(all_adapt$Cat %in% remove_cat)
all_adapt[36,]<-c("other",0,0)
all_adapt[36,2] <- sum(as.numeric(all_adapt$Freq[index]))
all_adapt[36,3] <- as.numeric(all_adapt[36,2])/(dim(gwas_top)[1])
index<-which(all_adapt$Cat %in% remove_cat)
all_adapt <- all_adapt[-index,]
all_adapt$Freq <- as.numeric(all_adapt$Freq)

all.equal(all_adapt[,1],meta_adapt[,1]) #TRUE

results<-matrix(nrow=21, ncol=3)
for (i in 1:21)
{
  results[i,1] <- meta_adapt$Cat[i]  
  mod<-fisher.test(matrix(c(meta_adapt$Freq[i],dim(gwas_initial)[1]-meta_adapt$Freq[i],all_adapt$Freq[i],dim(gwas_top)[1]-all_adapt$Freq[i]),ncol=2), alternative="greater")
  mod$p.value -> results[i,2]
  mod$estimate -> results[i,3]}

results<-as.data.frame(results)
names(results)<-c("cat","p","OR")
for (i in 2:3)
{results[,i]<-as.numeric(results[,i])}



