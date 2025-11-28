#read in list of significant contmeQTLs
all<-read.table("results_combined_all.txt",sep="\t",header=T) #n=5,120 contmeQTLs

########1. enrichment for CpGs affected by vSNPs, based on list from Ek et al. (450)#######
cpgs<-read.table("7195_CpG_vme_Ek_et_al.txt",sep="\t",header=T) #this is the list of CpGs identified in Ek et al.
length(unique(cpgs$CpG)) #n=7,195 CpGs
index<-which(all$CpG %in% cpgs$CpG) 
length(unique(all$CpG[index])) #n=411 CpGs out of 5,115 CpGs

#is this more than expected by chance?
#for this we need the full list
all_CpGs<-read.table("all_CpGs.txt",sep="\t",header=T) #this is the list of all tested CpGs (available in at least three of the six cohorts)
all_CpGs<-unique(all_CpGs) #n=155,640 CpGs

#take random subset of 5,115 CpGs and check how many overlap
x<-seq(1,155640,by=1)
overlap<-NA
for (i in 1 :1000)
{index<-sample(x,size=5115,rep=F)
overlap<-rbind(overlap,length(which(all_CpGs[index] %in% cpgs$CpG)))
print(i)}
overlap<-overlap[-1]
length(which(overlap>=411)) #n=0
#p_emp=1/1001=0.000999001 => overlap is higher than expected by chance

########2. enrichment for SNPs identified as vSNPs, based on list from Ek et al. (450)##########
snps<-read.table("7195_SNPs_vme_Ek_et_al.txt",sep="\t",header=T) ##this is the list of SNPs identified in Ek et al.
snps<-unique(snps) #n=6,456 unique SNPs
index<-which(all$SNPs %in% snps) 
length(unique(all$SNPs[index])) ##n=43 SNPs out of 5,466 SNPs

# is more than expected by chance?
#for this we need the full list
snpnames<-read.table("all_SNPs.txt",sep="\t",header=T) #this is the list of all tested SNPs 
snpnames<-unique(snpnames) #n=133,121 SNPs

#take random subset of 5,466 SNPs and check how many overlap
x<-seq(1,133121,by=1)
overlap<-NA
for (i in 1 :1000)
{index<-sample(x,size=5466,rep=F)
overlap<-rbind(overlap,length(which(all_SNPs[index] %in% snps)))
print(i)}
overlap<-overlap[-1]
length(which(overlap>=43)) #n=5
#p_emp=6/1001=0.005994006 => overlap is higher than expected by chance

