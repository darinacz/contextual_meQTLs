#are our CpGs enriched for MR results from Richardson et al.?
#https://academic.oup.com/hmg/article/27/18/3293/5034850?login=false#120839527
#systematic MR: meQTLs and complex traits
#meQTLs in blood in 1,000 samples
#complex traits from MRbase.org

#significant MR-results from Richardson et al.
mr<-read.delim("mr_results_richardson.txt",sep="\t",header=T)

#full list of results
meta<-read.table("meta.txt",sep="\t",header=T) #this is the full results list
all<-read.table("results_combined_all.txt",sep="\t",header=T) #this is the list of significant contmeQTLs

#how many of our tested CpGs are in the MR?
cpgs_all<-unique(meta$CpG) #n=155,640
length(which(cpgs_all %in% mr$CpG.site)) #n=517: 0.33 %
cpgs_initial<-unique(all$CpG) #n=5,115
length(which(cpgs_initial %in% mr$CpG.site)) #n=32: 0.6 %


#is this significantly more?
prop.test(x=c(32,517),n=c(5115,155640), alternative="greater") #p=0.00032

#with which complex traits are these associated?
#read in categories
cat<-read.delim("cat_traits.txt",sep="\t",header=T)
index<-(which(mr$CpG.site %in% cpgs_all)) #n=758
mr_all<-mr[index,]

cat_all<-merge(mr_all,cat, by.x="Complex.trait",by.y="trait")
table(cat_all$category)
index<-(which(cat_all$CpG.site %in% cpgs_initial))
table(cat_all$category[index])

#comparisons for each domain were run via Fisher's enrichment tests
all_cpgs<-cat_all
top_cpgs<-cat_all[index,]

#compare
all_cpg<-as.data.frame(table(all_cpgs$category))
all_cpg$Freq <- as.numeric(all_cpg$Freq)
all_cpg$freq<-all_cpg$Freq/(dim(all_cpgs)[1])
names(all_cpg)[1]<-"cat"
all_cpg[,1]<-as.character(all_cpg[,1])

top_cpg<-as.data.frame(table(top_cpgs$category))
top_cpg$Freq <- as.numeric(top_cpg$Freq)
top_cpg$freq<-top_cpg$Freq/(dim(top_cpgs)[1])
names(top_cpg)[1]<-'cat'

top_cpg<-as.data.frame(top_cpg)
top_cpg$cat<-as.character(top_cpg$cat)

#fill up the categories missing MR conmeQTLCpGs
top_cpg[10,]<-c("Cancer","0","0")
top_cpg[11,]<-c("Cognitive","0","0")
top_cpg[12,]<-c("Eating Disorders","0","0")
top_cpg[13,]<-c("Immunological","0","0")
top_cpg[14,]<-c("Kidney","0","0")
top_cpg[15,]<-c("Neurological","0","0")
top_cpg[16,]<-c("Pregnancy/Birth","0","0")
top_cpg[17,]<-c("Smoking","0","0")
top_cpg[18,]<-c("Urine","0","0")

top_cpg<-top_cpg[order(top_cpg$cat),]
top_cpg$Freq <- as.numeric(top_cpg$Freq)

all.equal(all_cpg[,1],top_cpg[,1]) #TRUE

#if relative frequency less than 1%, put these into the "other" category
which(top_cpg$freq<0.01) -> index
remove_cat<-top_cpg$cat[index]
top_adapt<-top_cpg
index<-which(top_adapt$cat %in% remove_cat)
top_adapt[19,]<-c("other",0,0)
top_adapt[19,2] <- sum(as.numeric(top_adapt$Freq[index]))
top_adapt[19,3] <- as.numeric(top_adapt[19,2])/(dim(top_cpgs)[1])
index<-which(top_adapt$cat %in% remove_cat)
top_adapt <- top_adapt[-index,]
top_adapt$Freq <- as.numeric(top_adapt$Freq)

all_adapt<-all_cpg
index<-which(all_adapt$cat %in% remove_cat)
all_adapt[19,]<-c("other",0,0)
all_adapt[19,2] <- sum(as.numeric(all_adapt$Freq[index]))
all_adapt[19,3] <- as.numeric(all_adapt[19,2])/(dim(all_cpgs)[1])
index<-which(all_adapt$cat %in% remove_cat)
all_adapt <- all_adapt[-index,]
all_adapt$Freq <- as.numeric(all_adapt$Freq)

all.equal(all_adapt[,1],top_adapt[,1]) #TRUE


#counts of events per group

events<-c(dim(top_cpgs)[1], dim(all_cpgs)[1] - dim(top_cpgs)[1]) #this is always the same, number of overall obs per group

results<-matrix(nrow=10, ncol=3)
for (i in 1:10)
{exposure<-c(top_adapt$Freq[i],all_adapt$Freq[i]-top_adapt$Freq[i]) #make disjoint groups, so trait association in contmeQTL CpGs and trait associations in non contmeQTLs CpGs
mod<-poisson.test(exposure, events, alternative="greater")
results[i,1] <- top_adapt$cat[i]
results[i,2] <- mod$p.value
results[i,3] <- mod$estimate}

results<-as.data.frame(results)
names(results)<-c("cat","p","rate_ratio")
for (i in 2:3)
{results[,i]<-as.numeric(results[,i])}






