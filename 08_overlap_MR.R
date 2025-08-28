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

#comparisons for each domain were run via prop.test


