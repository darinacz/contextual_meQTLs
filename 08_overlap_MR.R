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

results<-matrix(nrow=10, ncol=3)
for (i in 1:10)
{
  results[i,1] <- top_adapt$cat[i]  
  mod<-fisher.test(matrix(c(top_adapt$Freq[i],dim(top_cpgs)[1]-top_adapt$Freq[i],all_adapt$Freq[i],dim(all_cpgs)[1]-all_adapt$Freq[i]),ncol=2), alternative="greater")
  mod$p.value -> results[i,2]
  mod$estimate -> results[i,3]}

results<-as.data.frame(results)
names(results)<-c("cat","p","OR")
for (i in 2:3)
{results[,i]<-as.numeric(results[,i])}





index<-which(results$FDR <0.05)
results[index,] #none

#use color-coding to match with GWAS results
#Anthropometrics "#ACA985"
#Cardiovascular  "#CF948C"
#Endocrine "#9B9DBD"
#Gastrointestinal "#7B9DBD"
#Haematological "#F29066"
#Immunological "#D58EC4"
#Mental "#D51EC4"
#Metabolic "#E1C296"
#Musculoskeletal "#B7C371"
#Other "#B3B3B3"
#Reproductive "#EFCC6B",
#Respiratory  "#FBD63C"
#SES "#FBD08C"

#merge with results from GWAS enrichment and show both in one plot
res<-read.delim("../annotation/overlap_GWAS/results_GWAS_enrichment_fisher_120326.txt", header=T)
both<-merge(res,results,by.x="cat",by.y="cat",all.x=T,all.y=T)

names(both)<-c("category","p_SNPs","OR_SNPs","FDR_SNPs","p_CpGs","OR_CpGs","FDR_CpGs")

#convert to long format
df_long <- both %>%
  pivot_longer(cols = c(OR_SNPs, OR_CpGs),
               names_to = "group",
               values_to = "OR")

#remove category "other"
index<-which(df_long$category=="other")
df_long<-df_long[-index,]

index<-which(df_long$category=="Autoimmune immunological disease")
df_long$category[index]<-"Autoimmunological"

#order by ORs in SNPs
ref_group <- "OR_SNPs"


# shared = categories present in BOTH groups
#remove category other
index<-which(both$category=="other")
both<-both[-index,]
both$category[2] <- "Autoimmunological"

shared_cats <- both %>%
  filter(!is.na(OR_SNPs) & !is.na(OR_CpGs)) %>%
  pull(category)

unique_cats <- setdiff(both$category, shared_cats)


#shared order according to OR_SNPs
shared_order <- both %>%
  filter(category %in% shared_cats) %>%
  arrange(desc(OR_SNPs)) %>%  # descending OR_SNPs
  pull(category)

#unique order
unique_order <- both %>%
  filter(category %in% unique_cats) %>%
  mutate(max_OR = pmax(OR_SNPs, OR_CpGs, na.rm=TRUE)) %>%
  arrange(desc(max_OR)) %>%
  pull(category)

category_order <- c(shared_order, unique_order)

df_long$category <- factor(df_long$category, levels = category_order)

cat_colors <- c(
  "#66C2A5",  # teal-green
  "#DB9810",  # warm orange
  "#ACA985",  # muted olive
  "#CF948C",  # soft rose
  "#9B9DBD",  # muted periwinkle
  "#F29066",  # soft coral
  "#D58EC4",  # muted pink-purple
  "#E1F296",  # light lime
  "#B7C371",  # moss green
  "#B7D84C",  # yellow-green
  "#AB98C8",  # lavender
  "#E1D83B",  # muted yellow
  "#FBE63C",  # soft yellow
  "#EFCC6B",  # sandy yellow-orange
  "#C3B3B3",  # grey
  "#D59EA5",  # muted rose
  "#A1C296",  # soft green
  "#A29066",  # brown-beige
  "#8CB5D6",  # soft blue
  "#D6A07F"   # terracotta
)

df_long$group <- factor(df_long$group, levels = c("OR_SNPs", "OR_CpGs"))

custom_labels <- c(
  OR_CpGs = "CpG",
  OR_SNPs = "SNP"
)

ggplot(df_long,
       aes(x = OR, y = category, fill = category))+
  geom_col() +
  facet_wrap(~group, ncol = 2, labeller = labeller(group = custom_labels)) +
  scale_fill_manual(values = cat_colors) + 
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_bw(base_size = 11) +
  ylab(NULL) +
  theme(legend.position = "none")  + 
  #write categories significant in at least one omics in bold
  theme(axis.text.y = element_text(face = c('bold','bold',
                                            'bold', 'bold',
                                            'plain', 'plain',
                                            'plain', 'plain',
                                            'plain','bold',
                                            'bold', 'bold',
                                            'bold', 'bold',
                                            'plain','plain',
                                            'plain','plain',
                                            'plain','plain')))



ggplot(df_long,
       aes(x = OR, y = category, fill = category))+
  geom_col() +
  facet_wrap(~group, ncol = 2, labeller = labeller(group = custom_labels)) +
  scale_fill_manual(values = cat_colors) + 
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_bw(base_size = 11) +
  ylab(NULL) +
  theme(legend.position = "none")  + 
  #write categories significant in at least one omics in bold
  theme(axis.text.y = element_text(face = c('bold','plain',
                                            'plain', 'bold',
                                            'plain', 'plain',
                                            'plain', 'plain',
                                            'plain','bold',
                                            'bold', 'bold',
                                            'bold', 'bold',
                                            'plain','plain',
                                            'plain','plain',
                                            'plain','plain'))) 



