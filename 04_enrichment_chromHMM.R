#enrichment for chromHMM states
setwd("/Users/darina/Contextual_meQTLs/chromHMM/")

#this data was downloaded from http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/
states <- readRDS("chromhmm_states.Rds")
codes <- read.delim('epigenomes.tsv', header = T, stringsAsFactors = F)
hmm <- readRDS('/chromHMM_all_states.Rds')
hmm$group <- codes$GROUP[match(hmm$code, codes$EID)]

## blood data
table(codes$GROUP)
index<-which(codes$GROUP=="Blood & T-cell")
interesting.tissues <- codes$EID[index]
df_blood <- subset(hmm,  (code %in% interesting.tissues))
df_blood <- subset(df_blood,  !(seqnames %in% 'chrM'))
df_blood <- subset(df_blood,  !(seqnames %in% 'chrY'))
df_blood <- subset(df_blood,  !(seqnames %in% 'chrX'))

#full list of results
meta<-read.table("meta.txt",sep="\t",header=T) #this is the full results list
all<-read.table("results_combined_all.txt",sep="\t",header=T) #this is the list of significant contmeQTLs

####1.overlap of CpGs with chromHMM states#########
load("annot_epic.Rdata") #this is Illumina's EPIC v1 annotation
meta_CpGs <- unique(meta$CpG) #n=155,640 CpGs
meta_CpGs<-as.data.frame(meta_CpGs)
names(meta_CpGs)[1]<-"CpG"
meta_CpGs<-merge(meta_CpGs, annot_epic,by.x=1, by.y="Name")
meta_CpGs$end<-meta_CpGs$pos+1

#make genomic ranges dataset
#all
cpgs<-makeGRangesFromDataFrame(meta_CpGs,start.field ="pos",end.field="end",seqnames="chr",ignore.strand=T,keep.extra.columns=F)
#significant contmeQTLs
index<-which(meta_CpGs$CpG %in% all$CpG) #n=5,115
initial<-makeGRangesFromDataFrame(meta_CpGs[index,],start.field ="pos",end.field="end",seqnames="chr",ignore.strand=T,keep.extra.columns=F)
#all without initial
index<-which(meta_CpGs$CpG %in% all$CpG) #n=5,115
rest_initial<-makeGRangesFromDataFrame(meta_CpGs[-index,],start.field ="pos",end.field="end",seqnames="chr",ignore.strand=T,keep.extra.columns=F)

#test chromHMM states
results<-matrix(ncol=3,nrow=15)

###1. tssA###
x1_tssA = subset(df_blood, (name == "1_TssA"))
#initial vs. all
n_tssA_initial = sum(distanceToNearest(initial, x1_tssA)@elementMetadata@listData$distance == 0)
n_x1_tssA_no = sum(distanceToNearest(rest_initial, x1_tssA)@elementMetadata@listData$distance == 0)

a=n_tssA_initial 
b=length(initial) - a
c=n_x1_tssA_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"tssA"->results[1,1]
test$p.value->results[1,2]
test$estimate->results[1,3]

###2. TssAFlnk###
x2_tssAflnk = subset(df_blood, (name == "2_TssAFlnk"))
#initial vs. all
n_x2_tssAflnk_initial = sum(distanceToNearest(initial, x2_tssAflnk)@elementMetadata@listData$distance == 0)
n_x2_tssAflnk_no = sum(distanceToNearest(rest_initial, x2_tssAflnk)@elementMetadata@listData$distance == 0)

a=n_x2_tssAflnk_initial 
b=length(initial) - a
c=n_x2_tssAflnk_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"TssAFlnk"->results[2,1]
test$p.value->results[2,2]
test$estimate->results[2,3]

###3.TxFlnk###
x3_txflnk = subset(df_blood, (name == "3_TxFlnk"))

#initial vs. all
n_x3_txflnk_initial = sum(distanceToNearest(initial, x3_txflnk)@elementMetadata@listData$distance == 0)
n_x3_txflnk_no = sum(distanceToNearest(rest_initial, x3_txflnk)@elementMetadata@listData$distance == 0)

a=n_x3_txflnk_initial
b=length(initial) - a
c=n_x3_txflnk_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"TxFlnk"->results[3,1]
test$p.value->results[3,2]
test$estimate->results[3,3]

###4.Tx###
x4_tx= subset(df_blood, (name == "4_Tx"))

#initial vs. all
n_x4_tx_initial = sum(distanceToNearest(initial, x4_tx)@elementMetadata@listData$distance == 0)
n_x4_tx_no = sum(distanceToNearest(rest_initial, x4_tx)@elementMetadata@listData$distance == 0)

a=n_x4_tx_initial
b=length(initial) - a
c=n_x4_tx_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"Tx"->results[4,1]
test$p.value->results[4,2]
test$estimate->results[4,3]

###5.Tx_wk###
x5_txwk= subset(df_blood, (name == "5_TxWk"))

#initial vs. all
n_x5_txwk_initial = sum(distanceToNearest(initial, x5_txwk)@elementMetadata@listData$distance == 0)
n_x5_txwk_no = sum(distanceToNearest(rest_initial, x5_txwk)@elementMetadata@listData$distance == 0)

a=n_x5_txwk_initial
b=length(initial) - a
c=n_x5_txwk_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"Tx_wk"->results[5,1]
test$p.value->results[5,2]
test$estimate->results[5,3]


###6.EnhG###
x6_enhg= subset(df_blood, (name == "6_EnhG"))

#initial vs. all
n_x6_enhg_initial = sum(distanceToNearest(initial, x6_enhg)@elementMetadata@listData$distance == 0)
n_x6_enhg_no = sum(distanceToNearest(rest_initial, x6_enhg)@elementMetadata@listData$distance == 0)

a=n_x6_enhg_initial
b=length(initial) - a
c=n_x6_enhg_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"EnhG"->results[6,1]
test$p.value->results[6,2]
test$estimate->results[6,3]

###7.Enh###
x7_enh= subset(df_blood, (name == "7_Enh"))

#initial vs. all
n_x7_enh_initial = sum(distanceToNearest(initial, x7_enh)@elementMetadata@listData$distance == 0)
n_x7_enh_no = sum(distanceToNearest(rest_initial, x7_enh)@elementMetadata@listData$distance == 0)

a=n_x7_enh_initial
b=length(initial) - a
c=n_x7_enh_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"Enh"->results[7,1]
test$p.value->results[7,2]
test$estimate->results[7,3]


###8.ZNF###
x8_znf= subset(df_blood, (name == "8_ZNF/Rpts"))

#initial vs. all
n_x8_znf_initial = sum(distanceToNearest(initial, x8_znf)@elementMetadata@listData$distance == 0)
n_x8_znf_no = sum(distanceToNearest(rest_initial, x8_znf)@elementMetadata@listData$distance == 0)

a=n_x8_znf_initial
b=length(initial) - a
c=n_x8_znf_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"ZNF"->results[8,1]
test$p.value->results[8,2]
test$estimate->results[8,3]

###9.Het###
x9_het= subset(df_blood, (name == "9_Het"))

#initial vs. all
n_x9_het_initial = sum(distanceToNearest(initial, x9_het)@elementMetadata@listData$distance == 0)
n_x9_het_no = sum(distanceToNearest(rest_initial, x9_het)@elementMetadata@listData$distance == 0)

a=n_x9_het_initial
b=length(initial) - a
c=n_x9_het_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"Het"->results[9,1]
test$p.value->results[9,2]
test$estimate->results[9,3]

###10.TSSBiv###
x10_tssbiv= subset(df_blood, (name == "10_TssBiv"))


#initial vs. all
n_x10_tssbiv_initial = sum(distanceToNearest(initial, x10_tssbiv)@elementMetadata@listData$distance == 0)
n_x10_tssbiv_no = sum(distanceToNearest(rest_initial, x10_tssbiv)@elementMetadata@listData$distance == 0)

a=n_x10_tssbiv_initial
b=length(initial) - a
c=n_x10_tssbiv_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"TSSBiv"->results[10,1]
test$p.value->results[10,2]
test$estimate->results[10,3]

###11.TSSBivFlnk###
x11_bivflnk= subset(df_blood, (name == "11_BivFlnk"))


#initial vs. all
n_x11_bivflnk_initial = sum(distanceToNearest(initial, x11_bivflnk)@elementMetadata@listData$distance == 0)
n_x11_bivflnk_no = sum(distanceToNearest(rest_initial, x11_bivflnk)@elementMetadata@listData$distance == 0)

a=n_x11_bivflnk_initial
b=length(initial) - a
c=n_x11_bivflnk_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"TSSBivFlnk"->results[11,1]
test$p.value->results[11,2]
test$estimate->results[11,3]

###12.EnhBiv###
x12_enhbiv= subset(df_blood, (name == "12_EnhBiv"))


#initial vs. all
n_x12_enhbiv_initial = sum(distanceToNearest(initial, x12_enhbiv)@elementMetadata@listData$distance == 0)
n_x12_enhbiv_no = sum(distanceToNearest(rest_initial, x12_enhbiv)@elementMetadata@listData$distance == 0)

a=n_x12_enhbiv_initial
b=length(initial) - a
c=n_x12_enhbiv_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"EnhBiv"->results[12,1]
test$p.value->results[12,2]
test$estimate->results[12,3]

###13.ReprPC###
x13_reprpc= subset(df_blood, (name == "13_ReprPC"))


#initial vs. all
n_x13_reprpc_initial = sum(distanceToNearest(initial, x13_reprpc)@elementMetadata@listData$distance == 0)
n_x13_reprpc_no = sum(distanceToNearest(rest_initial, x13_reprpc)@elementMetadata@listData$distance == 0)

a=n_x13_reprpc_initial
b=length(initial) - a
c=n_x13_reprpc_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"ReprPC"->results[13,1]
test$p.value->results[13,2]
test$estimate->results[13,3]

###14.ReprPCWk###
x14_reprpcwk= subset(df_blood, (name == "14_ReprPCWk"))


#initial vs. all
n_x14_reprpcwk_initial = sum(distanceToNearest(initial, x14_reprpcwk)@elementMetadata@listData$distance == 0)
n_x14_reprpcwk_no = sum(distanceToNearest(rest_initial, x14_reprpcwk)@elementMetadata@listData$distance == 0)

a=n_x14_reprpcwk_initial
b=length(initial) - a
c=n_x14_reprpcwk_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"ReprPCWk"->results[14,1]
test$p.value->results[14,2]
test$estimate->results[14,3]

###15.Quies###
x15_quies= subset(df_blood, (name == "15_Quies"))

#initial vs. all
n_x15_quies_initial = sum(distanceToNearest(initial, x15_quies)@elementMetadata@listData$distance == 0)
n_x15_quies_no = sum(distanceToNearest(rest_initial, x15_quies)@elementMetadata@listData$distance == 0)

a=n_x15_quies_initial
b=length(initial) - a
c=n_x15_quies_no
d=length(rest_initial)-c

test<-fisher.test(matrix(c(a,b,c,d),ncol=2,byrow=T))

"Quies"->results[15,1]
test$p.value->results[15,2]
test$estimate->results[15,3]


class(results)
results<-as.data.frame(results)
names(results)<-c("ChromHMM_state","p_initial","OR_initial")

for (i in 2:5)
{results[,i]<-as.numeric(results[,i])}

#15 states: 0.05/15=0.003333333
index<-which(results$p_initial<0.05/15)
results$ChromHMM_state[index]
#"tssA"     
#"TssAFlnk" 
#"EnhG"     
#"Enh"      
#"ReprPCWk" (OR<1)
#"Quies"    (OR<1)

results_cpgs<-results
save(results_cpgs,file="CpGs_chromHMM.Rdata")

