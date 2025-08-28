#compute eQTMs 

library(biomaRt)
library(Biobase)
library(MatrixEQTL)

#prepare data for eQTM analysis 
#expression data
load("gex.Rdata") #gex, one column per ID, one row per gene

#we need start and end position for these genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_location<-getBM(attributes =c('hgnc_symbol','ensembl_gene_id','chromosome_name','start_position','end_position'),mart = ensembl,
                     values = rownames(gex))
index<-which(rownames(gex) %in% gene_location$ensembl_gene_id) #n=12,117
gex<-gex[index,]
index<-which(gene_location$ensembl_gene_id %in% rownames(gex)) #n=12,118
gene_locations<-gene_location[index,c(2,3,4,5)] #geneid chr s1 s2 
index<-which(gene_locations$chr=="GL000219.1" | gene_locations$chr=="KI270711.1" | gene_locations$chr=="KI270713.1" | gene_locations$chr=="MT")
gene_locations<-gene_locations[-index,]
gene_locations$chromosome_name<-paste("chr",gene_locations$chromosome_name,sep="")
gene_locations$chromosome_name[gene_locations$chromosome_name=="chrX"]<-"chr23"
index<-which(duplicated(gene_locations$ensembl_gene_id)) #ENSG00000280739
gene_locations<-gene_locations[-index,] #n=12,091 genes
index<-which(rownames(gex) %in% gene_locations$ensembl_gene_id) #n=12,091
gex<-gex[index,]
all.equal(rownames(gex),gene_locations$ensembl_gene_id) #FALSE
gex<-gex[order(rownames(gex)),]
gene_locations<-gene_locations[order(gene_locations$ensembl_gene_id),]
all.equal(rownames(gex),gene_locations$ensembl_gene_id) ##TRUE
names(gene_locations)<-c("gene","chr","s1","s2")
gene_locations$chr<-gsub("chr","",gene_locations$chr)

#methylation data
load("BetaSet_final.rda")
beta<-exprs(BetaSet)

#annotation of CpGs
load("/Users/darina/annot_epic.Rdata")
index<-which(annot_epic$Name %in% rownames(beta)) 
cpg_locations<-annot_epic[index,c("Name", "chr", "pos")]
beta<-beta[order(rownames(beta)),]
cpg_locations<-cpg_locations[order(rownames(cpg_locations)),]
all.equal(rownames(beta),rownames(cpg_locations)) #TRUE
names(cpg_locations)<-c("cpg","chr","pos")
cpg_locations$chr<-gsub("chr","",cpg_locations$chr)

#overlap of methylation and gex
index<-which(colnames(gex) %in% colnames(beta)) 
gex<-gex[,index]
index<-which(colnames(beta) %in% colnames(gex)) 
beta<-beta[,index]
beta<-beta[,order(colnames(beta))]
gex<-gex[,order(colnames(gex))]
all.equal(colnames(beta),colnames(gex)) #TRUE

#phenotypes/covariates
cov<-pData(BetaSet) 

#save data
save(gex,file="gex_for_MatrixEQTL.RData")
write.table(gex, file="gex_for_MatrixEQTL.txt")
save(beta,file="DNAm_for_MatrixEQTL.RData")
write.table(beta, file="DNAm_for_MatrixEQTL.txt")
save(cov,file="cov_for_MatrixEQTL.RData")
write.table(cov, file="cov_for_MatrixEQTL.txt")
write.table(gene_locations,file="genes_for_MatrixEQTL.txt")
write.table(cpg_locations,file="cpg_for_MatrixEQTL.txt")

#####run MatrixEQTL#####
base_dir <- getwd()

#set thresholds and parameters
cis_threshold <- 1 #cis-meqtlS cutoff=1, as we want to compare to results from KORA
cis_dist <- 5e5 #Distance for local (cis) gene-SNP pairs: cis window of 500kb
useModel = modelLINEAR; #linear model to model the effect of the methylation as additive linear and test for its significance using t-statistic.
errorCovariance = numeric() # Error covariance matrix

#set file names and paths
expression_file_name <- "gex_for_MatrixEQTL.txt"
methylation_file_name <- "DNAm_for_MatrixEQTL.txt";
covariates_file_name <- "cov_for_MatrixEQTL.txt";
gene_location_file_name <- "genes_for_MatrixEQTL.txt";
cpg_location_file_name <- "cpg_for_MatrixEQTL.txt";

methylation_file_path <- file.path(base_dir, methylation_file_name)
expression_file_path <- file.path(base_dir, expression_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
gene_location_file_path <- file.path(base_dir, gene_location_file_name)
cpg_location_file_path <- file.path(base_dir, cpg_location_file_name)

output_file_name_cis = "eqtm_cis.txt"

## Load methylation data
meth <- SlicedData$new()
meth$fileDelimiter = " " # the blank character
meth$fileOmitCharacters = "NA" # denote missing values
meth$fileSkipRows = 1 # one row of column labels
meth$fileSkipColumns = 1 # one column of row labels
meth$fileSliceSize = 2000 # read file in slices of 2,000 rows
meth$LoadFile(methylation_file_name)

## Load gene expression data
gene = SlicedData$new()
gene$fileDelimiter = " " # the blank character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1 # one row of column labels
gene$fileSkipColumns = 1 # one column of row labels
gene$fileSliceSize = 2000 # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

## Load covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = " "      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1      # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}

cpgpos = read.table(cpg_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

##run analysis
me_eqtm_become = Matrix_eQTL_main(
  snps = meth, # SlicedData object with methylation information. 
  gene = gene, # SlicedData object with gene expression information. 
  cvrt = cvrt, # SlicedData object with additional covariates.
  snpspos = cpgpos, # must have 3 columns - SNP name, chromosome, and position
  genepos = genepos, # data.frame with information about transcript locations, must have 4 columns - the name, chromosome, and positions of the left and right ends
  output_file_name.cis = output_file_name_cis, # local associations are saved to this file
  pvOutputThreshold.cis = cis_threshold, # cis threshold for reporting
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = FALSE,
  cisDist = cis_dist, 
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE, 
  noFDRsaveMemory = FALSE)

