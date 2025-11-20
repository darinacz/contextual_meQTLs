##use data across all 4 cohorts where we have raw data available in adults/adolscence 
setwd("/Users/darina/Contextual_meQTLs/UKBB/number_of_risk_alleles/")

library(Biobase)
library(dplyr)
library(tidyr)
library(readr)


#we will define the number of risk,i.e. most reactive, alleles, based on:
#for each SNP: compare 4 means for the specific DNAm
#1. A1/CA=0
#2. A1/CA=1
#3. A2/CA=0
#4. A2/CA=1
#homozygous individuals will be counted twice
#compare: abs(2-1) vs. abs (4-3)
#abs(2-1) > abs (4-3) => A1 as reactive allele
#abs(4-3) > abs (2-1) => A2 as reactive allele


#####1. get the SNPs and CpGs######
load("../PIRS/pgs_scoring.Rdata")  # this is the list of all SNPs after clumping plus the CpG-site they are regulating

####2. find the risk alleles####
#ALSPAC/LMU/BECOME/OPTIMA
#read in genotypes/DNAm data
#compute the respective means per cohort

#####ALSPAC#####
load("../../episcore/BetaSet_Alspac_final_616.rda")
beta<-exprs(Beta_Alspac_616)
index<-which(rownames(beta) %in% snp_cpg$CpG) 
beta_als<-beta[index,]
pt_als<-pData(Beta_Alspac_616)
rm(beta)
rm(Beta_Alspac_616)
geno<-read.plink("alspac_final_616.bed","alspac_final_616.bim", "alspac_final_616.fam")
genotypes = as(geno$genotypes, "numeric") 
index<-which(colnames(genotypes) %in% snp_cpg$SNP) 
geno_als<-genotypes[,index]

#these are the counts of copies of the second allele in the bim-file
snps_als<-read.table("alspac_final_616.bim", sep="\t",header=F)
index<-which(snps_als$V2 %in% snp_cpg$SNP)
snps_als<-snps_als[index,]
names(snps_als)<-c("CHR","SNP","CM","POS","A1","A2")

# --- Initialize results list
result_list <- vector("list", nrow(pair_list))

# --- Weighted mean 
weighted_mean <- function(values, weights) {
  if (all(is.na(values))) return(NA_real_)
  sum(values * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
}

for (i in seq_len(nrow(pair_list))) {
  
  snp <- pair_list$SNP[i]
  cpg <- pair_list$CpG[i]
  
  # --- Check SNP and CpG existence
  if (!(snp %in% colnames(geno_als)) | !(cpg %in% rownames(beta_als))) {
    result_list[[i]] <- tibble(
      SNP = snp,
      CpG = cpg,
      allele1 = NA_character_,
      allele2 = NA_character_,
      mean_phen0_A1 = NA_real_,
      mean_phen1_A1 = NA_real_,
      mean_phen0_A2 = NA_real_,
      mean_phen1_A2 = NA_real_,
      N_total = 0,
      delta_allele1 = NA_real_,
      delta_allele2 = NA_real_,
      allele_larger_change = NA_character_
    )
    next
  }
  
  # --- Get allele info from BIM
  bim_row <- snps_als %>% filter(SNP == snp)
  allele1 <- ifelse(nrow(bim_row) == 0, NA_character_, bim_row$A1)
  allele2 <- ifelse(nrow(bim_row) == 0, NA_character_, bim_row$A2)
  
  # --- Extract numeric vectors
  meth_vec <- as.numeric(beta_als[cpg, , drop = TRUE])
  geno_vec <- as.numeric(as.character(geno_als[, snp]))
  pheno_vec <- as.numeric(as.character(pt_als$status))
  
  df <- tibble(
    geno = geno_vec,
    meth = meth_vec,
    pheno = pheno_vec
  ) %>% filter(!is.na(meth))
  
  # Skip if no valid data
  if (nrow(df) == 0) {
    result_list[[i]] <- tibble(
      SNP = snp,
      CpG = cpg,
      allele1 = allele1,
      allele2 = allele2,
      mean_phen0_A1 = NA_real_,
      mean_phen1_A1 = NA_real_,
      mean_phen0_A2 = NA_real_,
      mean_phen1_A2 = NA_real_,
      N_total = 0,
      delta_allele1 = NA_real_,
      delta_allele2 = NA_real_,
      allele_larger_change = NA_character_
    )
    next
  }
  
  # --- Allele weights
  df <- df %>% mutate(
    w_A1 = ifelse(geno == 0, 2, 1),
    w_A2 = ifelse(geno == 2, 2, 1)
  )
  
  # --- Compute weighted means
  mean_A1 <- df %>% group_by(pheno) %>%
    summarise(mean_A1 = weighted_mean(meth, w_A1), .groups = "drop")
  mean_A2 <- df %>% group_by(pheno) %>%
    summarise(mean_A2 = weighted_mean(meth, w_A2), .groups = "drop")
  
  # --- Pivot to wide
  wide_df <- full_join(mean_A1, mean_A2, by = "pheno") %>%
    pivot_wider(
      names_from = pheno,
      values_from = c(mean_A1, mean_A2),
      names_glue = "mean_phen{pheno}_{.value}"
    )
  
  # --- Fix column names
  colnames(wide_df) <- gsub("mean_phen(\\d+)_mean_A1", "mean_phen\\1_A1", colnames(wide_df))
  colnames(wide_df) <- gsub("mean_phen(\\d+)_mean_A2", "mean_phen\\1_A2", colnames(wide_df))
  
  # --- Ensure all expected columns exist
  expected_cols <- c("mean_phen0_A1","mean_phen1_A1","mean_phen0_A2","mean_phen1_A2")
  for (col in expected_cols) if (!(col %in% names(wide_df))) wide_df[[col]] <- NA_real_
  
  # --- Add SNP, CpG, alleles, N
  wide_df <- wide_df %>%
    mutate(
      SNP = snp,
      CpG = cpg,
      allele1 = allele1,
      allele2 = allele2,
      N_total = nrow(df)
    ) %>%
    select(SNP, CpG, all_of(c("allele1","allele2")), any_of(expected_cols), N_total)
  
  # --- Compute Δ and larger allele change
  wide_df <- wide_df %>%
    mutate(
      delta_allele1 = abs(mean_phen1_A1 - mean_phen0_A1),
      delta_allele2 = abs(mean_phen1_A2 - mean_phen0_A2),
      allele_larger_change = case_when(
        delta_allele1 > delta_allele2 ~ allele1,
        delta_allele2 > delta_allele1 ~ allele2,
        delta_allele1 == delta_allele2 ~ "equal"
      )
    )
  
  result_list[[i]] <- wide_df
}

# --- Combine all results
results_all <- bind_rows(result_list)

results_als<-results_all
results_als<-as.data.frame(results_als)
names(results_als)<-c("SNP","CpG","A1","A2","mean_noCA_A1_als","mean_CA_A1_als","mean_noCA_A2_als","mean_CA_A2_als", "N_total_als","delta_A1_als", "delta_A2_als",
                      "allele_larger_change_als")
save(results_als,file="results_ALSPAC.RData")

#run the same for the other cohorts and merge all in one dataframe
#take weighted mean across all cohorts
# Define the cohort-specific columns
cohort_deltas1 <- c("delta_A1_als", "delta_A1_bec", "delta_A1_lmu", "delta_A1_opt")
cohort_deltas2 <- c("delta_A2_als", "delta_A2_bec", "delta_A2_lmu", "delta_A2_opt")
cohort_N <- c("N_total_als", "N_total_bec", "N_total_lmu", "N_total_opt")

results_weighted <- all %>%
  mutate(
    # Mask N_total where delta is NA
    N1_masked = as.data.frame(across(all_of(cohort_N))) * (!is.na(across(all_of(cohort_deltas1)))),
    N2_masked = as.data.frame(across(all_of(cohort_N))) * (!is.na(across(all_of(cohort_deltas2)))),
    
    
    # Weighted mean delta for allele 1
    delta_allele1_weighted = rowSums(across(all_of(cohort_deltas1)) * N1_masked, na.rm = TRUE) /
      rowSums(N1_masked, na.rm = TRUE),
    
    # Weighted mean delta for allele 2
    delta_allele2_weighted = rowSums(across(all_of(cohort_deltas2)) * N2_masked, na.rm = TRUE) /
      rowSums(N2_masked, na.rm = TRUE),
    
    # Determine which allele has the larger weighted delta
    allele_larger_change_weighted = case_when(
      delta_allele1_weighted > delta_allele2_weighted ~ all$A1,
      delta_allele2_weighted > delta_allele1_weighted ~ all$A2,
      delta_allele1_weighted == delta_allele2_weighted ~ "equal",
      TRUE ~ NA_character_
    ) 
  ) %>%
  select(-N1_masked, -N2_masked)

names(results_weighted)[c(37:39)]<-c("delta_A1_weighted","delta_A2_weighted","allele_larger_change_weighted")

save(results_weighted,file="results_weighted.RData")

#score:
write.table(results_weighted, "score_largest_diff.txt",sep="\t",quote=F,row.names=F)









##based on meta CpGs
meta<-read.table("meta.txt",sep="\t",header=T) #this is the full results list
meta<-unique(meta$CpG) 

load("BetaSet_Alspac_final_616.rda")
beta_als<-exprs(Beta_Alspac_616)
index<-which(rownames(beta_als) %in% meta) #n=4,394
beta_als<-beta_als[index,]
pt_als<-pData(Beta_Alspac_616)

load("BetaSet_LMU_final_551.rda")
beta_lmu<-exprs(Beta_LMU_551)
index<-which(rownames(beta_lmu) %in% meta) #n=4,943
beta_lmu<-beta_lmu[index,]
pt_lmu<-pData(Beta_LMU_551)

load("BetaSet_Become_final_249.rda")
beta_bec<-exprs(Beta_Become_249)
index<-which(rownames(beta_bec) %in% meta) #n=4,693
beta_bec<-beta_bec[index,]
pt_bec<-pData(Beta_Become_249)

load("BetaSet_Optima_final_96.rda")
beta_opt<-exprs(Beta_Optima_96)
index<-which(rownames(beta_opt) %in% meta) #n=4,693
beta_opt<-beta_opt[index,]
pt_opt<-pData(Beta_Optima_96)

#subset to CpGs available in all 4 cohorts
all.equal(rownames(beta_k2h),rownames(beta_opt)) #TRUE
all.equal(rownames(beta_k2h),rownames(beta_bec)) #TRUE
all.equal(rownames(beta_k2h),rownames(beta_lmu)) #TRUE
all.equal(rownames(beta_k2h),rownames(beta_als)) #TRUE

#combine all in one dataset 
beta_all<-cbind(beta_als,beta_lmu,beta_bec,beta_opt)
study<-rbind(cbind(colnames(beta_als),"ALSPAC"), cbind(colnames(beta_lmu),"LMU"),
             cbind(colnames(beta_bec),"BECOME"), cbind(colnames(beta_opt),"OPTIMA"))
all.equal(colnames(beta_all),study[,1]) #TRUE
study<-as.data.frame(study)
rownames(study)<-study$V1
names(study)<-c("ID","study")

CA<-rbind(cbind(rownames(pt_als),as.numeric(pt_als$CTQ_sexemotphys_01_modandsev)), cbind(rownames(pt_lmu),as.numeric(pt_lmu$CTQ_sexemotphys_01_modandsev)),
          cbind(rownames(pt_bec),as.numeric(pt_bec$CTQ_sexemotphys_01_modandsev)), cbind(rownames(pt_opt),as.numeric(pt_opt$CTQ_sexemotphys_01_modandsev)))
CA<-as.data.frame(CA)
names(CA)<-c("ID","CA")
all.equal(CA$ID, colnames(beta_all)) #TRUE

beta_resid<-beta_all
for (i in 1: dim(beta_resid)[1]) #residualize for cohort, to avoid cohort-specific effects
{beta_resid[i,]<- summary(lm(beta_all[i,]~as.factor(study[,2])))$residuals
print(i)}

#beta-values should be between 0 and 1 => reshift residualized beta-values to [0,1]
beta_resid2<-beta_resid
for (i in 1: dim(beta_resid2)[1])
{ if (min(beta_resid[i,]) < 0)
{beta_resid2[i,] <- beta_resid2[i,] + abs(min(beta_resid2[i,]))}}

#run LASSO regression
#predictors: CpGs
#predicted variable CA yes/no
#10-fold CV
#CV based on misclassification rate
misclass <- NULL
for (i in 1:100){
  cv <- cv.glmnet(t(beta_resid2),CA$CA, family="binomial", alpha=1, nfolds=10, type.measure="class")  
  misclass <- cbind(misclass, cv$cvm)
  print(i)}

#pick model with lowest cvm
rownames(misclass) <- cv$lambda
lambda.min <- as.numeric(names(which.min(rowMeans(misclass))))
save(misclass,file="cv_runs_4_cohorts.Rdata")
#run glmnet with the best lambda
fit <- glmnet(t(beta_resid2),CA$CA, family="binomial", alpha=1 ,lambda=lambda.min, type.measure="class", intercept=T)
fit
#Df %Dev  Lambda
#13  1.83 0.03552
coe[rownames(coef(fit, s = 'lambda.min'))[coef(fit, s = 'lambda.min')[,1]!= 0],]

#Intercept	-0.08882592
#cg04154027	0.19927924
#cg05809437	-0.01095841
#cg07201456	0.23763807
#cg08976526	-2.65959496 
#cg12094174	-0.30000245
#cg13404038	0.11755069
#cg13494126	-1.00065040
#cg16474684	-0.26653293 
#cg20252067	-0.18475204    
#cg23242489	 -0.42120220
#cg23513447	 -3.42134188
#cg24375209	  -0.45848415
#cg26288502	  2.59109143

