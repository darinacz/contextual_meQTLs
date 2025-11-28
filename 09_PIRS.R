##use data across all 4 cohorts where we have raw data available in adults/adolscence 

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
load("pgs_scoring.Rdata")  # this is the list of all SNPs after clumping plus the CpG-site they are regulating, i.e. interacting with CA on this site

####2. find the risk alleles####
#ALSPAC/LMU/BECOME/OPTIMA
#read in genotypes/DNAm data
#compute the respective means per cohort
#ALSPAC as example
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
results_als <- bind_rows(result_list)

results_als<-as.data.frame(results_als)
names(results_als)<-c("SNP","CpG","A1","A2","mean_noCA_A1_als","mean_CA_A1_als","mean_noCA_A2_als","mean_CA_A2_als", "N_total_als","delta_A1_als", "delta_A2_als",
                      "allele_larger_change_als")
save(results_als,file="results_ALSPAC.RData")

#run the same for the other cohorts and merge all into one dataframe
#take weighted mean across all cohorts
# Define the cohort-specific columns
cohort_deltas1 <- c("delta_A1_als", "delta_A1_bec", "delta_A1_bio", "delta_A1_opt")
cohort_deltas2 <- c("delta_A2_als", "delta_A2_bec", "delta_A2_bio", "delta_A2_opt")
cohort_N <- c("N_total_als", "N_total_bec", "N_total_bio", "N_total_opt")

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















