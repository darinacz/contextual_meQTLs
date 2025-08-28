#individual cohort results were split into chunks of 1 million combinations based on input list of CpG-SNP combinations
#to ensure that the same combinations are present in the same chunk across cohorts
#meta-analysis is carried out for chunks 1-10, 11-20 etc. as otherwise cpu memory is exceeded
#example here for chunks 1-10

######1-10########

#run METAL
/generic-metal/metal
SCHEME   STDERR
SEPARATOR  TAB

#set-up files
MARKER   CpG_SNP #this is the specific CpG-SNP combination
ALLELE   A1 A2
EFFECT   beta_int
STDERR   se_int
PVAL     p_int

PROCESS ../ALSPAC/alspac_GintE_1_10_final.txt
PROCESS ../BECOME/become_GintE_1_10_final.txt 
PROCESS ../KORA/kora_GintE_1_10_final.txt
PROCESS ../LMU/lmu_GintE_1_10_final.txt
PROCESS ../OPTIMA/optima_GintE_1_10_final.txt
PROCESS ../SHIP/ship_GintE_1_10_final.txt


#Execute random-effects meta-analysis and compute heterogeneity
OUTFILE meta_1_10_GintE_random .tbl
ANALYZE RANDOM

QUIT


