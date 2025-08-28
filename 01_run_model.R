#run contmeQTLs in one cohort
#BeCOME cohort as example

require(snpStats)
library(Biobase)
library(data.table)

#####load genotypes
genotype<-read.plink("become_final_249.bed",
                     "become_final_249.bim",
                     "become_final_249.fam")
genotypes = as(genotype$genotypes, "numeric") # Genotype as 0|1|2

#####load phenotypes and methylation data
load("BetaSet_Become_final_249.rda")
beta<-exprs(Beta_Become_249)
pheno<-pData(Beta_Become_249)

#check that ID-order is the same for all datasets
#all.equal(colnames(beta),rownames(pheno)) #TRUE
#all.equal(rownames(genotypes),rownames(pheno)) #TRUE


#####list of meQTLs to calculate, based on Min et al., https://doi.org/10.1038/s41588-021-00923-x 
comb<-fread("../../01_meQTLs/cpg_snp.txt",sep="\t",header=T)

####get A1 and A2
bim<-fread("become_final_249.bim",sep="\t",header=F)

#####calculate GxE model
#this is set-up as array job in SLURM, 1.000.000 combinations are run in one batch 
#cov: PCs, BCCs, sex, age

slurm_arrayid   <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i               <- as.numeric(slurm_arrayid)
begin           = (i-1)*1000000+1
end             = (i-1)*1000000+1000000

filetitle_int<- paste0("become_results_ginte_", i,".txt")
write.table(cbind("CpG", "SNP", "A1","A2","AIC", "adj_R2", "N", "beta_g", "se_g", "p_g", "beta_e","se_e","p_e",
                  "beta_int", "se_int", "p_int"),
            file=filetitle_int, sep="\t", quote=F, row.names=F, col.names=F)


for (j in begin:end)
{cpg<-which(rownames(beta)==comb$cpg[j])
snp<-which(colnames(genotypes)==comb$snp[j])
if (length(snp)!=1 | length(cpg)!=1)  #if SNP or CpG not available => stop
{res<-NA}
else 
  
{
  
  #######interactive model: GxE##########
  result_int =         lm(beta[cpg,] ~ pheno$age + pheno$sex+
                            pheno$C1 + pheno$C2 +
                            pheno$C3 + 
                            beta["cg05575921",] +
                            pheno$CD8T + pheno$CD4T +
                            pheno$NK + pheno$Bcell +
                            pheno$Mono + pheno$Neu +
                            genotypes[, snp] * pheno$CTQ_sexemotphys_01_modandsev) # Meth ~ cov + SNP x CTQ
  summary_result_int = summary(result_int)
  index<-which(bim$V2 ==colnames(genotypes)[snp])
  
  #interaction could be NA if low samplesize for specific CA/SNP combination
  #specific number here depends on number of parameters in the model
  if (dim(summary_result_int$coefficients)[1]==16)
  {res_int =                         cbind(rownames(beta)[cpg], #CpG
                                           colnames(genotypes)[snp], #SNP
                                           bim$V5[index], #A1
                                           bim$V6[index], #A2
                                           AIC(result_int), #AIC
                                           summary_result_int$adj.r.squared, # adjusted R2
                                           length(result_int$residuals), # N
                                           summary_result_int$coefficients[14,1], # beta_SNP
                                           summary_result_int$coefficients[14,2], # se_SNP
                                           summary_result_int$coefficients[14,4], # p_snp
                                           summary_result_int$coefficients[15,1], # beta_E
                                           summary_result_int$coefficients[15,2], # se_E
                                           summary_result_int$coefficients[15,4], # p_E
                                           summary_result_int$coefficients[16,1], # beta_int
                                           summary_result_int$coefficients[16,2], # se_int
                                           summary_result_int$coefficients[16,4])} # p_int
  if (dim(summary_result_int$coefficient)[1]!=16)
  {res_int =                         cbind(rownames(beta)[cpg], #CpG
                                           colnames(genotypes)[snp], #SNP
                                           bim$V5[index], #A1
                                           bim$V6[index], #A2
                                           AIC(result_int), #AIC
                                           summary_result_int$adj.r.squared, # adjusted R2
                                           length(result_int$residuals), # N
                                           summary_result_int$coefficients[14,1], # beta_SNP
                                           summary_result_int$coefficients[14,2], # se_SNP
                                           summary_result_int$coefficients[14,4], # p_snp
                                           summary_result_int$coefficients[15,1], # beta_E
                                           summary_result_int$coefficients[15,2], # se_E
                                           summary_result_int$coefficients[15,4], # p_E
                                           "NA", # beta_int
                                           "NA", # se_int
                                           "NA")} # p_int
  
  write.table(res_int, file=filetitle_int, sep="\t", quote=F, row.names=F, col.names=F, append=T)     
  print(j)}
}
