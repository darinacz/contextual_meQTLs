##train across all 4 cohorts where we have raw data available in adults/adolscence 
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

