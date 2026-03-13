library(car)
library(sandwich)
library(glmnet)

load("df_proinfl.Rdata") #this is the dataframe containing the UKBB-sample including the standardized OLINK markers, phenotypes and the PIRS

#we focus on the 31 pro-inflammatory markers
pro <- c(
  "il18","ccl13","ccl11","tnfsf14","ccl4","ccl19","cxcl5",
  "ccl7","il6","cxcl1","il17c","il17a","cxcl11","cxcl9",
  "il1a","osm","ccl2","cxcl8","cxcl10","cxcl6","s100a12",
  "ccl3","ccl23","tnf","il12b","cd40","ifng","lta","ccl20",
  "ccl25","ccl8"
) #n=31

#1. MANOVA: test is the overall pro-inflammatory profile is associated with PIRS
#we use a sandwich estimator as variances between CA and no CA are different

index<-which(colnames(df_final) %in% pro)
mlm_model <- lm(manova(as.matrix(df_final[,index]) ~ as.factor(df_final$sex) + df_final$PC1 + df_final$PC2 + df_final$PC3 + 
                         df_final$PC4 + df_final$PC5 + df_final$BMI+ df_final$age+ df_final$SCORE1_AVG * df_final$any_abuse))

#robust covariance
robust_cov <- vcovHC(mlm_model, type = "HC3") # heteroskedasticity-consistent
#significance test
Manova(mlm_model, test = "Wilks", vcov = robust_cov) 

#2. Are there significant differences between the CA and no CA groups?
#run single marker associations and compare effect sizes

table(df_final$any_abuse)
index<-which(df_final$any_abuse==0)
df0<-df_final[index,]

results<-matrix(ncol=4,nrow=92)
for (i in 81:172) #markers are in columns 81 -172
{
  mod<-summary(lm(df0[,i]~ as.factor(df0$sex) + df0$PC1 +df0$PC2 + df0$PC3 + df0$PC4 + df0$PC5 + df0$BMI+ df0$age+ df0$SCORE1_AVG)) 
  results[i-80,1]<-colnames(df_final)[i]
  results[i-80,2]<-mod$coefficients[10,1]
  results[i-80,3]<-mod$coefficients[10,2]
  results[i-80,4]<-mod$coefficients[10,4]}

results<-as.data.frame(results)
names(results)<-c("marker","beta","se","p")
for (i in 2:4)
{results[,i]<-as.numeric(results[,i])}

index<-which(results$marker %in% pro) #subset to the 31 markers of interest
results[index,]->res_PIRS_pro_noCA

table(df_final$any_abuse)
index<-which(df_final$any_abuse==1)
df1<-df_final[index,]

results<-matrix(ncol=4,nrow=92)
for (i in 81:172) #markers are in columns 81 -172
{
  mod<-summary(lm(df1[,i]~ as.factor(df1$sex) + df1$PC1 +df1$PC2 + df1$PC3 + df1$PC4 + df1$PC5 + df1$BMI+ df1$age+ df1$SCORE1_AVG)) 
  results[i-80,1]<-colnames(df_final)[i]
  results[i-80,2]<-mod$coefficients[10,1]
  results[i-80,3]<-mod$coefficients[10,2]
  results[i-80,4]<-mod$coefficients[10,4]}

results<-as.data.frame(results)
names(results)<-c("marker","beta","se","p")
for (i in 2:4)
{results[,i]<-as.numeric(results[,i])}


index<-which(results$marker %in% pro) #subset to the 31 markers of interest
results[index,]->res_PIRS_pro_CA

#combine all in one dataframe
all.equal(res_PIRS_noCA$marker,res_PIRS_CA$marker) #TRUE

all<-rbind(cbind(res_PIRS_noCA$marker,res_PIRS_noCA$beta,res_PIRS_noCA$p,rep("noCA",31)),
           cbind(res_PIRS_CA$marker,res_PIRS_CA$beta,res_PIRS_CA$p,rep("CA",31)))
all<-as.data.frame(all)
names(all)<-c("marker","beta","p","group")
all$beta<-as.numeric(all$beta)
all$p<-as.numeric(all$p)

wilcox.test(all$beta[all$group=="CA"],all$beta[all$group=="noCA"],alternative="greater")

#3.Define shift of the profile based on composite score
#31 pro markers as predictors
index<-which(colnames(df0) %in% pro)
X0 <- as.matrix(df0[, index]) 
index<-which(colnames(df1) %in% pro)
X1 <- as.matrix(df1[, index]) 

#PIRS
y0 <- df0$score
y1 <- df1$score

#fit, based on Ridge regression, so alpha=0
fit0 <- cv.glmnet(X0, y0, alpha = 0)
fit1 <- cv.glmnet(X1, y1, alpha = 0)


#build composite score based on min. CV-error
beta0 <- as.vector(coef(fit0, s = "lambda.min"))[-1]
beta1 <- as.vector(coef(fit1, s = "lambda.min"))[-1]

#add to data frame
index<-which(colnames(df0) %in% pro)
df0$compscore <- as.matrix(df0[,index]) %*% beta0
index<-which(colnames(df1) %in% pro)
df1$compscore <- as.matrix(df1[,index]) %*% beta1

#4. Test for association with diagnoses
summary(glm((df0$any_dx>=1)~df0$age + df0$sex + df0$PC1 + df0$PC2 + df0$PC3 + df0$PC4 + df0$PC5 + df0$BMI + df0$compscore,family="binomial")) 
summary(lm(df0$any_dx~df0$age + df0$sex + df0$PC1 + df0$PC2 + df0$PC3 + df0$PC4 + df0$PC5 + df0$BMI + df0$compscore))

summary(glm((df1$any_dx>=1)~df1$age + df1$sex + df1$PC1 + df1$PC2 + df1$PC3 + df1$PC4 + df1$PC5 + df1$BMI + df1$compscore,family="binomial")) 
summary(lm(df1$any_dx~df1$age + df1$sex + df1$PC1 + df1$PC2 + df1$PC3 + df1$PC4 + df1$PC5 + df1$BMI + df1$compscore))






