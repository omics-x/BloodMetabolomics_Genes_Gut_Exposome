######################################################################### 
### Code for Univariate association analysis of medication use and metabolomics 
#########################################################################
### Load R packages 
library(foreign)
library(dplyr)
library("readxl")

### Prepare metabolomics and medication use dataset 
### Load Medication use data from RSIII-2
medication<-read.spss("e5_MEDICATION_(12-MAR-2015).sav",to.data.frame=T)
df<-as.data.frame(medication)
df2 = as.data.frame(sapply(df, as.numeric))
rownames(df2)<-df$ergoid
### Load metabolomics data
load("KNN_imputed_metabolon_data_3rdMarch.RData")
metabo_data<-as.data.frame(imputeddata)
colnames(metabo_data)<-paste("metab",colnames(metabo_data),sep="_")
metabo_data$sampleID<-rownames(metabo_data)
### create two datasets
medication<-df2[metabo_data$sampleID,]
medication<-medication[,-c(1:4)]
## convert to 0 and 1 (medication use)
medication[medication==1]<-0
medication[medication==2]<-1
## Remove the medication with less than 1 percent users in the dataset
percent0 <- apply(medication,2, function (x) { sum(x)/length(x) })
med_qc<-colnames(medication)[percent0>0.01]
final_medication<-medication[,med_qc]
write.csv(final_medication,file="medication_data_all_1percent.csv",row.names=T)
write.csv(metabo_data[,-992],file="metabolomics_data_ML_1percent.csv",row.names=T)

######################################################################################
### Regression analysis for the association of medication use and metabolomics 
######################################################################################
medication<-read.csv("medication_data_all_1percent.csv")
metabolom<-read.csv("metabolomics_data_ML_1percent.csv")
colnames(medication)[1]<-"ergoid"
colnames(metabolom)[1]<-"ergoid"
covariates<-readRDS(file ="Study_RSI_IV_RSIII_2_covars.rds")
### Merge medication use and metabolomics
merge<-merge(metabolom,medication,by.x='ergoid',by.y='ergoid',x.all=T)
merge<-merge(merge,covariates[,c("ergoid","sex","Age_blood_collection")],by.x='ergoid',by.y='ergoid',x.all=T)

results <- data.frame(
  Metabolite=as.character(),
  Medication=as.character(),
  Beta=as.numeric(),
  Se=as.numeric(),
  p=as.numeric(),
  n=as.numeric(),
  stringsAsFactors=FALSE)

medication_list<-colnames(medication)[-1]
metabo<-colnames(metabolom)[-1]
for(i in 1:length(metabo)){
  merge$pheno <- scale(merge[,metabo[i]], scale=TRUE, center=TRUE)
  for (j in 1:length(medication_list)) {
    fit <- lm(merge$pheno ~ merge[,medication_list[j]] + Age_blood_collection + sex, data=merge)
    fit.sum <- summary(fit)
    tablerow <- data.frame(
      Metabolite = metabo[i],
      Medication = medication_list[j],
      Beta = fit.sum$coefficients[2,1],
      Se = fit.sum$coefficients[2,2],
      p = fit.sum$coefficients[2,4],
      n = dim(merge)[1],
      stringsAsFactors=FALSE)
    results <- rbind(results, tablerow)
  }
}

### Apply FDR correction
results$FDR<-p.adjust(results[,5], method = 'fdr', n = length(results[,5]))
key<-read.table("medication_use.txt",sep='\t',head=T)
colnames(key)[1]<-'Medication'
### Annotation of the general cognition results
anno_com <- read_excel("/home/sahmad/metabolomics/gut_liver_brain/DUKE-0304-19ML+_CLIENT DATA TABLE.XLSX",sheet = 2)
anno_com<-anno_com%>%as.data.frame(anno_com)%>%mutate(Metabolite = paste0("metab_",anno_com$CHEM_ID))
res_anno<-anno_com%>%select("Metabolite","CHEMICAL_NAME","SUPER_PATHWAY","SUB_PATHWAY")%>%left_join(results,.,by="Metabolite")
res_anno2<-merge(res_anno,key,by='Medication',all.x='True')

### Write the results as CSV file
write.csv(res_anno2, file="Association_analysis_medication.csv")
