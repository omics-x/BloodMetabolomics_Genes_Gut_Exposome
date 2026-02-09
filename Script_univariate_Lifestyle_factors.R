######################################################################### 
### Script for univariate association analysis of lifestyle factors with metabolomics 
#########################################################################
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library("haven")
library(dplyr)
library("readxl")
######################################################################### 
### Data preparation for univariate association analysis of lifestyle and clinical factors with metabolomics 
######################################################################### 

### Load Combined covariate information file for RS1-IV and RSIII-2
covariates<-readRDS(file ="Study_RSI_IV_RSIII_2_covars.rds")
oh1_5 <- read_sav("e5_intvw_Alcoholperday_05-11-2019.sav")
alc_e5<-as.data.frame(oh1_5)
alc_rsIII_2<-alc_e5[which(alc_e5$rs_cohort==3),]
### Load imputed metabolomics dataset 
load("KNN_imputed_metabolon_data_3rdMarch.RData")
metabo_data<-as.data.frame(imputeddata)
colnames(metabo_data)<-paste("metab",colnames(metabo_data),sep="_")
metabo_data$sampleID<-rownames(metabo_data)
### extract variables of interest
demographics<-covariates[covariates$ergoid%in%metabo_data$sampleID,]
demographics<-demographics[,c("ergoid","Age_blood_collection","sex","education","BMI","SysBP","DiasBP","Hypertension","Diabetes2","Smoke")]
demographics2<-merge(demographics,alc_rsIII_2[,c("ergoid","e5_Alc_tot")],by='ergoid',all.x="True")
demographics2<-demographics2[complete.cases(demographics2), ]
rownames(demographics2)<-demographics2$ergoid
### Keep metabolomics dataset for which lifestyle and clinical variable information is available 
metabo_data<-metabo_data[metabo_data$sampleID%in%demographics2$ergoid,]
demographics3<-demographics2[metabo_data$sampleID,]

demographics3$education_label<-NA
demographics3$education_label[which(demographics3$education=="primary_education")]<-0
demographics3$education_label[which(demographics3$education=="further_education")]<-1
demographics3$education_label[which(demographics3$education=="higher_education")]<-2
demographics3$education_label<-as.factor(demographics3$education_label)

demographics3$Smoke_label<-NA
demographics3$Smoke_label[which(demographics3$Smoke=="never")]<-0
demographics3$Smoke_label[which(demographics3$Smoke=="current")]<-2
demographics3$Smoke_label[which(demographics3$Smoke=="former")]<-1
demographics3$Smoke_label<-as.factor(demographics3$Smoke_label)

demographics3$Sex_label<-NA
demographics3$Sex_label[which(demographics3$sex=="male")]<-0
demographics3$Sex_label[which(demographics3$sex=="female")]<-1
demographics3$Sex_label<-as.factor(demographics3$Sex_label)

demographics3$Hypertension_label<-NA
demographics3$Hypertension_label[which(demographics3$Hypertension=="yes")]<-1
demographics3$Hypertension_label[which(demographics3$Hypertension=="no")]<-0
demographics3$Hypertension_label<-as.factor(demographics3$Hypertension_label)

demographics3$Alcohol_gm_d<-demographics3$e5_Alc_tot

### comorbidity contains information about clinical variables 
comorbidity<-demographics3[,c("Hypertension_label","SysBP","DiasBP","Diabetes2")]
### demographic file contains information about lifestyle factors 
age_sex_lifestyle<-demographics3[,c("education_label","Smoke_label","Age_blood_collection","Sex_label","Alcohol_gm_d","BMI")]

### Save the files 
write.csv(comorbidity,file="comorbidity_jun_2022.csv",row.names=T)
write.csv(age_sex_lifestyle,file="Demographic_jun_2022.csv",row.names=T)
write.csv(metabo_data[,-992],file="metabolomics_data_comorbidity_jun_2022.csv",row.names=T)

######################################################################### 
### Regression analysis 
#########################################################################

demographics<-read.csv("Demographic_jun_2022.csv")
metabolom<-read.csv("metabolomics_data_comorbidity_jun_2022.csv")
colnames(demographics)[1]<-"ergoid"
colnames(metabolom)[1]<-"ergoid"


merge<-merge(metabolom,demographics,by.x='ergoid',by.y='ergoid',x.all=T)

results <- data.frame(
  Metabolite=as.character(),
  Demographics=as.character(),
  Beta=as.numeric(),
  Se=as.numeric(),
  p=as.numeric(),
  n=as.numeric(),
  stringsAsFactors=FALSE)

demo_list<-colnames(demographics)[c(2,3,6,7)]
metabo<-colnames(metabolom)[-1]


for(i in 1:length(metabo)){
  merge$pheno <- scale(merge[,metabo[i]], scale=TRUE, center=TRUE)
  for (j in 1:length(demo_list)) {
    fit <- lm(merge$pheno ~ merge[,demo_list[j]] + Age_blood_collection + Sex_label, data=merge)
    fit.sum <- summary(fit)
    tablerow <- data.frame(
      Metabolite = metabo[i],
      Demographics = demo_list[j],
      Beta = fit.sum$coefficients[2,1],
      Se = fit.sum$coefficients[2,2],
      p = fit.sum$coefficients[2,4],
      n = dim(merge)[1],
      stringsAsFactors=FALSE)
    results <- rbind(results, tablerow)
  }
}
results$FDR<-p.adjust(results[,5], method = 'fdr', n = length(results[,5]))
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-anno_com%>%as.data.frame(anno_com)%>%mutate(Metabolite = paste0("metab_",anno_com$CHEM_ID))
res_anno<-anno_com%>%select("Metabolite","CHEMICAL_NAME","SUPER_PATHWAY","SUB_PATHWAY")%>%left_join(results,.,by="Metabolite")

### Save results of association of lifstyle factors with metabolomics 
write.csv(res_anno, file="Association_analysis_demographics.csv")

