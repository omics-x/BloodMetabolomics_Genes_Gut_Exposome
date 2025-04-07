.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library ("foreign")
library("Hmisc")
library("haven")
library("dplyr")
library("readxl")
## White matter lesions 
load("/home/sahmad/charge/wml_metabolomics/mri_data/dataset_WMH_dti_20200602.RData")
## brain volume data is same as we have already The variable Total_par is not in milliliters yet. If you divide it by a 1000 then you have total brain volume in ml. 
brain_size<- read_spss("/home/sahmad/metabolomics/tamo_project/association_mri/mri_dataset/RSS volumes archived 20160413 - 20160428.zsav")
brain_size<-as.data.frame(brain_size)
brain_size<-brain_size[,c(1:19)]
merg1<-merge(dataset_selection,brain_size,by.x="bigrfullname",by.y="BIGRfullname",all.x=T)

### Hippocampal volume 
hip<-read.csv("/home/sahmad/metabolomics/tamo_project/association_mri/mri_dataset/aseg.vol.table.csv")
hip_qc<-read.csv("/home/sahmad/metabolomics/tamo_project/association_mri/mri_dataset/autoQA_data_ergo.csv")
## select scanes with good quality
select<-hip_qc[which(hip_qc$t1_auto_qa_grad>=310),]
hip_good<-hip[which(hip$bigrfullname%in%select$bigrfullname),]
hip_good2<-hip_good[,c("bigrfullname","Left.Hippocampus","Right.Hippocampus","EstimatedTotalIntraCranialVol")]
hip_good2$total_Hippocampus<-hip_good2$Left.Hippocampus+hip_good2$Right.Hippocampus
MRIlong<-merge(merg1,hip_good2,by.x="bigrfullname",by.y="bigrfullname",all.x=T)
## Total brain volume is not in milliletter therefore
MRIlong$Total_par_ml<-MRIlong$Total_par/1000

###############################################################################################
### Metabolomics data
###############################################################################################
load("/home/sahmad/metabolomics/gut_liver_brain/quality_control_2021/KNN_imputed_metabolon_data_3rdMarch.RData")
imputeddata<-as.data.frame(imputeddata)
colnames(imputeddata)<-paste("metab",colnames(imputeddata),sep="_")
metabolites<-colnames(imputeddata)
imputeddata$ergoid<-rownames(imputeddata)

############### Blood collection data 
blood_date<- read.spss("/home/sahmad/charge/wml_metabolomics/nightingale/rs1_5_new_data/e5_(5)_BLDAFNAM_(10-jul-2017).sav",to.data.frame=T,as.factor=F)
blood_date<-blood_date[blood_date$cohort=="RS-III",]
RSIII<-merge(imputeddata,blood_date[,c("ergoid","e5_2686")],by.all="ergoid",all.x=T)

####### MRI merging
spss2date <- function(x) as.Date(x/86400, origin = "1582-10-14")
count_measurements <- function(data,ergoid){
  data$measurement_number <- rep(NA, length(data[[ergoid]]))
  for (i in 1:length(data[[ergoid]])) {
    if (i == 1) {
      data$measurement_number[i] <- 1
    }
    if (i > 1) {
      if (data[[ergoid]][i] == data[[ergoid]][i - 1]) {
        data$measurement_number[i] <- data$measurement_number[i - 1] + 1
      } else {
        data$measurement_number[i] <- 1
      }
    }
  }
  data
}

### chosing the closet scan function
select_closest_scan <- function(dataset, date1,date2){
  dataset$difference <- as.numeric(unlist(dataset[date2]-dataset[date1]))/365.25
  dataset$abs_difference<- abs(dataset$difference)
  dataset$select_scan <- NA
  for (i in unique(dataset$ergoid)){
    differences <- dataset[dataset$ergoid==i,"abs_difference"]
    dataset$select_scan[dataset$ergoid==i]<- as.numeric(differences==min(differences))
  }
  dataset
}


## 

MRIlong<-MRIlong[which(MRIlong$dementia_at_scan%nin%c(1) & MRIlong$stroke_at_scan%nin%c(1)),]
MRIlong$bigrfullname <- trimws(MRIlong$bigrfullname, "right")
# date of blood collection/Metabolomics data created above
determinant_date<-RSIII
determinant_date$blood_date <- spss2date(determinant_date$e5_2686)
merged <- merge(determinant_date,MRIlong,by="ergoid")
## remove NA
merged <- merged[!is.na(merged[["blood_date"]]),]
merged <- merged[!is.na(merged[["scandate"]]),]
## calculate
merged_nearest <- select_closest_scan(merged,"blood_date","scandate")
# now I only keep the scans that were closest to the bloodsample:
merc <- merged_nearest[merged_nearest$select_scan==1,]

write.table(merc,file="RS1_5_Metabolon_MRIdata_2Jun2021.txt", sep="\t", quote=F, row.names=T, col.names=T)

#######################################
####################################### Covariates files 
covariates<-readRDS(file ="/home/sahmad/metabolomics/gut_liver_brain/quality_control_2021/association_v2_age_corrected/covariat_files/Study_RSI_IV_RSIII_2_covars.rds")
m<-merge(merc,covariates,by.all="ergoid",all.x=T)
m<-m[m$sex=='male']
#########################################################
#########################################################
results <- data.frame(
            endo_pheno=as.character(),
            Metabolite=as.character(),
            Beta=as.numeric(),
            Se=as.numeric(),
            p=as.numeric(),
            n=as.numeric(),
            lower=as.numeric(),
            upper=as.numeric(),
            stringsAsFactors=FALSE)

m$log_Total_wml<-log(m$Total_wml)
endo_pheno<-c("total_Hippocampus","Total_par_ml","Total_wml")

for (j in 1:length(endo_pheno)) {
for ( i in 1:length(metabolites)){

m$metabo<-scale(m[,metabolites[i]],center = TRUE, scale = TRUE)
m$PD_model<-scale(log(m[,endo_pheno[j]]),center = TRUE, scale = TRUE)
test<- summary(lm (paste("PD_model ~ metabo + Age_blood_collection + sex + BMI + Lipilower + ICV_from_mask",sep=""), data=m))
#test<- summary(lm (paste("PD_model ~ metabo + Age_blood_collection + sex + BMI + Lipilower + Smoke + Diabetes2 + Hypertension + ICV_from_mask",sep=""), data=m))
tablerow <- data.frame(
            endo_pheno=endo_pheno[j],
            Metabolite=metabolites[i],
            Beta=test$coefficients[2,1],
            Se=test$coefficients[2,2],
            p=test$coefficients[2,4],
            n=dim(na.omit(m[,c("metabo","PD_model")]))[[1]],
            lower=test$coef[2,1] - qt(0.975, df = test$df[2]) * test$coef[2, 2],
            upper=test$coef[2,1] + qt(0.975, df = test$df[2]) * test$coef[2, 2],
            stringsAsFactors=FALSE)
results <- rbind(results, tablerow)
}
}
head (results[order(results$p),])
results$FDR<-p.adjust(results[,5], method = 'fdr', n = length(results[,4]))

write.table(results,file="Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",sep="\t",row.names=F,col.names=T,quote=F)
########################################################
# Sex stratified analysis 
########################################################
m<-merge(merc,covariates,by.all="ergoid",all.x=T)
m<-m[m$sex=='female',]
#########################################################
#########################################################
results <- data.frame(
            endo_pheno=as.character(),
            Metabolite=as.character(),
            Beta=as.numeric(),
            Se=as.numeric(),
            p=as.numeric(),
            n=as.numeric(),
            lower=as.numeric(),
            upper=as.numeric(),
            stringsAsFactors=FALSE)

m$log_Total_wml<-log(m$Total_wml)
endo_pheno<-c("total_Hippocampus","Total_par_ml","Total_wml")

for (j in 1:length(endo_pheno)) {
for ( i in 1:length(metabolites)){

m$metabo<-scale(m[,metabolites[i]],center = TRUE, scale = TRUE)
m$PD_model<-scale(log(m[,endo_pheno[j]]),center = TRUE, scale = TRUE)
test<- summary(lm (paste("PD_model ~ metabo + Age_blood_collection + BMI + Lipilower + ICV_from_mask",sep=""), data=m))
#test<- summary(lm (paste("PD_model ~ metabo + Age_blood_collection + sex + BMI + Lipilower + Smoke + Diabetes2 + Hypertension + ICV_from_mask",sep=""), data=m))
tablerow <- data.frame(
            endo_pheno=endo_pheno[j],
            Metabolite=metabolites[i],
            Beta=test$coefficients[2,1],
            Se=test$coefficients[2,2],
            p=test$coefficients[2,4],
            n=dim(na.omit(m[,c("metabo","PD_model")]))[[1]],
            lower=test$coef[2,1] - qt(0.975, df = test$df[2]) * test$coef[2, 2],
            upper=test$coef[2,1] + qt(0.975, df = test$df[2]) * test$coef[2, 2],
            stringsAsFactors=FALSE)
results <- rbind(results, tablerow)
}
}
head (results[order(results$p),])
results$FDR<-p.adjust(results[,5], method = 'fdr', n = length(results[,4]))
write.table(results,file="Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1_female.csv",sep="\t",row.names=F,col.names=T,quote=F)
