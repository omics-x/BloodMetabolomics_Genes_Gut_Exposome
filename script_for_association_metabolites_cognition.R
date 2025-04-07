library(foreign)
library(dplyr)
library("readxl")
### General cognition information
cognitionRS3<-read.spss("/home/sahmad/general_tasks/generalwork/some_test_analysis/costreamdata/data_analysis_jan06/e5_RS-I-5 RS-II-3 RS-III-2_cognition_complete 030217_analysis_file.sav",to.data.frame=T)
cognitionRS3<-cognitionRS3[,c("ergoid","FAC1_1","LDST5","STR1T5_adjusted","STR2T5_adjusted","STR3T5_adjusted","WFT5","WLTimm5","WLTdel5","WLTrecog5","PPB_sum5")]
### Rotterdam Study covariate information
covariates<-readRDS(file ="/home/sahmad/metabolomics/gut_liver_brain/quality_control_2021/association_v2_age_corrected/covariat_files/Study_RSI_IV_RSIII_2_covars.rds")
### Metabolomics data file 
load("/home/sahmad/metabolomics/gut_liver_brain/quality_control_2021/KNN_imputed_metabolon_data_3rdMarch.RData")
imputeddata<-as.data.frame(imputeddata)
colnames(imputeddata)<-paste("metab",colnames(imputeddata),sep="_")
imputeddata$ergoid<-rownames(imputeddata)
## merging the data 
merg1<-merge(cognitionRS3,covariates,by.x="ergoid",by.y="ergoid",all.x=T,sort=F)
merg2<-merge(imputeddata,merg1,by.x="ergoid",by.y="ergoid",all.x=T,sort=F)
merg2<-merg2[which(!is.na(merg2$FAC1_1)),]

#### Results matrix 
results <- data.frame(
            endo_pheno=as.character(),
            Metabolite=as.character(),
            Beta=as.numeric(),
            Se=as.numeric(),
            p=as.numeric(),
            n=as.numeric(),
            stringsAsFactors=FALSE)
            

#endo_pheno<-c("FAC1_1","LDST5","STR1T5_adjusted","STR2T5_adjusted","STR3T5_adjusted","WFT5","WLTimm5","WLTdel5","WLTrecog5","PPB_sum5")
endo_pheno<-c("FAC1_1")

metabolites<-colnames(imputeddata)[-dim(imputeddata)[2]]

for (j in 1:length(endo_pheno)) {
for ( i in 1:length(metabolites)){

merg2$metabo<-scale(merg2[,metabolites[i]], scale=TRUE, center=TRUE)
merg2$PD_model<-merg2[,endo_pheno[j]]
 test<- summary(lm (paste("PD_model ~ metabo + Age_blood_collection + sex + Lipilower + BMI",sep=""), data=merg2))
#test<- summary(lm (paste("PD_model ~ metabo + Age_blood_collection + sex + Lipilower + BMI + education",sep=""), data=merg2))
tablerow <- data.frame(
            endo_pheno=endo_pheno[j],
            Metabolite=metabolites[i],
            Beta=test$coefficients[2,1],
            Se=test$coefficients[2,2],
            p=test$coefficients[2,4],
            n=dim(na.omit(merg2[,c("metabo","PD_model")]))[[1]],
            stringsAsFactors=FALSE)
results <- rbind(results, tablerow)
}
}
head (results[order(results$p),])
results$FDR<-p.adjust(results[,5], method = 'fdr', n = length(results[,4]))

###### Annotation of the general cognition results
anno_com <- read_excel("/home/sahmad/metabolomics/gut_liver_brain/DUKE-0304-19ML+_CLIENT DATA TABLE.XLSX",sheet = 2)
anno_com<-anno_com%>%as.data.frame(anno_com)%>%mutate(Metabolite = paste0("metab_",anno_com$CHEM_ID))
res_anno<-anno_com%>%select("Metabolite","CHEMICAL_NAME","SUPER_PATHWAY","SUB_PATHWAY")%>%left_join(results,.,by="Metabolite")

#########################################################################################################
# forN<-merg2%>%select("metabo","Age_blood_collection","sex","Lipilower","BMI") %>% tidyr::drop_na()
# forN<-merg2%>%select("metabo","Age_blood_collection","sex","Lipilower","BMI","education") %>% tidyr::drop_na() 
# with education there are 898 samples and without education there are 896 samples 
########################################################################################################
write.table(res_anno,file="/home/sahmad/metabolomics/gut_liver_brain/quality_control_2021/association_v2_age_corrected/association_cognition/updated_jan2022_covars/Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated.csv",sep="\t",row.names=F,col.names=T,quote=F)

#########################################################################################################
# Sex stratified association of metabolites with general cognition 
########################################################################################################

library(foreign)
library(dplyr)
library("readxl")
### General cognition information
cognitionRS3<-read.spss("/home/sahmad/general_tasks/generalwork/some_test_analysis/costreamdata/data_analysis_jan06/e5_RS-I-5 RS-II-3 RS-III-2_cognition_complete 030217_analysis_file.sav",to.data.frame=T)
cognitionRS3<-cognitionRS3[,c("ergoid","FAC1_1","LDST5","STR1T5_adjusted","STR2T5_adjusted","STR3T5_adjusted","WFT5","WLTimm5","WLTdel5","WLTrecog5","PPB_sum5")]
### Rotterdam Study covariate information
covariates<-readRDS(file ="/home/sahmad/metabolomics/gut_liver_brain/quality_control_2021/association_v2_age_corrected/covariat_files/Study_RSI_IV_RSIII_2_covars.rds")
### Metabolomics data file 
load("/home/sahmad/metabolomics/gut_liver_brain/quality_control_2021/KNN_imputed_metabolon_data_3rdMarch.RData")
imputeddata<-as.data.frame(imputeddata)
colnames(imputeddata)<-paste("metab",colnames(imputeddata),sep="_")
imputeddata$ergoid<-rownames(imputeddata)
## merging the data 
merg1<-merge(cognitionRS3,covariates,by.x="ergoid",by.y="ergoid",all.x=T,sort=F)
merg2<-merge(imputeddata,merg1,by.x="ergoid",by.y="ergoid",all.x=T,sort=F)
merg2<-merg2[which(!is.na(merg2$FAC1_1)),]

# Repeat the association for both males and females seperately 
merg2_m<-merg2[which(merg2$sex=='female'),]
#merg2_m<-merg2[which(merg2$sex=='male'),]


#### Results matrix 
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
            

endo_pheno<-c("FAC1_1")

metabolites<-colnames(imputeddata)[-dim(imputeddata)[2]]

for (j in 1:length(endo_pheno)) {
for ( i in 1:length(metabolites)){

merg2_m$metabo<-scale(merg2_m[,metabolites[i]], scale=TRUE, center=TRUE)
merg2_m$PD_model<-merg2_m[,endo_pheno[j]]
#test<- summary(lm (paste("PD_model ~ metabo + Age_blood_collection + sex + Lipilower + BMI",sep=""), data=merg2))
test<- summary(lm (paste("PD_model ~ metabo + Age_blood_collection + Lipilower + BMI",sep=""), data=merg2_m))
tablerow <- data.frame(
            endo_pheno=endo_pheno[j],
            Metabolite=metabolites[i],
            Beta=test$coefficients[2,1],
            Se=test$coefficients[2,2],
            p=test$coefficients[2,4],
            n=dim(na.omit(merg2_m[,c("metabo","PD_model")]))[[1]],
            lower=test$coef[2,1] - qt(0.975, df = test$df[2]) * test$coef[2, 2],
            upper=test$coef[2,1] + qt(0.975, df = test$df[2]) * test$coef[2, 2],
            stringsAsFactors=FALSE)
results <- rbind(results, tablerow)
}
}
head (results[order(results$p),])
results$FDR<-p.adjust(results[,5], method = 'fdr', n = length(results[,4]))

###### Annotation of the general cognition results
anno_com <- read_excel("/home/sahmad/metabolomics/gut_liver_brain/DUKE-0304-19ML+_CLIENT DATA TABLE.XLSX",sheet = 2)
anno_com<-anno_com%>%as.data.frame(anno_com)%>%mutate(Metabolite = paste0("metab_",anno_com$CHEM_ID))
res_anno<-anno_com%>%select("Metabolite","CHEMICAL_NAME","SUPER_PATHWAY","SUB_PATHWAY")%>%left_join(results,.,by="Metabolite")

#########################################################################################################
# forN<-merg2%>%select("metabo","Age_blood_collection","sex","Lipilower","BMI") %>% tidyr::drop_na()
# forN<-merg2%>%select("metabo","Age_blood_collection","sex","Lipilower","BMI","education") %>% tidyr::drop_na() 
# with education there are 898 samples and without education there are 896 samples 
########################################################################################################
write.table(results,file="/home/sahmad/metabolomics/gut_liver_brain/quality_control_2021/association_v2_age_corrected/association_cognition/updated_jan2022_covars/sex_stratified_model//Gfactor_age_sex_BMI_lipid_low_female.csv",sep="\t",row.names=F,col.names=T,quote=F)
