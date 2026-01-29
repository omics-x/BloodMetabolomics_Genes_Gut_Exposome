### script to perform association of metabolites with MRI variables (Total brain volume, total hippocampal volume and Total white matter lesions)
### Results are provided in supplementary table 2
### Load packages 
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library("dplyr")
### Loading the MRI and metabolomics combined dataset for analysis 
mri_metabolomics_df<-read.table('RS1_5_Metabolon_MRIdata_2Jun2021.txt',head=T)
covariates<-readRDS(file ="Study_RSI_IV_RSIII_2_covars.rds")
m<-merge(mri_metabolomics_df,covariates,by.all="ergoid",all.x=T)
########################################################
########################################################
### Model 1: adjusted for Age_blood_collection + sex + BMI + Lipilowering medication + ICV_from_mask
########################################################
########################################################
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

# List of tested phenotype: Total hippocampal volume, total brain volume (total_par_ml), and total white matter lesions
endo_pheno<-c("total_Hippocampus","Total_par_ml","Total_wml")

for (j in 1:length(endo_pheno)) {
for ( i in 1:length(metabolites)){
# each metabolite is transformed to z-value before association
m$metabo<-scale(m[,metabolites[i]],center = TRUE, scale = TRUE)
# all MRI variables are log transformed and scaled, ICV_from_mask is intracranial volume variable 
m$mri_variable<-scale(log(m[,endo_pheno[j]]),center = TRUE, scale = TRUE)
test<- summary(lm (paste("mri_variable ~ metabo + Age_blood_collection + sex + BMI + Lipilower + ICV_from_mask",sep=""), data=m))
tablerow <- data.frame(
            endo_pheno=endo_pheno[j],
            Metabolite=metabolites[i],
            Beta=test$coefficients[2,1],
            Se=test$coefficients[2,2],
            p=test$coefficients[2,4],
            n=dim(na.omit(m[,c("metabo","mri_variable")]))[[1]],
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
########################################################
### Model 2: adjusted for Age_blood_collection + sex + BMI + Lipilower + Smoke + Diabetes2 + Hypertension + ICV_from_mask
########################################################
########################################################
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

# List of tested phenotype: Total hippocampal volume, total brain volume (total_par_ml), and total white matter lesions
endo_pheno<-c("total_Hippocampus","Total_par_ml","Total_wml")

for (j in 1:length(endo_pheno)) {
  for ( i in 1:length(metabolites)){
    # each metabolite is transformed to z-value before association
    m$metabo<-scale(m[,metabolites[i]],center = TRUE, scale = TRUE)
    # all MRI variables are log transformed and scaled, ICV_from_mask is intracranial volume variable 
    m$mri_variable<-scale(log(m[,endo_pheno[j]]),center = TRUE, scale = TRUE)
    test<- summary(lm (paste("mri_variable ~ metabo + Age_blood_collection + sex + BMI + Lipilower + Smoke + Diabetes2 + Hypertension + ICV_from_mask",sep=""), data=m))
    tablerow <- data.frame(
      endo_pheno=endo_pheno[j],
      Metabolite=metabolites[i],
      Beta=test$coefficients[2,1],
      Se=test$coefficients[2,2],
      p=test$coefficients[2,4],
      n=dim(na.omit(m[,c("metabo","mri_variable")]))[[1]],
      lower=test$coef[2,1] - qt(0.975, df = test$df[2]) * test$coef[2, 2],
      upper=test$coef[2,1] + qt(0.975, df = test$df[2]) * test$coef[2, 2],
      stringsAsFactors=FALSE)
    results <- rbind(results, tablerow)
  }
}
head (results[order(results$p),])
results$FDR<-p.adjust(results[,5], method = 'fdr', n = length(results[,4]))

write.table(results,file="Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m2.csv",sep="\t",row.names=F,col.names=T,quote=F)

########################################################
########################################################
## Association of metabolites with MRI variables in only participants which have less than equal to 1 year difference between MRI and blood collection
########################################################
########################################################

MRI_data<-read.table("RS1_5_Metabolon_MRIdata_2Jun2021.txt",head=T)
covariates<-readRDS(file ="Study_RSI_IV_RSIII_2_covars.rds")
covariates<-as.data.frame(covariates)
m<-left_join(MRI_data,covariates,by="ergoid")

# abs_difference variable defines the absolute time in years between MRI and blood collection for metabolomics 
m<-m %>% filter(abs_difference<=1)
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


endo_pheno<-c("total_Hippocampus","Total_par_ml","Total_wml")
metabolites<-grep("metab_",colnames(m),value = T)

for (j in 1:length(endo_pheno)) {
  for ( i in 1:length(metabolites)){
    
    m$metabo<-scale(m[,metabolites[i]],center = TRUE, scale = TRUE)
    m$mri_variable<-scale(log(m[,endo_pheno[j]]),center = TRUE, scale = TRUE)
    test<- summary(lm (paste("mri_variable ~ metabo + Age_blood_collection + sex.y + BMI + Lipilower + ICV_from_mask",sep=""), data=m))
    #test<- summary(lm (paste("mri_variable ~ metabo + Age_blood_collection + sex + BMI + Lipilower + Smoke + Diabetes2 + Hypertension + ICV_from_mask",sep=""), data=m))
    tablerow <- data.frame(
      endo_pheno=endo_pheno[j],
      Metabolite=metabolites[i],
      Beta=test$coefficients[2,1],
      Se=test$coefficients[2,2],
      p=test$coefficients[2,4],
      n=dim(na.omit(m[,c("metabo","mri_variable")]))[[1]],
      lower=test$coef[2,1] - qt(0.975, df = test$df[2]) * test$coef[2, 2],
      upper=test$coef[2,1] + qt(0.975, df = test$df[2]) * test$coef[2, 2],
      stringsAsFactors=FALSE)
    results <- rbind(results, tablerow)
  }
}
head (results[order(results$p),])
results$FDR<-p.adjust(results[,5], method = 'fdr', n = length(results[,4]))

write.table(results,file="Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1_YearLessthan1.csv",sep="\t",row.names=F,col.names=T,quote=F)

########################################################
########################################################
### Sex stratified analysis 
########################################################
########################################################
m<-merge(mri_metabolomics_df,covariates,by.all="ergoid",all.x=T)
m<-m[m$sex=='female',]
#m<-m[m$sex=='male',]
# Model 1 was run separately for females and males 
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


endo_pheno<-c("total_Hippocampus","Total_par_ml","Total_wml")

for (j in 1:length(endo_pheno)) {
for ( i in 1:length(metabolites)){

m$metabo<-scale(m[,metabolites[i]],center = TRUE, scale = TRUE)
m$mri_variable<-scale(log(m[,endo_pheno[j]]),center = TRUE, scale = TRUE)
test<- summary(lm (paste("mri_variable ~ metabo + Age_blood_collection + BMI + Lipilower + ICV_from_mask",sep=""), data=m))
tablerow <- data.frame(
            endo_pheno=endo_pheno[j],
            Metabolite=metabolites[i],
            Beta=test$coefficients[2,1],
            Se=test$coefficients[2,2],
            p=test$coefficients[2,4],
            n=dim(na.omit(m[,c("metabo","mri_variable")]))[[1]],
            lower=test$coef[2,1] - qt(0.975, df = test$df[2]) * test$coef[2, 2],
            upper=test$coef[2,1] + qt(0.975, df = test$df[2]) * test$coef[2, 2],
            stringsAsFactors=FALSE)
results <- rbind(results, tablerow)
}
}
head (results[order(results$p),])
results$FDR<-p.adjust(results[,5], method = 'fdr', n = length(results[,4]))
write.table(results,file="Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1_female.csv",sep="\t",row.names=F,col.names=T,quote=F)
#write.table(results,file="Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1_male.csv",sep="\t",row.names=F,col.names=T,quote=F)


###############################################################################################
## Sex interaction term in the model analysis: Supplementary Table 3 
###############################################################################################
MRI_data<-read.table("/Users/sahmad1/Documents/RS_projects/mri_associations/mri_new_covariat_file/RS1_5_Metabolon_MRIdata_2Jun2021.txt",head=T)
covariates<-readRDS(file ="/Users/sahmad1/Documents/RS_projects/RS1_replication/cognitionReplication_R01_dec21/data/Study_RSI_IV_RSIII_2_covars.rds")
covariates<-as.data.frame(covariates)
m<-left_join(MRI_data,covariates,by="ergoid")

#########################################################
#########################################################
results <- data.frame(
  endo_pheno = as.character(),
  Metabolite = as.character(),
  Beta = as.numeric(),     
  Se = as.numeric(),       
  p = as.numeric(),        
  n = as.numeric(),
  lower = as.numeric(),
  upper = as.numeric(),
  stringsAsFactors = FALSE
)

endo_pheno <- c("total_Hippocampus", "Total_par_ml", "Total_wml")
metabolites <- grep("metab_", colnames(m), value = TRUE)

# convert sex value into factor
m$sex.y <- factor(m$sex.y, levels = c("female", "male"))  

for (j in 1:length(endo_pheno)) {
  for (i in 1:length(metabolites)) {
    
    m$metabo <- scale(m[, metabolites[i]], center = TRUE, scale = TRUE)
    m$mri_variable <- scale(log(m[, endo_pheno[j]]), center = TRUE, scale = TRUE)
    
    model <- lm(mri_variable ~ metabo * sex.y + Age_blood_collection + BMI + Lipilower + ICV_from_mask, data = m)
    test <- summary(model)
    
    # Get interaction term row: metabo:sex.ymale
    interaction_row <- grep("metabo:sex\\.ymale", rownames(test$coefficients))
    
    # Check in case interaction term is missing due to singularity or coding error
    if (length(interaction_row) == 1) {
      tablerow <- data.frame(
        endo_pheno = endo_pheno[j],
        Metabolite = metabolites[i],
        Beta = test$coefficients[interaction_row, 1],
        Se = test$coefficients[interaction_row, 2],
        p = test$coefficients[interaction_row, 4],
        n = nrow(na.omit(m[, c("metabo", "mri_variable", "sex.y")])),
        lower = test$coefficients[interaction_row, 1] - qt(0.975, df = test$df[2]) * test$coefficients[interaction_row, 2],
        upper = test$coefficients[interaction_row, 1] + qt(0.975, df = test$df[2]) * test$coefficients[interaction_row, 2],
        stringsAsFactors = FALSE
      )
      results <- rbind(results, tablerow)
    }
  }
}
results <- results[order(results), ]
# Write to file
write.table(results,
            file = "MRI_Interaction_Metabolite_Sex.csv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
)



