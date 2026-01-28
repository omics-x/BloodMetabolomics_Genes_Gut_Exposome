### Supplementary Table 13: Mediation analysis 
### Load R packages 
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library(foreign)
library(dplyr)
library(mediation)
### Load General cognition information from phenotype files. General cognition (G-factor) is coded as FAC1_1 in the file 
cognitionRS3<-read.spss("e5_RS-I-5 RS-II-3 RS-III-2_cogmplete 030217_analysis_file.sav",to.data.frame=T)
cognitionRS3<-cognitionRS3[,c("ergoid","FAC1_1")]
### Load Rotterdam Study covariate information
covariates<-readRDS(file ="Study_RSI_IVcovars.rds")
### Upload Metabolomics data file
load("KNN_imputed_metabolon_data_3rdMarch.RData")
imputeddata<-as.data.frame(imputeddata)
colnames(imputeddata)<-paste("metab",colnames(imputeddata),sep="_")
imputeddata$ergoid<-rownames(imputeddata)
## Merging the data and selecting only ergothioneine (compound ID: metab_100002154) for further analysis 
merg1<-merge(cognitionRS3,covariates,by.x="ergoid",by.y="ergoid",all.x=T,sort=F)
merg2<-merge(imputeddata[,c("ergoid","metab_100002154")],merg1,by.x="ergoid",by.y="ergoid",all.x=T,sort=F)
merg2<-merg2[which(!is.na(merg2$FAC1_1)),]

#### Load medication information file for Ergo5 follow-up (RSIII-2 follow-up)
medication<-read.csv("medication_data_all_1percent.csv")
merg3<-merge(merg2,medication,by.x="ergoid",by.y="X",all.x=T,sort=F)

#### Define the variables 
endo_pheno<-c("FAC1_1")
metabo<-'metab_100002154'
merg3$metabolite <- scale(merg3[,metabo], scale=TRUE, center=TRUE)
merg3$metabolite<-as.numeric(merg3$metabolite)

### Mediation analysis for Antacid medications, adapted from https://library.virginia.edu/data/articles/introduction-to-mediation-analysis
### The e5_a02 is code for Antacid medications in Ergo5 (RSIII-2) dataset
model.0 <- lm(FAC1_1 ~ e5_a02, merg3)
model.M <- lm(metabolite ~ e5_a02, merg3)
model.Y <- lm(FAC1_1 ~ e5_a02 + metabolite, merg3)
results <- mediate(model.M, model.Y, treat='e5_a02', mediator='metabolite', boot=TRUE, sims=500)

summary(results)

#Causal Mediation Analysis for metab_100002154 with Antacids (e5_a02) code 
#Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#
#               Estimate 95% CI Lower 95% CI Upper p-value
#ACME             -0.074       -0.121        -0.04  <2e-16 ***
#ADE              -0.161       -0.304        -0.03   0.016 *
#Total Effect     -0.235       -0.368        -0.10  <2e-16 ***
#Prop. Mediated    0.315        0.155         0.71  <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 898
#Simulations: 500

#######################################################
#### Mediation analysis for Psychoanaleptics (e5_n06 )
#######################################################

model.0 <- lm(FAC1_1 ~ e5_n06, merg3)
model.M <- lm(metabolite ~ e5_n06, merg3)
model.Y <- lm(FAC1_1 ~ e5_n06 + metabolite, merg3)
results <- mediate(model.M, model.Y, treat='e5_n06', mediator='metabolite', boot=TRUE, sims=500)

summary(results)

# Causal Mediation Analysis
# 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# 
#                Estimate 95% CI Lower 95% CI Upper p-value
# ACME            -0.0464      -0.1062         0.01   0.108
# ADE             -0.1348      -0.3312         0.06   0.196
# Total Effect    -0.1812      -0.3749         0.02   0.084 .
# Prop. Mediated   0.2560      -0.4859         1.98   0.152
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Sample Size Used: 898
# 
# 
# Simulations: 500
#######################################################
#### Mediation analysis for Thyroid therapy (code:e5_h03)
#######################################################

model.0 <- lm(FAC1_1 ~ e5_h03, merg3)
model.M <- lm(metabolite ~ e5_h03, merg3)
model.Y <- lm(FAC1_1 ~ e5_h03 + metabolite, merg3)
results <- mediate(model.M, model.Y, treat='e5_h03', mediator='metabolite', boot=TRUE, sims=500)

summary(results)
# Causal Mediation Analysis
# 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# 
#                Estimate 95% CI Lower 95% CI Upper p-value
# ACME            -0.0703      -0.1336        -0.02  <2e-16 ***
# ADE              0.0444      -0.2016         0.29    0.72
# Total Effect    -0.0258      -0.2761         0.22    0.81
# Prop. Mediated   2.7192      -6.5213         5.43    0.81
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Sample Size Used: 898
# 
# 
# Simulations: 500
