####################################################################################
### Source data files generated from raw files for extended figures 
####################################################################################
####################################################################################
### Extended Figure 1: Correlation plot for the metabolites associated with cognition and MRI 
####################################################################################
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library("corrplot")
library("readxl")
## Load results of association of MRI and general cognition with metabolites 
mri<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_MRI_M1_RSIII_2.csv",sep='\t')
cog<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2.csv",sep='\t')
## Unique metabolites associated with cognition and the three MRI markers
metabolites_cog<-cog$Metabolite[which(cog$FDR<0.05)]
metabolites_mri<-mri$Metabolite[which(mri$FDR<0.05)]
all_metabolites<-unique(c(metabolites_cog,metabolites_mri))
### Load metabolomics data for RSIII-2
load("/Users/sahmad1/Documents/RS_projects/KNN_imputed_metabolon_data_3rdMarch.RData")
imputeddata<-as.data.frame(imputeddata)
colnames(imputeddata)<-paste("metab",colnames(imputeddata),sep="_")
metabolites<-colnames(imputeddata)
imputeddata$ergoid<-rownames(imputeddata)
### Annotation file
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)

## Subset the full metabolite matrix and annotate metabolite names
subset<-imputeddata[,all_metabolites]
sub_set_anno<-anno_com[anno_com$Name%in%colnames(subset),]
subset<-subset[,sub_set_anno$Name]
## Use the CHEMICAL_NAME_new column, as it is already corrected for 1-carboxyethyl metabolites
colnames(subset)<-sub_set_anno$CHEMICAL_NAME_new
### Plotting and saving the PDF
cor_matrix<-cor(subset, method = c("spearman"))
cor_df <- as.data.frame(cor_matrix)
cor_df$Metabolites <- rownames(cor_df)
cor_df <- cor_df[, c(ncol(cor_df), 1:(ncol(cor_df)-1))]
write_xlsx(
  list(Extended_Figure1 = cor_df),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Extended_Figure1.xlsx"
)


#################################################### 
## Extended Figure 2a: UpSet plot
#################################################### 
## Load results from association analyses of general cognition and MRI measures
mri<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_MRI_M1_RSIII_2.csv",sep='\t')
cog<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2.csv",sep='\t')

## Unique metabolites associated with cognition and the three MRI markers
metabolites_cog<-cog$Metabolite[which(cog$p<0.05)]
metabolites_mri<-mri$Metabolite[which(mri$p<0.05)]
all_metabolites<-unique(c(metabolites_cog,metabolites_mri))
### Create three data frames for the three MRI variables
tbv<-mri[mri$endo_pheno=="Total_par_ml",]
hcv<-mri[mri$endo_pheno=="total_Hippocampus",]
wml<-mri[mri$endo_pheno=="Total_wml",]
## Load the annotation file 
anno_com <- read_excel("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)

## Merge annotation file with MRI result files 
tbv<-merge(tbv,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
hcv<-merge(hcv,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
wml<-merge(wml,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)

## Extract only metabolites with p<0.05
all_sig<- rbind(
  G_factor = cog$Metabolite[which(cog$p<0.05)], 
  Total_BV = tbv$Metabolite[which(tbv$p<0.05)], 
  Total_HCV =  hcv$Metabolite[which(hcv$p<0.05)],
  Total_WML = wml$Metabolite[which(wml$p<0.05)]
)

## Convert to dataframe for saving in excel 
all_sig_df <- data.frame(
  Outcome = rownames(all_sig),
  Metabolite = as.vector(all_sig)
)

#################################################### 
## Extended Figure 2b-d: Scatter plots for overlaping metabolites between cognition and MRI
#################################################### 
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library(ggplot2)
library(jcolors)
library(gridExtra)

## Load results from association analyses of general cognition and MRI measures
mri<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_MRI_M1_RSIII_2.csv",sep='\t')
cog<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2.csv",sep='\t')

### Unique metabolites associated with cognition and the three MRI markers
metabolites_cog<-cog$Metabolite[which(cog$p<0.05)]
metabolites_mri<-mri$Metabolite[which(mri$p<0.05)]
all_metabolites<-unique(c(metabolites_cog,metabolites_mri))
### Create three data frames for the three MRI variables
tbv<-mri[mri$endo_pheno=="Total_par_ml",]
hcv<-mri[mri$endo_pheno=="total_Hippocampus",]
wml<-mri[mri$endo_pheno=="Total_wml",]

### Load annotation file MRI data

anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)
## Merge annotation file with MRI result files 
tbv<-merge(tbv,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
hcv<-merge(hcv,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
wml<-merge(wml,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)


cog_set<-cog[cog$Metabolite%in%all_metabolites,]
tbv_set<-tbv[tbv$Metabolite%in%all_metabolites,]
hcv_set<-hcv[hcv$Metabolite%in%all_metabolites,]
wml_set<-wml[wml$Metabolite%in%all_metabolites,]

### Merging the MRI variables with cognition
# cognition results combined with total brain volume
cog_tbv<-merge(cog_set[,c("Beta","Metabolite","CHEMICAL_NAME","SUPER_PATHWAY","FDR")],tbv_set[,c("Beta","CHEMICAL_NAME","FDR")],by='CHEMICAL_NAME')
cog_tbv$SUPER_PATHWAY[is.na(cog_tbv$SUPER_PATHWAY)]<-"Unknowns"

# added to get ID and names 
colnames(cog_tbv)<-c("CHEMICAL_NAME","Beta_cognition","Metabolite","SUPER_PATHWAY","FDR_cognition","Beta_TBV","FDR_TBV")
# consistent correct naming 3-methyl catechol sulfate (2) to 3-methylcatechol sulfate
cog_tbv$CHEMICAL_NAME[cog_tbv$CHEMICAL_NAME=='3-methyl catechol sulfate (2)']<-'3-methylcatechol sulfate'
cog_tbv$significant_metabolites<-NA
cog_tbv$significant_metabolites<-ifelse(cog_tbv$FDR_cognition<0.05 | cog_tbv$FDR_TBV<0.05,cog_tbv$CHEMICAL_NAME,"")

## Cognition results combined with hippocampal volume
cog_hcv <- merge(cog_set[,c("Beta","Metabolite","CHEMICAL_NAME","SUPER_PATHWAY","FDR")],hcv_set[,c("Beta","CHEMICAL_NAME","FDR")],by = "CHEMICAL_NAME")
cog_hcv$SUPER_PATHWAY[is.na(cog_hcv$SUPER_PATHWAY)] <- "Unknowns"
colnames(cog_hcv) <- c("CHEMICAL_NAME","Beta_cognition","Metabolite","SUPER_PATHWAY","FDR_cognition","Beta_HCV","FDR_HCV")
cog_hcv$CHEMICAL_NAME[cog_hcv$CHEMICAL_NAME == "3-methyl catechol sulfate (2)"] <- "3-methylcatechol sulfate"
cog_hcv$significant_metabolites <- ifelse(cog_hcv$FDR_cognition < 0.05 | cog_hcv$FDR_HCV < 0.05,cog_hcv$CHEMICAL_NAME,"")

## Cognition results combined with white matter lesions
cog_wml <- merge(cog_set[,c("Beta","Metabolite","CHEMICAL_NAME","SUPER_PATHWAY","FDR")],wml_set[,c("Beta","CHEMICAL_NAME","FDR")],by = "CHEMICAL_NAME")
cog_wml$SUPER_PATHWAY[is.na(cog_wml$SUPER_PATHWAY)] <- "Unknowns"
colnames(cog_wml) <- c("CHEMICAL_NAME","Beta_cognition","Metabolite","SUPER_PATHWAY","FDR_cognition","Beta_WML","FDR_WML")

cog_wml$CHEMICAL_NAME[cog_wml$CHEMICAL_NAME == "3-methyl catechol sulfate (2)"] <- "3-methylcatechol sulfate"
cog_wml$significant_metabolites <- ifelse( cog_wml$FDR_cognition < 0.05 | cog_wml$FDR_WML < 0.05,cog_wml$CHEMICAL_NAME,"")

cog_tbv <- cog_tbv %>%
  select(Metabolite, CHEMICAL_NAME, SUPER_PATHWAY,
         Beta_cognition, FDR_cognition,
         Beta_TBV, FDR_TBV,
         significant_metabolites)
cog_hcv <- cog_hcv %>%
  select(Metabolite, CHEMICAL_NAME, SUPER_PATHWAY,
         Beta_cognition, FDR_cognition,
         Beta_HCV, FDR_HCV,
         significant_metabolites)
cog_wml <- cog_wml %>%
  select(Metabolite, CHEMICAL_NAME, SUPER_PATHWAY,
         Beta_cognition, FDR_cognition,
         Beta_WML, FDR_WML,
         significant_metabolites)

write_xlsx(
  list(
    Extended_Figure2a = all_sig_df,
    Extended_Figure2b = cog_tbv,
    Extended_Figure2c = cog_hcv,
    Extended_Figure2d = cog_wml
  ),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Extended_Figure2.xlsx"
)
################################################################################
### Extended Figure 3: Relationship of metabolites strongly influenced by the gut microbiome (EV ≥ 5%) with metabolites nominally associated (P < 0.05) with general cognition 
################################################################################
### Load the files 
cog<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2.csv",sep='\t')
female<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2_Sex_stratified_Female.csv",sep='\t')
male<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2_Sex_stratified_male.csv",sep='\t')
### Load annotation file 
library("readxl")
anno_com <- read_excel("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)

male_anno<-merge(male,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
female_anno<-merge(female,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
### Load the files 
mic<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/EV_microbiota_results.csv")
mic$fdr<-p.adjust(mic$spearman_p,method='fdr')
mic$EV_Microbiota<-ifelse((mic$fdr<0.05 & mic$explained_variance_score>0),mic$explained_variance_score*100,0)
mic_sub<-mic[mic$EV_Microbiota>=5,]

merg_all<-merge(mic_sub,cog,by.x="X",by.y='Metabolite',all.x=T)
merg_female<-merge(mic_sub,female_anno,by.x="X",by.y='Metabolite',all.x=T)
merg_male<-merge(mic_sub,male_anno,by.x="X",by.y='Metabolite',all.x=T)

merg_all<-merg_all[,c("X","CHEMICAL_NAME","EV_Microbiota","Beta","p","FDR","lower","upper")]
merg_female<-merg_female[,c("X","CHEMICAL_NAME","EV_Microbiota","Beta","p","FDR","lower","upper")]
merg_male<-merg_male[,c("X","CHEMICAL_NAME","EV_Microbiota","Beta","p","FDR","lower","upper")]

### Replace corrected names 
name_map <- c(
  "1-carboxyethyltyrosine" = "N-lactoyltyrosine",
  "1-carboxyethylphenylalanine" = "N-lactoylphenylalanine",
  "1-carboxyethylisoleucine" = "N-lactoylisoleucine",
  "3-methyl catechol sulfate (2)" ="3-methylcatechol sulfate",
  "3-methyl catechol sulfate (1)" ="3-methylcatechol sulfate (1)"
)

merg_all$CHEMICAL_NAME[merg_all$CHEMICAL_NAME %in% names(name_map)] <- 
  name_map[merg_all$CHEMICAL_NAME[merg_all$CHEMICAL_NAME %in% names(name_map)]]

merg_female$CHEMICAL_NAME[merg_female$CHEMICAL_NAME %in% names(name_map)] <- 
  name_map[merg_female$CHEMICAL_NAME[merg_female$CHEMICAL_NAME %in% names(name_map)]]

merg_male$CHEMICAL_NAME[merg_male$CHEMICAL_NAME %in% names(name_map)] <- 
  name_map[merg_male$CHEMICAL_NAME[merg_male$CHEMICAL_NAME %in% names(name_map)]]

write_xlsx(
  list(
    Extended_Figure3a = merg_all,
    Extended_Figure3b = merg_male,
    Extended_Figure3c = merg_female
  ),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Extended_Figure3.xlsx"
)

######################################################################### 
### Extended Figure 4a: Heatmap of cognition associated metabolites with clinical factors 
######################################################################### 

cog<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2.csv",head=T,sep='\t')
metabolites<-cog$Metabolite[which(cog$FDR<0.05)]
### Load results
res_anno<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_medication_RSIII_2.csv",head=T)
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
res_anno$medication_names<-str_trim(res_anno$name)
### Plot files 
res_anno<-res_anno[which(res_anno$Metabolite%in%metabolites),]
# Name correction
name_map <- c(
  "1-carboxyethyltyrosine" = "N-lactoyltyrosine",
  "1-carboxyethylphenylalanine" = "N-lactoylphenylalanine",
  "1-carboxyethylisoleucine" = "N-lactoylisoleucine",
  "3-methyl catechol sulfate (2)" ="3-methylcatechol sulfate",
  "3-methyl catechol sulfate (1)" ="3-methylcatechol sulfate (1)"
)
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME %in% names(name_map)] <-name_map[res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME %in% names(name_map)]]
extFig4a<-res_anno[,c("Metabolite","CHEMICAL_NAME","medication_names","Beta","Se","zvalue","FDR")]

######################################################################### 
### Extended Figure 4b: Heatmap of cognition associated metabolites with lifestyle factors 
#########################################################################

cog<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2.csv",head=T,sep='\t')
metabolites<-cog$Metabolite[which(cog$FDR<0.05)]
### Load results
res_anno<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_demographics_RSIII_2.csv",head=T)
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
res_anno<-res_anno[which(res_anno$Metabolite%in%metabolites),]
# Name correction
name_map <- c(
  "1-carboxyethyltyrosine" = "N-lactoyltyrosine",
  "1-carboxyethylphenylalanine" = "N-lactoylphenylalanine",
  "1-carboxyethylisoleucine" = "N-lactoylisoleucine",
  "3-methyl catechol sulfate (2)" ="3-methylcatechol sulfate",
  "3-methyl catechol sulfate (1)" ="3-methylcatechol sulfate (1)"
)
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME %in% names(name_map)] <-name_map[res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME %in% names(name_map)]]
extFig4b<-res_anno[,c("Metabolite","CHEMICAL_NAME","Demographics","Beta","Se","zvalue","FDR")]

######################################################################### 
### Extended Figure 4c: Heatmap of cognition associated metabolites with clinical factors 
######################################################################### 
cog<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2.csv",head=T,sep='\t')
metabolites<-cog$Metabolite[which(cog$FDR<0.05)]
### Load results
res_anno<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_clinical_factors_RSIII_2.csv",head=T)
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
### Prepare files for heatmap
res_anno<-res_anno[which(res_anno$Metabolite%in%metabolites),]
# Name correction
name_map <- c(
  "1-carboxyethyltyrosine" = "N-lactoyltyrosine",
  "1-carboxyethylphenylalanine" = "N-lactoylphenylalanine",
  "1-carboxyethylisoleucine" = "N-lactoylisoleucine",
  "3-methyl catechol sulfate (2)" ="3-methylcatechol sulfate",
  "3-methyl catechol sulfate (1)" ="3-methylcatechol sulfate (1)"
)
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME %in% names(name_map)] <-name_map[res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME %in% names(name_map)]]
extFig4c<-res_anno[,c("Metabolite","CHEMICAL_NAME","Demographics","Beta","Se","zvalue","FDR")]

write_xlsx(
  list(
    Extended_Figure4a = extFig4a,
    Extended_Figure4b = extFig4b,
    Extended_Figure4c = extFig4c
  ),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Extended_Figure4.xlsx"
)


######################################################################### 
### Extended Figure 5a: Heatmap of MRI associated metabolites with medication use 
#########################################################################
library(gplots)

mri<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_MRI_M1_RSIII_2.csv",head=T,sep='\t')
metabolites<-mri$Metabolite[which(mri$FDR<0.05)]
### Load association results of medication use with metabolomics 
res_anno<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_medication_RSIII_2.csv",head=T)
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
res_anno$medication_names<-str_trim(res_anno$name)
## P-value and regression coefficient matrix 
res_anno<-res_anno[which(res_anno$Metabolite%in%metabolites),]
## Correct the name of metabolite
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyltyrosine"
extFig5a<-res_anno[,c("Metabolite","CHEMICAL_NAME","medication_names","Beta","Se","zvalue","FDR")]

######################################################################### 
### Extended Figure 5b: Heatmap of MRI associated metabolites with lifestyle factors  
#########################################################################
mri<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_MRI_M1_RSIII_2.csv",head=T,sep='\t')
metabolites<-mri$Metabolite[which(mri$FDR<0.05)]
### Load results
res_anno<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_demographics_RSIII_2.csv",head=T)
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
### Prepare files 
res_anno<-res_anno[which(res_anno$Metabolite%in%metabolites),]
## Correct the name of metabolite
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyltyrosine"
extFig5b<-res_anno[,c("Metabolite","CHEMICAL_NAME","Demographics","Beta","Se","zvalue","FDR")]

######################################################################### 
### Extended Figure5c: Heatmap of MRI associated metabolites with clinical factors  
#########################################################################
library(stringr)
library(gplots)

mri<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_MRI_M1_RSIII_2.csv",head=T,sep='\t')
metabolites<-mri$Metabolite[which(mri$FDR<0.05)]
### Load results
res_anno<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_clinical_factors_RSIII_2.csv",head=T)
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
### Prepare heatmap files 
res_anno<-res_anno[which(res_anno$Metabolite%in%metabolites),]
## Correct the name of metabolite
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyltyrosine"
extFig5c<-res_anno[,c("Metabolite","CHEMICAL_NAME","Demographics","Beta","Se","zvalue","FDR")]


write_xlsx(
  list(
    Extended_Figure5a = extFig5a,
    Extended_Figure5b = extFig5b,
    Extended_Figure5c = extFig5c
  ),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Extended_Figure5.xlsx"
)



#########################################################################
### Extended Figure 6a: Relationship of general cognition and brain imaging associated metabolites with gut microbiota
#########################################################################

cog<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_Cognition_M1_RSIII_2.csv",sep='\t')
metabolites<-cog$Metabolite[which(cog$FDR<0.05)]

### Load results of association of gut-microbiome ASVs with metabolomics
result<-read.table("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_analysis_metabolites_gut_microbiota_RSIII_2_full_model_Nov1_clr.txt",head=T)
### Load Annotation file 
anno_com <- read_excel("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)
res_anno<-merge(result,anno_com[,c("Name","CHEMICAL_NAME","SUPER_PATHWAY","SUB_PATHWAY")],by.x="Pheno",by.y="Name",all.x=T)

res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyltyrosine"
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "3-methyl catechol sulfate (2)"]<- "3-methylcatechol sulfate"
### FDR correction 
res_anno$FDR<-p.adjust(res_anno[,5], method = 'fdr', n = length(res_anno[,4]))
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
res_anno<-res_anno[grep("genus",res_anno$Taxon),]
### Generate files for heatmap 
res_anno<-res_anno[which(res_anno$Pheno%in%metabolites),]

extFig6a<-res_anno[,c("Pheno","CHEMICAL_NAME","Taxon","Beta","Se","zvalue","FDR")]

#########################################################################
### Extended Figure 6b
#########################################################################
### Load results of association of MRI markers with metabolomics
mri<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_regression_MRI_M1_RSIII_2.csv",sep='\t')
metabolites<-mri$Metabolite[which(mri$FDR<0.05)]

### Load results of association of gut-microbiome ASVs with metabolomics
result<-read.table("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/Association_analysis_metabolites_gut_microbiota_RSIII_2_full_model_Nov1_clr.txt",head=T)
anno_com <- read_excel("/Users/sahmad1/Downloads/SOURCE_FILES/inputfilessourcedata/DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)
res_anno<-merge(result,anno_com[,c("Name","CHEMICAL_NAME","SUPER_PATHWAY","SUB_PATHWAY")],by.x="Pheno",by.y="Name",all.x=T)

res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyltyrosine"
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "3-methyl catechol sulfate (2)"]<- "3-methylcatechol sulfate"
### FDR correction 
res_anno$FDR<-p.adjust(res_anno[,5], method = 'fdr', n = length(res_anno[,4]))
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
res_anno<-res_anno[grep("genus",res_anno$Taxon),]
### Plot files 
res_anno<-res_anno[which(res_anno$Pheno%in%metabolites),]
extFig6b<-res_anno[,c("Pheno","CHEMICAL_NAME","Taxon","Beta","Se","zvalue","FDR")]

write_xlsx(
  list(
    Extended_Figure6a = extFig6a,
    Extended_Figure6b = extFig6b
  ),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Extended_Figure6.xlsx"
)
