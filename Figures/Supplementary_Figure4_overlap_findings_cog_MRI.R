#################################################### 
## Code for Supplementary Figure 4: Concordance of metabolite signatures between cognition and MRI phenotypes
#################################################### 
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library(UpSetR)
library("readxl")

#################################################### 
## Supplementary Figure 4A: UpSet plot
#################################################### 
## Load results from association analyses of general cognition and MRI measures
mri<-read.csv("Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",sep='\t')
cog<-read.csv("Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated_95CI.csv",sep='\t')

## Unique metabolites associated with cognition and the three MRI markers
metabolites_cog<-cog$Metabolite[which(cog$p<0.05)]
metabolites_mri<-mri$Metabolite[which(mri$p<0.05)]
all_metabolites<-unique(c(metabolites_cog,metabolites_mri))
### Create three data frames for the three MRI variables
tbv<-mri[mri$endo_pheno=="Total_par_ml",]
hcv<-mri[mri$endo_pheno=="total_Hippocampus",]
wml<-mri[mri$endo_pheno=="Total_wml",]
## Load the annotation file 
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)

## Merge annotation file with MRI result files 
tbv<-merge(tbv,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
hcv<-merge(hcv,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
wml<-merge(wml,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)

## Extract only metabolites with p<0.05
all_sig<- list(
  G_factor = cog$Metabolite[which(cog$p<0.05)], 
  Total_BV = tbv$Metabolite[which(tbv$p<0.05)], 
  Total_HCV =  hcv$Metabolite[which(hcv$p<0.05)],
  Total_WML = wml$Metabolite[which(wml$p<0.05)]
)

### Plot the UpSet plot
pdf("UpSet_plot_overlap_associations_MRI_cognition.pdf")
upset(fromList(all_sig),order.by="freq")
dev.off()

#################################################### 
## Supplementary Figure 4B-D: Scatter plots for overlaping metabolites between cognition and MRI
#################################################### 
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library(ggplot2)
library(jcolors)
library(gridExtra)

## Load results from association analyses of general cognition and MRI measures
mri<-read.csv("Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",sep='\t')
cog<-read.csv("Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated_95CI.csv",sep='\t')

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
cog_tbv<-merge(cog_set[,c("Beta","CHEMICAL_NAME","SUPER_PATHWAY","FDR")],tbv_set[,c("Beta","CHEMICAL_NAME","FDR")],by='CHEMICAL_NAME')
cog_tbv$SUPER_PATHWAY[is.na(cog_tbv$SUPER_PATHWAY)]<-"Unknowns"
cog_tbv$name<-NA
cog_tbv$name<-ifelse(cog_tbv$FDR.x<0.05 | cog_tbv$FDR.y<0.05,cog_tbv$CHEMICAL_NAME,"")

## Cognition results combined with hippocampal volume
cog_hcv<-merge(cog_set[,c("Beta","CHEMICAL_NAME","SUPER_PATHWAY","FDR")],hcv_set[,c("Beta","CHEMICAL_NAME","FDR")],by='CHEMICAL_NAME')
cog_hcv$SUPER_PATHWAY[is.na(cog_hcv$SUPER_PATHWAY)]<-"Unknowns"
cog_hcv$name<-NA
cog_hcv$name<-ifelse(cog_hcv$FDR.x<0.05 | cog_hcv$FDR.y<0.05,cog_hcv$CHEMICAL_NAME,"")

## Cognition results combined with white matter lesions
cog_wml<-merge(cog_set[,c("Beta","CHEMICAL_NAME","SUPER_PATHWAY","FDR")],wml_set[,c("Beta","CHEMICAL_NAME","FDR")],by='CHEMICAL_NAME')
cog_wml$SUPER_PATHWAY[is.na(cog_wml$SUPER_PATHWAY)]<-"Unknowns"
cog_wml$name<-NA
cog_wml$name<-ifelse(cog_wml$FDR.x<0.05 | cog_wml$FDR.y<0.05,cog_wml$CHEMICAL_NAME,"")

### Plot scatter plots for the overlap between cognition and the three MRI phenotypes and save as a PDF
pdf(file="Correlation_plot_cognition_TBV_Sup.Fig4B.pdf",width=10,height=8)
# Calculate correlation
correlation <- cor(cog_tbv$Beta.x, cog_tbv$Beta.y)
# Calculate p-value using cor.test()
cor_test <- cor.test(cog_tbv$Beta.x, cog_tbv$Beta.y)
p_value <- cor_test$p.value
# Convert correlation and p-value to strings
cor_str <- sprintf("r = %.2f", correlation)
p_value_str <- sprintf("p = %.2e", p_value)  # Display p-value in scientific notation
# Plot with annotations
ggplot(cog_tbv, aes(x = Beta.x, y = Beta.y, color = factor(SUPER_PATHWAY))) +
  geom_point() +
  geom_text(label = cog_tbv$name, show.legend = FALSE, size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_light() +
  labs(title = "(A)", x = "Regression_Coefficient_Cognition", y = "Regression_Coefficient_TBV") +
  scale_color_jcolors(palette = "pal8") +
  annotate("text", x = min(cog_tbv$Beta.x), y = max(cog_tbv$Beta.y),
           label = paste(cor_str, p_value_str, sep = ", "), hjust = 0, vjust = 1)
dev.off()

pdf(file="Correlation_plot_cognition_HCV_Sup.Fig4C.pdf",width=10,height=8)                                                                                                 
# Calculate correlation
correlation <- cor(cog_hcv$Beta.x, cog_hcv$Beta.y)
# Calculate p-value using cor.test()
cor_test <- cor.test(cog_hcv$Beta.x, cog_hcv$Beta.y)
p_value <- cor_test$p.value
# Convert correlation and p-value to strings
cor_str <- sprintf("r = %.2f", correlation)
p_value_str <- sprintf("p = %.2e", p_value)  # Display p-value in scientific notation
# Plot with annotations
ggplot(cog_hcv, aes(x = Beta.x, y = Beta.y, color = factor(SUPER_PATHWAY))) +
  geom_point() +
  geom_text(label = cog_hcv$name, show.legend = FALSE, size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_light() +
  labs(title = "(A)", x = "Regression_Coefficient_Cognition", y = "Regression_Coefficient_HCV") +
  scale_color_jcolors(palette = "pal8") +
  annotate("text", x = min(cog_hcv$Beta.x), y = max(cog_hcv$Beta.y),
           label = paste(cor_str, p_value_str, sep = ", "), hjust = 0, vjust = 1)

dev.off()

pdf(file="Correlation_plot_cognition_WML_Sup.Fig4D.pdf",width=10,height=8)                                                                                               
# Calculate correlation
correlation <- cor(cog_wml$Beta.x, cog_wml$Beta.y)
# Calculate p-value using cor.test()
cor_test <- cor.test(cog_wml$Beta.x, cog_wml$Beta.y)
p_value <- cor_test$p.value
# Convert correlation and p-value to strings
cor_str <- sprintf("r = %.2f", correlation)
p_value_str <- sprintf("p = %.2e", p_value)  # Display p-value in scientific notation
# Plot with annotations
ggplot(cog_wml, aes(x = Beta.x, y = Beta.y, color = factor(SUPER_PATHWAY))) +
  geom_point() +
  geom_text(label = cog_wml$name, show.legend = FALSE, size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_light() +
  labs(title = "(A)", x = "Regression_Coefficient_Cognition", y = "Regression_Coefficient_WML") +
  scale_color_jcolors(palette = "pal8") +
  annotate("text", x = min(cog_wml$Beta.x), y = max(cog_wml$Beta.y),
           label = paste(cor_str, p_value_str, sep = ", "), hjust = 0, vjust = 1)

dev.off()