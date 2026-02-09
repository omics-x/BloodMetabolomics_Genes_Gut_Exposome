#################################################### 
## Code for Main Figure 1: Forest plot for the association of metabolites with general cognition and MRI markers
#################################################### 
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
## Load R packages 
library("readxl")
library(tidyverse)
library(gridExtra)
library(dplyr)

## Load Result files 
## Load results of Association of metabolites with MRI markers (M1)
mri<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",sep='\t')
## Load results of Association of metabolites with general cognition (M1)
cog<-read.csv("/Users/sahmad1/Downloads/SOURCE_FILES/Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated_95CI.csv",sep='\t')

## Unique metabolites with FDR < 0.05
metabolites_cog<-cog$Metabolite[which(cog$FDR<0.05)]
metabolites_mri<-mri$Metabolite[which(mri$FDR<0.05)]
all_metabolites<-unique(c(metabolites_cog,metabolites_mri))


### seperate the MRI variables (Total_par_ml = Total Brain volume)
tbv<-mri[mri$endo_pheno=="Total_par_ml",]
hcv<-mri[mri$endo_pheno=="total_Hippocampus",]
wml<-mri[mri$endo_pheno=="Total_wml",]

### Load annotation file for metabolomics data RSIII-2
anno_com <- read_excel("/Users/sahmad1/Downloads/SOURCE_FILES/DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)

### combine the result files with annotation (new chemical name corrects the N-lactoyl- metabolite names)
tbv<-merge(tbv,anno_com[,c("Name","CHEMICAL_NAME_new")],by.x="Metabolite",by.y="Name",all.x=T)
hcv<-merge(hcv,anno_com[,c("Name","CHEMICAL_NAME_new")],by.x="Metabolite",by.y="Name",all.x=T)
wml<-merge(wml,anno_com[,c("Name","CHEMICAL_NAME_new")],by.x="Metabolite",by.y="Name",all.x=T)
cog<-merge(cog[,c("Metabolite","Beta","Se","p","lower","upper","FDR")],anno_com[,c("Name","CHEMICAL_NAME_new")],by.x="Metabolite",by.y="Name",all.x=T)

cog_set<-cog[cog$Metabolite%in%all_metabolites,]
tbv_set<-tbv[tbv$Metabolite%in%all_metabolites,]
hcv_set<-hcv[hcv$Metabolite%in%all_metabolites,]
wml_set<-wml[wml$Metabolite%in%all_metabolites,]

#### Color coding 
metabolites_cog<-cog$Metabolite[which(cog$FDR<0.05)]
metabolites_tbv<-tbv_set$Metabolite[which(tbv_set$FDR<0.05)]
metabolites_hcv<-hcv_set$Metabolite[which(hcv_set$FDR<0.05)]
metabolites_wml<-wml_set$Metabolite[which(wml_set$FDR<0.05)]

cog_set$color<-NA
cog_set$color[cog_set$Metabolite%in%metabolites_cog]<-"red"
cog_set$color[!cog_set$Metabolite%in%metabolites_cog]<-"black"

tbv_set$color<-NA
tbv_set$color[tbv_set$Metabolite%in%metabolites_tbv]<-"red"
tbv_set$color[!tbv_set$Metabolite%in%metabolites_tbv]<-"black"

hcv_set$color<-NA
hcv_set$color[hcv_set$Metabolite%in%metabolites_hcv]<-"red"
hcv_set$color[!hcv_set$Metabolite%in%metabolites_hcv]<-"black"

wml_set$color<-NA
wml_set$color[wml_set$Metabolite%in%metabolites_wml]<-"red"
wml_set$color[!wml_set$Metabolite%in%metabolites_wml]<-"black"

### preserve the order of metabolites in the final pdf, as we changed the name of metabolite 1-carboxyethytyrosin ==> N-lactoyl tyrosine
order_metabolites <- c("X - 26107","X - 25420","X - 24418","X - 11849","X - 11847","X - 11787","uridine","theophylline","sphingomyelin (d18:2/24:2)","sphingomyelin (d18:2/18:1)",
                       "sphingomyelin (d18:1/20:1, d18:2/20:0)","S-adenosylhomocysteine (SAH)","paraxanthine","o-cresol sulfate","glycerophosphorylcholine (GPC)",
                       "glutamine conjugate of C6H10O2 (2)","glutamine conjugate of C6H10O2 (1)","ergothioneine","cyclo(leu-pro)","caffeine","argininate","6-bromotryptophan",
                       "4-vinylguaiacol sulfate","4-vinylcatechol sulfate","3-methyl catechol sulfate (2)","3-hydroxysebacate","3-hydroxyoleoylcarnitine","3-hydroxyhexanoylcarnitine (1)",
                       "3-hydroxy-2-methylpyridine sulfate","3-acetylphenol sulfate","2'-deoxyuridine","2-naphthol sulfate","1,3,7-trimethylurate","1,3-dimethylurate",
                       "N-lactoyl tyrosine","(S)-3-hydroxybutyrylcarnitine")
cog_set$CHEMICAL_NAME <- factor(cog_set$CHEMICAL_NAME_new, levels = rev(order_metabolites))
tbv_set$CHEMICAL_NAME <- factor(tbv_set$CHEMICAL_NAME_new, levels = rev(order_metabolites))
hcv_set$CHEMICAL_NAME <- factor(hcv_set$CHEMICAL_NAME_new, levels = rev(order_metabolites))
wml_set$CHEMICAL_NAME <- factor(wml_set$CHEMICAL_NAME_new, levels = rev(order_metabolites))

### Plotting and saving results of general cognition, and MRI markers 
pdf(file="Forest_plot_overlap_cognition_MRI_metabolites_associations_FDR_0.05_color.pdf",width=18,height=8)
### Variance explained 
plot1<-ggplot(cog_set, aes(y = CHEMICAL_NAME, x = Beta)) +
geom_point(shape = 16, size = 3, colour = factor(cog_set$color)) +  
geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) +
geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
xlab("Beta (95% CI)-Cognition") + 
ylab(" ") +  theme_light() 

plot2<-ggplot(tbv_set, aes(y = CHEMICAL_NAME, x = Beta)) +
geom_point(shape = 16, size = 3, colour = factor(tbv_set$color)) +  
geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) +
geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
xlab("Beta (95% CI)-TBV") + 
ylab(" ") +  theme_light() 

plot3<-ggplot(hcv_set, aes(y = CHEMICAL_NAME, x = Beta)) +
geom_point(shape = 16, size = 3, colour = factor(hcv_set$color)) +  
geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) +
geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
xlab("Beta (95% CI)-HCV") + 
ylab(" ") +  theme_light()

plot4<-
ggplot(wml_set, aes(y = CHEMICAL_NAME, x = Beta)) +
geom_point(shape = 16, size = 3, colour = factor(wml_set$color)) +  
geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) +
geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
xlab("Beta (95% CI)-WML") + 
ylab(" ") +  theme_light() 
grid.arrange(plot1, plot2, plot3, plot4,ncol=4)
dev.off()


# In pdf, we corrected the annotation for 1-carboxyethyltyrosine to N-lactoyl tyrosine