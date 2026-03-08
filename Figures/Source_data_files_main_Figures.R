#################################################### 
## Code for Source Data for Main Figures 
#################################################### 
#################################################### 
## Figure 1
#################################################### 

.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
## Load R packages 
library("readxl")
library(tidyverse)
library(dplyr)
library(writexl)
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

### preserve the order of metabolites in the final pdf, as we changed the name of metabolite 1-carboxyethytyrosin ==> N-lactoyltyrosine
order_metabolites <- c("X - 26107","X - 25420","X - 24418","X - 11849","X - 11847","X - 11787","uridine","theophylline","sphingomyelin (d18:2/24:2)","sphingomyelin (d18:2/18:1)",
                       "sphingomyelin (d18:1/20:1, d18:2/20:0)","S-adenosylhomocysteine (SAH)","paraxanthine","o-cresol sulfate","glycerophosphorylcholine (GPC)",
                       "glutamine conjugate of C6H10O2 (2)","glutamine conjugate of C6H10O2 (1)","ergothioneine","cyclo(leu-pro)","caffeine","argininate","6-bromotryptophan",
                       "4-vinylguaiacol sulfate","4-vinylcatechol sulfate","3-methylcatechol sulfate","3-hydroxysebacate","3-hydroxyoleoylcarnitine","3-hydroxyhexanoylcarnitine (1)",
                       "3-hydroxy-2-methylpyridine sulfate","3-acetylphenol sulfate","2'-deoxyuridine","2-naphthol sulfate","1,3,7-trimethylurate","1,3-dimethylurate",
                       "N-lactoyltyrosine","(S)-3-hydroxybutyrylcarnitine")
cog_set$CHEMICAL_NAME <- factor(cog_set$CHEMICAL_NAME_new, levels = rev(order_metabolites))
tbv_set$CHEMICAL_NAME <- factor(tbv_set$CHEMICAL_NAME_new, levels = rev(order_metabolites))
hcv_set$CHEMICAL_NAME <- factor(hcv_set$CHEMICAL_NAME_new, levels = rev(order_metabolites))
wml_set$CHEMICAL_NAME <- factor(wml_set$CHEMICAL_NAME_new, levels = rev(order_metabolites))


# Write the source data file 
write_xlsx(
  list(
    Figure1a_cognition = cog_set,
    Figure1b_TBV = tbv_set,
    Figure1c_HCV = hcv_set,
    Figure1d_WML = wml_set
  ),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Figure1.xlsx"
)

#################################################### 
## Figure 2
#################################################### 

.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library("readxl")
library("dplyr")
## Load AD association results of RSI-IV
ad<-read.csv("AD_association_metabolites_age_sex_antilipid_BMI_model1_excStroke.csv",sep='\t')
## Load general cognition association results of RSIII-2
cog<-read.csv("Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated_95CI.csv",sep='\t')
## General cognition associated metabolites with p<0.05
metabolites_cog<-cog$Metabolite[which(cog$p<0.05)]
## Extract only rows with metabolites p<0.05 for association with cognition
cog_set<-cog[cog$Metabolite%in%metabolites_cog,]
ad_set<-ad[ad$metab_name%in%metabolites_cog,]
## combine both cognition and AD Beta matrix
cog_ad<-merge(cog_set[,c("Beta","Metabolite","CHEMICAL_NAME","SUPER_PATHWAY","FDR")],ad_set[,c("Beta","CHEMICAL_NAME","FDR")],by='CHEMICAL_NAME')
cog_ad$SUPER_PATHWAY[is.na(cog_ad$SUPER_PATHWAY)]<-"Unknowns"
colnames(cog_ad)<-c("CHEMICAL_NAME","Beta_cognition","Metabolite","SUPER_PATHWAY","FDR_cognition","Beta_incAD","FDR_incAD")
# consistent correct naming 3-methyl catechol sulfate (2) to 3-methylcatechol sulfate
cog_ad$CHEMICAL_NAME[cog_ad$CHEMICAL_NAME=='3-methyl catechol sulfate (2)']<-'3-methylcatechol sulfate'
# retain names of only top metabolites 
cog_ad$significant_metabolites<-NA
cog_ad$significant_metabolites<-ifelse(cog_ad$FDR_cognition<0.05 | cog_ad$FDR_incAD<0.05,cog_ad$CHEMICAL_NAME,"")
cog_ad<-cog_ad %>% select(Metabolite,CHEMICAL_NAME,SUPER_PATHWAY,Beta_cognition,FDR_cognition,Beta_incAD,FDR_incAD,significant_metabolites)
write_xlsx(
  list(
    Figure2 = cog_ad
  ),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Figure2.xlsx"
)


######################################################
### Figure 3
######################################################
### Figure 3A
cog<-read.csv("Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated_95CI.csv",sep='\t')
sig<-cog$Metabolite[cog$FDR<0.05]
####### Load annotation file 
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$names<-paste0("metab_",anno_com$CHEM_ID)
anno_com<-anno_com[,c("names","CHEMICAL_NAME")]
######
select<-anno_com[anno_com$names%in%sig,]
######
mic<-read.csv("EV_microbiota_results.csv")
med<-read.csv("EV_medication_results.csv")
demo<-read.csv("EV_lifestyle_results.csv")
como<-read.csv("EV_comborbidity_results.csv")
gene<-read.csv("EV_gene_results.csv")

mic$fdr<-p.adjust(mic$spearman_p,method='fdr')
med$fdr<-p.adjust(med$spearman_p,method='fdr')
demo$fdr<-p.adjust(demo$spearman_p,method='fdr')
como$fdr<-p.adjust(como$spearman_p,method='fdr')
gene$fdr<-p.adjust(gene$spearman_p,method='fdr')


mic$EV_Microbiota<-ifelse((mic$fdr<0.05 & mic$explained_variance_score>0),mic$explained_variance_score*100,0)
med$EV_Medication<-ifelse((med$fdr<0.05 & med$explained_variance_score>0),med$explained_variance_score*100,0)
demo$EV_Lifestyle<-ifelse((demo$fdr<0.05 & demo$explained_variance_score>0),demo$explained_variance_score*100,0)
como$EV_Clinical<-ifelse((como$fdr<0.05 & como$explained_variance_score>0),como$explained_variance_score*100,0)
gene$EV_Genetics<-ifelse((gene$fdr<0.05 & gene$explained_variance_score>0),gene$explained_variance_score*100,0)


merg1<-merge(select,med[,c(1,10)],by.x="names",by.y="X")
merg2<-merge(merg1,demo[,c(1,10)],by.x="names",by.y="X")
merg3<-merge(merg2,gene[,c(1,10)],by.x="names",by.y="X")
merg4<-merge(merg3,mic[,c(1,10)],by.x="names",by.y="X")
merg5<-merge(merg4,como[,c(1,10)],by.x="names",by.y="X")

final<-merg5[,c(1, 2,5,6,4,7,3)]
# consistent correct naming 3-methyl catechol sulfate (2) to 3-methylcatechol sulfate
final$CHEMICAL_NAME[final$CHEMICAL_NAME=='3-methyl catechol sulfate (2)']<-'3-methylcatechol sulfate'

Figure3A<-final[order(final$EV_Genetics,decreasing = TRUE),]

### Figure 3B
mri<-read.csv("Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",sep='\t')
sig<-mri$Metabolite[mri$FDR<0.05]
### Load annotation file 
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$names<-paste0("metab_",anno_com$CHEM_ID)
anno_com<-anno_com[,c("names","CHEMICAL_NAME")]
### Select only significant metabolites 
select<-anno_com[anno_com$names%in%sig,]
### Load the result files of explained variance by gut microbiota, medication, lifestyle, clinical factors and genes 
mic<-read.csv("EV_microbiota_results.csv")
med<-read.csv("EV_medication_results.csv")
demo<-read.csv("EV_lifestyle_results.csv")
como<-read.csv("EV_comborbidity_results.csv")
gene<-read.csv("EV_gene_results.csv")


mic$fdr<-p.adjust(mic$spearman_p,method='fdr')
med$fdr<-p.adjust(med$spearman_p,method='fdr')
demo$fdr<-p.adjust(demo$spearman_p,method='fdr')
como$fdr<-p.adjust(como$spearman_p,method='fdr')
gene$fdr<-p.adjust(gene$spearman_p,method='fdr')

mic$EV_Microbiota<-ifelse((mic$fdr<0.05 & mic$explained_variance_score>0),mic$explained_variance_score*100,0)
med$EV_Medication<-ifelse((med$fdr<0.05 & med$explained_variance_score>0),med$explained_variance_score*100,0)
demo$EV_Lifestyle<-ifelse((demo$fdr<0.05 & demo$explained_variance_score>0),demo$explained_variance_score*100,0)
como$EV_Clinical<-ifelse((como$fdr<0.05 & como$explained_variance_score>0),como$explained_variance_score*100,0)
gene$EV_Genetics<-ifelse((gene$fdr<0.05 & gene$explained_variance_score>0),gene$explained_variance_score*100,0)


merg1<-merge(select,med[,c(1,10)],by.x="names",by.y="X")
merg2<-merge(merg1,demo[,c(1,10)],by.x="names",by.y="X")
merg3<-merge(merg2,gene[,c(1,10)],by.x="names",by.y="X")
merg4<-merge(merg3,mic[,c(1,10)],by.x="names",by.y="X")
merg5<-merge(merg4,como[,c(1,10)],by.x="names",by.y="X")

final<-merg5[,c(1, 2,5,6,4,7,3)]
final$CHEMICAL_NAME[final$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyltyrosine"
Figure3B<-final[order(final$EV_Genetics,decreasing = TRUE),]

write_xlsx(
  list(
    Figure3a = Figure3A,
    Figure3b = Figure3B
  ),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Figure3.xlsx"
)


#################################################### 
## Figure 4
#################################################### 
mic<-read.csv("EV_microbiota_results.csv")
med<-read.csv("EV_medication_results.csv")
demo<-read.csv("EV_lifestyle_results.csv")
como<-read.csv("EV_comborbidity_results.csv")
gene<-read.csv("EV_gene_results.csv")

mic$fdr<-p.adjust(mic$spearman_p,method='fdr')
med$fdr<-p.adjust(med$spearman_p,method='fdr')
demo$fdr<-p.adjust(demo$spearman_p,method='fdr')
como$fdr<-p.adjust(como$spearman_p,method='fdr')
gene$fdr<-p.adjust(gene$spearman_p,method='fdr')


mic_sig<-mic[mic$explained_variance_score>0 & mic$fdr<0.05,]
med_sig<-med[med$explained_variance_score>0 & med$fdr<0.05,]
demo_sig<-demo[demo$explained_variance_score>0 & demo$fdr<0.05,]
como_sig<-como[como$explained_variance_score>0 & como$fdr<0.05,]
gene_sig<-gene[gene$explained_variance_score>0 & gene$fdr<0.05,]

mic_sig$Features<-rep("1_Microbiota",length(mic_sig$explained_variance_score))
med_sig$Features<-rep("2_Medication",length(med_sig$explained_variance_score))
como_sig$Features<-rep("3_Clinical",length(como_sig$explained_variance_score))
demo_sig$Features<-rep("4_Demographics",length(demo_sig$explained_variance_score))
gene_sig$Features<-rep("5_Genetics",length(gene_sig$explained_variance_score))

combined<-rbind(mic_sig[,c("Features","explained_variance_score","X")],med_sig[,c("Features","explained_variance_score","X")], demo_sig[,c("Features","explained_variance_score","X")],como_sig[,c("Features","explained_variance_score","X")], gene_sig[,c("Features","explained_variance_score","X")])
combined$Percent_variance_explained<-combined$explained_variance_score*100

## Load general cognition and MRI association results to tag significant metabolites in the Beeswarm plot
mri<-read.csv("Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",sep='\t')
cog<-read.csv("Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated_95CI.csv",sep='\t')

## Unique metabolites 
metabolites_cog<-cog$Metabolite[which(cog$FDR<0.05)]
metabolites_mri<-mri$Metabolite[which(mri$FDR<0.05)]
all_metabolites<-unique(c(metabolites_cog,metabolites_mri))
################################### 
combined$color<-NA
combined$color<-ifelse(combined$X%in%all_metabolites,'red','grey')
colnames(combined)[colnames(combined) == "X"] <- "Metabolites"

#################################################### 
## Main Figure 4B: UpSet plot to show the number of metabolites explained by different determinants
#################################################### 
mic<-read.csv("EV_microbiota_results.csv")
med<-read.csv("EV_medication_results.csv")
demo<-read.csv("EV_lifestyle_results.csv")
como<-read.csv("EV_comborbidity_results.csv")
gene<-read.csv("EV_gene_results.csv")


mic$fdr<-p.adjust(mic$spearman_p,method='fdr')
med$fdr<-p.adjust(med$spearman_p,method='fdr')
demo$fdr<-p.adjust(demo$spearman_p,method='fdr')
como$fdr<-p.adjust(como$spearman_p,method='fdr')
gene$fdr<-p.adjust(gene$spearman_p,method='fdr')


mic_sig<-mic[mic$explained_variance_score>0 & mic$fdr<0.05,]
med_sig<-med[med$explained_variance_score>0 & med$fdr<0.05,]
demo_sig<-demo[demo$explained_variance_score>0 & demo$fdr<0.05,]
como_sig<-como[como$explained_variance_score>0 & como$fdr<0.05,]
gene_sig<-gene[gene$explained_variance_score>0 & gene$fdr<0.05,]

mic_sig$Features<-rep("1_Microbiota",length(mic_sig$explained_variance_score))
med_sig$Features<-rep("2_Medication",length(med_sig$explained_variance_score))
como_sig$Features<-rep("3_Clinical",length(como_sig$explained_variance_score))
demo_sig$Features<-rep("4_Lifestyle",length(demo_sig$explained_variance_score))
gene_sig$Features<-rep("5_Genetics",length(gene_sig$explained_variance_score))

all_sig <- rbind(
  mic_sig,
  med_sig,
  como_sig,
  demo_sig,
  gene_sig
)
all_sig <- all_sig[, c("Features","X","explained_variance_score","spearman_p","fdr")]
colnames(all_sig)[colnames(all_sig) == "X"] <- "Metabolites"

write_xlsx(
  list(
    Figure4a = combined,
    Figure4b = all_sig
  ),
  "/Users/sahmad1/Downloads/SOURCE_FILES/source_data_main_figures/Source_data_Figure4.xlsx"
)

