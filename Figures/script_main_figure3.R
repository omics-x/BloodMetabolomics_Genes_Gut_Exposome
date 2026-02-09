#################################################### 
## Code for Main Figure 3: Barplots for Explained Variance (EV) for metabolites associated with cognition and MRI markers 
#################################################### 
### Load R Packages 
### install.packages("stringi")
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library("readxl")


######################################################
### Figure 3A: Explained variance (EV) of all metabolites associated with general cognition
######################################################
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

final<-merg5[,c(2,5,6,4,7,3)]

final<-final[order(final$EV_Genetics,decreasing = TRUE),]


Df <- final %>% 
  gather(keys, values, EV_Genetics:EV_Medication)
# To define the order of determinants
Df$keys_factor = factor(Df$keys, levels=c('EV_Genetics','EV_Lifestyle','EV_Medication','EV_Clinical','EV_Microbiota'))
# Order of metabolites in the final figure 
order_metabolites <- c("X - 25420","X - 11849","X - 11847","uridine","o-cresol sulfate","ergothioneine","4-vinylguaiacol sulfate",
  "4-vinylcatechol sulfate","3-methyl catechol sulfate (2)","3-hydroxy-2-methylpyridine sulfate","3-acetylphenol sulfate",
  "2-naphthol sulfate","X - 24418","2'-deoxyuridine")

Df$CHEMICAL_NAME <- factor(Df$CHEMICAL_NAME, levels = rev(order_metabolites))

pdf(file="EV_metabolites_associated_with_General_cognition.pdf",width=12,height=5)
colors <- c("#999999", "#E69F00", "#0072B2", "#D55E00", "#CC79A7")
Df %>% mutate(CHEMICAL_NAME = CHEMICAL_NAME, desc(values)) %>%     
  ggplot(aes(CHEMICAL_NAME, values)) +
  geom_col(aes(fill = keys_factor), width = 0.5) +
  scale_fill_manual(values = colors) +
  facet_wrap(~ keys_factor, ncol = 5) +
  coord_flip() + theme_light() 
dev.off()
    
######################################################
### Figure 3B. Load MRI association results Model 1 
######################################################
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

final<-merg5[,c(2,5,6,4,7,3)]

final<-final[order(final$EV_Genetics,decreasing = TRUE),]
Df <- final %>% 
  gather(keys, values, EV_Genetics:EV_Medication)
# To define the order of determinants
Df$keys_factor = factor(Df$keys, levels=c('EV_Genetics','EV_Lifestyle','EV_Medication','EV_Clinical','EV_Microbiota'))

Df$CHEMICAL_NAME[Df$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyl tyrosine"

order_metabolites_mri <- c(
  "X - 26107",
  "X - 11787",
  "theophylline",
  "sphingomyelin (d18:1/20:1, d18:2/20:0)*",
  "S-adenosylhomocysteine (SAH)",
  "paraxanthine",
  "cyclo(leu-pro)",
  "argininate*",
  "6-bromotryptophan",
  "3-hydroxyoleoylcarnitine",
  "3-hydroxyhexanoylcarnitine (1)",
  "1,3,7-trimethylurate",
  "1,3-dimethylurate",
  "(S)-3-hydroxybutyrylcarnitine",
  "glutamine conjugate of C6H10O2 (1)*",
  "caffeine",
  "glycerophosphorylcholine (GPC)",
  "glutamine conjugate of C6H10O2 (2)*",
  "sphingomyelin (d18:2/18:1)*",
  "3-hydroxysebacate",
  "sphingomyelin (d18:2/24:2)*",
  "N-lactoyl tyrosine"
)
Df$CHEMICAL_NAME <- factor(Df$CHEMICAL_NAME, levels = rev(order_metabolites_mri))

pdf(file="EV_metabolites_associated_with_MRI.pdf",width=12,height=5)
colors <- c("#999999", "#E69F00", "#0072B2", "#D55E00", "#CC79A7")
Df %>% mutate(CHEMICAL_NAME = CHEMICAL_NAME, desc(values)) %>%     
  ggplot(aes(CHEMICAL_NAME, values)) +
  geom_col(aes(fill = keys_factor), width = 0.5) +
  scale_fill_manual(values = colors) +
  facet_wrap(~ keys_factor, ncol = 5) +
  coord_flip() + theme_light() 
dev.off()

# In pdf, we corrected the annotation for 1-carboxyethyltyrosine to N-lactoyl tyrosine