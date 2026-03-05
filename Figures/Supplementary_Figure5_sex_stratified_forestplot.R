################################################################################
### Supplementary Figure 5: Relationship of metabolites strongly influenced by the gut microbiome (EV ≥ 5%) with metabolites nominally associated (P < 0.05) with general cognition 
################################################################################
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")

library(ggplot2)
library(dplyr)
library(gridExtra)
library(forcats)

### Load the files 
cog<-read.csv("Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated_95CI.csv",sep='\t')
female<-read.csv("Gfactor_age_sex_BMI_lipid_low_female.csv",sep='\t')
male<-read.csv("Gfactor_age_sex_BMI_lipid_low_male.csv",sep='\t')
### Load annotation file 
library("readxl")
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)

male_anno<-merge(male,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
female_anno<-merge(female,anno_com[,c("Name","CHEMICAL_NAME")],by.x="Metabolite",by.y="Name",all.x=T)
### Load the files 
mic<-read.csv("EV_microbiota_results.csv")
mic$fdr<-p.adjust(mic$spearman_p,method='fdr')
mic$EV_Microbiota<-ifelse((mic$fdr<0.05 & mic$explained_variance_score>0),mic$explained_variance_score*100,0)
mic_sub<-mic[mic$EV_Microbiota>=5,]

merg_all<-merge(mic_sub,cog,by.x="X",by.y='Metabolite',all.x=T)
merg_female<-merge(mic_sub,female_anno,by.x="X",by.y='Metabolite',all.x=T)
merg_male<-merge(mic_sub,male_anno,by.x="X",by.y='Metabolite',all.x=T)

merg_all<-merg_all[,c("CHEMICAL_NAME","EV_Microbiota","Beta","p","FDR","lower","upper")]
merg_female<-merg_female[,c("CHEMICAL_NAME","EV_Microbiota","Beta","p","FDR","lower","upper")]
merg_male<-merg_male[,c("CHEMICAL_NAME","EV_Microbiota","Beta","p","FDR","lower","upper")]

### Replace corrected names 
name_map <- c(
  "1-carboxyethyltyrosine" = "N-lactoyl tyrosine",
  "1-carboxyethylphenylalanine" = "N-lactoyl phenylalanine",
  "1-carboxyethylisoleucine" = "N-lactoyl isoleucine"
)

merg_all$CHEMICAL_NAME[merg_all$CHEMICAL_NAME %in% names(name_map)] <- 
  name_map[merg_all$CHEMICAL_NAME[merg_all$CHEMICAL_NAME %in% names(name_map)]]

merg_female$CHEMICAL_NAME[merg_female$CHEMICAL_NAME %in% names(name_map)] <- 
  name_map[merg_female$CHEMICAL_NAME[merg_female$CHEMICAL_NAME %in% names(name_map)]]

merg_male$CHEMICAL_NAME[merg_male$CHEMICAL_NAME %in% names(name_map)] <- 
  name_map[merg_male$CHEMICAL_NAME[merg_male$CHEMICAL_NAME %in% names(name_map)]]

pdf(file="Supplementary_Figure5_sex_stratified_analysis.pdf",
    width=20, height=15)

plot1<-merg_all %>% mutate(CHEMICAL_NAME = fct_reorder(CHEMICAL_NAME, desc(EV_Microbiota))) %>%     
  ggplot(aes(CHEMICAL_NAME, EV_Microbiota)) +
  geom_bar(stat = "identity", fill = "#CC79A7", width = 0.5) +
  coord_flip() + theme_light() 

plot2<-merg_all %>% mutate(CHEMICAL_NAME = fct_reorder(CHEMICAL_NAME, desc(EV_Microbiota))) %>%
  ggplot(aes(y = CHEMICAL_NAME, x = Beta)) +
  geom_point(shape = 16, size = 3) +  
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
  xlab("Beta (95% CI)-all") + 
  ylab(" ") +  theme_light() 

plot3<-merg_female %>% mutate(CHEMICAL_NAME = fct_reorder(CHEMICAL_NAME, desc(EV_Microbiota))) %>%
  ggplot(aes(y = CHEMICAL_NAME, x = Beta)) +
  geom_point(shape = 16, size = 3) +  
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
  xlab("Beta (95% CI)-female") + 
  ylab(" ") +  theme_light() 

plot4<-merg_male %>% mutate(CHEMICAL_NAME = fct_reorder(CHEMICAL_NAME, desc(EV_Microbiota))) %>%
  ggplot(aes(y = CHEMICAL_NAME, x = Beta)) +
  geom_point(shape = 16, size = 3) +  
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
  xlab("Beta (95% CI)-male") + 
  ylab(" ") +  theme_light() 


grid.arrange(plot2, plot4, plot3, plot1,ncol=4)
dev.off()
