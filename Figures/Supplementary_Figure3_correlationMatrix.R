####################################################################################
### Supplementary Figure 3: Correlation plot for the metabolites associated with cognition and MRI 
####################################################################################
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library("corrplot")
library("readxl")
## Load results of association of MRI and general cognition with metabolites 
mri<-read.csv("Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",sep='\t')
cog<-read.csv("Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated_95CI.csv",sep='\t')
## Unique metabolites associated with cognition and the three MRI markers
metabolites_cog<-cog$Metabolite[which(cog$FDR<0.05)]
metabolites_mri<-mri$Metabolite[which(mri$FDR<0.05)]
all_metabolites<-unique(c(metabolites_cog,metabolites_mri))
### Load metabolomics data for RSIII-2
load("KNN_imputed_metabolon_data_3rdMarch.RData")
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
pdf(file="SupplementaryFig3_correlationMatrix.pdf",width=18,height=8)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor_matrix, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", number.cex=0.4, # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 0.75,  #Text label color and rotation
         diag=FALSE 
)
dev.off()       
