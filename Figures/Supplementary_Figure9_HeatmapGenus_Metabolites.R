#########################################################################
### Code for Supplementary Figure 9: Relationship of general cognition and brain imaging associated metabolites with gut microbiota
#########################################################################
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library("readxl")
library(gplots)

#########################################################################
### Supplementary Figure 9A
#########################################################################

cog<-read.csv("Gfactor_Association_metabolites_age_sex_antilipid_BMI_M1_annotated_95CI.csv",sep='\t')
metabolites<-cog$Metabolite[which(cog$FDR<0.05)]

### Load results of association of gut-microbiome ASVs with metabolomics
result<-read.table("Association_analysis_ful_model_metabolon_gut_microbiota_Nov1_clr.txt",head=T)
### Load Annotation file 
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)
res_anno<-merge(result,anno_com[,c("Name","CHEMICAL_NAME","SUPER_PATHWAY","SUB_PATHWAY")],by.x="Pheno",by.y="Name",all.x=T)

### FDR correction 
res_anno$FDR<-p.adjust(res_anno[,5], method = 'fdr', n = length(res_anno[,4]))
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
res_anno<-res_anno[grep("genus",res_anno$Taxon),]
### Generate files for heatmap 
res_anno<-res_anno[which(res_anno$Pheno%in%metabolites),]
pvalue<-res_anno[,c("CHEMICAL_NAME","Taxon","FDR")]
effect<-res_anno[,c("CHEMICAL_NAME","Taxon","zvalue")]
effect1 <- xtabs( zvalue ~ as.character(CHEMICAL_NAME) + as.character(Taxon) , data=effect)
pvalue1 <- xtabs( FDR ~ as.character(CHEMICAL_NAME) + as.character(Taxon) , data=pvalue)
effect2<-as.data.frame.matrix(effect1)
pvalue2<-as.data.frame.matrix(pvalue1)
cors <- data.matrix (effect2)
pcors <- data.matrix(pvalue2)
vecCol<-numeric()
for(i in 1:dim(pcors)[2]){if(min(pcors[,i])>0.05){vecCol<-c(vecCol,i)}}
cors<-cors[,-vecCol]
pcors<-pcors[,-vecCol]
### color code
data_colors<-cors
data_colors[is.nan(data_colors)] = 0
data_colors[is.na(data_colors)] = 0
max<-max(c(abs(min(data_colors)),abs(max(data_colors))))
palette.breaks <- seq(-max, max, length.out=300)
color.palette  <- colorRampPalette(c("darkblue", "white", "red"))(length(palette.breaks) - 1)
### Significance code inside heatmap
pcors<--log10(pcors)
pcors[pcors>=1.3] = "*"
pcors[pcors>=0] = " "
pcors[is.na(pcors)] = "NA"
pcors<-pcors[ match(rownames(cors), rownames(pcors)), match(colnames(cors), colnames(pcors))]

pdf("Heatmaps_metabolite_cognition_cluster_pathways2.pdf", height=8, width=8)
heatmap.2(
  cors,
  Rowv = TRUE,
  Colv = TRUE,
  distfun = function(x) dist(x,method = 'euclidian'),
  hclustfun = hclust,
  colsep=1:ncol(cors),
  rowsep=1:nrow(cors),
  sepcolor="grey",
  cellnote= pcors, 
  notecex=1.2,
  notecol="black",
  lhei = c(1,5),
  lwid =  c(1,5),
  keysize=0.6,
  scale      = "none",  
  trace      = "none",  
  key        = TRUE, 
  density.info = "none",
  col    = color.palette,
  breaks= palette.breaks,
  margins =  c(14,14), 
  srtCol =65,
  cexRow=1, 
  cexCol=0.9)
dev.off()

#########################################################################
### Supplementary Figure 9B
#########################################################################
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library("readxl")
library(gplots)

### Load results of association of MRI markers with metabolomics
mri<-read.csv("Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",sep='\t')
metabolites<-mri$Metabolite[which(mri$FDR<0.05)]

### Load results of association of gut-microbiome ASVs with metabolomics
result<-read.table("Association_analysis_ful_model_metabolon_gut_microbiota_Nov1_clr.txt",head=T)
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)
res_anno<-merge(result,anno_com[,c("Name","CHEMICAL_NAME","SUPER_PATHWAY","SUB_PATHWAY")],by.x="Pheno",by.y="Name",all.x=T)

res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyl tyrosine"
### FDR correction 
res_anno$FDR<-p.adjust(res_anno[,5], method = 'fdr', n = length(res_anno[,4]))
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
res_anno<-res_anno[grep("genus",res_anno$Taxon),]
### Plot files 
res_anno<-res_anno[which(res_anno$Pheno%in%metabolites),]
pvalue<-res_anno[,c("CHEMICAL_NAME","Taxon","FDR")]
effect<-res_anno[,c("CHEMICAL_NAME","Taxon","zvalue")]
effect1 <- xtabs( zvalue ~ as.character(CHEMICAL_NAME) + as.character(Taxon) , data=effect)
pvalue1 <- xtabs( FDR ~ as.character(CHEMICAL_NAME) + as.character(Taxon) , data=pvalue)
effect2<-as.data.frame.matrix(effect1)
pvalue2<-as.data.frame.matrix(pvalue1)
cors <- data.matrix (effect2)
pcors <- data.matrix(pvalue2)
vecCol<-numeric()
for(i in 1:dim(pcors)[2]){if(min(pcors[,i])>0.05){vecCol<-c(vecCol,i)}}
cors<-cors[,-vecCol]
pcors<-pcors[,-vecCol]
### color codes 
data_colors<-cors
data_colors[is.nan(data_colors)] = 0
data_colors[is.na(data_colors)] = 0
max<-max(c(abs(min(data_colors)),abs(max(data_colors))))
palette.breaks <- seq(-max, max, length.out=300)
color.palette  <- colorRampPalette(c("darkblue", "white", "red"))(length(palette.breaks) - 1)
### annotation of heatmap for significance 
pcors<--log10(pcors)
pcors[pcors>=1.3] = "*"
pcors[pcors>=0] = " "
pcors[is.na(pcors)] = "NA"
pcors<-pcors[ match(rownames(cors), rownames(pcors)), match(colnames(cors), colnames(pcors))]

pdf("Heatmaps_metabolite_MRI_microbiota_association.pdf", height=8, width=8)
heatmap.2(
  cors,
  Rowv = TRUE,
  Colv = TRUE,
  distfun = function(x) dist(x,method = 'euclidian'),
  hclustfun = hclust,
  colsep=1:ncol(cors),
  rowsep=1:nrow(cors),
  sepcolor="grey",
  cellnote= pcors, 
  notecex=1.2,
  notecol="black",
  lhei = c(1,5),
  lwid =  c(1,5),
  keysize=0.6,
  scale      = "none",  
  trace      = "none",  
  key        = TRUE, 
  density.info = "none",
  col    = color.palette,
  breaks= palette.breaks,
  margins =  c(16,16), 
  srtCol =65,
  cexRow=1, 
  cexCol=0.9)
dev.off()