######################################################################### 
### Supplementary Figure 8A: Heatmap of MRI associated metabolites with medication use 
#########################################################################
library(gplots)

mri<-read.csv("Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",head=T,sep='\t')
metabolites<-mri$Metabolite[which(mri$FDR<0.05)]
### Load association results of medication use with metabolomics 
res_anno<-read.csv("Association_analysis_medication.csv",head=T)
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
res_anno$medication_names<-str_trim(res_anno$name)
## P-value and regression coefficient matrix 
res_anno<-res_anno[which(res_anno$Metabolite%in%metabolites),]
## Correct the name of metabolite
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyl tyrosine"
## P-value and Beta matrix 
pvalue<-res_anno[,c("CHEMICAL_NAME","medication_names","FDR")]
effect<-res_anno[,c("CHEMICAL_NAME","medication_names","zvalue")]
effect1 <- xtabs( zvalue ~ as.character(CHEMICAL_NAME) + as.character(medication_names) , data=effect)
pvalue1 <- xtabs( FDR ~ as.character(CHEMICAL_NAME) + as.character(medication_names) , data=pvalue)
effect2<-as.data.frame.matrix(effect1)
pvalue2<-as.data.frame.matrix(pvalue1)

### prepare data for heatmap
cors <- data.matrix (effect2)
pcors <- data.matrix(pvalue2)
vecCol<-numeric()
for(i in 1:dim(pcors)[2]){if(min(pcors[,i])>0.05){vecCol<-c(vecCol,i)}}
cors<-cors[,-vecCol]
pcors<-pcors[,-vecCol]
data_colors<-cors
data_colors[is.nan(data_colors)] = 0
data_colors[is.na(data_colors)] = 0
max<-max(c(abs(min(data_colors)),abs(max(data_colors))))
palette.breaks <- seq(-max, max, length.out=300)
color.palette  <- colorRampPalette(c("darkblue", "white", "red"))(length(palette.breaks) - 1)
pcors<--log10(pcors)
pcors[pcors>=1.3] = "*"
pcors[pcors>=0] = " "
pcors[is.na(pcors)] = "NA"
pcors<-pcors[ match(rownames(cors), rownames(pcors)), match(colnames(cors), colnames(pcors))]
pdf("Heatmaps_metabolite_MRI_overlap_medication.pdf", height=8, width=6)
heatmap.2(
  cors,
  Rowv = TRUE,
  Colv = TRUE,
  distfun = function(x) dist(x,method = 'euclidian'),
  hclustfun = hclust,
  cellnote= pcors, 
  notecex=0.8,
  notecol="black",
  lhei = c(1,5),
  lwid =  c(2,5),
  
  scale      = "none",  
  trace      = "none",  
  key        = TRUE, 
  density.info = "none",
  col    = color.palette,
  breaks= palette.breaks,
  margins =  c(10,12), 
  srtCol =65,
  
  cexRow=0.8, 
  cexCol=0.8)
dev.off()

######################################################################### 
### Supplementary Figure 8B: Heatmap of MRI associated metabolites with lifestyle factors  
#########################################################################
library(stringr)
library(gplots)

mri<-read.csv("Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",head=T,sep='\t')
metabolites<-mri$Metabolite[which(mri$FDR<0.05)]
### Load results
res_anno<-read.csv("Association_analysis_demographics.csv",head=T)
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
### Prepare files 
res_anno<-res_anno[which(res_anno$Metabolite%in%metabolites),]
## Correct the name of metabolite
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyl tyrosine"
## P-value and Beta matrix 
pvalue<-res_anno[,c("CHEMICAL_NAME","Demographics","FDR")]
effect<-res_anno[,c("CHEMICAL_NAME","Demographics","zvalue")]
effect1 <- xtabs( zvalue ~ as.character(CHEMICAL_NAME) + as.character(Demographics) , data=effect)
pvalue1 <- xtabs( FDR ~ as.character(CHEMICAL_NAME) + as.character(Demographics) , data=pvalue)
effect2<-as.data.frame.matrix(effect1)
pvalue2<-as.data.frame.matrix(pvalue1)
colnames(effect2)<-c("Alcohol_intake","BMI","Education","Smoking")
colnames(pvalue2)<-c("Alcohol_intake","BMI","Education","Smoking")
### Plot healtmap 
cors <- data.matrix (effect2)
pcors <- data.matrix(pvalue2)
### color of z-values
data_colors<-cors
data_colors[is.nan(data_colors)] = 0
data_colors[is.na(data_colors)] = 0
max<-max(c(abs(min(data_colors)),abs(max(data_colors))))
palette.breaks <- seq(-max, max, length.out=300)
color.palette  <- colorRampPalette(c("darkblue", "white", "red"))(length(palette.breaks) - 1)
pcors<--log10(pcors)
pcors[pcors>=1.3] = "*"
pcors[pcors>=0] = " "
pcors[is.na(pcors)] = "NA"
pcors<-pcors[ match(rownames(cors), rownames(pcors)), match(colnames(cors), colnames(pcors))]
pdf("Heatmaps_metabolite_MRI_overlap_lifestyle_factors.pdf", height=8, width=6)
heatmap.2(
  cors,
  Rowv = TRUE,
  Colv = TRUE,
  distfun = function(x) dist(x,method = 'euclidian'),
  hclustfun = hclust,
  cellnote= pcors, 
  notecex=0.8,
  notecol="black",
  lhei = c(1,5),
  lwid =  c(2,5),
  scale      = "none",  
  trace      = "none",  
  key        = TRUE, 
  density.info = "none",
  col    = color.palette,
  breaks= palette.breaks,
  margins =  c(10,12), 
  srtCol =65,
  
  cexRow=0.8, 
  cexCol=0.8)
dev.off()

######################################################################### 
### Supplementary Figure 8C: Heatmap of MRI associated metabolites with clinical factors  
#########################################################################
library(stringr)
library(gplots)

mri<-read.csv("Association_Metabolon_fullmodel_excludingStroke_AD_RS1_5_m1.csv",head=T,sep='\t')
metabolites<-mri$Metabolite[which(mri$FDR<0.05)]
### Load results
res_anno<-read.csv("Association_analysis_clinical.csv",head=T)
res_anno$zvalue<-(res_anno$Beta/res_anno$Se)
### Prepare heatmap files 
res_anno<-res_anno[which(res_anno$Metabolite%in%metabolites),]
## Correct the name of metabolite
res_anno$CHEMICAL_NAME[res_anno$CHEMICAL_NAME == "1-carboxyethyltyrosine"] <- "N-lactoyl tyrosine"

pvalue<-res_anno[,c("CHEMICAL_NAME","Demographics","FDR")]
effect<-res_anno[,c("CHEMICAL_NAME","Demographics","zvalue")]
effect1 <- xtabs( zvalue ~ as.character(CHEMICAL_NAME) + as.character(Demographics) , data=effect)
pvalue1 <- xtabs( FDR ~ as.character(CHEMICAL_NAME) + as.character(Demographics) , data=pvalue)
effect2<-as.data.frame.matrix(effect1)
pvalue2<-as.data.frame.matrix(pvalue1)
colnames(effect2)<-c("Diabetes","DiasBP","Hypertension","SysBP")
colnames(pvalue2)<-c("Diabetes","DiasBP","Hypertension","SysBP")
cors <- data.matrix (effect2)
pcors <- data.matrix(pvalue2)
### Color code based on correlation
data_colors<-cors
data_colors[is.nan(data_colors)] = 0
data_colors[is.na(data_colors)] = 0
max<-max(c(abs(min(data_colors)),abs(max(data_colors))))
palette.breaks <- seq(-max, max, length.out=300)
color.palette  <- colorRampPalette(c("darkblue", "white", "red"))(length(palette.breaks) - 1)
### P-value annotation
pcors<--log10(pcors)
pcors[pcors>=1.3] = "*"
pcors[pcors>=0] = " "
pcors[is.na(pcors)] = "NA"
pcors<-pcors[ match(rownames(cors), rownames(pcors)), match(colnames(cors), colnames(pcors))]

pdf("Heatmaps_metabolite_MRI_overlap_clinical_factors.pdf", height=8, width=6)
heatmap.2(
  cors,
  Rowv = TRUE,
  Colv = TRUE,
  distfun = function(x) dist(x,method = 'euclidian'),
  hclustfun = hclust,
  cellnote= pcors, 
  notecex=0.8,
  notecol="black",
  lhei = c(1,5),
  lwid =  c(2,5),
  
  scale      = "none",  
  trace      = "none",  
  key        = TRUE, 
  density.info = "none",
  col    = color.palette,
  breaks= palette.breaks,
  margins =  c(10,12), 
  srtCol =65,
  cexRow=0.8, 
  cexCol=0.8)
dev.off()
