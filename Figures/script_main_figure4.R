#################################################### 
## Code for Main Figure 4: EV plots comparing different class of determinants of metabolites 
#################################################### 
## Load R packages 
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library(beeswarm)
library(UpSetR)

#################################################### 
## Main Figure 4A: Beeswarm plot to compare EV between different determinants
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

pdf("Beeswarm plot_EV_all_features_Figure4A.pdf")

beeswarm(Percent_variance_explained ~ Features, data=combined, pch=20, corral = "random", cex=0.8,labels = FALSE,xlab='', pwcol=combined$color,method='swarm')
text(x = c(1:5), y = par("usr")[3] - 0.45, xpd = NA, adj=1.2, labels=c("Gut-Microbiota","Medication","Clinical","Lifestyle","Genetics"), srt =60, cex = 0.7)
boxplot(Percent_variance_explained ~ Features, data = combined, add = T, col="#0000ff22",names = NULL)

dev.off()


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

all_sig<- list(
  Gut_Microbiome = mic_sig$X, 
  Medication = med_sig$X, 
  Clinical = como_sig$X,
  Lifestyle = demo_sig$X,
  Genetics = gene_sig$X
  )

pdf("upset_diagram_Figure4B.pdf")
upset(fromList(all_sig),order.by="freq")
dev.off()