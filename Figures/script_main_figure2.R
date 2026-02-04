#################################################### 
## Code for Main Figure 2: Correlation plot between regression coefficients of association of generation cognition (RSIII-2) and Alzeimer's disease (RSI-IV)
#################################################### 
## Load R library 
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
library("readxl")
library(ggplot2)
library(jcolors)
library(gridExtra)
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
cog_ad<-merge(cog_set[,c("Beta","CHEMICAL_NAME","SUPER_PATHWAY","FDR")],ad_set[,c("Beta","CHEMICAL_NAME","FDR")],by='CHEMICAL_NAME')
cog_ad$SUPER_PATHWAY[is.na(cog_ad$SUPER_PATHWAY)]<-"Unknowns"
cog_ad$name<-NA
cog_ad$name<-ifelse(cog_ad$FDR.x<0.05 | cog_ad$FDR.y<0.05,cog_ad$CHEMICAL_NAME,"")

pdf(file="Correlation_plot_cognition_AD_regression_coefficients.pdf",width=10,height=8)
# Calculate correlation
correlation <- cor(cog_ad$Beta.x, cog_ad$Beta.y)
# Calculate p-value using cor.test()
cor_test <- cor.test(cog_ad$Beta.x, cog_ad$Beta.y)
p_value <- cor_test$p.value
# Convert correlation and p-value to strings
cor_str <- sprintf("r = %.2f", correlation)
p_value_str <- sprintf("p = %.2e", p_value)  # Display p-value in scientific notation
# Plot with annotations
ggplot(cog_ad, aes(x = Beta.x, y = Beta.y, color = factor(SUPER_PATHWAY))) +
  geom_point() +
  geom_text(label = cog_ad$name, show.legend = FALSE, size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_light() +
  labs(title = "(A)", x = "Regression Coefficient Cognition", y = "Regression Coefficient Alzheimer's disease") +
  scale_color_jcolors(palette = "pal8") +
  annotate("text", x = min(cog_ad$Beta.x), y = max(cog_ad$Beta.y),
           label = paste(cor_str, p_value_str, sep = ", "), hjust = 0, vjust = 1)
dev.off()