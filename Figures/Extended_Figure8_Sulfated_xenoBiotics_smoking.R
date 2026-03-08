#################################################
# Code for Extended Figure 8. The data contains personal information from participants therefore not provided as Extended source data for this figure 
#################################################

library(foreign)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(jcolors)
library(stringr)
library(cowplot) 


#################################################
# Load cognition data: FAC1_1 is G_factor variable 
#################################################
cognitionRS3 <- read.spss("/Users/sahmad1/Documents/RS_projects/association_cognition/e5_RS-I-5 RS-II-3 RS-III-2_cognition_complete 030217_analysis_file.sav",to.data.frame = TRUE)
cognitionRS3 <- cognitionRS3 %>%
  select(ergoid,
         FAC1_1)
#################################################
# Load covariate files 
#################################################

covariates <- readRDS(
  "/Users/sahmad1/Documents/RS_projects/Study_RSI_IV_RSIII_2_covars.rds"
)

#################################################
# Load metabolomics data
#################################################

load("/Users/sahmad1/Documents/RS_projects/KNN_imputed_metabolon_data_3rdMarch.RData")

imputeddata <- as.data.frame(imputeddata)
colnames(imputeddata) <- paste0("metab_", colnames(imputeddata))
imputeddata$ergoid <- rownames(imputeddata)

#################################################
# Metabolites with readable names
#################################################

metabolites <- c(
  metab_100006126 = "4-vinylguaiacol sulfate",
  metab_100001806 = "o-cresol sulfate",
  metab_100020975 = "3-hydroxy-2-methylpyridine sulfate",
  metab_100020515 = "2-naphthol sulfate",
  metab_100021208 = "4-vinylcatechol sulfate",
  metab_100004110 = "3-methylcatechol sulfate",
  metab_100004326 = "3-acetylphenol sulfate"
)

#################################################
# Merge datasets
#################################################

data <- imputeddata %>%
  select(ergoid, names(metabolites)) %>%
  merge(cognitionRS3, by = "ergoid") %>%
  merge(covariates, by = "ergoid") %>%
  filter(!is.na(FAC1_1), !is.na(Smoke))

data$G_factor <- data$FAC1_1

#################################################
# Scale metabolites
#################################################

data[names(metabolites)] <- scale(data[names(metabolites)],scale = T, center = T)

#################################################
# Plotting function
#################################################

new_theme <- theme_cowplot(font_size = 8) + 
  theme(
    axis.line         = element_line(linewidth = 0.3),
    axis.ticks        = element_line(linewidth = 0.3),
    legend.title      = element_text(size = 7, face = "bold"),
    legend.text       = element_text(size = 6),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 8),
    plot.title        = element_text(hjust = 0.5, size = 9, face = "bold")
  )

#################################################
# Updated Plotting Function
#################################################

make_plots <- function(df, metab_id, metab_name){
  
  # Scatter plot
  p1 <- ggplot(df, aes_string(x = metab_id, y = "G_factor", color = "factor(Smoke)")) +
    geom_point(alpha = 0.4, size = 0.8, stroke = 0) + 
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, fill = "lightgrey", alpha = 0.2) +
    scale_color_brewer(palette = "Set1", name = "Smoking Status") +
    labs(x = metab_name, y = "G factor") + 
    new_theme +
    theme(legend.position = "none") 
  
  # Boxplot
  p2 <- ggplot(df, aes_string(x = "factor(Smoke)", y = metab_id, fill = "factor(Smoke)")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, linewidth = 0.3) +
    geom_jitter(width = 0.15, color = "black", size = 0.2, alpha = 0.3) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Smoking Status", y = metab_name) +
    new_theme +
    theme(legend.position = "none")
  
  # Density plot
  p3 <- ggplot(df, aes_string(x = metab_id, color = "factor(Smoke)")) +
    geom_density(linewidth = 0.7) + 
    scale_color_brewer(palette = "Set1") +
    labs(x = metab_name, y = "Density") +
    new_theme +
    theme(legend.position = "bottom")
  
  list(p1, p2, p3)
}

#################################################
# Generate plot and save pdf
#################################################

plot_list <- list()
for(i in seq_along(metabolites)){
  metab_id <- names(metabolites)[i]
  metab_name <- metabolites[i]
  plot_list <- c(plot_list, make_plots(data, metab_id, metab_name))
}

pdf("Metabolomics_smoking_stratified_extFig8.pdf", width = 8.27, height = 11.69)

grid.arrange(
  grobs = plot_list, 
  ncol = 3, 
  padding = unit(1, "line"),
  top = "Association between Metabolites, Smoking, and Cognition"
)

dev.off()