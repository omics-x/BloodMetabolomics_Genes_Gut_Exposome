####################################################################################
### Supplementary Figure 6: Correlation plot between explained variance (EV) by different tested features
####################################################################################
.libPaths("/home/sahmad/R/x86_64-pc-linux-gnu-library/4.1")
### install.packages("GGally")
library(GGally)
library(ggplot2)

### Load EV results from various feature sets; the demo file contains lifestyle factors, and the como file contains clinical factors
mic<-read.csv("EV_microbiota_results.csv")
med<-read.csv("EV_medication_results.csv")
demo<-read.csv("EV_lifestyle_results.csv")
como<-read.csv("EV_comborbidity_results.csv")
gene<-read.csv("EV_gene_results.csv")

## FDR correction for spearman p-value
mic$fdr<-p.adjust(mic$spearman_p,method='fdr')
med$fdr<-p.adjust(med$spearman_p,method='fdr')
demo$fdr<-p.adjust(demo$spearman_p,method='fdr')
como$fdr<-p.adjust(como$spearman_p,method='fdr')
gene$fdr<-p.adjust(gene$spearman_p,method='fdr')

## If FDR < 0.05 and EV > 0, convert EV to a percentage; otherwise, set to zero
mic$EV_Microbiota<-ifelse((mic$fdr<0.05 & mic$explained_variance_score>0),mic$explained_variance_score*100,0)
med$EV_Medication<-ifelse((med$fdr<0.05 & med$explained_variance_score>0),med$explained_variance_score*100,0)
demo$EV_Lifestyle<-ifelse((demo$fdr<0.05 & demo$explained_variance_score>0),demo$explained_variance_score*100,0)
como$EV_Clinical<-ifelse((como$fdr<0.05 & como$explained_variance_score>0),como$explained_variance_score*100,0)
gene$EV_Genetics<-ifelse((gene$fdr<0.05 & gene$explained_variance_score>0),gene$explained_variance_score*100,0)


## Merge files 
merg1<-merge(med[,c(1,10)],demo[,c(1,10)],by.x="X",by.y="X")
merg2<-merge(merg1,gene[,c(1,10)],by.x="X",by.y="X")
merg3<-merge(merg2,mic[,c(1,10)],by.x="X",by.y="X")
merg4<-merge(merg3,como[,c(1,10)],by.x="X",by.y="X")

colnames(merg4) <- c(
  "Feature",
  "EV Medication",
  "EV Lifestyle",
  "EV Genetics",
  "EV Microbiota",
  "EV Clinical"
)


## Plot the save the PDF
# pdf("Supplementary_Figure6_CorrelationEV_Features.pdf")
# ggpairs(merg4, columns = 2:6)
# dev.off()
library(GGally)
library(ggplot2)
library(grid)

# Custom function for Pearson correlation + p-value + n
panel_cor_np <- function(data, mapping, digits = 3, ...) {
  
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  # Remove missing values pairwise
  complete <- complete.cases(x, y)
  x <- x[complete]
  y <- y[complete]
  
  n <- length(x)
  
  test <- cor.test(x, y, method = "pearson")
  
  r <- round(test$estimate, digits)
  p <- signif(test$p.value, digits)
  
  label <- paste0(
    "r = ", r,
    "\n",
    "p = ", p,
    "\n",
    "n = ", n
  )
  
  ggplot(data = data, mapping = mapping) +
    annotate(
      "text",
      x = mean(range(x, na.rm = TRUE)),
      y = mean(range(y, na.rm = TRUE)),
      label = label,
      size = 5,
      color = "#1B1B1B"
    ) +
    theme_void()
}

# Save figure
pdf("Supplementary_Figure6_CorrelationEV_Features.pdf", width = 8, height = 8)

ggpairs(
  merg4,
  columns = 2:6,
  
  lower = list(
    continuous = wrap(
      "points",
      alpha = 0.55,
      size = 1.6,
      color = "#2C7FB8"
    )
  ),
  
  upper = list(
    continuous = panel_cor_np
  ),
  
  diag = list(
    continuous = wrap(
      "densityDiag",
      fill = "#1B9E77",
      alpha = 0.65,
      color = "white",
      size = 0.3
    )
  )
) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    panel.spacing = unit(0.8, "lines")
  )

dev.off()