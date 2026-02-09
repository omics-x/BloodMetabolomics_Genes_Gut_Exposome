###################################################################
### Code for multivariate prediction of MRI variables using elastic net regularization
###################################################################
### Load R packages, R version 4.5.1 (2025-06-13)
library(foreign)
library(dplyr)
library(glmnet)
library(caret)
library(ggplot2)
library(Metrics)

### The combined Metabolomics and MRI data in RSIII-2
MRI_data<-read.table("/Users/sahmad1/Documents/RS_projects/mri_associations/mri_new_covariat_file/RS1_5_Metabolon_MRIdata_2Jun2021.txt",head=T)
### Log transform the WML variable
MRI_data$log_Total_wml<-log(MRI_data$Total_wml)
### Remove participants with missing information on WML
metabolites_wml_data<-MRI_data %>% filter(!is.na(log_Total_wml))

set.seed(123)
### Define outcome (Y) and metabolite predictors (X)
Y <- metabolites_wml_data$log_Total_wml
X <- as.matrix(metabolites_wml_data[grep("^metab_", names(metabolites_wml_data))])
data_caret <- data.frame(Y = Y, X)

# Create training (80%) and testing (20%) split
train_index <- createDataPartition(data_caret$Y, p = 0.8, list = FALSE)
train_data <- data_caret[train_index, ]
test_data <- data_caret[-train_index, ]

# 10 fold cross-validation 
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  verboseIter = TRUE
)

# Training the model using elastic net model with automatic tuning
enet_model <- train(
  Y ~ .,
  data = train_data,
  method = "glmnet",
  trControl = ctrl,
  tuneLength = 10,  
  metric = "RMSE"
)

# Evaluate on test set
test_pred <- predict(enet_model, newdata = test_data)
#test_rsq <- 1 - sum((test_data$Y - test_pred)^2) / sum((test_data$Y - mean(test_data$Y))^2)
test_rsq <- R2(pred = test_pred, obs = test_data$Y)


# Output results
cat("Best alpha:", enet_model$bestTune$alpha, "\n")
cat("Best lambda:", enet_model$bestTune$lambda, "\n")
cat("Test R-squared:", round(test_rsq, 4), "\n")

######################################################
### Extract metabolites that contributed to prediction of MRI variable: WML
######################################################
coef_enet <- coef(enet_model$finalModel, 
                  s = enet_model$bestTune$lambda)
# Convert sparse matrix to data frame
coef_df <- data.frame(
  feature = rownames(coef_enet),
  coefficient = as.numeric(coef_enet),
  row.names = NULL
)
# Keep only non-zero predictors (exclude intercept)
nonzero_coef <- subset(coef_df, coefficient != 0 & feature != "(Intercept)")
nonzero_coef <- nonzero_coef[order(abs(nonzero_coef$coefficient), decreasing = TRUE), ]
nonzero_coef$contributing_rank<-seq(1:length(nonzero_coef$coefficient))
# Annotate the metabolites 
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)
top_features_anno_WML<-merge(nonzero_coef,anno_com[,c("Name","CHEMICAL_NAME_new")],by.x="feature",by.y="Name",all.x=T,sort = F)
# Save metabolites with non-zero elastic net coefficients
write.csv(top_features_anno_WML,file='Elastic_net_model_prediction_WML.csv')

######################################################
### Metabolomics data for total brain volume (TBV)
######################################################
metabolites_totalBrainvol_data<-MRI_data %>% filter(!is.na(Total_par_ml))
set.seed(123)
### Define outcome (Y) and metabolite predictors (X)
Y <- metabolites_totalBrainvol_data$Total_par_ml
X <- as.matrix(metabolites_totalBrainvol_data[grep("^metab_", names(metabolites_totalBrainvol_data))])
data_caret <- data.frame(Y = Y, X)

# Create training and testing split
train_index <- createDataPartition(data_caret$Y, p = 0.8, list = FALSE)
train_data <- data_caret[train_index, ]
test_data <- data_caret[-train_index, ]

# 10 fold cross-validation 
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  verboseIter = TRUE
)

# training the model using elastic net model with automatic tuning
enet_model <- train(
  Y ~ .,
  data = train_data,
  method = "glmnet",
  trControl = ctrl,
  tuneLength = 10,  
  metric = "RMSE"
)

# Evaluate on test set
test_pred <- predict(enet_model, newdata = test_data)
#test_rsq <- 1 - sum((test_data$Y - test_pred)^2) / sum((test_data$Y - mean(test_data$Y))^2)
test_rsq <- R2(pred = test_pred, obs = test_data$Y)

# Output results
cat("Best alpha:", enet_model$bestTune$alpha, "\n")
cat("Best lambda:", enet_model$bestTune$lambda, "\n")
cat("Test R-squared:", round(test_rsq, 4), "\n")

######################################################
### Extract metabolites that contributed to prediction of MRI variable: Total brain volume
######################################################
coef_enet <- coef(enet_model$finalModel, 
                  s = enet_model$bestTune$lambda)
# Convert sparse matrix to data frame
coef_df <- data.frame(
  feature = rownames(coef_enet),
  coefficient = as.numeric(coef_enet),
  row.names = NULL
)
# Keep only non-zero predictors (exclude intercept)
nonzero_coef <- subset(coef_df, coefficient != 0 & feature != "(Intercept)")
nonzero_coef <- nonzero_coef[order(abs(nonzero_coef$coefficient), decreasing = TRUE), ]
nonzero_coef$contributing_rank<-seq(1:length(nonzero_coef$coefficient))
# Annotate the metabolites 
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)
top_features_anno_TBV<-merge(nonzero_coef,anno_com[,c("Name","CHEMICAL_NAME_new")],by.x="feature",by.y="Name",all.x=T,sort = F)
# Save metabolites with non-zero elastic net coefficients
write.csv(top_features_anno_TBV,file='Elastic_net_model_prediction_TBV.csv')



######################################################
### Metabolomics data for total Hippocampal volume 
######################################################
metabolites_totalHCV_data<-MRI_data %>% filter(!is.na(total_Hippocampus))

set.seed(123)
### Define outcome (Y) and metabolite predictors (X)
Y <- metabolites_totalHCV_data$total_Hippocampus
X <- as.matrix(metabolites_totalHCV_data[grep("^metab_", names(metabolites_totalHCV_data))])
data_caret <- data.frame(Y = Y, X)

# Create training and testing split
train_index <- createDataPartition(data_caret$Y, p = 0.8, list = FALSE)
train_data <- data_caret[train_index, ]
test_data <- data_caret[-train_index, ]

# 10 fold cross-validation 
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  verboseIter = TRUE
)

# training the model using elastic net model with automatic tuning
enet_model <- train(
  Y ~ .,
  data = train_data,
  method = "glmnet",
  trControl = ctrl,
  tuneLength = 10,  
  metric = "RMSE"
)

# Evaluate on test set
test_pred <- predict(enet_model, newdata = test_data)
#test_rsq <- 1 - sum((test_data$Y - test_pred)^2) / sum((test_data$Y - mean(test_data$Y))^2)
test_rsq <- R2(pred = test_pred, obs = test_data$Y)


# Output results
cat("Best alpha:", enet_model$bestTune$alpha, "\n")
cat("Best lambda:", enet_model$bestTune$lambda, "\n")
cat("Test R-squared:", round(test_rsq, 4), "\n")

######################################################
### Extract metabolites that contributed to prediction of MRI variable: HCV
######################################################
coef_enet <- coef(enet_model$finalModel, 
                  s = enet_model$bestTune$lambda)
# Convert sparse matrix to data frame
coef_df <- data.frame(
  feature = rownames(coef_enet),
  coefficient = as.numeric(coef_enet),
  row.names = NULL
)
# Keep only non-zero predictors (exclude intercept)
nonzero_coef <- subset(coef_df, coefficient != 0 & feature != "(Intercept)")
nonzero_coef <- nonzero_coef[order(abs(nonzero_coef$coefficient), decreasing = TRUE), ]
nonzero_coef$contributing_rank<-seq(1:length(nonzero_coef$coefficient))
# Annotate the metabolites 
anno_com <- read_excel("DUKE-0304-19ML_annotation_fileRSIII_2.XLSX",sheet = 1)
anno_com<-as.data.frame(anno_com)
anno_com$Name<-paste0("metab_",anno_com$CHEM_ID)
top_features_anno_HCV<-merge(nonzero_coef,anno_com[,c("Name","CHEMICAL_NAME_new")],by.x="feature",by.y="Name",all.x=T,sort = F)
# Save metabolites with non-zero elastic net coefficients
write.csv(top_features_anno_HCV,file='Elastic_net_model_prediction_HCV.csv')


