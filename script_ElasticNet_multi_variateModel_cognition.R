###################################################################
### Code for multivariate prediction of general cognition using elastic net regularization
###################################################################
### Load R packages, R version 4.5.1 (2025-06-13)
library(foreign)
library(dplyr)
library(glmnet)
library(caret)
library(ggplot2)
library(Metrics)
### the phenotype file for general cognition. Variable FAC1_1 represents the G-factor (general cognition score)
cognitionRS3<-read.spss("/Users/sahmad1/Documents/RS_projects/association_cognition/e5_RS-I-5 RS-II-3 RS-III-2_cognition_complete 030217_analysis_file.sav",to.data.frame=T)
cognitionRS3<-cognitionRS3[,c("ergoid","FAC1_1")]
### Load a combined Rotterdam Study covariate information file 
covariates<-readRDS(file ="/Users/sahmad1/Documents/RS_projects/Study_RSI_IV_RSIII_2_covars.rds")
### Load Metabolomics data file for RSIII-2
load("/Users/sahmad1/Documents/RS_projects/KNN_imputed_metabolon_data_3rdMarch.RData")
imputeddata<-as.data.frame(imputeddata)
colnames(imputeddata)<-paste("metab",colnames(imputeddata),sep="_")
imputeddata$ergoid<-rownames(imputeddata)
### Merge cognition, covariates, and metabolomics datasets by participant ID (ergoid)
merg1<-merge(cognitionRS3,covariates,by.x="ergoid",by.y="ergoid",all.x=T,sort=F)
merg2<-merge(imputeddata,merg1,by.x="ergoid",by.y="ergoid",all.x=T,sort=F)
### Remove participants with missing outcome (FAC1_1)
merg2<-merg2[which(!is.na(merg2$FAC1_1)),]
set.seed(123)
### Define outcome (Y) and metabolite predictors (X)
Y <- merg2$FAC1_1
X <- as.matrix(merg2[grep("^metab_", names(merg2))])
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
# Evaluate the model on test set
test_pred <- predict(enet_model, newdata = test_data)
#test_rsq <- 1 - sum((test_data$Y - test_pred)^2) / sum((test_data$Y - mean(test_data$Y))^2)
test_rsq <- R2(pred = test_pred, obs = test_data$Y)

# Output results
cat("Best alpha:", enet_model$bestTune$alpha, "\n")
cat("Best lambda:", enet_model$bestTune$lambda, "\n")
cat("Test R-squared:", round(test_rsq, 4), "\n")

######################################################
### Extract metabolites that contributed to prediction of general cognition
######################################################

# Extract coefficients from the final model
best_alpha <- enet_model$bestTune$alpha
best_lambda <- enet_model$bestTune$lambda
# Extract coefficients using coef() from glmnet
coef_mat <- coef(enet_model$finalModel, s = best_lambda)
#Convert to data frame
coef_df <- as.data.frame(as.matrix(coef_mat))
coef_df$feature <- rownames(coef_df)
colnames(coef_df)[1] <- "coefficient"
#Filter non-zero coefficients and exclude intercept
nonzero_coefs <- subset(coef_df, coefficient != 0 & feature != "(Intercept)")
# Sort by absolute value (importance)
top_features <- nonzero_coefs[order(abs(nonzero_coefs$coefficient), decreasing = TRUE), ]
# Save metabolites with non-zero elastic net coefficients
write.csv(top_features,file='Elastic_net_model_prediction_G_factor.csv')

######################################################
### Plot predicted versus actual outcomes with regression line and 95% CI
######################################################

plot_data <- data.frame(Predicted = test_pred, Actual = test_data$Y)
r2_val <- round(R2(pred = plot_data$Predicted, obs = plot_data$Actual), 3)
rmse_val <- round(rmse(actual = plot_data$Actual, predicted = plot_data$Predicted), 3)
ggplot(plot_data, aes(x = Actual, y = Predicted)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "lightblue", alpha = 0.4) +
  annotate("text", 
           x = min(plot_data$Actual), 
           y = max(plot_data$Predicted), 
           label = paste0("R² = ", r2_val, "\nRMSE = ", rmse_val),
           hjust = 0, vjust = 1, size = 4.5, color = "black") +
  labs(
    title = "Predicted vs Actual with 95% CI",
    x = "Actual Y (General Cognition)",
    y = "Predicted Y (General Cognition)"
  ) +
  theme_minimal()
