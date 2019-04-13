#! /usr/bin/Rscript

# Set up working libraries
library(withr)
library(ggplot2)
library(randomForest)
library(dplyr)

# Parse input arguments from command line
args = commandArgs(trailingOnly=TRUE)
training_file <- args[1]
test_file <- args[2]
output_fp <- args[3]
varimp_fp <- args[4]
seed_num <- 123

# Read in training set of annotated variants
training_df <- read.table(training_file, header=TRUE, sep="\t")

# Remove columns with the same-value
training_df <- training_df %>% select_if(~ length(unique(.)) > 1)
training_df <- na.roughfix(training_df)

# Train random forest model
set.seed(seed_num)
model <- randomForest(HepG2_delta_logit_psi ~ . -ID, data=training_df, ntree=2000)

# Record mean decrease in node impurity for each feature used by the model
varImp <- model$importance
write.table(varImp, varimp_fp, append=TRUE, sep="\t", col.names=FALSE)

# Read in test set of annotated variants
test_df <- read.table(test_file, header=TRUE, sep="\t")
test_df <- na.roughfix(test_df)

# Make predictions on variants in test-set using trained model
test_pred <- predict(model, newdata=test_df)
prediction_summary <- data.frame(test_df$ID, test_pred)
write.table(prediction_summary, output_fp, sep="\t", row.names=FALSE, col.names=FALSE)
