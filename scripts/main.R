library(glmnet)
library(pROC)
library(ROCR)

rm(list=ls())

source("functions/manualStratify.R")
source("functions/computeWeights.R")



kFoldsEval = 10
log_lambda_seq <- seq(log(1e-4),log(1e4),length.out=100)
lambda_seq <- exp(log_lambda_seq)

# data

data_dir <- "C:/Work/Projects/MultipleSclerosis/Results/2015-10-06/2015-10-06 17.57.03/"
study_name <- c("combined_relapse_fu_any_01")

# read and transform the data

dataset <- read.csv(paste(data_dir, study_name,"_data_for_model.csv", sep=""), 
                    header=TRUE, sep=",", check.names=FALSE)
# first column is index
dataset[,1] <- NULL
y <- dataset[, 1]
X <- dataset[, 2:ncol(dataset)]
X <- data.matrix(X)
n_data <- length(y)


# parameters for the best glmnet model

selected_alpha <- 0.85
selected_lambda <- 0.01462
# selected_lambda <- 0.01

# fit the glmnet and glm models

folds <- manualStratify(y, kFoldsEval)

predprobs_alldata_glmnet <- matrix(data=-1, nrow=n_data, ncol=1)
predprobs_alldata_glm <- matrix(data=-1, nrow=n_data, ncol=1)

timeStamp <- as.character(Sys.time())
timeStamp <- gsub(":", ".", timeStamp)  # replace ":" by "."
result_dir <- paste("./Results/", timeStamp, "/", sep = '')
dir.create(result_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")

cat("Fold ")
selected_var_names <- NULL
for (iFold in 1:length(folds))
{
  cat(paste(iFold, "..", sep=""))
  
  train_ids <- folds[[iFold]]
  test_ids <- which(!((1:n_data) %in% train_ids))
  X_train <- X[train_ids,]
  X_test <- X[-train_ids,]
  y_train <- y[train_ids]
  y_test <- y[-train_ids]
  
  weight_vec <- computeWeights(y_train)
  
  fit_glmnet <- glmnet(X_train, y_train, family="binomial", 
                       alpha=selected_alpha, lambda=lambda_seq)
  predprobs_test_glmnet <- 
    predict(fit_glmnet, newx=X_test, s=selected_lambda, type="response")
  
  data_train <- data.frame(cbind(y_train,X_train))
  colnames(data_train)[1] <- "y"
  fit_glm <- glm(y ~ ., family="binomial", data=data_train)
  data_test <- data.frame(cbind(y_test, X_test))
  colnames(data_test)[1] <- "y"
  predprobs_test_glm <- 
    predict(fit_glm, newdata=data_test, type="response")
  
  # 
  
  coefs_test <- predict(fit_glmnet, newx=X_test, s=selected_lambda, type="coefficients")
  selected_var_names <- c(selected_var_names, rownames(coefs_test)[coefs_test[,1]!=0])
  
  # 
  predprobs_alldata_glmnet[test_ids] <- predprobs_test_glmnet
  predprobs_alldata_glm[test_ids] <- predprobs_test_glm
  # lambda_allfolds[iFold, 1] <- fit$lambda.min
}
cat("\n")

# test the AUC

pred_glmnet <- prediction(predictions=predprobs_alldata_glmnet, labels=y)
perf_glmnet <- performance(pred_glmnet, measure = "tpr", x.measure = "fpr") 
png(filename=paste(result_dir, study_name, "_roc_alpha", 
                   selected_alpha, "glmnet.png", sep=""))
plot(perf_glmnet, col=rainbow(10))
dev.off()

pred_glm <- prediction(predictions=predprobs_alldata_glm, labels=y)
perf_glm <- performance(pred_glm, measure = "tpr", x.measure = "fpr") 
png(filename=paste(result_dir, study_name, "_roc_alpha", 
                   selected_alpha, "glm.png", sep=""))
plot(perf_glm, col=rainbow(10))
dev.off()

roc_obj_glmnet <- roc(response=as.vector(y), predictor=as.vector(predprobs_alldata_glmnet))
auc_glmnet <- roc_obj_glmnet$auc
cat(paste("AUC glmnet: ", auc_glmnet, "\n", sep=""))

ci_obj_glmnet <- ci.auc(roc_obj_glmnet)

file_ciauc_glmnet <- file(paste(result_dir, "ci_auc_glmnet.csv", sep=""), "w")
writeLines(paste(c(ci_obj_glmnet[1], ci_obj_glmnet[2], ci_obj_glmnet[3]),collapse=","), file_ciauc_glmnet)
close(file_ciauc_glmnet)


roc_obj_glm <- roc(response=as.vector(y), predictor=as.vector(predprobs_alldata_glm))
auc_glm <- roc_obj_glm$auc
cat(paste("AUC glm: ", auc_glm, "\n", sep=""))

ci_obj_glm <- ci.auc(roc_obj_glm)

file_ciauc_glm <- file(paste(result_dir, "ci_auc_glm.csv", sep=""), "w")
writeLines(paste(c(ci_obj_glm[1], ci_obj_glm[2], ci_obj_glm[3]),collapse=","), file_ciauc_glm)
close(file_ciauc_glm)



hh <- roc.test(roc_obj_glm, roc_obj_glmnet)
file_roc_compare <- file(paste(result_dir, "roc_compare.txt", sep=""), "w")
for (i in 1:length(hh))
{
  if ((names(hh)[i] == "roc1") | (names(hh)[i] == "roc2"))
    next
  
  writeLines(paste(names(hh)[i], ": ",sep=""), file_roc_compare)
  if (is.null(names(hh[[i]])))
  {
    writeLines(paste(hh[[i]], "\n", sep=""), file_roc_compare)
    # writeLines(paste(hh[[i]],"\n", sep=""), file_roc_compare)
  }
    
  else
  {
    if (names(hh)[i] == "p.value")
      writeLines(paste(hh[[i]][1], "\n", sep=""), file_roc_compare)
    else
    {
      for (j in 1:length(hh[[i]]))
      {
        writeLines(paste(names(hh[[i]][j]), ": ", hh[[i]][j], sep=""), file_roc_compare)
      }
      writeLines("\n", file_roc_compare)
    }
  }
}

writeLines("Are the two ROCs paired?\n", file_roc_compare)
writeLines(paste(are.paired(roc_obj_glm, roc_obj_glmnet), "\n"), file_roc_compare)
close(file_roc_compare)


# use selected features for non-regularised LR

selected_var_names <- unique(selected_var_names)
selected_var_names <- selected_var_names[which(selected_var_names!="(Intercept)")]
write.table(selected_var_names, sep=",", 
            file=paste(result_dir, "selected_vars.csv", sep=""), col.names=NA)

# get the confidence intervals and p-values for the coefficients

fit_glm <- glm(y ~ ., family="binomial", data=dataset[, colnames(dataset) %in% c("y", selected_var_names)])
# p_values <- coef(summary(fit_glm))[,4]
ci_coefs <- confint(fit_glm)

coef_info <- cbind(coef(summary(fit_glm)), ci_coefs)
odds_ratios <- exp(coef_info[,1])
odds_ratios_low <- exp(ci_coefs[,1])
odds_ratios_high <- exp(ci_coefs[,2])
coef_info <- cbind(coef_info, odds_ratios, odds_ratios_low, odds_ratios_high)
colnames(coef_info)[ncol(coef_info)-2] <- "odds"
colnames(coef_info)[ncol(coef_info)-1] <- "odds_2.5%"
colnames(coef_info)[ncol(coef_info)] <- "odds_97.5%"
write.table(round(coef_info, digits=3), sep=",", 
            file=paste(result_dir, "coef_info.csv", sep=""), col.names=NA)

# 

