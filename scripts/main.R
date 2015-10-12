library(glmnet)
library(pROC)
library(ROCR)

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

# fit the glmnet model

folds <- manualStratify(y, kFoldsEval)

predprobs_alldata <- matrix(data=-1, nrow=n_data, ncol=1)
lambda_allfolds <- matrix(data=-1, nrow=kFoldsEval, ncol=1)

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
  
  fit <- glmnet(X_train, y_train, family="binomial", 
                alpha=selected_alpha, lambda=lambda_seq)
  predprobs_test <- predict(fit, newx=X_test, s=selected_lambda, type="response")
  
  # 
  
  coefs_test <- predict(fit, newx=X_test, s=selected_lambda, type="coefficients")
  selected_var_names <- c(selected_var_names, rownames(coefs_test)[coefs_test[,1]!=0])
  
  # 
  predprobs_alldata[test_ids] <- predprobs_test
  # lambda_allfolds[iFold, 1] <- fit$lambda.min
}
cat("\n")

# test the AUC

pred <- prediction(predictions=predprobs_alldata, labels=y)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
png(filename=paste(result_dir, study_name, "_roc_alpha", 
                   selected_alpha, ".png", sep=""))
plot(perf, col=rainbow(10))
dev.off()

roc_obj <- roc(response=as.vector(y), predictor=as.vector(predprobs_alldata))
auc <- roc_obj$auc
cat(paste("AUC: ", auc, "\n", sep=""))

ci_obj <- ci.auc(roc_obj)

file_ciauc <- file(paste(result_dir, "ci_auc.csv", sep=""), "w")
writeLines(paste(c(ci_obj[1], ci_obj[2], ci_obj[3]),collapse=","), file_ciauc)
close(file_ciauc)


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
write.table(coef_info, sep=",", 
            file=paste(result_dir, "coef_info.csv", sep=""), col.names=NA)

# 


