library(glmnet)
library(pROC)
library(ROCR)
library(dplyr)
library(MASS)

rm(list=ls())

source("functions/manualStratify.R")
source("functions/computeWeights.R")


bTestLambda0 <- T
kFoldsEval = 5
# log_lambda_seq <- seq(log(1e-20),log(1e4),length.out=100)
# lambda_seq <- exp(log_lambda_seq)

# data




cohortNames <- c("B2B")
# cohortNames <- c("Cmp", "BConti", "B2B", "B2Fir", "B2Sec")
outcomeNames <- c("relapse_fu_any_01", "edssprog", "edssconf3",
                  "relapse_or_prog", "relapse_and_prog", "relapse_or_conf")

bTopVarsOnly <- T

myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}


#
# !!
# The input modelling data needs be using the same number of top variables.
# Otherwise the 'best' hyper-parameters found cannot be used below. 
# 
#
if (bTopVarsOnly) 
{
  nTopVars <- 10 
  rootDataDir <- "F:/Jie/MS/03_Result/2016-08-11/2016-08-11 12.09.42/"
} else 
{
  nTopVars <- NULL
  rootDataDir <- "F:/Jie/MS/03_Result/2016-08-11/2016-08-11 12.00.36/"
}



timeStamp <- as.character(Sys.time())
timeStamp <- gsub(":", ".", timeStamp)  # replace ":" by "."
result_dir <- paste("./Results/", timeStamp, "/", sep = '')
dir.create(result_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")

set.seed(1)


for (cohortName in cohortNames)
{
  cat(paste0(cohortName, ": "))
  resultDirPerCohort <- paste0(result_dir, cohortName, "/")
  dir.create(resultDirPerCohort, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  for (outcomeName in outcomeNames)
  {
    cat(paste0(outcomeName, ","))
    resultDirPerOutcome <- paste0(resultDirPerCohort, outcomeName, "/")
    dir.create(resultDirPerOutcome, showWarnings = TRUE, recursive = TRUE, mode = "0777")
    
    
    dataDir <- paste0(rootDataDir, "1/", cohortName, "/", outcomeName, "/")
    
    # use top m features only
    
    if (bTopVarsOnly)
    {
      topVarsDir <- paste0(rootDataDir, "1/", cohortName, '/', outcomeName , '/')
      
      avRank <- read.table(
        paste0(topVarsDir, "av_ranking_", cohortName, ".csv"), sep=',', header = T, stringsAsFactors = F
      )
    
      topVarNames <- rownames(avRank)[order(avRank$x, decreasing=F)][1:nTopVars]
      
      if (any(grepl("Intercept", topVarNames)))
        stop("Error! 'Intercept' among the top variables!")
    }
   
    
    # read and transform the data
    
    dataset <- 
      tbl_df(
        read.csv(paste0(dataDir, cohortName, "_data_for_model.csv"), 
                 header=TRUE, sep=",", check.names=FALSE)
        ) %>%
      {
        # check
        print(paste0("The input modelling uses data of ", ncol(.)-1, 
                     " predictors. Is that really what you want?"))
        if (bTopVarsOnly)
        {
          if ((ncol(.)-1) != nTopVars)
            stop(paste0("Error! The input modelling result must use the same", 
                        " dimensional data as nTopVars. Otherwise the selected",
                        " beset hyper-parameters cannot be used."))
        }
        .
        
      } %>%
      {
        if (bTopVarsOnly)
          dplyr::select(., one_of(c("y", topVarNames)))
        else
          .
      } %>%
      as.data.frame()
    
    
    selectedVarNames <- setdiff(colnames(dataset), 'y')
    write.table(selectedVarNames, sep=",", 
                file=paste(resultDirPerOutcome, "selectedVarNames.csv", sep=""), col.names=NA)
    # 
    y <- dataset[, 1]
    X <- dataset[, selectedVarNames]
    X <- data.matrix(X)
    n_data <- length(y)
    
    
    # parameters for the best glmnet model
    
    paramInfo <- read.csv(paste0(dataDir, cohortName, "_params.csv"))
    
    # fit the glmnet and glm models
    
    # folds <- manualStratify(y, kFoldsEval)
    
    predprobs_alldata_glmnet <- matrix(data=-1, nrow=n_data, ncol=1)
    predprobs_alldata_glm <- matrix(data=-1, nrow=n_data, ncol=1)
    
    # cat("Fold ")
    train_ids_lst <- readRDS(paste0(dataDir, 'trainIds4EvalFolds.RDS'))
    for (iFold in 1:kFoldsEval)
    {
      cat(paste(iFold, "..", sep=""))
      
      # train_ids <- folds[[iFold]]
      train_ids <- train_ids_lst[[iFold]]
      test_ids <- which(!((1:n_data) %in% train_ids))
      X_train <- X[train_ids,]
      X_test <- X[-train_ids,]
      y_train <- y[train_ids]
      y_test <- y[-train_ids]
      
      weight_vec <- computeWeights(y_train)
      
      selected_alpha <- paramInfo$alpha[iFold]
      if(bTestLambda0==T){
        selected_lambda<- 0
      }
      
      selected_lambda <- paramInfo$lambda[iFold]
      
      fit_glmnet <- glmnet(X_train, y_train, family="binomial", 
                           alpha=selected_alpha, lambda=selected_lambda
                           , weights = weight_vec)
      if(bTestLambda0==T){
        selected_lambda<- 0
      }
      
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
      
#       coefs_test <- predict(fit_glmnet, s=selected_lambda, type="coefficients")
#       selected_var_names <- c(selected_var_names, rownames(coefs_test)[coefs_test[,1]!=0])
      
      # 
      predprobs_alldata_glmnet[test_ids,1] <- predprobs_test_glmnet
      predprobs_alldata_glm[test_ids,1] <- predprobs_test_glm
      # lambda_allfolds[iFold, 1] <- fit$lambda.min
    }
    cat("\n")
    
       
    # test the AUC
    pred_glmnet <- prediction(predictions=predprobs_alldata_glmnet[,1], labels=y)
    perf_glmnet <- performance(pred_glmnet, measure = "tpr", x.measure = "fpr") 
    png(filename=paste(resultDirPerOutcome, cohortName, "_", outcomeName, "_roc_alpha", 
                       selected_alpha, "glmnet.png", sep=""))
    plot(perf_glmnet, col=rainbow(10))
    dev.off()
    
    pred_glm <- prediction(predictions=predprobs_alldata_glm[,1], labels=y)
    perf_glm <- performance(pred_glm, measure = "tpr", x.measure = "fpr") 
    png(filename=paste(resultDirPerOutcome, cohortName, "_", outcomeName, "_roc_alpha", 
                       selected_alpha, "glm.png", sep=""))
    plot(perf_glm, col=rainbow(10))
    dev.off()
    
    roc_obj_glmnet <- roc(response=as.vector(y), predictor=as.vector(predprobs_alldata_glmnet[,1]))
    auc_glmnet <- roc_obj_glmnet$auc
    cat(paste("AUC glmnet: ", auc_glmnet, "\n", sep=""))
    
    ci_obj_glmnet <- ci.auc(roc_obj_glmnet)
    
    file_ciauc_glmnet <- file(paste(resultDirPerOutcome, "ci_auc_glmnet.csv", sep=""), "w")
    writeLines(paste(c(ci_obj_glmnet[1], ci_obj_glmnet[2], ci_obj_glmnet[3]),collapse=","), file_ciauc_glmnet)
    close(file_ciauc_glmnet)
    
    
    roc_obj_glm <- roc(response=as.vector(y), predictor=as.vector(predprobs_alldata_glm[,1]))
    auc_glm <- roc_obj_glm$auc
    cat(paste("AUC glm: ", auc_glm, "\n", sep=""))
    
    ci_obj_glm <- ci.auc(roc_obj_glm)
    
    file_ciauc_glm <- file(paste(resultDirPerOutcome, "ci_auc_glm.csv", sep=""), "w")
    writeLines(paste(c(ci_obj_glm[1], ci_obj_glm[2], ci_obj_glm[3]),collapse=","), file_ciauc_glm)
    close(file_ciauc_glm)
    
    
    
    hh <- roc.test(roc_obj_glm, roc_obj_glmnet)
    file_roc_compare <- file(paste(resultDirPerOutcome, "roc_compare.txt", sep=""), "w")
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
    
    
    
    
    # get the confidence intervals and p-values for the coefficients
    
    fit_glm_lst <- myTryCatch(glm(y ~ ., family="binomial", data=dataset[, colnames(dataset) %in% c("y", selectedVarNames)]))
    cat(fit_glm_lst$warning$message, '\n')
    fit_glm <- fit_glm_lst$value
    p_values <- coef(summary(fit_glm))[,4]
    ci_coefs_lst <- myTryCatch(confint(fit_glm))
    cat(ci_coefs_lst$error$message, '\n\n\n')
    ci_coefs <- ci_coefs_lst$value
    summary_fit_glm <- coef(summary(fit_glm))
    if(!is.null(ci_coefs_lst$error$message)) next
    
#     vars2Remove <- is.na(ci_coefs[, "2.5 %"])
#     ci_coefs <- ci_coefs[!vars2Remove, ]
    coef_info <- dplyr::left_join(
        tbl_df(as.data.frame(ci_coefs)) %>% mutate(rownameCol = rownames(ci_coefs)),
        tbl_df(as.data.frame(summary_fit_glm)) %>% mutate(rownameCol = rownames(summary_fit_glm)),
        by = "rownameCol"
      ) %>%
      {
        rownames(.) <- .$rownameCol
        .
      } %>%
      dplyr::select(-rownameCol)
    coef_info <- coef_info[, c("Estimate","Std. Error","z value","Pr(>|z|)","2.5 %","97.5 %")]
    # coef_info <- cbind(summary_fit_glm, ci_coefs)
    odds_ratios <- exp(coef_info[,1])
    odds_ratios_low <- exp(ci_coefs[,1])
    odds_ratios_high <- exp(ci_coefs[,2])
    coef_info <- cbind(coef_info, odds_ratios, odds_ratios_low, odds_ratios_high)
    colnames(coef_info)[ncol(coef_info)-2] <- "odds"
    colnames(coef_info)[ncol(coef_info)-1] <- "odds_2.5%"
    colnames(coef_info)[ncol(coef_info)] <- "odds_97.5%"
    write.table(round(coef_info, digits=3), sep=",", 
                file=paste(resultDirPerOutcome, "coef_info.csv", sep=""), col.names=NA)
  }
  cat("\n")
}



# 


