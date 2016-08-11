selectAlphaLambda_BuiltInCV <- function(alphaVals, X_train_val, y_train_val, 
                                        weight_vec, bClassWeights, kFoldsVal, 
                                        bParallel, iFold, n_alphas,
                                        # the next are output variables
                                        coefs_allalphas_folds, ranks_allalphas_folds)
{
  
    # generate an initial lambda sequence for this alpha
    initial_lambda<-glmnet(x=X_train_val, y=y_train_val
                           , family="binomial", alpha=alpha
                           , standardize=F)$lambda  # calculating the initial lambda
    
    cv.fit=cv.glmnet(X_train_val, y_train_val, family="binomial",
                     type.measure="auc", alpha=alpha,
                     weights=weight_vec, nfolds=kFoldsVal,
                     foldid=fold_ids,
                     parallel=bParallel)
    
    cv_result_auc <- cv.fit$cvm[which(cv.fit$lambda == cv.fit$lambda.min)]

    coefficients <-
      coef(cv.fit, s=cv.fit$lambda.min)[2:nrow(coef(cv.fit, s=cv.fit$lambda.min)),]
    coefs_allalphas_folds[, iAlpha+(iFold-1)*n_alphas] <- coefficients
    colnames(coefs_allalphas_folds)[iAlpha+(iFold-1)*n_alphas] <-
      paste("alpha_", alphaVals[iAlpha], "_fold_", iFold, sep="")
    
    ranks_allalphas_folds[, iAlpha+(iFold-1)*n_alphas] <-
      getRanking(coefficients)
    colnames(ranks_allalphas_folds)[iAlpha+(iFold-1)*n_alphas] <-
      paste("alpha_", alphaVals[iAlpha], "_fold_", iFold, sep="")

  
  # find the best alpha
  best_auc <- 0
  for (iAlpha in 1:length(alphaVals))
  {
    if (cv_results_allalphas[[iAlpha]][[2]] > best_auc)
    {
      best_auc <- cv_results_allalphas[[iAlpha]][[2]]
      best_alpha <- alphaVals[iAlpha]
      best_lambda <- cv_results_allalphas[[iAlpha]][[1]]$lambda.min
    }
  }
  
  result <- list(alpha=best_alpha, lambda=best_lambda, 
                 coefs_allalphas_folds=coefs_allalphas_folds, 
                 ranks_allalphas_folds=ranks_allalphas_folds,
                 initial_lambda_fromCV_lst=initial_lambda_fromCV_lst,
                 initial_lambda_lst=initial_lambda_lst)
  return (result)
}



run_evalFold_par <- function(iFold)
{
  cat(paste(iFold, "..", sep=""))
  
  train_val_ids <- folds[[iFold]]
  trainIds4EvalFolds[[iFold]] <- train_val_ids
#   if(iFold==length(folds)){
#     saveRDS(trainIds4EvalFolds, paste0(resultDirPerOutcome, 'trainIds4EvalFolds.RDS'))
#     # break
#   }
  test_ids <- which(!((1:n_data) %in% train_val_ids))
  X_train_val <- X[train_val_ids,]
  X_test <- X[-train_val_ids,]
  y_train_val <- y[train_val_ids]
  y_test <- y[-train_val_ids]
  
  # compute weights
  
  weight_vec <- computeWeights(y_train_val)
  
  cat("a0")
  best_alpha <- selected_alpha_lambda$alpha
  best_lambda <- selected_alpha_lambda$lambda
  coefs_allalphas_folds <- selected_alpha_lambda$coefs_allalphas_folds
  ranks_allalphas_folds <- selected_alpha_lambda$ranks_allalphas_folds
  initial_lambda_fromCV_lst <- selected_alpha_lambda$initial_lambda_fromCV_lst
  initial_lambda_lst <- selected_alpha_lambda$initial_lambda_lst
  init_lambda_fromCV_allEvalFolds_lst[[iFold]] <- initial_lambda_fromCV_lst
  init_lambda_allEvalFolds_lst[[iFold]] <- initial_lambda_lst
  
  if(iFold == kFoldsEval){
    saveRDS(init_lambda_fromCV_allEvalFolds_lst
            , paste0(resultDirPerOutcome, "init_lambda_seq_fromCV_allEvalFolds.RDS")
    )
    saveRDS(init_lambda_allEvalFolds_lst
            , paste0(resultDirPerOutcome, "init_lambda_seq_allEvalFolds.RDS")
    )
    
  }
  # train with the selected params
  lambdaSeq4BestAlpha <- initial_lambda_fromCV_lst[[which(alphaVals==best_alpha)]]
  if (bClassWeights){
    fit_glmnet <- glmnet(X_train_val,y_train_val, family="binomial", 
                         weights=weight_vec,
                         alpha=best_alpha, lambda=lambdaSeq4BestAlpha)
    
  }else{
    
    fit_glmnet <- glmnet(X_train_val,y_train_val, family="binomial", 
                         alpha=best_alpha, lambda=lambdaSeq4BestAlpha)
    
  }
  cat("a")
  
  # test immediately on the training
  
  predprobs_train <- 
    predict(fit_glmnet, newx = X_train_val, type="response", s=best_lambda)
  rocValues_train <- 
    roc(response=as.vector(y_train_val), 
        predictor=as.vector(predprobs_train),
        direction="<")
  auc_train_allfolds[iFold, 1] <- rocValues_train$auc
  
  cat("b")
  
  
  # test
  
  predprobs_test <- 
    predict(fit_glmnet, newx = X_test, type="response", s=best_lambda)
  
  cat("c")
  
  # keep the prediction probs
  
  predprobs_alldata[test_ids, 1] <- 
    predprobs_test
  # keep the lables of the test data of eval fold i
  predprobs_alldata[test_ids, 2] <- 
    y_test
  # 
  cat("d")
  params_allfolds[iFold, 1] <- best_alpha
  params_allfolds[iFold, 2] <- best_lambda
  rownames(params_allfolds)[iFold] <- paste("fold_", iFold, sep="")
}


eachRun_par <- function(iRun){
  seed <- repeat_evalFold[iRun, 'seed']
  iEvalFold <- repeat_evalFold[iRun, 'evalFold']
#   saveRDS(repeats, paste0(dirThisOutcome, 'index4AllRepeatsRun.RDS'))
  set.seed(seed)
  evalFoldsIds <- manualStratify(y, n.evalFolds)
  
  samp_idx <- evalFoldsIds[[iEvalFold]]
  samp_data <- data[samp_idx,]
  
  X_samp <- as.matrix(samp_data[,-match('y', names(samp_data))])
  y_samp <- samp_data[, match('y', names(samp_data))]
  
  weight_vec <- computeWeights(y_samp)
  
  initLmd <- glmnet(X_samp, y_samp, family='binomial', alpha=avgAlpha, weights=weight_vec)$lambda
  myLmdSeq <- c(avgLmd, min(initLmd), max(initLmd))
  fit_glmnet <- glmnet(X_samp, y_samp, family="binomial", 
                       alpha=avgAlpha, lambda=myLmdSeq
                       , weights = weight_vec)
  coefThisLmd <- fit_glmnet$beta[, 1]
  temp <- list(coefThisLmd=coefThisLmd, myLmdSeq=myLmdSeq)
  return(temp)
  
}

run_glmnet_repeatTimes <- function(iExp){
  outcome <- Exp.df[iExp, "outcome"]
  targetVars <- Exp.df[iExp, "targetVars"]
  dirThisOutcome <- paste0(subEnetDir, "1/", cohort, '/', outcome, '/')
  paraAllEval <- read.table(paste0(dirThisOutcome, paraFileNm)
                            , sep=','
                            , header = T
                            , stringsAsFactors = F
  )

  data <- read.table(paste0(dirThisOutcome, modelDtFileNm)
                     , sep=','
                     , header = T
                     , stringsAsFactors = F
  )
  y <- data$y
  # stratify the model into 100 parts to the repeat running
  set.seed(seedLst[iRepeat])
  evalFoldsIds <- manualStratify(y, n.evalFolds)
  
  
  sfInit(parallel=TRUE, cpus=num_pros, type='SOCK')
  sfSource("functions/manualStratify.R")
  sfSource("functions/computeWeights.R")
  sfSource("functions/funs_getCoefCi.R")
  
  sfExport('data', 'avgAlpha', 'avgLmd', 'repeats'
  )
  sfClusterEval(library("glmnet"))
  sfClusterEval(library("ROCR"))
  sfClusterEval(library("plyr"))
  sfClusterEval(library("dplyr"))
  temp <- sfClusterApplyLB(1:nrow(repeat_evalFold), eachRun_par)
  sfStop()
  coefAllRepeatsRun <- ldply(lapply(temp, function(X)X$coefThisLmd), quickdf)
  myLmdSeqAllRepeatsRun <- ldply(lapply(temp, function(X)X$myLmdSeq), quickdf)
  names(myLmdSeqAllRepeatsRun) <- c("lmdUsedThisRun", "minLmdOfInitSeq", "maxLmdOfInitSeq")
  # avgCoef <- apply(coefAllRepeatsRun, 2, mean)
  coef4TargetVar <- coefAllRepeatsRun[, match(targetVars, names(coefAllRepeatsRun))]
  
  write.table(coefAllRepeatsRun, paste0(result_dir, cohort, '_', outcome, '_Top10VarsCoefAllRepeats.csv'))
  write.table(myLmdSeqAllRepeatsRun, paste0(result_dir, cohort, '_', outcome, '_myLmdSeqAllRepeatsRun.csv'))
  
  # avgCoef4TargetVar <- avgCoef[match(targetVars, names(avgCoef))]
  return(coef4TargetVar)
  
}

getCoefCi <- function(iExp){
  vct <- coefAllRepeats4TargetVarsAddExp[iExp, -match(c("outcome", 'targetVars'), names(coefAllRepeats4TargetVarsAddExp))]
  qtl <- quantile(vct, probs=seq(0, 1, 0.025))
  qtl_2.5 <- as.numeric(qtl[match('2.5%', names(qtl))])
  qtl_97.5 <- as.numeric(qtl[match('97.5%', names(qtl))])
  ci <- c(Exp.df[iExp,], Coef_2.5=qtl_2.5, Coef_97.5=qtl_97.5)
  return(ci)
}
