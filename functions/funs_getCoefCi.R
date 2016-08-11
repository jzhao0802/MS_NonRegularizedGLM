eachRun_par <- function(iRepeat){
  samp_idx <- repeats[[iRepeat]]
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
  if(usedMethod=='mean'){
    avgPara <- apply(paraAllEval[, -1], 2, mean)
    
  }else if(usedMethod=="median"){
    avgPara <- apply(paraAllEval[, -1], 2, median)
    
  }
  avgAlpha <- as.numeric(avgPara[match("alpha", names(avgPara))])
  avgLmd <- as.numeric(avgPara[match("lambda", names(avgPara))])
  
  data <- read.table(paste0(dirThisOutcome, modelDtFileNm)
                     , sep=','
                     , header = T
                     , stringsAsFactors = F
  )
  y <- data$y
  # stratify the model into 100 parts to the repeat running
  set.seed(1)
  repeats <- manualStratify(y, n.repeat)
  saveRDS(repeats, paste0(dirThisOutcome, 'index4AllRepeatsRun.RDS'))
  
  
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
  temp <- sfClusterApplyLB(1:n.repeat, eachRun_par)
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
