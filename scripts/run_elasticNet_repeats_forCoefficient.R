rm(list=ls())
library(glmnet)
library(plyr)
library(dplyr)
library(snowfall)
source("functions/manualStratify.R")
source("functions/computeWeights.R")

subEnetDir <- "F:/Jie/MS/03_Result/2016-08-08/2016-08-08 09.24.44/"
cohort <- 'B2B'
n.repeat <- 100
num_pros <- 39

paraFileNm <- paste0(cohort, "_params.csv")
modelDtFileNm <- paste0(cohort, "_data_for_model.csv")

Exp.df <- data.frame(outcome=c("relapse_fu_any_01"
                               , "relapse_or_prog"
                               , "relapse_or_conf")
                     , targetVars=c("pre2_edss_score__gt4"
                                    , "pre3_edss_score__gt4"
                                    , "pre3_edss_score__gt4")
                     )

timeStamp <- as.character(Sys.time())
timeStamp <- gsub(":", ".", timeStamp)  # replace ":" by "."
result_dir <- paste("./Results/", timeStamp, "/", sep = '')
dir.create(result_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")

eachRun_par <- function(iRepeat){
  samp_idx <- repeats[[iRepeat]]
  samp_data <- data[samp_idx,]
  
  X_samp <- as.matrix(samp_data[,-match('y', names(samp_data))])
  y_samp <- samp_data[, match('y', names(samp_data))]

  weight_vec <- computeWeights(y_samp)
  
  initLmd <- glmnet(X_samp, y_samp, family='binomial', alpha=avgAlpha, weights=weight_vec)$lambda
  initLmd <- c(avgLmd, initLmd[1:5])
  fit_glmnet <- glmnet(X_samp, y_samp, family="binomial", 
                       alpha=avgAlpha, lambda=initLmd
                       , weights = weight_vec)
  coefThisLmd <- fit_glmnet$beta[, 1]
  return(coefThisLmd)
  
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
  avgPara <- apply(paraAllEval[, -1], 2, median)
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
  
  sfExport('data', 'avgAlpha', 'avgLmd', 'repeats'
           )
  sfClusterEval(library("glmnet"))
  sfClusterEval(library("ROCR"))
  sfClusterEval(library("plyr"))
  sfClusterEval(library("dplyr"))
  temp <- sfClusterApplyLB(1:n.repeat, eachRun_par)
  sfStop()
  coefAllRepeatsRun <- ldply(temp, quickdf)
  # avgCoef <- apply(coefAllRepeatsRun, 2, mean)
  coef4TargetVar <- coefAllRepeatsRun[, match(targetVars, names(coefAllRepeatsRun))]
  
  write.table(coefAllRepeatsRun, paste0(result_dir, cohort, '_', outcome, '_Top10VarsCoefAllRepeats.csv'))
  # avgCoef4TargetVar <- avgCoef[match(targetVars, names(avgCoef))]
  return(coef4TargetVar)
  
}

coefAllRepeats4TargetVars <- ldply(lapply(1:nrow(Exp.df), run_glmnet_repeatTimes), quickdf)

coefAllRepeats4TargetVarsAddExp <- cbind(Exp.df, coefAllRepeats4TargetVars)

avgCoef4TargetVarsAddExp <- cbind(Exp.df, avgCoef=apply(coefAllRepeats4TargetVars, 1, mean))
write.table(coefAllRepeats4TargetVarsAddExp
            , paste0(result_dir, '/coefAllRepeats4TargetVarsAddExp.csv')
            , sep=','
            , row.names = F
            )

write.table(avgCoef4TargetVarsAddExp
            , paste0(result_dir, '/avgCoef4TargetVarsAddExp.csv')
            , sep=','
            , row.names = F
)



