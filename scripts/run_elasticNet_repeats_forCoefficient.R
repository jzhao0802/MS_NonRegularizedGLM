rm(list=ls())
library(glmnet)
library(plyr)
library(dplyr)
library(snowfall)
source("functions/manualStratify.R")
source("functions/computeWeights.R")
source("functions/funs_getCoefCi.R")

subEnetDir <- "F:/Jie/MS/03_Result/2016-08-08/2016-08-08 09.24.44/"
cohort <- 'B2B'
n.repeat <- 5
num_pros <- 35
alpha <- 0.1
n.evalFolds <- 5
set.seed(1)
seedList <- sample(1e10, 100)
repeat_evalFold <- expand.grid(seed=seedList, evalFold=1:n.evalFolds)

paraFileNm <- paste0(cohort, "_params.csv")
modelDtFileNm <- paste0(cohort, "_data_for_model.csv")

Exp.df <- data.frame(outcome=c("relapse_fu_any_01"
                               , "relapse_or_prog"
                               , "relapse_or_conf")
                     , targetVars=c("pre2_edss_score__gt4"
                                    , "pre3_edss_score__gt4"
                                    , "pre3_edss_score__gt4")
                     )
# usedMethod <- "median"

timeStamp <- as.character(Sys.time())
timeStamp <- gsub(":", ".", timeStamp)  # replace ":" by "."
result_dir <- paste("./Results/", timeStamp, "/", sep = '')
dir.create(result_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")

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

# calculate the confidence Interval of the 100 coefficients
coefAllRepeats4TargetVarsAddExp <- read.table(paste0(result_dir, '/coefAllRepeats4TargetVarsAddExp.csv')
                              , sep=','
                              , header = T
                              , stringsAsFactors = F
                              )

coefCi <- ldply(lapply(1:nrow(Exp.df), getCoefCi), quickdf)

coefAndOrCi <- coefCi %>%
  mutate(OR_2.5=exp(coefCi$Coef_2.5)
         , OR_97.5=exp(coefCi$Coef_97.5)
         )

write.table(coefAndOrCi
            , paste0(result_dir, 'coefAndOrCi.csv')
            , row.names = F
            , sep=',')

