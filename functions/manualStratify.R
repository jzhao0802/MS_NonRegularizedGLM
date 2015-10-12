library(caret)

manualStratify <- function(y, kFolds)
{
  # to prevent unstratified data subsets (due to the data size)
  # use a check
  bStratValid <- FALSE
  nPosTot <- sum(y == 1)
  nNegTot <- sum(y == 0)
  pAvPos <- nPosTot / (nPosTot + nNegTot)
  pAvNeg <- nNegTot / (nPosTot + nNegTot)
  pPosThreshL <- pAvPos * 0.8
  pNegThreshL <- pAvNeg * 0.8
  pPosThreshH <- pAvPos * 1.25
  pNegThreshH <- pAvNeg * 1.25
  nTrialsLimit <- 100
  iTrialsLimit <- 1
  while (bStratValid == FALSE)
  {
    folds <- createFolds(y, k=kFolds, returnTrain=TRUE)
    # check
    bLocalEvidence <- TRUE
    for (iFold in 1:kFolds)
    {
      rest_ids <- folds[[iFold]]
      test_ids <- (1:length(y))[-rest_ids]
      nPos <- sum(y[test_ids] == 1)
      nNeg <- sum(y[test_ids] == 0)
      pPos <- nPos / (nPos + nNeg)
      pNeg <- 1 - pPos
      
      bLocalEvidence <- bLocalEvidence & 
        ((pPos >= pPosThreshL) & (pNeg >= pNegThreshL) & (pPos <= pPosThreshH) & (pNeg <= pNegThreshH))
    }
    
    if (bLocalEvidence)
    {
      bStratValid <- TRUE
    }
    
    # 
    iTrialsLimit = iTrialsLimit + 1
    if (iTrialsLimit > nTrialsLimit)
    {
      break
    }
  }
  
  if (iTrialsLimit > nTrialsLimit)
  {
    cat(paste("Not able to stratify. \n"))
  }
  
  return (folds)
}