

computeWeights_DiscreteClasses <- function(label_vec)
{
  itsLevels <- levels(factor(label_vec))
  counts <- vector()
  
  for (iLevel in 1:length(itsLevels))
  {
    sizeThisClass <- sum(label_vec == itsLevels[iLevel])
    counts <- c(counts, sizeThisClass)
  }
  
  weight_vec <- rep(0.0, length(label_vec))
  for (iLevel in 1:length(itsLevels))
  {
    weight_vec[which(label_vec == itsLevels[iLevel])] <- 1/(counts[iLevel]+1e-6) / (sum(1/counts) + +1e-6)
  }
  
  return (weight_vec)
}

computeWeights <- function(label_vec)
{
  weight_vec <- computeWeights_DiscreteClasses(label_vec)
  
  return (weight_vec)
}