library(caret)
library(compiler)

DivideIntoFolds <- function(IDs, kFolds)
{
  if (length(IDs) < kFolds)
    stop("Error! Number of data is smaller than kFolds!")
  
  avNDataPerFold <- length(IDs) / kFolds
  pointer <- 1
  frac <- 0
  
  IDsAllFolds <- list()
  
  for (iFold in 1:kFolds)
  {
    if (iFold != kFolds)
    {
      frac <- frac + avNDataPerFold - floor(avNDataPerFold)
      if (frac >= 1)
      {
        IDsThisFold <- IDs[pointer:(pointer+floor(avNDataPerFold))]
        frac <- frac-1
      } else
        IDsThisFold <- IDs[pointer:(pointer+floor(avNDataPerFold)-1)]
      
      pointer <- pointer + length(IDsThisFold)
      
    } else 
    {
      IDsThisFold <- IDs[pointer:length(IDs)]
    }
    
    IDsAllFolds[[iFold]] <- IDsThisFold
  }
  
  return (IDsAllFolds)
}

manualStratify <- cmpfun(function(y, k_folds, seed=NULL)
{
  if (!(all(levels(y) %in% c(0,1))))
  {
    stop(paste("Error! Invalid y-values for stratification. ", 
               "stratifySmallSample currently only supports one type ",
               "of y values: c(0, 1). "))
  }
  
  # get the positive and negative
  
  pos_indices <- which(y==1)
  if (length(pos_indices) < k_folds)
    stop("Error! Too few positives. StratifyEasyDifficultPositives failed.")
  neg_indices <- which(y==0)
  if (length(neg_indices) < k_folds)
    stop("Error! Too few negatives. StratifyEasyDifficultPositives failed.")
  if (!is.null(seed))
    set.seed(seed=seed)
  pos_indices <- sample(pos_indices)
  neg_indices <- sample(neg_indices)
  
  # 
  
  pos_ids_allfolds <- DivideIntoFolds(pos_indices, k_folds)
  neg_ids_allfolds <- DivideIntoFolds(neg_indices, k_folds)
  
  folds <- list()
  
  for (i_fold in 1:k_folds)
  {
    IDsInThisFold <- c(pos_ids_allfolds[[i_fold]], 
                       neg_ids_allfolds[[i_fold]])
    folds[[i_fold]] <- 
      (1:length(y))[which(!((1:length(y)) %in% IDsInThisFold))]
  }
  
  return (folds)
}, option=list(optimize=3))

stratifyFoldIDs <- function(y, k_folds, seed=NULL)
{
  ids <- 1:k_folds
  if (!is.null(seed))
    set.seed(seed=seed)
  ids_every_pos <- sample(rep(ids, length.out=sum(y==1)))
  ids_every_neg <- sample(rep(ids, length.out=sum(y!=1)))
  
  ids_every_datum <- rep(-1, length(y))
  ids_every_datum[which(y==1)] <- ids_every_pos
  ids_every_datum[which(y!=1)] <- ids_every_neg
  
  return (ids_every_datum)
}