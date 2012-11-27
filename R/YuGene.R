YuGene <-
function(data.prop, progressBar = TRUE){
  YuGene= data.frame(matrix(nrow = nrow(data.prop), ncol = ncol(data.prop)) )
  rownames(YuGene) <- rownames(data.prop)
  colnames(YuGene) <- colnames(data.prop)
  if(progressBar == TRUE) pb <- txtProgressBar(style=3);
  
  for (col in 1:ncol(data.prop)) {                # for each sample (col)
    sample <- data.prop[,col]
    o <- order(sample, decreasing=TRUE) # highest to lowest expression
    m <- cbind(1:length(sample),sample)  
    byOrder <- m[o,] # sort from highest expression to lowest
    sampleTotal <- sum(sample)
    cs <- cumsum(byOrder[,2])/sampleTotal  # get cummulative sum of expression
    byOrder[,2] <- cs  # replace with cumulative values
    byIndex <- byOrder[order(byOrder[,1]),] # re-sort by index
    YuGene[,col] <- byIndex[,2] # replace with cum prop values
    if(progressBar == TRUE) setTxtProgressBar(pb,col/ncol(data.prop)) # show progress
  }
  return(1-as.matrix(YuGene,nrow = nrow(data.prop), ncol = ncol(data.prop)))
}
