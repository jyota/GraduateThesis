# Perform upper quartile normalization without additional scaling for Fisher's exact test.

upperQuartileNormalize <- function(data) {
    sumtot <- rowSums(data)
    counts0 <- which(sumtot == 0)
    if (length(counts0) > 0) {
      dataAdjusted <- data[-counts0,]
    } else {
      dataAdjusted <- data
    }
    q3 <- apply(dataAdjusted, 2, quantile, probs = 0.75)
    d <- mean(q3)/q3
    dataNormalized <- t(t(data)*d)
  
  return(na.omit(dataNormalized))  
}
