.sizeFactor <- function(n, useTotal=FALSE) {
  if (useTotal) {
    temp <- log(colSums(n))
    temp <- temp - mean(temp)
    s_size <- exp(temp)
  } else {
    n <- pmax(n,1e-5)
    log_n <- log(n)
    pseudo <- rowMeans(log_n)
    ratio <- log_n-pseudo
    s_size <- exp(apply(ratio,2,median)) }
  return(s_size)
}