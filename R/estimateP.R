.estimateP <- function(meth, unmeth, size_t, size_c) {
  temp_t <- t(t(meth)/size_t)
  temp_c <- t(t(unmeth)/size_c)
  temp_n <- temp_t+temp_c
  p <- rowSums(temp_t)/rowSums(temp_n)
  return(p)
  # which(is.na(p))
  # which(is.infinite(p))  
}