.sizeFactor2 <- function(n) {
  temp <- log(colSums(n))
  temp <- temp - mean(temp)
  s_size <- exp(temp)
  return(s_size)
}