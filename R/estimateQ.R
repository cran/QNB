.estimateQ <- function(meth,unmeth,size_t,size_c,p,useAll=FALSE){
  if (useAll) {
    temp_t <- t(t(meth)/size_t)
    temp_c <- t(t(unmeth)/size_c)
    temp_n <- temp_t+temp_c
    q <- rowSums(temp_n)/length(size_t)
  } else {
    temp_c <- t(t(unmeth)/size_c)
    qc <- rowMeans(temp_c)
    q <- qc/(1-p)
  }
  
  q[is.na(q)] <- 0
  return(q)
  # which(is.na(q))
  # which(is.infinite(q))  
}