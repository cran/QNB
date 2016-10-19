.calculateZ <- function(q,p,size,e){
  temp <- p*q/length(size)
  
  norow <- length(q)
  nocol <- length(size)
  temp2 <- matrix(1,nrow=norow, ncol=nocol )
  temp3 <- t(t(temp2)/size)/e
  
  z <- rowSums(temp3)*temp
  
  #  z[is.infinite(z)] <- 0
  #  z[is.na(z)] <- 0
  return(z)
}