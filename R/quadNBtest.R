.quadNBtest <- function(t1,t,n1,n2,mu1_t,mu2_t,mu1_c,mu2_c,size1_t,size2_t,size1_c,size2_c){
  nrows <- length(t)
  pval4 <- rep(1,nrows)
  
  t2 <- t-t1
  
  for (irow in 1:nrows) {
    
    trip <- t[irow]
    
    if (trip<1) {p <- NA} else {
      
      trip_t1 <- 0:t[irow]
      trip_t2 <- t[irow] - trip_t1
      trip_c1 <- n1[irow] - trip_t1
      trip_c2 <- n2[irow] - trip_t2
      
      
      p1 <- dnbinom(x=trip_t1, size=size1_t[irow], mu=mu1_t[irow], log = TRUE)
      p2 <- dnbinom(x=trip_t2, size=size2_t[irow], mu=mu2_t[irow], log = TRUE)
      p3 <- dnbinom(x=trip_c1, size=size1_c[irow], mu=mu1_c[irow], log = TRUE)
      p4 <- dnbinom(x=trip_c2, size=size2_c[irow], mu=mu2_c[irow], log = TRUE)
      
      #p4
      p <- p1+p2+p3+p4
      p <- p-max(p)
      p <- exp(p)/sum(exp(p))
      p <- min(1,2*min(sum(p[1:(t1[irow]+1)]),sum(p[(t1[irow]+1):(t[irow]+1)])))- p[(t1[irow]+1)]/2  
      pval4[irow] <- p
    }
    
  }
  
  mu <- (mu1_t+mu2_t+mu1_c+mu2_c)/4
  
  res <- data.frame(pval4)
  return(res)
}