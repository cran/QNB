.getBaseMean <- function(meth1,meth2,unmeth1,unmeth2,s_t1,s_t2,s_c1,s_c2){
  
  # estimate probability of methylation under a condition
  p1 <- .estimateP(meth1, unmeth1, s_t1, s_c1)
  p2 <- .estimateP(meth2, unmeth2, s_t2, s_c2)
  p0 <- .estimateP(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2))
  
  # estimate the abundance of feature
  q0 <- .estimateQ(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2),p0)
  q1 <- .estimateQ(meth1,unmeth1,s_t1,s_c1,p0)
  q2 <- .estimateQ(meth2,unmeth2,s_t2,s_c2,p0)
  
  # estimate size e
  e1 <- .estimateE(meth1,unmeth1,s_t1,s_c1,q0)
  e2 <- .estimateE(meth2,unmeth2,s_t2,s_c2,q0)
  
  res <- list(p0=p0,p1=p1,p2=p2,q0=q0,q1=q1,q2=q2,e1=e1,e2=e2)
  return(res)
}