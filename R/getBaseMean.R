.getBaseMean <- function(meth1,meth2,unmeth1,unmeth2){
  s <- .sizeFactor2(cbind(meth1,meth2,unmeth1,unmeth2))
  s_t1 <- s[1:length(meth1[1,])]
  s_t2 <- s[(length(meth1[1,])+1):(length(cbind(meth1,meth2)[1,]))]
  s_c1 <- s[(length(cbind(meth1,meth2)[1,])+1):(length(cbind(meth1,meth2,unmeth1)[1,]))]
  s_c2 <- s[(length(cbind(meth1,meth2,unmeth1)[1,])+1):(length(cbind(meth1,meth2,unmeth1,unmeth2)[1,]))]
  
  # estimate probability of methylation under a condition
  p1 <- .estimateP(meth1, unmeth1, s_t1, s_c1)
  p2 <- .estimateP(meth2, unmeth2, s_t2, s_c2)
  p0 <- .estimateP(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2))
  
  # estimate the abundance of feature
  q0 <- .estimateQ(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2),p0,useAll=TRUE)
  q1 <- .estimateQ(meth1,unmeth1,s_t1,s_c1,p0,useAll=TRUE)
  q2 <- .estimateQ(meth2,unmeth2,s_t2,s_c2,p0,useAll=TRUE)
  
  # estimate size e
  e1 <- .estimateE(meth1,unmeth1,s_t1,s_c1,q0)
  e2 <- .estimateE(meth2,unmeth2,s_t2,s_c2,q0)
  
  res <- list(s_t1,s_t2,s_c1,s_c2,p0,p1,p2,q0,q1,q2,e1,e2)
  return(res)
}