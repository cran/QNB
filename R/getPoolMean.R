.getPoolMean <- function(meth,unmeth)
{
  s <- .sizeFactor2(cbind(meth,unmeth))
  s_t <- s[1:length(meth[1,])]
  s_c <- s[(length(meth[1,])+1):(length(cbind(meth,unmeth)[1,]))]
  # estimate probability of methylation under a condition
  p <- .estimateP(meth,unmeth, s_t,s_c)
  
  # estimate the abundance of feature
  q <- .estimateQ(meth,unmeth,s_t,s_c,p,useAll=TRUE)
  
  # estimate size e
  e <- .estimateE(meth,unmeth,s_t,s_c,q)
  
  res <- list(s_t,s_c,p,q,e)
  return(res)
}