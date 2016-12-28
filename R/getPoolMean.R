.getPoolMean <- function(meth,unmeth,s_t,s_c)
{
  # estimate probability of methylation under a condition
  p <- .estimateP(meth,unmeth, s_t,s_c)
  
  # estimate the abundance of feature
  q <- .estimateQ(meth,unmeth,s_t,s_c,p,useAll=TRUE)
  
  # estimate size e
  e <- .estimateE(meth,unmeth,s_t,s_c,q)
  
  res <- list(p=p,q=q,e=e)
  return(res)
}