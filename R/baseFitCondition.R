.baseFitCondition <- function(meth,unmeth,p,q,e,s_t,s_c){
  w_t <-.calculateW(meth,s_t,e)
  w_c <-.calculateW(unmeth,s_c,e)
  
  # locfit 
  fit_t <- .locfitW(p,q,w_t)
  fit_c <- .locfitW(p,q,w_c)
  
  res <- list(fit_t,fit_c)
  return(res)
}