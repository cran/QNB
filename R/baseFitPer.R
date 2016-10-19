.baseFitPer <- function(meth1,meth2,unmeth1,unmeth2,p1,p2,q,e1,e2,s_t1,s_t2,s_c1,s_c2){
  w_t1 <-.calculateW(meth1,s_t1,e1)
  w_t2 <-.calculateW(meth2,s_t2,e2)
  w_c1 <-.calculateW(unmeth1,s_c1,e1)
  w_c2 <-.calculateW(unmeth2,s_c2,e2)
  
  # locfit 
  fit_t1 <- .locfitW(p1,q,w_t1)
  fit_t2 <- .locfitW(p2,q,w_t2)
  fit_c1 <- .locfitW(p1,q,w_c1)
  fit_c2 <- .locfitW(p2,q,w_c2)
  
  res <- list(fit_t1,fit_t2,fit_c1,fit_c2)
  return(res)
}