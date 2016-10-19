.calculateW <- function(meth,size_t,e_t){
  temp_t <- t(t(meth)/size_t)
  q <- temp_t/e_t
  #  q[is.na(q)] <- 0
  resi <- q-rowMeans(q)
  w <- rowSums(resi^2)/(length(size_t)-1)
  w <- pmax(w,1e-8)
  return(w)
}