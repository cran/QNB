.estimateE <- function(meth,unmeth,size_t,size_c,q){
  temp_t <- t(t(meth)/size_t)
  temp_c <- t(t(unmeth)/size_c)
  temp_n <- temp_t+temp_c
  e <- rowSums(temp_n)/length(size_t)/q
}