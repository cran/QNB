.fittedW <- function(p,q,fit){ 
  l <- log(q+1)
  c <- c(rep(1,length(p)))
  data=data.frame(cbind(p,l,c))
  w_fit <- predict(fit,data)
  return(w_fit)
}