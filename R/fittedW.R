.fittedW <- function(p,q,fit){ 
  l <- log(q+1)
  data=data.frame(cbind(p,l))
  w_fit <- predict(fit,data)
  return(w_fit)
}