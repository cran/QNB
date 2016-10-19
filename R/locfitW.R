.locfitW <- function(p,q,w) {
  l <- log(q+1)
  data=data.frame(cbind(p,l,w))
  ID <- which(rowSums(is.na(data))>0)
  #  data <- data[-ID,]
  fit=locfit(w~lp(p,l),data=data,family="gamma")
  return(fit)
}