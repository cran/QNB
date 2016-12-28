.locfitW <- function(p,q,w) {
  l <- log(q+1)
  c=c(rep(1,length(p)))
  data=data.frame(cbind(p,l,c,w))
  ID <- which(rowSums(is.na(data))>0)
  #  data <- data[-ID,]
  fit=locfit(w~lp(p,l,c),data=data,family="gamma")
  return(fit)
}