.plotDispersion <-function(fit_t1,fit_c1,fit_t2,fit_c2,path) {
  p <- rep(seq(from=0,to=1,by=0.01),10)
  q <- rep(seq(from=0,to=10,by=0.1),10)
  p <- matrix(p,nrow = 101, ncol = 101,byrow=TRUE)
  q <- matrix(q,nrow = 101, ncol = 101,byrow=FALSE)
  
  p <- p[1:10201]
  q <- q[1:10201]
  
  w1 <- .fittedW(p,q,fit_t1)
  w1 <- matrix(log(w1),nrow = 101, ncol = 101,byrow=FALSE)
  w2 <-  .fittedW(p,q,fit_c1)
  w2 <- matrix(log(w2),nrow = 101, ncol = 101,byrow=FALSE)
  
  w3 <- .fittedW(p,q,fit_t2)
  w3 <- matrix(log(w3),nrow = 101, ncol = 101,byrow=FALSE)
  w4 <- .fittedW(p,q,fit_c2)
  w4 <- matrix(log(w4),nrow = 101, ncol = 101,byrow=FALSE)
  
  pdf(path,height=8,width=7)
  #pdf("C:/Users/S41-70/Documents/dispersion.pdf",height=4,width=7)
  par(mfrow=c(2,2))
  
   
  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=10,by=0.1),w2)
  title(main="log(w)",sub="untreated Input",ylab="log(q+1)", xlab="p")
  
  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=10,by=0.1),w1)
  title(main="log(w)",sub="untreated IP",ylab="log(q+1)", xlab="p") 
  
  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=10,by=0.1),w4)
  title(main="log(w)",sub="treated Input",ylab="log(q+1)", xlab="p") 
  
  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=10,by=0.1),w3)
  title(main="log(w)",sub="treated IP",ylab="log(q+1)", xlab="p") 
  dev.off()
}