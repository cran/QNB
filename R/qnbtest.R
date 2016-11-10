qnbtest <-
function(meth1,meth2,unmeth1,unmeth2,
                    mode="auto",
                    plot.dispersion=TRUE,
                    output.dir = NA) {
  print("Estimating dispersion for each RNA methylation site, this will take a while ...")
  if(mode=="per-condition"){
    mean <- .getBaseMean(meth1,meth2,unmeth1,unmeth2)
    s_t1 <- mean[[1]]
    s_t2 <- mean[[2]]
    s_c1 <- mean[[3]]
    s_c2 <- mean[[4]]
    p0 <- mean[[5]]
    p1 <- mean[[6]]
    p2 <- mean[[7]]
    q0 <- mean[[8]]
    q1 <- mean[[9]]
    q2 <- mean[[10]]
    e1 <- mean[[11]]
    e2 <- mean[[12]]
    res <- .baseFitPer(meth1,meth2,unmeth1,unmeth2,p1,p2,q0,e1,e2,s_t1,s_t2,s_c1,s_c2)
    fit_t1<-res[[1]]
    fit_t2<-res[[2]]
    fit_c1<-res[[3]]
    fit_c2<-res[[4]]
  }else if(mode=="pooled"){
    meth<-rbind(meth1,meth2)
    unmeth<-rbind(unmeth1,unmeth2)
    mean <- .getPoolMean(meth,unmeth)
    s_t1 <- mean[[1]]
    s_t2 <- mean[[1]]
    s_c1 <- mean[[2]]
    s_c2 <- mean[[2]]
    p1 <- mean[[3]]
    p2 <- p1
    p0 <- p1
    q0 <- mean[[4]]
    q1 <- q0
    q2 <- q0
    e1 <- mean[[5]]
    e2 <- e1
    res <- .baseFitCondition(meth,unmeth,p0,q0,e1,s_t1,s_c1)
    fit_t1<-res[[1]]
    fit_t2<-res[[1]]
    fit_c1<-res[[2]]
    fit_c2<-res[[2]]
  }else if(mode=="blind"){
    meth<-cbind(meth1,meth2)
    unmeth<-cbind(unmeth1,unmeth2)
    mean <- .getPoolMean(meth,unmeth)
    s_t1 <- mean[[1]]
    s_t2 <- mean[[1]]
    s_c1 <- mean[[2]]
    s_c2 <- mean[[2]]
    p1 <- mean[[3]]
    p2 <- p1
    p0 <- p1
    q0 <- mean[[4]]
    q1 <- q0
    q2 <- q0
    e1 <- mean[[5]]
    e2 <- e1
    res <- .poolFit(meth,unmeth,p0,q0,e1,s_t1,s_c1)
    fit_t1<-res[[1]]
    fit_t2<-res[[1]]
    fit_c1<-res[[2]]
    fit_c2<-res[[2]]
  }else if(mode=="auto"){
    rep1 <- ncol(meth1)
    rep2 <- ncol(meth2)
    if((rep1==rep2&&rep1>1)||(rep1!=rep2&&min(rep1,rep2)>1)){
      mean <- .getBaseMean(meth1,meth2,unmeth1,unmeth2)
      s_t1 <- mean[[1]]
      s_t2 <- mean[[2]]
      s_c1 <- mean[[3]]
      s_c2 <- mean[[4]]
      p0 <- mean[[5]]
      p1 <- mean[[6]]
      p2 <- mean[[7]]
      q0 <- mean[[8]]
      q1 <- mean[[9]]
      q2 <- mean[[10]]
      e1 <- mean[[11]]
      e2 <- mean[[12]]
      res <- .baseFitPer(meth1,meth2,unmeth1,unmeth2,p1,p2,q0,e1,e2,s_t1,s_t2,s_c1,s_c2)
      fit_t1<-res[[1]]
      fit_t2<-res[[2]]
      fit_c1<-res[[3]]
      fit_c2<-res[[4]]
    }else if((rep1==rep2&&rep1==1)||(rep1!=rep2&&min(rep1,rep2)<2)){
      meth<-cbind(meth1,meth2)
      unmeth<-cbind(unmeth1,unmeth2)
      mean <- .getPoolMean(meth,unmeth)
      s_t1 <- mean[[1]]
      s_t2 <- mean[[1]]
      s_c1 <- mean[[2]]
      s_c2 <- mean[[2]]
      p1 <- mean[[3]]
      p2 <- p1
      p0 <- p1
      q0 <- mean[[4]]
      q1 <- q0
      q2 <- q0
      e1 <- mean[[5]]
      e2 <- e1
      res <- .poolFit(meth,unmeth,p0,q0,e1,s_t1,s_c1)
      fit_t1<-res[[1]]
      fit_t2<-res[[1]]
      fit_c1<-res[[2]]
      fit_c2<-res[[2]]
    }
  }
  
  
  #   s <- .sizeFactor2(cbind(meth1,meth2,unmeth1,unmeth2))
  #   s_t1 <- s[1:length(meth1[1,])]
  #   s_t2 <- s[(length(meth1[1,])+1):(length(cbind(meth1,meth2)[1,]))]
  #   s_c1 <- s[(length(cbind(meth1,meth2)[1,])+1):(length(cbind(meth1,meth2,unmeth1)[1,]))]
  #   s_c2 <- s[(length(cbind(meth1,meth2,unmeth1)[1,])+1):(length(cbind(meth1,meth2,unmeth1,unmeth2)[1,]))]
  #   
  #   # estimate probability of methylation under a condition
  #   p1 <- .estimateP(meth1, unmeth1, s_t1, s_c1)
  #   p2 <- .estimateP(meth2, unmeth2, s_t2, s_c2)
  #   p0 <- .estimateP(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2))
  #   
  #   # estimate the abundance of feature
  #   q0 <- .estimateQ2(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2),p0,useAll=TRUE)
  #   q1 <- .estimateQ2(meth1,unmeth1,s_t1,s_c1,p0,useAll=TRUE)
  #   q2 <- .estimateQ2(meth2,unmeth2,s_t2,s_c2,p0,useAll=TRUE)
  #   
  #   # estimate size e
  #     e1 <- .estimateE(meth1,unmeth1,s_t1,s_c1,q0)
  #     e2 <- .estimateE(meth2,unmeth2,s_t2,s_c2,q0)
  #   
  #   # calculate methylation reads count variance on common scale, condition 1
  #   w_t1 <-.calculateW(meth1,s_t1,e1)
  #   w_t2 <-.calculateW(meth2,s_t2,e2)
  #   w_c1 <-.calculateW(unmeth1,s_c1,e1)
  #   w_c2 <-.calculateW(unmeth2,s_c2,e2)
  #   
  #   # locfit 
  #   fit_t1 <- .locfitW(p1,q0,w_t1)
  #   fit_t2 <- .locfitW(p2,q0,w_t2)
  #   fit_c1 <- .locfitW(p1,q0,w_c1)
  #   fit_c2 <- .locfitW(p2,q0,w_c2)
  
  if (is.na(output.dir)) {
    output.dir <- getwd()
  }
  
  path <- paste(output.dir,"dispersion.pdf",sep = '/')
  if(plot.dispersion){
    .plotDispersion(fit_t1,fit_c1,fit_t2,fit_c2,path)
  }
  
  
  # calculate z
  z_t1 <- .calculateZ(q0,p1,s_t1,e1)
  z_t2 <- .calculateZ(q0,p2,s_t2,e2)
  z_c1 <- .calculateZ(q0,(1-p1),s_c1,e1)
  z_c2 <- .calculateZ(q0,(1-p2),s_c2,e2)
  
  # get estimate w
  w_fit_t1 <- .fittedW(p0,q0,fit_t1)
  w_fit_t2 <- .fittedW(p0,q0,fit_t2)
  w_fit_c1 <- .fittedW(p0,q0,fit_c1)
  w_fit_c2 <- .fittedW(p0,q0,fit_c1)
  
  # get estimate of upi
  ups_t1 <- pmax(w_fit_t1 - z_t1, 1e-8)
  ups_t2 <- pmax(w_fit_t2 - z_t2, 1e-8)
  ups_c1 <- pmax(w_fit_c1 - z_c1, 1e-8)
  ups_c2 <- pmax(w_fit_c2 - z_c2, 1e-8)
  
  # get all means
  mu_t1 <- (e1*q0*p0)%*%t(as.numeric(s_t1))
  mu_t2 <- (e2*q0*p0)%*%t(as.numeric(s_t2))
  mu_c1 <- (e1*q0*(1-p0))%*%t(as.numeric(s_c1))
  mu_c2 <- (e2*q0*(1-p0))%*%t(as.numeric(s_c2))
  
  # get all variance
  raw_t1 <- (e1%*%t(s_t1))^2*ups_t1
  raw_t2 <- (e2%*%t(s_t2))^2*ups_t2
  raw_c1 <- (e1%*%t(s_c1))^2*ups_c1
  raw_c2 <- (e2%*%t(s_c2))^2*ups_c2
  
  # put mu together
  mu1_t <- rowSums(mu_t1)
  mu2_t <- rowSums(mu_t2)
  mu1_c <- rowSums(mu_c1)
  mu2_c <- rowSums(mu_c2)
  
  # put size together
  size1_t <- (mu1_t^2)/rowSums(raw_t1)
  size2_t <- (mu2_t^2)/rowSums(raw_t2)
  size1_c <- (mu1_c^2)/rowSums(raw_c1)
  size2_c <- (mu2_c^2)/rowSums(raw_c2)
  
  # observation together
  t1 <- rowSums(meth1)
  t2 <- rowSums(meth2)
  c1 <- rowSums(unmeth1)
  c2 <- rowSums(unmeth2)
  t <- t1 + t2
  n1 <- t1 + c1
  n2 <- t2 + c2
  
  raw <- (rowSums(raw_t1)+rowSums(raw_t2)+rowSums(raw_c1)+rowSums(raw_c2))/4
  # go to test
  res <- .quadNBtest(t1,t,n1,n2,mu1_t,mu2_t,mu1_c,mu2_c,size1_t,size2_t,size1_c,size2_c)
  
  # add fc
  fc <- log2(p1/p2)
  m1 <- rowSums(t(t(meth1)/s_t1))
  m2 <- rowSums(t(t(meth2)/s_t2))
  u1 <- rowSums(t(t(unmeth1)/s_c1))
  u2 <- rowSums(t(t(unmeth2)/s_c2))
  mfc <- log2(m1)-log2(m2)
  ufc <- log2(u1)-log2(u2)
  
  padj <- p.adjust( res[,1], method="BH" )
  res <- data.frame(res,fc,q0,padj)
    #res <- res2[,c(1,3:7)]
  #path=getwd()
  path <- paste(output.dir,"dif_meth.xls",sep = '/')
  #path=paste(path,"dif_meth.xls",sep="/")
  write.table(res,path,sep='\t',col.names=TRUE,row.names =FALSE)
  return(res)}
