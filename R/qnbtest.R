qnbtest <-
function(control_ip,treated_ip,control_input,treated_input,
                    size.factor=NA,
                    mode="auto",
                    plot.dispersion=TRUE,
                    output.dir = NA) {
  print("Estimating dispersion for each RNA methylation site, this will take a while ...")
  if(mode=="per-condition"){
    if(anyNA(size.factor)){
      s <- .sizeFactor2(cbind(control_ip,treated_ip,control_input,treated_input))
      s_t1 <- s[1:length(control_ip[1,])]
      s_t2 <- s[(length(control_ip[1,])+1):(length(cbind(control_ip,treated_ip)[1,]))]
      s_c1 <- s[(length(cbind(control_ip,treated_ip)[1,])+1):(length(cbind(control_ip,treated_ip,control_input)[1,]))]
      s_c2 <- s[(length(cbind(control_ip,treated_ip,control_input)[1,])+1):(length(cbind(control_ip,treated_ip,control_input,treated_input)[1,]))]
    }else{
      s_t1 <- size.factor$control_ip
      s_t2 <- size.factor$treated_ip
      s_c1 <- size.factor$control_input
      s_c2 <- size.factor$treated_ip
    }
    mean <- .getBaseMean(control_ip,treated_ip,control_input,treated_input,s_t1,s_t2,s_c1,s_c2)
    p0 <- mean[[1]]
    p1 <- mean[[2]]
    p2 <- mean[[3]]
    q0 <- mean[[4]]
    q1 <- mean[[5]]
    q2 <- mean[[6]]
    e1 <- mean[[7]]
    e2 <- mean[[8]]
    res <- .baseFitPer(control_ip,treated_ip,control_input,treated_input,p1,p2,q0,e1,e2,s_t1,s_t2,s_c1,s_c2)
    fit_t1<-res[[1]]
    fit_t2<-res[[2]]
    fit_c1<-res[[3]]
    fit_c2<-res[[4]]
   
  }else if(mode=="pooled"){
    meth<-rbind(control_ip,treated_ip)
    unmeth<-rbind(control_input,treated_input)
    if(anyNA(size.factor)){
      s <- .sizeFactor2(cbind(meth,unmeth))
      s_t <- s[1:length(meth[1,])]
      s_c <- s[(length(meth[1,])+1):(length(cbind(meth,unmeth)[1,]))]
      s_t1 <- s_t
      s_t2 <- s_t
      s_c1 <- s_c
      s_c2 <- s_c
    }else{
      s_t1 <- size.factor$control_ip
      s_t2 <- size.factor$treated_ip
      s_c1 <- size.factor$control_input
      s_c2 <- size.factor$treated_ip
      s_t <- c(s_t1,s_t2)
      s_c <- c(s_c1,s_c2)
    }
    mean <- .getPoolMean(meth,unmeth,s_t,s_c)
    p1 <- mean[[1]]
    p2 <- p1
    p0 <- p1
    q0 <- mean[[2]]
    q1 <- q0
    q2 <- q0
    e1 <- mean[[3]]
    e2 <- e1
    res <- .baseFitCondition(meth,unmeth,p0,q0,e1,s_t,s_c)
    fit_t1<-res[[1]]
    fit_t2<-res[[1]]
    fit_c1<-res[[2]]
    fit_c2<-res[[2]]
    
  }else if(mode=="blind"){
    meth<-cbind(control_ip,treated_ip)
    unmeth<-cbind(control_input,treated_input)
    if(anyNA(size.factor)){
      s <- .sizeFactor2(cbind(meth,unmeth))
      s_t <- s[1:length(meth[1,])]
      s_c <- s[(length(meth[1,])+1):(length(cbind(meth,unmeth)[1,]))]
      s_t1 <- s_t
      s_t2 <- s_t
      s_c1 <- s_c
      s_c2 <- s_c
    }else{
      s_t1 <- size.factor$control_ip
      s_t2 <- size.factor$treated_ip
      s_c1 <- size.factor$control_input
      s_c2 <- size.factor$treated_ip
      s_t <- c(s_t1,s_t2)
      s_c <- c(s_c1,s_c2)
    }
    mean <- .getPoolMean(meth,unmeth,s_t,s_c)
    p1 <- mean[[1]]
    p2 <- p1
    p0 <- p1
    q0 <- mean[[2]]
    q1 <- q0
    q2 <- q0
    e1 <- mean[[3]]
    e2 <- e1
    res <- .poolFit(meth,unmeth,p0,q0,e1,s_t,s_c)
    fit_t1<-res[[1]]
    fit_t2<-res[[1]]
    fit_c1<-res[[2]]
    fit_c2<-res[[2]]
  }else if(mode=="auto"){
    rep1 <- ncol(control_ip)
    rep2 <- ncol(treated_ip)
    if((rep1==rep2&&rep1>1)||(rep1!=rep2&&min(rep1,rep2)>1)){
      if(anyNA(size.factor)){
        s <- .sizeFactor2(cbind(control_ip,treated_ip,control_input,treated_input))
        s_t1 <- s[1:length(control_ip[1,])]
        s_t2 <- s[(length(control_ip[1,])+1):(length(cbind(control_ip,treated_ip)[1,]))]
        s_c1 <- s[(length(cbind(control_ip,treated_ip)[1,])+1):(length(cbind(control_ip,treated_ip,control_input)[1,]))]
        s_c2 <- s[(length(cbind(control_ip,treated_ip,control_input)[1,])+1):(length(cbind(control_ip,treated_ip,control_input,treated_input)[1,]))]
      }else{
        s_t1 <- size.factor$control_ip
        s_t2 <- size.factor$treated_ip
        s_c1 <- size.factor$control_input
        s_c2 <- size.factor$treated_ip
      }
      mean <- .getBaseMean(control_ip,treated_ip,control_input,treated_input,s_t1,s_t2,s_c1,s_c2)
      p0 <- mean[[1]]
      p1 <- mean[[2]]
      p2 <- mean[[3]]
      q0 <- mean[[4]]
      q1 <- mean[[5]]
      q2 <- mean[[6]]
      e1 <- mean[[7]]
      e2 <- mean[[8]]
      res <- res <- .baseFitPer(control_ip,treated_ip,control_input,treated_input,p1,p2,q0,e1,e2,s_t1,s_t2,s_c1,s_c2)
      fit_t1<-res[[1]]
      fit_t2<-res[[2]]
      fit_c1<-res[[3]]
      fit_c2<-res[[4]]

    }else if((rep1==rep2&&rep1==1)||(rep1!=rep2&&min(rep1,rep2)<2)){
      meth=cbind(control_ip,treated_ip)
      unmeth=cbind(control_input,treated_input)
      
      if(anyNA(size.factor)){
        s <- .sizeFactor2(cbind(meth,unmeth))
        s_t <- s[1:length(meth[1,])]
        s_c <- s[(length(meth[1,])+1):(length(cbind(meth,unmeth)[1,]))]
        s_t1 <- s_t
        s_t2 <- s_t
        s_c1 <- s_c
        s_c2 <- s_c
      }else{
        s_t1 <- size.factor$control_ip
        s_t2 <- size.factor$treated_ip
        s_c1 <- size.factor$control_input
        s_c2 <- size.factor$treated_ip
        s_t <- c(s_t1,s_t2)
        s_c <- c(s_c1,s_c2)
      }
      mean <- .getPoolMean(meth,unmeth,s_t,s_c)
      p1 <- mean[[1]]
      p2 <- p1
      p0 <- p1
      q0 <- mean[[2]]
      q1 <- q0
      q2 <- q0
      e1 <- mean[[3]]
      e2 <- e1
      
      res <- .poolFit(meth,unmeth,p0,q0,e1,s_t,s_c)
      fit_t1<-res[[1]]
      fit_t2<-res[[1]]
      fit_c1<-res[[2]]
      fit_c2<-res[[2]]
    
    }
  }
  

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
  t1 <- rowSums(control_ip)
  t2 <- rowSums(treated_ip)
  c1 <- rowSums(control_input)
  c2 <- rowSums(treated_input)
  t <- t1 + t2
  n1 <- t1 + c1
  n2 <- t2 + c2
  
  raw <- (rowSums(raw_t1)+rowSums(raw_t2)+rowSums(raw_c1)+rowSums(raw_c2))/4
  # go to test
  res <- .quadNBtest(t1,t,n1,n2,mu1_t,mu2_t,mu1_c,mu2_c,size1_t,size2_t,size1_c,size2_c)
  
  # add fc
  p1 <- .estimateP(control_ip,control_input,s_t1,s_c1)
  p2 <- .estimateP(treated_ip,treated_input,s_t2,s_c2)
  log2.RR <- log2(p2/p1)
  
  p.treated <- p2
  p.control <- p1
  
  log2.OR <- log2((rowSums(t(t(treated_ip)/s_t2))/rowSums(t(t(treated_input)/s_c2)))/(rowSums(t(t(control_ip)/s_t1))/rowSums(t(t(control_input)/s_c1))))
  m1 <- rowSums(t(t(control_ip)/s_t1))
  m2 <- rowSums(t(t(treated_ip)/s_t2))
  u1 <- rowSums(t(t(control_input)/s_c1))
  u2 <- rowSums(t(t(treated_input)/s_c2))
  mfc <- log2(m1)-log2(m2)
  ufc <- log2(u1)-log2(u2)
  
  padj <- p.adjust( res[,1], method="BH" )
  res <- data.frame(p.treated,p.control,log2.RR,log2.OR,res[,1],q0,padj)
  colnames(res) <- c("p.treated","p.control","log2.RR","log2.OR","pvalue","q","padj")
    #res <- res2[,c(1,3:7)]
  #path=getwd()
  path <- paste(output.dir,"dif_meth.xls",sep = '/')
  #path=paste(path,"dif_meth.xls",sep="/")
  write.table(res,path,sep='\t',col.names=TRUE,row.names =FALSE)
  return(res)}
