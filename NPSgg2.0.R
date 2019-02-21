#probing estimating Sgg non parametrically
loadlib=function(){
  library(foreach)
  library(doParallel) #parallel processing backend for Windows & unix
  library(smoothmest) #contains double exponential random generation 
  library(stringr)
  library(pracma)#
  library(plyr)
  library(dfoptim)
  library(mvtnorm)
  library(lcmix)
  library(MCMCpack)
}

loadlib()

set.cores = function(cr.want){
  #returns a number of cores  based on cr.want
  # does not make a cluster with more cores than detected on machine
  # if cr.want is a fraction, returns cluster with cores = trunc(cr.want*maxcore)
  # If cr.want is a negagive number returns cores = maxcore less cr.want
  crmx=detectCores()
  if ((cr.want<1)&( cr.want>0)) {crs=cr.want*crmx}
  else {crs=cr.want}
  if (cr.want>=crmx) {crs=crmx}
  if (cr.want<0) {crs=crmx+cores}
  if (cr.want < 1){crs=1}   
  return(crs)
}


cores=4
crs=set.cores(cores)
cl=makeCluster(crs)
registerDoParallel(cl,cores=crs)

balsmpl=function(B,sz){
  #Returns a BXsz matrix of indices for a balanced bootstrap.  
  #All indicies occur = number of times in matrix  
  #Balanced bootstraps are more efficient than regular bootstraps
  #Use each row of the matrix to select a bootstrapped sample
    idxes=sample(rep(1:sz,times=B))
  out=matrix(idxes,nrow=B)
  return=out
} 

NParHest <- function(X,method="KL"){
  #method="KL" uses H.KL
  #method="SFG" use H.SFG
  switch(method,
         KL = {out <- H.KL(X)},
         SFG = {out <- H.SFG(X)},
         stop("Unknown method in NParHest")
         )
  return(out)
}
  
  

H.SFG <- function(X){
  #non parametric H estimate as described by
  #Santamaria-Bonfil, G., N. Fernandez, and C. Gershenson. 2016. 
  #Measuring the Complexity of Continuous Distributions. Entropy 18.
  ss=length(X)
  del <- (max(X)-min(X))/ss #Bin width
  brks=seq(from=min(X),to=max(X),by=del)
  q=hist(x = X,breaks = brks,plot = F)
  p=q$counts/ss # a vector of probability estimates
  P=p[p>0] #Stripping out cells with 0 probability
  f=P/del  #converting probability masses to probability densities
  out <- -sum(del*f*log(f)) #calculating H as in equation 8 & 9
  return(out)
}

H.KL <- function(X,k=1,condOnY=NULL){
  # Kozachenko-Leonenko estimator of entropy as given by eq 1 of 
  # Bekrrett, Samworth & Yuan 2016 (arXiv:1606.00304v2)
  X=as.matrix(X)
  if (is.null(condOnY)) {
    Y <- X
  } else Y <- as.matrix(condOnY)
  if (dim(X)[2] != dim(Y)[2]) stop("dimenstions not equal for X and conditioning empirical dist")
  N <- dim(X)[1]
  d <- dim(X)[2]
  Vd <- pi^(d/2)/gamma(1 + (d/2))
  offst <- log((N-1)*(pi^(d/2))) - lgamma(1 + (d/2)) - digamma(k)
  Kdist <- function(dvec,k){sort(dvec[dvec>0])[k]} #finds kth distance greater than 0
  Dmat <- pracma::pdist2(X = X,Y=Y)
  rho <- plyr::aaply(.data = Dmat,.margins = 1,.fun = Kdist,k=k) # rho = vec of kth distances
  H <- offst + sum(log(rho))/N
  return(H)
}

H.wKL.loop <- function(X,k=0){
  # weighted Kozachenko-Leonenko estimator of entropy as given on page 3 of 
  # Bekrrett, Samworth & Yuan 2016 (arXiv:1606.00304v3)
  # initial implementation using for loops. Retained to check vector based calculations 
  X=as.matrix(X)
  Y <- X
  if (dim(X)[2] != dim(Y)[2]) stop("dimenstions not equal for X and conditioning empirical dist")
  wkdn <- Find.wkdN(X,k)
  w=wkdn$w # w is the weight
  k <- wkdn$k # is neighborhood size i.e. kth nearest neighbor
  d <- wkdn$d # d is demension of an observation 
  N <- wkdn$N  # N is the sample size
  lVd <- log((pi^(d/2)))-lgamma(1 + (d/2)) # log volume of a unit d dimensional ball
  offst <- log(N-1)+lVd
  lRHOd <- matrix(NA,nrow=N,ncol=k) #Matrix for nearest neighbor distances 1 - k
  PSI <- rep(NA,k) #vector of digamma of 1-k
  Dmat <- pracma::pdist2(X = X,Y=Y)
  Kdist <- function(dvec,k,dmin){sort(dvec[dvec>dmin])[k]} #finds kth distance greater than 0
  dmin <- max(diag(Dmat)) # dmin used instead of 0 because of round off errors
  for (j in 1:k) {
    PSI[j] <- digamma(j)
    lRHOd[ ,j] <-  plyr::aaply(.data = Dmat,.margins = 1,.fun = Kdist,k=j,dmin=dmin) # rho = vec of kth distances
  }
  lRHOd <- log(lRHOd^(d))
  H.w <- 0
  for (i in 1:N){  #these summation calculate the estimate as given on page 3 of BS&Y 2016
    for (j in 1:k){
      H.w <- H.w + w[j]*(-PSI[j]+offst + lRHOd[i,j])
    }
  }
  H.w <- H.w/N
  return(H.w)
}


Hse.wKL <- function(X,k=0, runprlll=F,progbr="none"){
  # weighted Kozachenko-Leonenko estimator of entropy as given on page 3 of 
  # Berrett, Samworth & Yuan 2016 (arXiv:1606.00304v3)
  # standard error estimate given on page 7
  # runprlll is usefull for large data sets.  Parallel processing backend needs to be registered.
  # with 4 cores parallel processing is faster with sample size >300
  # progbr ="text" can keep you from going crazy with large data sets.  will slow down small data sets
  # X is a matrix of multivariate observtions with rows being observations
  # k is the max nearest neighbor distance to use in weighting.  
  # if k==0, then k will be automatically calculated as floor(N^.3333)
  # function returns a list of the weighted H estimate and its standard error
  #   out <- list(Hw = Hw, Hwse = Hwse)
  X=as.matrix(X)
  wkdn <- Find.wkdN(X,k)
  w=as.matrix(wkdn$w) # w is the weight vector
  k <- wkdn$k # is max neighborhood size i.e. kth nearest neighbor
  d <- wkdn$d # d is demension of an observation 
  N <- wkdn$N  # N is the sample size
  lVd <- log((pi^(d/2)))-lgamma(1 + (d/2)) # log volume of a unit d dimensional ball
  offst <- log(N-1)+lVd
  Kdist <- function(dvec,k,dmin){sort(dvec[dvec>dmin])[k]} #finds kth distance greater than 0
  # dmin used instead of 0 because of round off errors
  PSI <- digamma(1:k) # vector of digamma of 1:k. digamma is a vectorized function
  H.wi <- function(xrow, Xmat){ # function to apply returns the H contribution of an observation
    Dveci <- pracma::pdist2(X = Xmat,Y=xrow) # vector distances to an observation
    dmin <- min(Dveci) #dimin is used because roundoff error in distance calculation sometimes >0
    lRHOdi <-  log(Kdist(dvec=Dveci,k = 1:k,dmin=dmin)^d) # rho = vec of kth distances
    logxi <- (-PSI+offst + lRHOdi) # xi defined on page 3 just below equation for Hwn
    Hwi <- logxi%*%w
    H2wi <- (logxi^2)%*%w #see page 7 second momment
    return(c(Hwi,H2wi))
  }
  HVwi <- plyr::aaply(.data = X,.margins = 1,.fun = H.wi,Xmat=X,.parallel=runprlll,.progress=progbr) # contributions of each observation
  HVw <- plyr::aaply(.data = HVwi,.margins = 2,.fun = mean) # averaged contributions
  Hw <- HVw[1]
  VH <- (HVw[2]-HVw[1]^2)/N
  if (VH > 0) {Hwse <- sqrt(VH)
  } else Hwse=0
  out <- list(Hw = Hw, Hwse = Hwse)
  return(out)
}


Find.wkdN <- function(X,k) {
  # email communication from Berrett.  These equations work up through d=7
  # \[\begin{align}
  # & g1=\sum\limits_{j=1}^{k}{\Gamma }(j+2/d)/\Gamma (j) \\ 
  # & g2=\sum\limits_{j=1}^{k}{\Gamma }{{(j+2/d)}^{2}}/\Gamma {{(j)}^{2}} \\ 
  # & {{w}_{j}}=\left( g2-g1\Gamma (j+2/d)/\Gamma (j) \right)/(kg2-g{{1}^{2}}) \\ 
  # \end{align}\] 

  X=as.matrix(X)
  if (dim(X)[1]<dim(X)[2]) { #checking that observations are rows
    X <- t(X)
  }
  numobs <- dim(X)[1]
  if (k==0) {  k <- floor(numobs^(1/3))}
  d <- dim(X)[2]
  g1 <-0
  g2 <- 0
  w <- rep(NA,k)
  for (j in 1:k) {
    g1 <- g1 + (gamma(j + 2/d)/gamma(j))
    g2 <- g2 + (gamma(j + 2/d)/gamma(j))^2
  }
  for (j in 1:k){
    w[j] <- (g2 - (g1*gamma(j + 2/d)/gamma(j)))
  }
  w <- w/(k*g2 -g1^2)
  out <- list(w=w,k=k,d=d,N=numobs)
  return(out)
}

mu=rep(10,7)
sigma=diag(rep(1,7))
X=rmvnorm(n =300,mean = mu,sigma = sigma)



bestk=rep(NA,10)
for (k in 1:10) {
  bestk[k] <- H.wKL(X,k=k,runprlll = T) # change this function call
}
H.true <- H.MVnorm(par.lst = list(mu=mu,sigma=sigma))
H.true
bestk
ke=Find.wkdN(X,k=0)$k 
ke
bestk[ke]/H.true
max(bestk,na.rm = T)/H.true

btd4.50 =matrix(rep(NA,500),ncol=2)
btd4.25 =matrix(rep(NA,500),ncol=2)
btd7.25 =matrix(rep(NA,100000),ncol=2)
btd7.10 = matrix(rep(NA,1000),ncol=1)
X=rmvnorm(n = 10,mean = mu,sigma = sigma)

ke=Find.wkdN(X,k=0)$k
pt1=proc.time()
for (i in 1:1000){
  X=rmvnorm(n = 10,mean = mu,sigma = sigma)
    btd7.10[i] <- H.wKL(X,k=0,runprlll = F)
}
pt2=proc.time()
pt2-pt1

fesH <- function(X,k){
  out=list(H=H.KL(X=X,k=k,condOnY = X),CO=X)
  return(out)
}

fevH <- function(X,k,Y){
  out=list(H=H.KL(X=X,k=k,condOnY = Y),CO=Y)
  return(out)
}



KK_H.KL=function(X,fesH,fevH,B,M,k){
  X=as.matrix(X)
  Xlen = dim(as.matrix(X))[1]
  est0=fesH(X,k=k)
  theta0=est0$CO
  H00=est0$H
  HB=matrix(NA,ncol=7,nrow=B,dimnames = list(NULL,c("H00","H01","H10","H11","H12","H21","H22")))
  Bidx=balsmpl(B,sz=Xlen)
  for (i in 1:B) {
    Xi=as.data.frame(X[Bidx[i,],])
    esti=fesH(Xi,k=k)
    thetai=esti$CO
    H11=esti$H
    H10=fevH(Xi,k=k,theta0)$H
    H01=fevH(X,k=k,thetai)$H
    if (M == 0) {H12=NA; H21 = NA; H22 = NA}
    else {
      if (M == 1) Midx=matrix(sample(1:Xlen,size = Xlen,replace = T),nrow=1)
      else Midx = balsmpl(M,sz=Xlen)
      H12 <- H21 <- H22 <- 0
      for (j in 1:M){
        Xib=as.data.frame(Xi[Midx[j, ], ])
        estib=fesH(Xib,k=k)
        thetaib=estib$CO
        H22=fevH(Xib,k=k,thetaib)$H + H22
        H21=fevH(Xib,k=k,thetai)$H + H21
        H12=fevH(Xi,k=k,thetaib)$H + H12
      }
      H22 = H22/M
      H21 = H21/M
      H12 = H12/M
    }  
    rowi=c(H00,H01,H10,H11,H12,H21,H22)
    HB[i,]=rowi
  }
  b1b=with(as.data.frame(HB),mean(H11-H10+H00-H01))
  b2b=2*b1b-with(as.data.frame(B),mean(H22-H21+H11-H12))
  NPH <- H00
  NPH.bc <- H00 - b1b
  NPH.bc2 <- H00 - b2b
  ests=list(NPH=NPH,NPH.bc=NPH.bc, NPH.bc2=NPH.bc2)
  out=list(ests=ests,X=X,B=B,M=M,k=k)
  return(out)
}

pt1=proc.time()
tst=KK_H.KL(X=X,fesH=fesH,fevH=fevH,B=100,M=100,k=1)
pt2=proc.time()

pt2-pt1

tst$ests
tst$ests
BtStrp <- function(X,fn,B){
  #Simple function for balanced univariate bootstrap
  #returns estimate, bootstrap bias estimate, bootstrap bias corrected estimate
  # and mean of bootstrapped estimates
  ss <- length(X)
  thetai <- rep(NA,B)
  IDXes <- balsmpl(B = B, sz = ss)
  for (i in 1:B) {
    thetai[i] <- fn(X[IDXes[i,]])
  }
  est <- fn(X)
  bias <- mean(thetai) - est
  est.bc <-  est - bias
  theta.bar <- mean(thetai)
  out <- list(est=est, bias=bias, est.bc=est.bc, theta.bar=theta.bar)
}

# BcNPH <- function(X,Bin) {
#   Bt <- BtStrp(X,NParHest,Bin)
#   out <- list(npH=Bt$est, npH.bc=Bt$est.bc)
#   return(out)
# }

BcNPH <- function(X,Bin) {
  Bt <- DblBtStrp(X,NParHest,Bin,Bin)
  out <- list(npH=Bt$est1, npH.bc=Bt$est2, npH.bc2=Bt$est3)
  return(out)
}

DblBtStrp <- function(X,fn,B1,B2){
  #Simple double bootstrp function
  # Based on MA Martin 1990. Standford U. Dept. of Statistics Tech report 347
  # https://statistics.stanford.edu/sites/default/files/EFS%20NSF%20347.pdf
  ss <- length(X)
  theta1 <- rep(fn(X),B1)
  theta2 <- rep(NA,B1)
  theta3 <- rep(NA,B1)
  IDXes <- balsmpl(B = B1, sz = ss)
    
 for (i in 1:B1){
   bt <- BtStrp(X[IDXes[i,]],fn,B2)
   theta2[i] <- bt$est
   theta3[i] <- bt$theta.bar
 } 
  est1 <- mean(theta1)
  est2 <- 2*est1 - mean(theta2)
  est3 <- 3*est1 - 3*mean(theta2) + mean(theta3)
  out <- list(est1=est1,est2=est2,est3=est3)
  return(out)
}

r.lnorm <- function(n,par.lst) {
  arglst=c(list(n=n),par.lst)
  out=do.call(rlnorm,arglst)
  return(out)
}

MLE.lnorm <- function(X){
  ss=length(X)
  x=log(X)
  meanlog <- sum(x)/ss
  sdlog <- sum((x-meanlog)^2)/ss
  MLEs <- list(meanlog=meanlog, sdlog=sdlog)
  return(MLEs)
}

H.lnorm <- function(par.lst=list(meanlog=1, sdlog=1)) {
  meanlog=par.lst$meanlog
  sdlog=par.lst$sdlog
  H <- meanlog + log(2*pi*exp(1)*sdlog^2)/2
  return(H)
}

r.mvlnorm <- function(n,par.lst) {
  #par.lst = list(meanlog,sigmalog)
  prlst <- list(mean=par.lst$meanlog, sigma=par.lst$sigmalog)
  arglst=c(list(n=n),prlst)
  out=do.call(rmvnorm,arglst)
  out=exp(out)
  return(out)
}

MLE.mvlnorm <- function(X){
  lX=log(X)
  meanlog <- plyr::aaply(.data=lX,.fun=mean,.margins=2)
  sigmalog <- cov(lX)
  MLEs <- list(meanlog=meanlog, sigmalog=sigmalog)
  return(MLEs)
}

H.mvlnorm <- function(par.lst=list(meanlog=rep(0,3), sigmalog=diag(rep(1,3)))) {
  cat(par.lst$sigmalog,"\n")
  meanlog=par.lst$meanlog
  sigmalog <- as.matrix(par.lst$sigmalog)

  d=length(meanlog)
  H <- sum(meanlog) + (d/2)*(1+log(2*pi)) + log(det(sigmalog))/2
  return(H)
}

r.dirchlet <- function(n,par.lst){
  #requires MCMCpack loaded
  alpha=par.lst$alpha
  X=rdirichlet(n=n,alpha=alpha)
  return(X)
}
MLE.dirchlet <- function(X) {}
H.dirchlet <- function(par.lst=list(alpha=c(1:4))) {
  alpha <- par.lst$alpha
  K <- length(alpha)
  a0 <- sum(alpha)
  B <- function(alpha){prod(gamma(alpha))/gamma(sum(alpha))}
  H <- log(B(alpha))+(a0-K)*digamma(a0) - sum((alpha-1)*digamma(alpha))
  return(H)
}

r.norm <- function(n,par.lst=list(mean=0, sd=1)){
  arglst=c(list(n=n),par.lst)
  out=do.call(rnorm,arglst)
  return(out)
}

MLE.norm <- function(X){
  ss=length(X)
  mean <- sum(X)/ss
  sd <- sqrt(sum((X-mean)^2)/ss)
  MLEs <- list(mean=mean, sd=sd)
  return(MLEs)
}


H.norm <- function(par.lst=list(mean=0, sd=1)){
  mean=par.lst$mean
  sd <- par.lst$sd
  H=log(sd*sqrt(2*pi*exp(1)))
  return(H)
}

H.MVnorm <- function(par.lst=list(mu=c(0,0),sigma=matrix(c(1,0,0,1),nrow=2))) {
  sigma=par.lst$sigma
  d=dim(sigma)[1]
  out <- log(det(2*pi*exp(1)*sigma))/2
  return(out)
}

r.exp <- function(n,par.lst=list(lambda)){
  arglst=c(list(n=n),rate=par.lst$lambda)
  out=do.call(rexp,arglst)
  return(out)
}

MLE.exp <- function(X){
  lambda <- 1/mean(X)  
  MLEs <- list(lambda=lambda)
  return(MLEs)
}


H.exp <- function(par.lst=list(lambda=1)){
  lambda <- par.lst$lambda
  H <- 1 - log(lambda)
  return(H)
}


r.Laplace <- function(n, par.lst=list(mu=0, lambda=1)){
  out <- do.call(smoothmest::rdoublex,c(list(n=n),par.lst))
  return(out)
}

MLE.Laplace <- function(X){
  N <- length(X)
  mu <- median(X)
  lambda <- sum(abs(X-mu))/N
  par.lst <- list(mu=mu, lambda=lambda)
  return(par.lst)
}

H.Laplace <- function(par.lst=list(mu=0, sd,lambda=1)){
  mu=par.lst$mu
  lambda=par.lst$lambda
  H=log(2*lambda*exp(1))
  return(H)
}

r.normal_mixture <- function(n, par.lst=list(mu=c(0,0),sd=c(1,1), mixP=c(0.5,0.5))){
  arglst=list(n=n,mu=par.lst$mu,sigma=par.lst$sd,lambda=par.lst$mixP)
  out=do.call(mixtools::rnormmix,arglst)
  return(out)
}

d.normal_mixture <- function(x,par.lst=list(mu=c(0,0),sd=c(1,1), mixP=c(0.5,0.5))){
  mu <- par.lst$mu
  sd <- par.lst$sd
  mixP <- par.lst$mixP
  d1 <- dnorm(x = x,mean = mu[1],sd = sd[1])
  d2 <- dnorm(x = x,mean = mu[2],sd = sd[2])
  d <- mixP[1]*d1 + mixP[2]*d2
  return(d)
}

H.normal_mixture <- function(par.lst=list(mu=c(0,0),sd=c(1,1), mixP=c(0.5,0.5))){
  mu <- par.lst$mu
  sd <- par.lst$sd
  mixP <- par.lst$mixP
  h <- function(x){
    out <- d.normal_mixture(x,par.lst)*log(d.normal_mixture(x,par.lst))
    return(out)
  }
  H <- -integrate(f = h,lower = min(mu)-10*max(sd),upper = max(mu)+10*max(sd))$value
  return(H)
}



MLE.normal_mixture <- function(X){
  lambda <- 0.5 #guess
  mu <- c(-1,1) #guess
  sigsqrd <- c(1,1) #guess
  value <- mixtools::normalmixEM2comp(x=X,lambda=lambda,mu=mu, sigsqrd = sigsqrd)
  par.lst <- list(mu=value$mu, sd=value$sigma,mixP=value$lambda)
  return(par.lst) # does not return any of the fit statistics
}

r.multinom <-function(n,par.lst){
  # par.lst=list(size,prob)
  # each row is a random observation vector 
  # note that the output is the transpose of rmultinom
    arglst=c(list(n=n),par.lst)
    out=t(do.call(rmultinom,arglst))
    return(out)
}

AllSame <- function(vec,tol=1e-8){
  # Fast way to check that all elements of a vector are the same to within
  # a specified numerical tolerance
  var(vec)<tol}

MLE.multinom <- function(X){
  #X is matrix of discrete MV observations, rows are observations
  # Uses package plyr
  #https://onlinecourses.science.psu.edu/stat504/node/45
  # the value is a list of (size, prob) where prob is an estimated vector of probabilities
  # output list can be used as par.lst in H.multinom/H.multinom.loop
  AllSame <- function(vec,tol=1e-8){var(vec)<tol} #fast testing of equality
  if (! AllSame(plyr::aaply(.data=X,.margins=1,.fun=sum))) stop("observations not all same size")
  prob <- plyr::aaply(.data = X,.margins = 2,.fun = mean)
  prob <- prob/sum(prob)
  size <- sum(X[1,])
  out <- list(size=size, prob=prob)
  return(out)
}


H.multinom.loop <- function(par.lst=list(size=15,prob=c(.1,.3,.6))) {
  # https://stats.stackexchange.com/questions/207893/entropy-of-the-multinomial-distribution/208075#208075
  # function will tolerate 0s in prob vector
  n <- par.lst$size
  p <- par.lst$prob
  k <- length(p)
  piece1 <- -lfactorial(n)
  piece2 <- -n*sum(p*log(p),na.rm = T) #0*log(0) is NaN should be 0 so na.rm
  piece3 <-0
  for (j in 1: 1:k){
    for (xi in 0:n) {
      piece3 <- piece3 + dbinom(x = xi,size = n,prob = p[j])*lfactorial(xi)
    }
  }
  H <- piece1+piece2+piece3
  out <- list(H=H)
  return(out)  
} 



H.multinom <- function(par.lst=list(size=15,prob=c(.1,.3,.6))) {
  # https://stats.stackexchange.com/questions/207893/entropy-of-the-multinomial-distribution/208075#208075
  # H.multinom is a vectorized version of H.multinom.loop.  H.multinom returns exactly
  # the same values as H.multinom.loop, but more quickly.  However the logic of the loop
  # function is easier to follow.
  # function will tolerate 0s in prob vector
  n <- par.lst$size
  p <- par.lst$prob
  k <- length(p)
  piece1 <- -lfactorial(n)
  piece2 <- -n*sum(p*log(p),na.rm = T)  #0*log(0) is NaN should be 0 so na.rm
  piece3 <-0
  xi <- 0:n
  ncxi <- choose(n,xi)
  lfxi <- lfactorial(xi)
  nclfxi <- matrix(rep(ncxi,k),ncol=k,byrow=F)*matrix(rep(lfxi,k),ncol=k,byrow=F)
  PJ <- matrix(rep(p,n+1),ncol=k,byrow=T)
  XI <- matrix(rep(xi,k),ncol=k,byrow=F)
  piece3=sum(nclfxi*(PJ^XI)*((1-PJ)^(n-XI)))
  H <- piece1+piece2+piece3
  out <- list(H=H)
  return(out)  
} 



H4dist <- function(n,par.lst,Bin,distname){
  r.dist <- get(paste("r",distname,sep="."))
  MLE.dist <- get(paste("MLE",distname,sep="."))
  H.dist <- get(paste("H",distname,sep="."))
  X <- r.dist(n=n,par.lst=par.lst)
  trH <- H.dist(par.lst=par.lst)
  mlH <- H.dist(par.lst=MLE.dist(X))
  bcnph <- BcNPH(X = X,Bin=Bin)
  npH <- bcnph$npH
  npH.bc <- bcnph$npH.bc
  npH.bc2 <- bcnph$npH.bc2
  nrmHbnd <- H.norm(par.lst=list(mean=0,sd=sd(X)))
  out <- c(trH,mlH,npH,npH.bc,npH.bc2,nrmHbnd)
  names(out) <- c("trH","mlH","npH","npH.bc","npH.bc2","nrmHbnd")
  return(out)
}

BigBt4dist <- function(nvec,Bout,Bin,distname,par.lst,talk=T){
  if (talk) cat("Working on distribution:",distname,"\n")
  out <- list()
  exvec <- c(paste(c("MLE","r","H"),distname,sep="."),"BcNPH","NParHest","H4dist","DblBtStrp","balsmpl","BtStrp","H.norm","H.KL","d.normal_mixture")
  for (i in 1:length(nvec)){
    if (talk) cat("Working on i=",i,"element of nvec=",nvec[i],"\n")
    n=nvec[i]
    Btmat <- foreach(k=1:Bout,.combine=rbind, .inorder = F,.export = exvec) %dopar% {
      btH <- H4dist(n=n,par.lst=par.lst,Bin=Bin,distname=distname)
      return(btH)
    }
    BigBt <- list(Btmat = Btmat, n=n, Bout=Bout,Bin=Bin, distname=distname, par.lst=par.lst)
    Btname <- paste(paste("BBt",distname,sep="."),paste('n',n,sep=""),sep=".")
    tmplst <- list(BigBt)
    names(tmplst) <- Btname
    out <- c(out,tmplst)
  }
  return(out)
}

H4dist.KK <- function(n,par.lst,B,M,k=1,distname,fesH,fevH){
  r.dist <- get(paste("r",distname,sep="."))
  MLE.dist <- get(paste("MLE",distname,sep="."))
  H.dist <- get(paste("H",distname,sep="."))
  X <- r.dist(n=n,par.lst=par.lst)
  trH <- H.dist(par.lst=par.lst)
  mlH <- H.dist(par.lst=MLE.dist(X))
  bcnph <- KK_H.KL(X = X,fesH = fesH,fevH = fevH,B=B,M=M,k=k)
  npH <- bcnph$ests$NPH
  npH.bc <- bcnph$ests$NPH.bc
  npH.bc2 <- bcnph$ests$NPH.bc2
  nrmHbnd <- H.norm(par.lst=list(mean=0,sd=sd(X)))
  out <- c(trH,mlH,npH,npH.bc,npH.bc2,nrmHbnd)
  names(out) <- c("trH","mlH","npH","npH.bc","npH.bc2","nrmHbnd")
  return(out)
}

BigBt4dist.KK <- function(nvec,Bout,Bin,distname,par.lst,talk=T){
  if (talk) cat("Working on distribution:",distname,"\n")
  out <- list()
  exvec <- c(paste(c("MLE","r","H"),distname,sep="."),"BcNPH","NParHest","H4dist.KK",
             "KK_H.KL","balsmpl","H.norm","H.KL","d.normal_mixture","fesH","fevH")
  for (i in 1:length(nvec)){
    if (talk) cat("Working on i=",i,"element of nvec=",nvec[i],"\n")
    n=nvec[i]
    Btmat <- foreach(k=1:Bout,.combine=rbind, .inorder = F,.export = exvec) %dopar% {
                btH <- H4dist.KK(n=n,par.lst=par.lst,B=Bin,M=10,distname=distname,fesH=fesH,fevH = fevH,k=1)
                return(btH)
             }
    BigBt <- list(Btmat = Btmat, n=n, Bout=Bout,Bin=Bin, distname=distname, par.lst=par.lst)
    Btname <- paste(paste("BBt",distname,sep="."),paste('n',n,sep=""),sep=".")
    tmplst <- list(BigBt)
    names(tmplst) <- Btname
    out <- c(out,tmplst)
  }
  return(out)
}



stopCluster(cl)
cores=4
crs=set.cores(cores)
cl=makeCluster(crs)
registerDoParallel(cl)

distname <- "norm"
#nvec <- c(25,50,200,800,1600)
nvec <- c(25,50,200,800)

par.lst=list(mean=0, sd=1)
pt1=proc.time()
Bt.lst.norm <- BigBt4dist(nvec=nvec,Bout=100,Bin=10,distname=distname,par.lst=par.lst)
pt2=proc.time()
pt2-pt1

distname <- "exp"
par.lst=list(lambda=1)
nvec <- c(25,50,200,800)
Bt.lst.exp <- BigBt4dist(nvec=nvec,Bout=100,Bin=10,distname=distname,par.lst=par.lst)

distname <- "lnorm"
par.lst=list(meanlog=1, sdlog=1)
nvec <- c(25,50,200,800)
Bt.lst.lnorm <- BigBt4dist(nvec=nvec,Bout=100,Bin=10,distname=distname,par.lst=par.lst)


distname <- "laplace"
par.lst <- list(mu=0, lambda=1)
nvec <- c(25,50,200,800)

Bt.lst.laplace <- BigBt4dist(nvec=nvec,Bout=100,Bin=10,distname=distname,par.lst=par.lst)



distname <- "normal_mixture"
par.lst=list(mu=c(-1,1),sd=c(.25,.25),mixP=c(.5,.5))
nvec <- c(25,50,200,800)
Bt.lst.normal_mixture <- BigBt4dist(nvec=nvec,Bout=100,Bin=1,distname=distname,par.lst=par.lst)

distname <- "normal_mixture"
par.lst=list(mu=c(-1,1),sd=c(.25,.25),mixP=c(.5,.5))
nvec <- c(25,75,225)
pt1=proc.time()
Bt.lst.normal_mixture.KK <- BigBt4dist.KK(nvec=nvec,Bout=100,Bin=100,distname=distname,par.lst=par.lst)
pt2=proc.time()

nvec <- c(25)
pt1=proc.time()
Bt.lst.normal_mixture.KK.B1000 <- BigBt4dist.KK(nvec=nvec,Bout=100,Bin=1000,distname=distname,par.lst=par.lst)
pt2=proc.time()
pt2-pt1



tst=H4dist.KK(n=nvec[1],par.lst,B=2,M,k=1,distname,fesH,fevH)

distname <- "normal_mixture"
par.lst=list(mu=c(-1,1),sd=c(1,1),mixP=c(.5,.5))
nvec <- c(50,1600)

Bt.lst.normal_mixture.2 <- BigBt4dist(nvec=nvec,Bout=100,Bin=100,distname=distname,par.lst=par.lst)
PltBt(Btmat =  Bt.lst.normal_mixture.2[[1]], distname = distname, n = 50,par.lst = par.lst)
PltBt(Btmat =  Bt.lst.normal_mixture.2[[2]], distname = distname, n = 1600,par.lst = par.lst)

stopCluster(cl)

######
cores=4
crs=set.cores(cores)
cl=makeCluster(crs)
registerDoParallel(cl)

exvec=c("MLE.lnorm","r.lnorm","H.lnorm","H.norm","BcNPH","NParHest")
par.lst=list(meanlog=1,sdlog=1)
n=50
Bin=100
Bout=100
distname="lnorm"
pt1=proc.time()
BigBt.50 <- foreach(k=1:Bout,.combine=rbind, .inorder = F,.export = exvec) %dopar% {
  btH <- H4dist(n=n,par.lst=par.lst,Bin=Bin,distname=distname)
  return(btH)
}
pt2=proc.time()
stopCluster(cl)

time.50=pt2-pt1

cores=4
crs=set.cores(cores)
cl=makeCluster(crs)
registerDoParallel(cl)

exvec=c("MLE.lnorm","r.lnorm","H.lnorm","H.norm","BcNPH","NParHest")
par.lst=list(meanlog=1,sdlog=1)
n=200
Bin=100
Bout=100
distname="lnorm"
pt1=proc.time()
BigBt.200 <- foreach(k=1:Bout,.combine=rbind, .inorder = F,.export = exvec) %dopar% {
  btH <- H4dist(n=n,par.lst=par.lst,Bin=Bin,distname=distname)
  return(btH)
}
pt2=proc.time()
stopCluster(cl)
time.200=pt2-pt1

cores=4
crs=set.cores(cores)
cl=makeCluster(crs)
registerDoParallel(cl)

exvec=c("MLE.lnorm","r.lnorm","H.lnorm","H.norm","BcNPH","NParHest")
par.lst=list(meanlog=1,sdlog=1)
n=800
Bin=100
Bout=100
distname="lnorm"
pt1=proc.time()
BigBt.800 <- foreach(k=1:Bout,.combine=rbind, .inorder = F,.export = exvec) %dopar% {
  btH <- H4dist(n=n,par.lst=par.lst,Bin=Bin,distname=distname)
  return(btH)
}
pt2=proc.time()
stopCluster(cl)
time.800=pt2-pt1

cores=4
crs=set.cores(cores)
cl=makeCluster(crs)
registerDoParallel(cl)

exvec=c("MLE.lnorm","r.lnorm","H.lnorm","H.norm","BcNPH","NParHest")
par.lst=list(meanlog=1,sdlog=1)
n=1600
Bin=100
Bout=100
distname="lnorm"
pt1=proc.time()
BigBt.1600 <- foreach(k=1:Bout,.combine=rbind, .inorder = F,.export = exvec) %dopar% {
  btH <- H4dist(n=n,par.lst=par.lst,Bin=Bin,distname=distname)
  return(btH)
}
pt2=proc.time()
stopCluster(cl)
time.1600=pt2-pt1

pdf(file = "FirstHsim.pdf")
boxplot(x=list(BigBt.50[,"trH"],BigBt.50[,"mlH"],BigBt.50[,3],BigBt.50[,4],BigBt.50[,5],BigBt.50[,6]),names = colnames(BigBt.50),ylim=c(0,4))
title(main="lnorm(1,1), n=50")
boxplot(x=list(BigBt.200[,"trH"],BigBt.200[,"mlH"],BigBt.200[,3],BigBt.200[,4],BigBt.200[,5],BigBt.200[,6]),names = colnames(BigBt.200),ylim=c(0,4))
title(main="lnorm(1,1), n=200")
boxplot(x=list(BigBt.800[,"trH"],BigBt.800[,"mlH"],BigBt.800[,3],BigBt.800[,4],BigBt.800[,5],BigBt.800[,6]),names = colnames(BigBt.800),ylim=c(0,4))
title(main="lnorm(1,1), n=800")
boxplot(x=list(BigBt.1600[,"trH"],BigBt.1600[,"mlH"],BigBt.1600[,3],BigBt.1600[,4],BigBt.1600[,5],BigBt.1600[,6]),names = colnames(BigBt.1600),ylim=c(0,4))
title(main="lnorm(1,1), n=1600")
dev.off()

PltBt <- function(BBt.obj,ylim=NULL) {
  Btmat <- BBt.obj$Btmat
  n <- BBt.obj$n
  distname <- BBt.obj$distname
  par.lst <- BBt.obj$par.lst
  boxplot(x=list(Btmat[,"trH"],Btmat[,"mlH"],Btmat[,3],Btmat[,4],Btmat[,5],Btmat[,6]),names = colnames(Btmat),ylim=ylim)
  title(main=str_c(distname,"(",str_c(names(par.lst),par.lst,sep="=",collapse = ", "),") n=",n))
}

pltbt <- function(BBt.obj,ylim=NULL) {
  Btmat <- BBt.obj$Btmat
  Sgghat <- -Btmat
  cnames=c("true Sgg","ML Sgg","non-par\nSgg","normal\nbound")
  n <- BBt.obj$n
  distname <- BBt.obj$distname
  par.lst <- BBt.obj$par.lst
  boxplot(x=list(Sgghat[,"trH"],Sgghat[,"mlH"],Sgghat[,3],Sgghat[,6]),names = cnames,ylim=ylim)
  title(main=str_c(distname,"(",str_c(names(par.lst),par.lst,sep="=",collapse = ", "),") n=",n))
}


lst.lst <- list(Bt.lst.lnorm,Bt.lst.laplace,Bt.lst.exp,Bt.lst.normal_mixture)
pdf(file="First.H.KL.sim.pdf")
for (i in 1:length(lst.lst)) {
  bt.ob <- lst.lst[[i]]
  for (j in 1:4) {
    PltBt(bt.ob[[j]])
  }
}
dev.off()

lst.lst <- list(Bt.lst.norm,Bt.lst.lnorm,Bt.lst.laplace,Bt.lst.exp,Bt.lst.normal_mixture)
tiff(filename = "Sgg.norm.tiff")
  plot.new()
  par(mfrow=c(1,2))
  bt.ob <- lst.lst[[5]]
  ymin <- -max(bt.ob[[1]]$Btmat[,-c(4,5)],bt.ob[[4]]$Btmat[,-c(4,5)])
  ymax <- -min(bt.ob[[1]]$Btmat[,-c(4,5)],bt.ob[[4]]$Btmat[,-c(4,5)])
  for (j in c(1,4)) {
    pltbt(bt.ob[[j]],ylim=c(ymin,ymax))
  }

dev.off()





PltBt(Bt.lst.laplace[[5]],distname=distname,n=1600,par.lst=par.lst)

PltBt(Bt.lst.normal_mixture[[5]],distname=distname,n=1600,par.lst=par.lst)

PltBt(Bt.lst.laplace[[1]])

mu=0
sd=2
ss=100
X=rnorm(n = ss,mean = mu,sd=sd) # data simulated from a Gaussian
e=exp(1)

sdo=sd(X) #estimated sd

H.tr <- log(sd^2*2*pi*e)/2
H.tr
H.sd.est <- log(sdo^2*2*pi*e)/2
H.sd.est
H.NP.est <- NParHest(X)
H.NP.est
H.NP.bt.bc.est <- BtStrp(X,fn=NParHest,B=1000)
H.NP.bt.bc.est
H.NP.dbbt.bc.est <- DblBtStrp(X,fn=NParHest,B1 = 1000, B2 = 1000)
H.NP.dbbt.bc.est


stopCluster(cl)
cores=4
crs=set.cores(cores)
cl=makeCluster(crs)
registerDoParallel(cl)

distname <- "norm"
#nvec <- c(25,50,200,800,1600)
nvec <- c(25)

par.lst=list(mean=0, sd=1)
pt1=proc.time()
tst.norm.Bin50 <- BigBt4dist(nvec=nvec,Bout=100,Bin=50,distname=distname,par.lst=par.lst)
pt2=proc.time()
pt2-pt1
