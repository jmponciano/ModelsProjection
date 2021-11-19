
# Figure 10
# Simulation for part (A), page 56, section 3.3.1.4 "Misspecified models"

library(smoothmest)
library(entropy)
library(MASS)

muo <- 350
sig <- 100# 15, 33, 50
b   <- sig/sqrt(2) # so that Var(samples) = 2*(b^2) from Laplace distrib

nmax <- 2003
x1 <- runif(n=nmax,min=-3, max=3)
x2 <- runif(n=nmax, min=10, max=20)
x3 <- runif(n=nmax, min=40, max=44)

beta0 <- muo
beta1 <- 2.5
beta2 <- 0.75
beta3 <- 5

beta.true <- c(beta0, beta1)
beta.alt  <- c(beta0, beta1, beta2, beta3)

p1 <- length(beta.true)
p2 <- length(beta.alt)

nvec <- seq(6,nmax,by=10)
len <- length(nvec)
B    <- 10000
pct.err.null.aic <- rep(0,len)
pct.err.null.bic <- rep(0,len)

Delta.crit <- 2

for(i in 1:len){
  
  is.f2.best.aic <- rep(0,B);
  is.f2.best.bic <- rep(0,B);
  
  ni <- nvec[i]
  
  x1.ni <- x1[1:ni]
  x2.ni <- x2[1:ni]
  x3.ni <- x3[1:ni]
  
  # Design matrices
  Xnull <- cbind(rep(1,ni), x1.ni)
  Xalt  <- cbind(rep(1,ni), x1.ni, x2.ni, x3.ni)		
  
  mean.true <- Xnull%*%beta.true
  #mean.alt  <- Xalt%*%beta.alt
  
  for(j in 1:B){		
    
    
    
    ############ Simulating closer to the null and testing
    Y   <- rdoublex(n=ni, mu=mean.true, lambda=b)
    
    beta.hat1   <- ginv(t(Xnull)%*%Xnull)%*%t(Xnull)%*%Y
    pred1       <- Xnull%*%beta.hat1
    sigsq.hat1  <- (1/ni)*t(Y-pred1)%*%(Y-pred1)
    
    beta.hat2   <- ginv(t(Xalt)%*%Xalt)%*%t(Xalt)%*%Y
    pred2       <- Xalt%*%beta.hat2
    sigsq.hat2  <- (1/ni)*t(Y-pred2)%*%(Y-pred2)
    
    m2log.L1   <- ni*(log(2*pi*sigsq.hat1) + 1)
    m2log.L2   <- ni*(log(2*pi*sigsq.hat2) + 1)		
    
    #AIC testing
    AIC1    <- m2log.L1 + 2*p1
    AIC2    <- m2log.L2 + 2*p2;
    
    is.f2.best.aic[j] <- AIC1-AIC2 >=Delta.crit; # If so, that would be saying model 2 is best!! 
    
    #BIC testing
    BIC1     <- m2log.L1 + p1*log(ni); # There are 0 parameters ESTIMATED in this case: null is fixed	
    BIC2    <- m2log.L2+ p2*log(ni);
    is.f2.best.bic[j] <- BIC1-BIC2 >=Delta.crit; # If so, that would be saying model 2 is best!!
    
  }
  pct.err.null.aic[i]  <- sum(is.f2.best.aic)/B
  pct.err.null.bic[i]  <- sum(is.f2.best.bic)/B
  
}

pct.err.null.aic.A <- pct.err.null.aic
pct.err.null.bic.A <- pct.err.null.bic



#### Simulation of part B: f1 and f2 are overlappin but f1 has a part closer to truth

muo <- 350
sig <- 100# 15, 33, 50
b   <- sig/sqrt(2) # so that Var(samples) = 2*(b^2) from Laplace distrib

nmax <- 2003
x1 <- runif(n=nmax,min=-3, max=3)
x2 <- runif(n=nmax, min=10, max=20)
x3 <- runif(n=nmax, min=40, max=44)

beta0 <- muo
beta1 <- 2.5
beta2 <- 1.75
beta3 <- 5

beta.true <- c(beta1,beta2)
beta.alt  <- c(beta0, beta2, beta3)

p1 <- length(beta.true)
p2 <- length(beta.alt)

nvec <- seq(6,nmax,by=10)
len <- length(nvec)
B    <- 10000
pct.err.null.aic <- rep(0,len)
pct.err.null.bic <- rep(0,len)

Delta.crit <- 2

for(i in 1:len){
  
  is.f2.best.aic <- rep(0,B);
  is.f2.best.bic <- rep(0,B);
  
  ni <- nvec[i]
  
  x1.ni <- x1[1:ni]
  x2.ni <- x2[1:ni]
  x3.ni <- x3[1:ni]
  
  # Design matrices
  Xnull <- cbind(x1.ni, x2.ni)
  Xalt  <- cbind(rep(1,ni), x2.ni, x3.ni)		
  mean.true <- Xnull%*%beta.true
  
  for(j in 1:B){		
    
    
    
    ############ Simulating closer to the null and testing
    Y   <- rdoublex(n=ni, mu=mean.true, lambda=b)
    
    beta.hat1   <- ginv(t(Xnull)%*%Xnull)%*%t(Xnull)%*%Y
    pred1       <- Xnull%*%beta.hat1
    sigsq.hat1  <- (1/ni)*t(Y-pred1)%*%(Y-pred1)
    
    beta.hat2   <- ginv(t(Xalt)%*%Xalt)%*%t(Xalt)%*%Y
    pred2       <- Xalt%*%beta.hat2
    sigsq.hat2  <- (1/ni)*t(Y-pred2)%*%(Y-pred2)
    
    m2log.L1   <- ni*(log(2*pi*sigsq.hat1) + 1)
    m2log.L2   <- ni*(log(2*pi*sigsq.hat2) + 1)		
    
    #AIC testing
    AIC1    <- m2log.L1 + 2*p1
    AIC2    <- m2log.L2 + 2*p2;
    
    is.f2.best.aic[j] <- AIC1-AIC2 >=Delta.crit; # If so, model 2 is best!! 
    
    #BIC testing
    BIC1     <- m2log.L1 + p1*log(ni); 
    BIC2    <- m2log.L2 + p2*log(ni);
    is.f2.best.bic[j] <- BIC1-BIC2 >=Delta.crit; # If so, model 2 is best!!
    
  }
  pct.err.null.aic[i]  <- sum(is.f2.best.aic)/B
  pct.err.null.bic[i]  <- sum(is.f2.best.bic)/B
  
}

figpath1 <- "~/Dropbox/BrianMarkSubhash/BriansMisSpecPaper/PublishedFigures&Code"
figname  <- "Figure10.eps" #OLD NAME: "M&Q-withBMapprox.tiff"
figpath <- paste0(figpath1,"/",figname)
postscript(figpath, horizontal=FALSE, paper="letter", width=8, height=8)
par(mfrow=c(2,1),mar=c(4,4,1,1), oma=c(3,3,1,1), mgp=c(2,1,0), cex.axis=1.25)

nmax <- 2003
nvec <- seq(6,nmax,by=10)
plot(nvec,pct.err.null.aic.A, type="l", lty=1, col="red", lwd=3, ylim=c(0,0.15),xlab="",ylab="",axes=F);
points(nvec,pct.err.null.bic.A,type="l", lty=2, col="blue", lwd=3);
abline(h=0.05)
legend("topright", legend=c(expression(paste(Delta," AIC")), expression(paste(Delta," SIC")) ), lty=c(1,2), lwd=c(3,3), col=c("red", "blue"),bty="n", cex=1.5)
axis(1, cex.axis=1.5, outer=FALSE) 
axis(2, cex.axis=1.5, outer=FALSE) 
box(bty="l")
mtext("A",side=3,adj=0.08,cex=2)

nmax <- 2003
nvec <- seq(6,nmax,by=10)
plot(nvec,pct.err.null.aic, type="l", lty=1, col="red", lwd=3, ylim=c(0,0.15),xlab="",ylab="",axes=F);
points(nvec,pct.err.null.bic,type="l", lty=2, col="blue", lwd=3);
abline(h=0.05)
legend("topright", legend=c(expression(paste(Delta," AIC")), expression(paste(Delta," SIC")) ), lty=c(1,2), lwd=c(3,3), col=c("red", "blue"),bty="n", cex=1.5)
axis(1, cex.axis=1.5, outer=FALSE) 
axis(2, cex.axis=1.5, outer=FALSE) 
box(bty="l")
mtext("B",side=3,adj=0.08,cex=2)

mtext(text=expression(paste(M[1], ", with k=2")), side=2, cex=1.5, outer=TRUE, at=0.5)
mtext(text="Sample size", side=1, cex=1.5, outer=TRUE, at=0.5)

dev.off()

save.image(file="Figure10.RData")

