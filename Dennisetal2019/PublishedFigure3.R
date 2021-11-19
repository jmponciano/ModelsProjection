#truth
p1   <- 0.75 # P(x=1)
omp1 <- 1-p1 # P(X=0)

#false
p2   <- 0.50 # P(X=1)
omp2 <- 1-p2 # P(X=0)

# KL(g,f)
K    <- p1*log(p1/p2) + omp1*log(omp1/omp2) 

# Set the maximum sample size
nmax <- 100

# Compute E(log(f1(X)/f2(X))) and its variance under truth
E1 <- K
sigma.sq <- p1*(log(p1/p2))^2 + omp1*(log(omp1/omp2))^2 - E1^2
sigma    <- sqrt(sigma.sq)


# Now, pick 'C' so that n.star = 20 approx.
cs        <- 2:20
nstar.vec <- log(cs)/K
print(nstar.vec)
cbind(cs,nstar.vec)

#  Then we let C = 14
C  <- 14
c1.vec <-  (1/sqrt(1:nmax))*log(C) - sqrt(1:nmax)*K
c2.vec <- -((1/sqrt(1:nmax))*log(C) + sqrt(1:nmax)*K)
M.vec  <- pnorm(q=c2.vec/sigma)
Q.vec  <- 1-pnorm(q=c1.vec/sigma)

#plot(1:nmax, Q.vec, type="l", col="blue", lwd=2, ylim=c(0,1))
#plot(1:nmax, M.vec, type="l", col="red", lwd=2)

# Simulating a single realization of the accumulation of evidence
#lnL.ratio.vec <- rep(0,nmax)
#BM.approx  <- rep(0,nmax)
#xs.obs      <- rbinom(n=nmax, size=1,prob=p1)
#lnL.ratio <- 0
#
#for(n in 1:nmax){
#	
#	xn   <- xs.obs[n]
#	
# 	lnL1   <- log((p1^xn)*(omp1^(1-xn)))
#	lnL2   <- log((p2^xn)*(omp2^(1-xn)))
#	
#	lnL.ratio <- lnL.ratio + (lnL1-lnL2)
#	lnL.ratio.vec[n] <- lnL.ratio
#	BM.approx[n]   <- (1/sqrt(n))*lnL.ratio - sqrt(n)*K
# 	
#}
#par(mfrow=c(2,1))
#plot(1:nmax, lnL.ratio.vec, type="b")
#plot(1:nmax, BM.approx, type="b", pch=16)

# Simulating many realizations of the accumulation of evidence
# and computing the proportion of times, for each n, the
# BM process goes past the thresholds c2 and c1


nsims <- 100000

Q.countmat <- matrix(0,nrow=nsims, ncol=nmax)
M.countmat <- matrix(0,nrow=nsims, ncol=nmax)

for(i in 1:nsims){
  Q.count   <- rep(0,nmax)
  M.count   <- rep(0,nmax)
  
  L.ratio.vec <- rep(0,nmax)
  BM.approx  <- rep(0,nmax)
  xs.obs      <- rbinom(n=nmax, size=1,prob=p1)
  L.ratio <- 1
  
  for(n in 1:nmax){
    
    xn   <- xs.obs[n]
    
    L1   <- (p1^xn)*(omp1^(1-xn))
    L2   <- (p2^xn)*(omp2^(1-xn))
    # Continuity correction calculation in the probabilities
    #L1a  <- pnorm(q=(xn-1/2), mean=p1,sd=p1*(1-p1))
    #L1b  <- pnorm(q=(xn+1/2), mean=p1,sd=p1*(1-p1))
    #L1   <- L1b-L1a 
    
    #L2a  <- pnorm(q=(xn-1/2), mean=p2, sd=p2*(1-p2))
    #L2b  <- pnorm(q=(xn+1/2), mean=p2, sd=p2*(1-p2))
    #L2   <- L2b-L2a 
    
    L.ratio <- L.ratio*(L1/L2)
    lnL     <- log(L.ratio)
    L.ratio.vec[n] <- L.ratio
    BM.approxn     <- (1/sqrt(n))*lnL - sqrt(n)*K
    BM.approx[n]   <- BM.approxn
    
    Q.count[n]  <- BM.approxn > (c1.vec[n])
    M.count[n]  <- BM.approxn < (c2.vec[n])
    
  }
  
  Q.countmat[i,] <- Q.count
  M.countmat[i,] <- M.count
  
}

Q.n <- apply(Q.countmat,2,sum)/nsims
M.n <- apply(M.countmat,2,sum)/nsims



figpath1 <- "~/Dropbox/BrianMarkSubhash/BriansMisSpecPaper/FrontiersInLatex/PublishedFigures&Code"
figname  <- "Figure3.eps" #OLD NAME: "M&Q-withBMapprox.tiff"
figpath <- paste0(figpath1,"/",figname)

#tiff(figpath, width=6,height=7, units="in", res=600, compression="lzw", type="cairo")

postscript(figpath, horizontal=FALSE, paper="letter")

par(mfrow=c(2,1), oma=c(2,2,2,2), mar=c(3,5,2,2))
plot(1:nmax, Q.vec, type="l", lty="dashed", col="black", ylim=c(0,1),
     ylab=expression(V[1]), lwd=2, axes=F,cex.lab=1.25)
points(1:nmax, Q.n, type="l", lty="solid", col="black", lwd=1)
legend(x=72,y=0.88, legend=c("Simulated", "B.M. approx."), border="white", bty="n", lty=c("solid","dashed"), lwd=c(1,2), col=c("black","black"), cex=1.05)
mtext("A",side=3,adj=0.01,cex=2)
axis(side=1)
axis(side=2)

plot(1:nmax, M.vec, type="l", lty="dashed", col="black", lwd=2, ylim=c(0,max(M.n)), ylab=expression(M[1]), axes=F,cex.lab=1.25)
points(1:nmax, M.n,type="l", lty="solid", col="black", lwd=1)
legend("topright", legend=c("Simulated", "B.M. approx."), border="white", bty="n", lty=c("solid","dashed"), lwd=c(1,2), col=c("black","black"), cex=1.05)
mtext("B",side=3,adj=0.01,cex=2)
axis(side=1)
axis(side=2)

# Single x-axis
mtext(text="Sample size", side=1, outer=TRUE, cex=1.25)


dev.off()



