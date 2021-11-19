
# Figure 6 in Dennis et al 2019

#model1
p1   = 0.75 # P(x=1)
omp1 = 1-p1 # P(X=0)

#model2
p2   = 0.50 # P(X=1)
omp2 = 1-p2 # P(X=0)

# KL(f1,f2)  model1 true
K    = p1*log(p1/p2) + omp1*log(omp1/omp2)

#misspecified
preal  =  .65
ompreal  =  1-preal
DK  =  preal*log(p1/p2) + ompreal*log(omp1/omp2)

#plot(preal,DK,type="l")
#abline(h=0,lty=2)


# Set the maximum sample size
nmax = 2000

# Compute Eg(log(f1(X)/f2(X))) and its variance under truth
Eg = abs(DK)
sigma.sq = preal*(log(p1/p2))^2 + ompreal*(log(omp1/omp2))^2 - Eg^2
sigma  = sqrt(sigma.sq)


# Now, pick 'k' so that n.star = 125 approx.
ks        = 2:20
nstar.vec = log(ks)/Eg
cbind(ks,nstar.vec)

#  Then we let k = 14
k  = 14
c1.vec <-  (1/sqrt(1:nmax))*log(k) - sqrt(1:nmax)*Eg
c2.vec <- -((1/sqrt(1:nmax))*log(k) + sqrt(1:nmax)*Eg)
M.vec  <- pnorm(q=c2.vec/sigma)
Q.vec  <- 1-pnorm(q=c1.vec/sigma)

#plot(1:nmax, Q.vec, type="l", col="blue", lwd=2, ylim=c(0,1))
#windows()
#plot(1:nmax, M.vec, type="l", col="red", lwd=2)

# Simulating a single realization of the accumulation of evidence
L.ratio.vec = rep(0,nmax)
BM.approx  = rep(0,nmax)
xs.obs  = rbinom(n=nmax, size=1,prob=preal)
L.ratio = 1

for(n in 1:nmax){
  
  xn   = xs.obs[n]
  
  L1   = xn*p1 + (1-xn)*omp1
  L2   = xn*p2 + (1-xn)*omp2
  
  L.ratio = L.ratio*(L1/L2)
  lnL     = log(L.ratio)
  L.ratio.vec[n] = L.ratio
  BM.approx[n]   = (1/sqrt(n))*lnL - sqrt(n)*DK
  
}

#windows()
#par(mfrow=c(2,1))
#plot(1:nmax, L.ratio.vec, type="b")
#plot(1:nmax, BM.approx, type="b", pch=16)

# Simulating many realizations of the accumulation of evidence
# and computing the proportion of times, for each n, the
# BM process goes past the thresholds c2 and c1


nsims = 50000

Q.countmat = matrix(0,nrow=nsims, ncol=nmax)
M.countmat = matrix(0,nrow=nsims, ncol=nmax)

for(i in 1:nsims){
  Q.count   = rep(0,nmax)
  M.count   = rep(0,nmax)
  
  L.ratio.vec = rep(0,nmax)
  BM.approx  = rep(0,nmax)
  xs.obs      = rbinom(n=nmax, size=1,prob=preal)
  L.ratio = 1
  
  for(n in 1:nmax){
    
    xn   = xs.obs[n]
    
    L1   = xn*p1 + (1-xn)*omp1
    L2   = xn*p2 + (1-xn)*omp2
    
    L.ratio = L.ratio*(L1/L2)
    lnL     = log(L.ratio)
    L.ratio.vec[n] = L.ratio
    BM.approxn     = (1/sqrt(n))*lnL - sqrt(n)*Eg
    BM.approx[n]   = BM.approxn
    
    Q.count[n]  = BM.approxn > c1.vec[n]
    M.count[n]  = BM.approxn < c2.vec[n]
    
  }
  
  Q.countmat[i,] = Q.count
  M.countmat[i,] = M.count
  
}

Q.n = apply(Q.countmat,2,sum)/nsims
M.n = apply(M.countmat,2,sum)/nsims


figname  = "Figure6.eps"

#tiff(figname, width=7,height=6, units="in", res=600, compression="lzw", type="cairo")
postscript(figname, horizontal=TRUE, paper="letter")

par(mfrow=c(1,2), oma=c(2,2,1,1), mar=c(4,4,2,1))
plot(1:nmax, Q.n, type="l", col="black", lwd=2, ylim=c(0,1),
     ylab="P(Rejecting H1|H1 is closer to truth)", axes=F,cex.lab=1.25)
points(1:nmax, Q.vec, type="l", col="grey", lwd=2)
legend("topleft", legend=c("Simulated", "B.M. approx."), border="white", bty="n", lty=c(1,1), lwd=c(1,2), col=c("black","black"), cex=1.05)
mtext("A", side=3,adj=0.01,cex=2)
axis(side=1)
axis(side=2)

plot(1:nmax, M.n, type="l", col="black", lwd=2, ylim=c(0,max(M.n)),
     ylab="P(Misleading evidence for H2)",axes=F, cex.lab=1.25)
points(1:nmax, M.vec, type="l", col="grey",lwd=2)
legend("topright", legend=c("Simulated", "B.M. approx."), border="white", bty="n", lty=c(1,1), lwd=c(1,2), col=c("black","black"), cex=1.05)
mtext("B", side=3,adj=0.01,cex=2)
axis(side=1)
axis(side=2)

# Single x-axis
mtext(text="Sample size", side=1, outer=TRUE, cex=1.25)
dev.off()



