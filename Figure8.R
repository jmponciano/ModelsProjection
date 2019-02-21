# Testing the non-parametric Sgg estimator from 
# Berrett, Samworth & Yuan 2016 (arXiv:1606.00304v3)
#  This program contains only 3 functions and a test of Berret et al's estimator

library(pracma)#
library(plyr)
library(mvtnorm)

# Analytical Neg entropy for a Multivariate Normal (MVN) distribution
H.MVnorm <- function(par.lst=list(mu=c(0,0),sigma=matrix(c(1,0,0,1),nrow=2))) {
  sigma=par.lst$sigma
  d=dim(sigma)[1]
  out <- log(((2*pi*exp(1))^d)*det(sigma))/2
  return(out)
}

# Auxiliary function for the function Hse.wKL which computes Berrett et al's estimator:
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


#########  PROCEDURAL PROGRAM SECTION ###########

# 1. Generate data from a Multivariate Normal distribution, for which we know the
#    analytical Sgg form.

mu <- rep(10,7)
sigma <- diag(rep(1,7))
X <- rmvnorm(n =300,mean = mu,sigma = sigma)

# 2. Compute the NP Sgg estimate from Berrett et al and
NPSgg.hat <- Hse.wKL(X=X, k=0)$Hw 

# 3. Repeat
B <- 2000
nvec <- c(10,25,50,75,150)
nsamps <- length(nvec)
NPSgg.hats <- matrix(0,nrow=B,ncol=nsamps)

for(i in 1:nsamps){
	
	n <- nvec[i]
	
	for(j in 1:B){
	
		X.star   <- rmvnorm(n = n,mean = mu,sigma = sigma) 
		NPSgg.hats[j,i] <- Hse.wKL(X=X.star, k=0)$Hw  
	}
}



# 4. Compute the true Parametric Sgg for this MVN distribution and compare to 
#    Monte Carlo simulations
H.true <- H.MVnorm(par.lst = list(mu=mu,sigma=sigma))

rel.NPSgg.hats <- NPSgg.hats/H.true
save.image("NonParametricSggBootTest.RData")


postscript("Figure8.eps",width=22,height=22,horizontal=FALSE)
#tiff("Figure8.tiff",width=22,height=22,units="cm",res=600,compression="lzw",type="cairo") # w/l 22/22
par(oma=c(2.5,2.5,1.5,1.5),mar=c(3.5,3.5,1.5,1), mgp=c(2,0.75,0))
boxplot(rel.NPSgg.hats,axes=FALSE,outline=FALSE, width=rep(0.75,nsamps))
abline(h=1,col="black",lty=2, lwd=2)
axis(side=1, at=1:nsamps,labels=nvec)
mtext(text="Sample size", side=1, adj=0.5,outer=TRUE, cex=2)
axis(side=2, at=c(0.5,0.75,1,1.25,1.5), labels=c(0.5,0.75,1,1.25,1.5))
mtext(text="Relative estimate", side=2, adj=0.5,outer=TRUE, cex=2)
dev.off()





