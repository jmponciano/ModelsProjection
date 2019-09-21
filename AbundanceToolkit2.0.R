#    Function to retrieve an initial parameter guess
#    Input: vecor of count frequencies Y = [y1 = # of quadrats with 0 individs. present , y2 = #  of
#    quadrats with 1 individ. present, y3 = # of quadrats with 2 individs present,...,yk = # of quadrats with k
#    individuals present]
sample.mean <- function(y){
	
	k    <- length(y)
	km1  <- k-1
	ntot <- sum(y)
	xvec <- 0:km1
	lambda.hat <- sum(xvec*y)/ntot 
	return(lambda.hat)
}


# 1. Poisson negative log-likelihood
#    Input:  loglam = log of the value for the parameter lambda, Y = vector of count frequencies
poiss.negloglike <- function(loglam,y){
	
	lam <- exp(loglam)
	k   <- length(y)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1:km1] <- dpois(x[1:km1],lambda=lam)
	p[k] <- 1-sum(p[1:km1]);if(p[k]==0){p[k]<- .Machine$double.xmin}
	any.pi0 <- sum(p==0|(p<0))
	if(any.pi0>0){return(.Machine$double.xmax/1e100)}else{
		negloglike <-  -sum(y*log(p)) # constant of proportionality doesn't matter for minimization
		
		return(negloglike)
		}
}

# 2 Zero-Inflated Poisson a la Pielou
# guess = c(logitpi,loglam)
zipois.negloglike <- function(guess,y){
	
	pi1 <- 1/(1+exp(-guess[1]))
	lam <- exp(guess[2])
	k   <- length(y)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1] <- pi1 + (1-pi1)*exp(-lam)
	p[2:km1] <- (1-pi1)*dpois(x[2:km1],lambda=lam)
	p[k] <- 1-sum(p[1:km1]);if(p[k]==0){p[k]<- .Machine$double.xmin}
	#print(p)
	any.pi0 <- sum(p==0|(p<0))
	if(any.pi0>0){return(.Machine$double.xmax/1e100)}else{
		negloglike <-  -sum(y*log(p)) # constant of proportionality doesn't matter for minimization
		
		return(negloglike)
		}
	}

# 2.2 One-Inflated Poisson distribution
# guess = c(logitpi,loglam)
onein.pois.negloglike <- function(guess,y){
	
	pi1 <- 1/(1+exp(-guess[1]))
	lam <- exp(guess[2])
	k   <- length(y)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1] <- (1-pi1)*exp(-lam)
	p[2] <- pi1+ (1-pi1)*dpois(x[2],lambda=lam)
	p[3:km1] <- (1-pi1)*dpois(x[3:km1],lambda=lam)
	p[k] <- 1-sum(p[1:km1]);if(p[k]==0){p[k]<- .Machine$double.xmin}
	#print(p)
	any.pi0 <- sum(p==0|(p<0))
	if(any.pi0>0){return(.Machine$double.xmax/1e100)}else{
		negloglike <-  -sum(y*log(p)) # constant of proportionality doesn't matter for minimization
		
		return(negloglike)
		}
	}



#2.4 Poisson for x=0,1 and Neg Bin for x>=2
Pois.NB.nloglike <- function(guess,y){
	
    lam   <-  exp(guess[1])
	P     <- 1/(1+exp(-guess[2]))
	kappa <- exp(guess[3])
	kappap1 <- kappa+1
	Q     <- 1-P
	k     <- length(y)
	km1   <- k-1
	x     <- 0:km1
	p     <- rep(0,k)
	PXleq1<- sum(dpois(x[1:2], lambda=lam))
	PXgeq2<- 1- sum(dnbinom(x[1:2],size=kappa,prob=P))
	TotP <- PXleq1+PXgeq2
	p[1:2] <- dpois(x[1:2],lambda=lam)/TotP
	dnbinprobs <- dnbinom(x[3:km1],size=kappa,prob=P)/TotP
	dnbinprobs[is.na(dnbinprobs)] <- 1
	p[3:km1] <- dnbinprobs
	p[k]   <- 1-sum(p[1:km1]);if(p[k]==0){p[k]<- .Machine$double.xmin}
	any.pi0 <- sum(p==0|(p<0))
	if(any.pi0>0){return(.Machine$double.xmax/1e100)}else{
		negloglike <-  -sum(y*log(p)) # constant of proportionality doesn't matter for minimization
		
		return(negloglike)
		}
}

#2.6 NegBin for x=0,1 and Poisson for x>=2
NB.Pois.nloglike <- function(guess,y){
	
    lam   <-  exp(guess[1])
	P     <- 1/(1+exp(-guess[2]))
	kappa <- exp(guess[3])
	Q     <- 1-P
	k     <- length(y)
	km1   <- k-1
	x     <- 0:km1
	p     <- rep(0,k)
	PXleq1<- sum(dnbinom(x[1:2],size=kappa,prob=P))
	PXgeq2<- 1- sum(dpois(x[1:2], lambda=lam))
	TotP <- PXleq1+PXgeq2
	p[1:2] <- dnbinom(x[1:2],size=kappa,prob=P)/TotP
	dpoiprobs <- dpois(x[3:km1],lambda=lam)/TotP
	dpoiprobs[is.na(dpoiprobs)] <- 1
	p[3:km1] <- dpoiprobs
	p[k]   <- 1-sum(p[1:km1]);if(p[k]==0){p[k]<- .Machine$double.xmin}
	any.pi0 <- sum(p==0|(p<0))
	if(any.pi0>0){return(.Machine$double.xmax/1e100)}else{
		negloglike <-  -sum(y*log(p)) # constant of proportionality doesn't matter for minimization
		
		return(negloglike)
		}
}



# 3. Negative Binomial neg-log-likelihood
# Input: guess = c(log.alpha, log.kappa), y=vector of count frequencies
negbin.nloglike<- function(guess,y){
	
	guess <- exp(guess)
	alpha <- guess[1]
	kappa <- guess[2]
	P     <- alpha/(1+alpha)
	k     <- length(y)
	km1   <- k-1
	x     <- 0:km1
	p     <- rep(0,k)
	dnbinprobs <- dnbinom(x[1:km1],size=kappa,prob=P)
	dnbinprobs[is.na(dnbinprobs)] <- 1
	p[1:km1] <- dnbinprobs
	p[k]   <- 1-sum(p[1:km1]);if(p[k]==0){p[k]<- .Machine$double.xmin}
	any.pi0 <- sum(p==0|(p<0))
	if(any.pi0>0){return(.Machine$double.xmax/1e100)}else{
		negloglike <-  -sum(y*log(p)) # constant of proportionality doesn't matter for minimization
		
		return(negloglike)
		}
	
	}

# 4 Zero-Inflated Neg-Bin a la Pielou 1969, p. 88
# guess = c(logitpi,logalpha,logkappa)
zinbinom.negloglike <- function(guess,y){
	
	pi1   <- 1/(1+exp(-guess[1]))
	P      <- 1/(1+exp(-guess[2]))
	kappa <- exp(guess[3])
	#P     <- alpha/(1+alpha)
	k     <- length(y)
	km1   <- k-1
	x     <- 0:km1
	p <- rep(0,k)
	p[1] <- pi1 + (1-pi1)*P^(kappa) #P^(kappa) = P(X=0) for the Neg Bin model
	dnbinprobs <- dnbinom(x[2:km1],size=kappa,prob=P)
	dnbinprobs[is.na(dnbinprobs)] <- 1
	p[2:km1] <- (1-pi1)*dnbinprobs
	p[k] <- 1-sum(p[1:km1]);if(p[k]==0){p[k]<- .Machine$double.xmin}
	any.pi0 <- sum(p==0|(p<0))
	test <- any.pi0>0
	test[is.na(test)] <- 1
	if(any.pi0>0){return(.Machine$double.xmax/1e100)}else{
		negloglike <-  -sum(y*log(p)) # constant of proportionality doesn't matter for minimization
		
		return(negloglike)
		}
	}



# 5.2 One-Inflated Neg-Bin a la Pielou 1969, p. 88
# guess = c(logitpi,logalpha,logkappa)
onein.nbinom.negloglike <- function(guess,y){
	
	pi1   <- 1/(1+exp(-guess[1]))
	P      <- 1/(1+exp(-guess[2]))
	kappa <- exp(guess[3])
	#P     <- alpha/(1+alpha)
	k     <- length(y)
	km1   <- k-1
	x     <- 0:km1
	p <- rep(0,k)
	p[1] <- (1-pi1)*P^(kappa) #P^(kappa) = P(X=0) for the Neg Bin model
	dnbinprobs <- dnbinom(x[2:km1],size=kappa,prob=P)
	dnbinprobs[is.na(dnbinprobs)] <- 1
	p[2:km1] <- (1-pi1)*dnbinprobs
	p[2] <- pi1 + p[2]
	p[k] <- 1-sum(p[1:km1]);if(p[k]==0){p[k]<- .Machine$double.xmin}
	any.pi0 <- sum(p==0|(p<0))
	test <- any.pi0>0
	test[is.na(test)] <- 1
	if(any.pi0>0){return(.Machine$double.xmax/1e100)}else{
		negloglike <-  -sum(y*log(p)) # constant of proportionality doesn't matter for minimization
		
		return(negloglike)
		}
	}

# 5.5 Hurdle model: 0 term happens with prob. 'p', then for x>0, probs are scaled Neg.Binom.without the 0 term. 
hurdnbinom.negloglike <- function(guess,y){
	
	pi1   <- 1/(1+exp(-guess[1]))
	P      <- 1/(1+exp(-guess[2]))
	kappa <- exp(guess[3])
	#P     <- alpha/(1+alpha)
	k     <- length(y)
	km1   <- k-1
	x     <- 0:km1
	p <- rep(0,k)
	p[1] <- pi1  
	dnbinprobs <- dnbinom(x[2:km1],size=kappa,prob=P) #P^(kappa) = P(X=0) for the Neg Bin model
	dnbinprobs[is.na(dnbinprobs)] <- 1
	restofps <- ((1-pi1)/(1-P^(kappa)))*dnbinprobs
	restofps[is.na(restofps)] <- 1
	p[2:km1] <- restofps
	sumps.km1 <- sum(p[1:km1])
	if(is.na(sumps.km1)){sumps.km1 <- 0}
	p[k] <- 1-sumps.km1;if(p[k]<=0){p[k]<- .Machine$double.xmin}
	any.pi0 <- sum(p==0|(p<0))
	if(is.na(any.pi0)){any.pi0 <-1}
	if(any.pi0>0){return(.Machine$double.xmax/1e100)}else{
		negloglike <-  -sum(y*log(p)) # constant of proportionality doesn't matter for minimization
		
		return(negloglike)
		}
	}


# 6 Simulating from Zero-Inflated Poisson
zipoiss.sim <- function(theta,nn){
	
	pi1 <- theta[1]
	lam <- theta[2]
		
	# First flip a coin w.p. pi1:
	Uvec   <- runif(nn)
	outvec	 <- rep(0,nn)
	index.fornbin <- which((Uvec > pi1),arr.ind=TRUE)
	howmany.nbin <- length(index.fornbin)
	outvec[index.fornbin] <- rpois(n=howmany.nbin,lambda=lam)
	return(outvec)	
	}


# 7 Computing by hand the pmf of the zipois
dzipois<- function(theta,x){
	
	pi1   <- theta[1]
	lam   <- theta[2]	
	k <- length(x)
	logps      <- rep(0,k)
	logps[1]   <- log(pi1 + (1-pi1)*(exp(-lam)))  
	logps[2:k] <- log((1-pi1)*dpois(x=x[2:k],lambda=lam))
	return(exp(logps))	
	}




# Testing simulator of zinbinom:
#samp.size <- 10000
#test.samp <- zipoiss.sim(theta=c(0.5123541,2.771765),nn=samp.size)
#rel.freqs <- table(test.samp)/samp.size
#xx <- 0:(length(rel.freqs)-1)
#pred.relfreqs <- dzipois(theta=c(0.5123541,2.771765),x=xx)
#plot(xx,rel.freqs,pch=1)
#points(xx,pred.relfreqs,pch=16)


# 8. Function to compute the poisson ML estimate, maximized likelihood score, Expected frequencies, G-squared 
#    statistic and the p-value of the Goodness of fit LRT

Poisson.ml <- function(guess,yvec,my.method="BFGS"){
	
	max.out  <- optim(par=guess, fn=poiss.negloglike, method="Brent",lower=0.0001,upper=2*exp(guess),y=yvec)
	proplogL.hat <- -max.out$value
	lam.hat  <- exp(max.out$par[1])
	n        <- sum(yvec)
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	k   <- length(yvec)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1:km1] <- dpois(x[1:km1],lambda=lam.hat)
	p[k] <- 1-sum(p[1:km1]);if(p[k]==0|p[k]==-.Machine$double.xmin){p[k]<- .Machine$double.xmin}	
	Expect <- n*p
	y.nonzeros <- yvec
	y.nonzeros[y.nonzeros==0] <- 1 # protect against log(0) in G-squared
	E.nonzeros <- Expect
	E.nonzeros[E.nonzeros==0] <- 1
	Gsq <- 2*sum(yvec*log(y.nonzeros/E.nonzeros))
	if(k>2){pval<- 1-pchisq(Gsq,k-1-1)}else{pval<-0};
	bic <- -2*logL.hat+ length(guess)*log(n)
	
	return(list(Expected.freqs = Expect,lambda.hat = lam.hat, logLhat = logL.hat, Gsq=Gsq,pval=pval,bic=bic))
		
}


# 9. ML estimation and overall GOF LRTest for the ZIP (same results as in 8.)
zipois.ml<-function(guess,yvec,my.method="BFGS"){

	max.out  <- optim(par=guess, fn=zipois.negloglike, method="Nelder-Mead",y=yvec)
	proplogL.hat <- -max.out$value
	pi1.hat  <- 1/(1+exp(-max.out$par[1]))
	lam.hat  <- exp(max.out$par[2])
	n        <- sum(yvec)
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	k   <- length(yvec)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1] <- pi1.hat + (1-pi1.hat)*exp(-lam.hat)
	p[2:km1] <- (1-pi1.hat)*dpois(x[2:km1],lambda=lam.hat)
	p[k] <- 1-sum(p[1:km1]);if(p[k]==0|p[k]==-.Machine$double.xmin){p[k]<- .Machine$double.xmin}
	Expect <- n*p
	y.nonzeros <- yvec
	y.nonzeros[y.nonzeros==0] <- 1 # protect against log(0) in G-squared
	E.nonzeros <- Expect
	E.nonzeros[E.nonzeros==0] <- 1
	Gsq <- 2*sum(yvec*log(y.nonzeros/E.nonzeros))
	if(k>3){pval<- 1-pchisq(Gsq,k-1-2)}else{pval<-0};
	bic <- -2*logL.hat+ length(guess)*log(sum(yvec))
	
	return(list(Expected.freqs = Expect,pi1.hat = pi1.hat, lambda.hat = lam.hat, logLhat = logL.hat, Gsq=Gsq,pval=pval,bic=bic))	
}

# 9.5 ML estimation and overall GOF LRTest for the OIP (One Inflated Poisson)
oipois.ml<-function(guess,yvec,my.method="BFGS"){

	max.out  <- optim(par=guess, fn=onein.pois.negloglike, method="Nelder-Mead",y=yvec)
	proplogL.hat <- -max.out$value
	pi1.hat  <- 1/(1+exp(-max.out$par[1]))
	lam.hat  <- exp(max.out$par[2])
	n        <- sum(yvec)
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	k   <- length(yvec)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1] <- (1-pi1.hat)*exp(-lam.hat)
	p[2] <- pi1.hat+ (1-pi1.hat)*dpois(x[2],lambda=lam.hat)
	p[3:km1] <- (1-pi1.hat)*dpois(x[3:km1],lambda=lam.hat)
	p[k] <- 1-sum(p[1:km1]);if(p[k]==0|p[k]==-.Machine$double.xmin){p[k]<- .Machine$double.xmin}
	Expect <- n*p
	y.nonzeros <- yvec
	y.nonzeros[y.nonzeros==0] <- 1 # protect against log(0) in G-squared
	E.nonzeros <- Expect
	E.nonzeros[E.nonzeros==0] <- 1
	Gsq <- 2*sum(yvec*log(y.nonzeros/E.nonzeros))
	if(k>3){pval<- 1-pchisq(Gsq,k-1-2)}else{pval<-0};
	bic <- -2*logL.hat+ length(guess)*log(sum(yvec))
	
	return(list(Expected.freqs = Expect,pi1.hat = pi1.hat, lambda.hat = lam.hat, logLhat = logL.hat, Gsq=Gsq,pval=pval,bic=bic))	
}


# 10. Function to compute the Negative-Binomial ML estimate, maximized likelihood score, Expected frequencies
#    statistic and the p-value of the Goodness of fit LRT
Negbin.ml <- function(guess,yvec){
	
	max.out <- optim(par=guess, fn=	negbin.nloglike, method="Nelder-Mead", y=yvec)
	proplogL.hat <- -max.out$value
	mls          <- exp(max.out$par[1:2])
	alpha.hat  <- mls[1]
	kappa.hat  <- mls[2]
	P.hat      <- alpha.hat/(1+alpha.hat)
	n        <- sum(yvec)
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	k   <- length(yvec)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1:km1] <- dnbinom(x[1:km1],size=kappa.hat,prob=P.hat)
	p[k] <- 1-sum(p[1:km1]);if(p[k]==0|p[k]==-.Machine$double.xmin){p[k]<- .Machine$double.xmin}	
	Expect <- n*p
	y.nonzeros <- yvec
	y.nonzeros[y.nonzeros==0] <- 1 # protect against log(0) in G-squared
	E.nonzeros <- Expect
	E.nonzeros[E.nonzeros==0] <- 1
	Gsq <- 2*sum(yvec*log(y.nonzeros/E.nonzeros))
	if(k>3){pval<- 1-pchisq(Gsq,k-1-2)}else{pval<-0};
	bic <- -2*logL.hat+ length(guess)*log(sum(yvec))
	
	return(list(Expected.freqs = Expect,alpha.hat = alpha.hat, kappa.hat = kappa.hat, logLhat = logL.hat, Gsq=Gsq,pval=pval,bic=bic))
			
}



# 11. same as 10, but with the ZINegBinom
ZINegbin.ml <- function(guess,yvec){
	
	max.out <- optim(par=guess, fn=	zinbinom.negloglike, method="Nelder-Mead", y=yvec)
	proplogL.hat <- -max.out$value
	mls.a      <- 1/(1+exp(-max.out$par[1:2]))   
   pi1.hat    <- mls.a[1]
	P.hat      <- mls.a[2]
	kappa.hat  <- exp(max.out$par[3])
	#P.hat      <- alpha.hat/(1+alpha.hat)
	n        <- sum(yvec)
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	k   <- length(yvec)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1] <- pi1.hat + (1-pi1.hat)*P.hat^(kappa.hat) #P^(kappa) = P(X=0) for the Neg Bin model
	p[2:km1] <- (1-pi1.hat)*dnbinom(x[2:km1],size=kappa.hat,prob=P.hat)
	p[k] <- 1-sum(p[1:km1]);if(p[k]<0){p[k]<- .Machine$double.xmin}
	Expect <- n*p
	y.nonzeros <- yvec
	y.nonzeros[y.nonzeros==0] <- 1 # protect against log(0) in G-squared
	E.nonzeros <- Expect
	E.nonzeros[E.nonzeros==0] <- 1
	Gsq <- 2*sum(yvec*log(y.nonzeros/E.nonzeros))
	if(k>4){pval<- 1-pchisq(Gsq,k-1-3)}else{pval<-0}
	bic <- -2*logL.hat+ length(guess)*log(sum(yvec))
	
	return(list(Expected.freqs = Expect,pi1.hat=pi1.hat,P.hat=P.hat, kappa.hat = kappa.hat, logLhat = logL.hat, Gsq=Gsq,pval=pval, bic=bic))
			
}


# 11.5 . same as 10, but with the oninegBinom
OINegbin.ml <- function(guess,yvec){
	
	max.out <- optim(par=guess, fn=	onein.nbinom.negloglike, method="Nelder-Mead", y=yvec)
	proplogL.hat <- -max.out$value
	mls.a      <- 1/(1+exp(-max.out$par[1:2]))   
   pi1.hat    <- mls.a[1]
	P.hat      <- mls.a[2]
	kappa.hat  <- exp(max.out$par[3])
	#P.hat      <- alpha.hat/(1+alpha.hat)
	n        <- sum(yvec)
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	k   <- length(yvec)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1] <- (1-pi1.hat)*P.hat^(kappa.hat) #P^(kappa) = P(X=0) for the Neg Bin model
	p[2:km1] <- (1-pi1.hat)*dnbinom(x[2:km1],size=kappa.hat,prob=P.hat)
	p[2] <- pi1.hat + p[2] #one inflation
	p[k] <- 1-sum(p[1:km1]);if(p[k]<0){p[k]<- .Machine$double.xmin}
	Expect <- n*p
	y.nonzeros <- yvec
	y.nonzeros[y.nonzeros==0] <- 1 # protect against log(0) in G-squared
	E.nonzeros <- Expect
	E.nonzeros[E.nonzeros==0] <- 1
	Gsq <- 2*sum(yvec*log(y.nonzeros/E.nonzeros))
	if(k>4){pval<- 1-pchisq(Gsq,k-1-3)}else{pval<-0}
	bic <- -2*logL.hat+ length(guess)*log(sum(yvec))
	
	return(list(Expected.freqs = Expect,pi1.hat=pi1.hat,P.hat=P.hat, kappa.hat = kappa.hat, logLhat = logL.hat, Gsq=Gsq,pval=pval, bic=bic))
			
}





# hurdnbinom.negloglike
# 12. same as 10, but with the ZINegBinom

Hurdnegbin.ml <- function(guess,yvec){
	
	max.out <- optim(par=guess, fn=	hurdnbinom.negloglike, method="Nelder-Mead", y=yvec)
	proplogL.hat <- -max.out$value
	mls.a      <- 1/(1+exp(-max.out$par[1:2]))   
   pi1.hat    <- mls.a[1]
	P.hat      <- mls.a[2]
	kappa.hat  <- exp(max.out$par[3])
	#P.hat      <- alpha.hat/(1+alpha.hat)
	n        <- sum(yvec)
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	k   <- length(yvec)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1] <- pi1.hat #P^(kappa) = P(X=0) for the Neg Bin model
	p[2:km1] <- ((1-pi1.hat)/(1-P.hat^(kappa.hat)))*dnbinom(x[2:km1],size=kappa.hat,prob=P.hat)
	p[k] <- 1-sum(p[1:km1]);if(p[k]<=0){p[k]<- .Machine$double.xmin}
	Expect <- n*p
	y.nonzeros <- yvec
	y.nonzeros[y.nonzeros==0] <- 1 # protect against log(0) in G-squared
	y.nonzeros <- y.nonzeros*(1+.Machine$double.eps)
	E.nonzeros <- Expect
	E.nonzeros[E.nonzeros==0] <- 1
	E.nonzeros <- E.nonzeros*(1+.Machine$double.eps)
	Gsq <- 2*sum(yvec*log(y.nonzeros/E.nonzeros))
	if(k>4){pval<- 1-pchisq(Gsq,k-1-3)}else{pval<-0};
	bic <- -2*logL.hat+ length(guess)*log(sum(yvec))
	
	return(list(Expected.freqs = Expect,pi1.hat=pi1.hat,P.hat=P.hat, kappa.hat = kappa.hat, logLhat = logL.hat, Gsq=Gsq,pval=pval, bic=bic))
			
}


#12.4 Poisson for x=0,1 and Neg Bin for x>=2
Pois.NB.ml <- function(guess,yvec,my.method="BFGS"){
	n       <- sum(yvec)
	#Pois.NB.nloglike(guess=guess,y=yvec)
	max.out <- optim(par=guess, fn=Pois.NB.nloglike, method="Nelder-Mead", y=yvec)
	proplogL.hat <- -max.out$value
    lam   <-  exp(max.out$par[1])
	P     <- 1/(1+exp(-max.out$par[2]))
	kappa <- exp(max.out$par[3])
	kappap1 <- kappa+1
	Q     <- 1-P
	k     <- length(yvec)
	km1   <- k-1
	x     <- 0:km1
	p     <- rep(0,k)
	PXleq1<- sum(dpois(x[1:2], lambda=lam))
	PXgeq2<- 1-sum(dnbinom(x[1:2],size=kappa,prob=P))
	TotP <- PXleq1+PXgeq2
	p[1:2] <- dpois(x[1:2],lambda=lam)/TotP
	dnbinprobs <- dnbinom(x[3:km1],size=kappa,prob=P)/TotP
	dnbinprobs[is.na(dnbinprobs)] <- 1
	p[3:km1] <- dnbinprobs
	p[k]   <- 1-sum(p[1:km1]);if(p[k]==0|p[k]==-.Machine$double.xmin){p[k]<- .Machine$double.xmin}
	Expect <- n*p
	y.nonzeros <- yvec
	y.nonzeros[y.nonzeros==0] <- 1 # protect against log(0) in G-squared
	E.nonzeros <- Expect
	E.nonzeros[E.nonzeros==0] <- 1
	Gsq <- 2*sum(yvec*log(y.nonzeros/E.nonzeros))
	if(k>3){pval<- 1-pchisq(Gsq,k-1-2)}else{pval<-0};
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	bic <- -2*logL.hat+ length(guess)*log(sum(yvec))
	
	return(list(Expected.freqs = Expect,lam.hat=lam, kappa.hat = kappa, P.hat=P, logLhat = logL.hat, Gsq=Gsq,pval=pval,bic=bic))


}

#12.6 NegBin for x=0,1 and Poisson for x>=2
NB.Pois.ml <- function(guess,yvec,my.method="BFGS"){

	n       <- sum(yvec)
	max.out <- optim(par=guess, fn=NB.Pois.nloglike, method="Nelder-Mead", y=yvec)
	proplogL.hat <- -max.out$value
    lam   <-  exp(max.out$par[1])
	P     <- 1/(1+exp(-max.out$par[2]))
	kappa <- exp(max.out$par[3])
	Q     <- 1-P
	k     <- length(yvec)
	km1   <- k-1
	x     <- 0:km1
	p     <- rep(0,k)
	PXleq1<- sum(dnbinom(x[1:2],size=kappa,prob=P))
	PXgeq2<- 1- sum(dpois(x[1:2], lambda=lam))
	TotP <- PXleq1+PXgeq2
	p[1:2] <- dnbinom(x[1:2],size=kappa,prob=P)/TotP
	dpoiprobs <- dpois(x[3:km1],lambda=lam)/TotP
	dpoiprobs[is.na(dpoiprobs)] <- 1
	p[3:km1] <- dpoiprobs
	p[k]   <- 1-sum(p[1:km1]);if(p[k]==0|p[k]==-.Machine$double.xmin){p[k]<- .Machine$double.xmin}
	Expect <- n*p
	y.nonzeros <- yvec
	y.nonzeros[y.nonzeros==0] <- 1 # protect against log(0) in G-squared
	E.nonzeros <- Expect
	E.nonzeros[E.nonzeros==0] <- 1
	Gsq <- 2*sum(yvec*log(y.nonzeros/E.nonzeros))
	if(k>3){pval<- 1-pchisq(Gsq,k-1-2)}else{pval<-0};
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	bic <- -2*logL.hat+ length(guess)*log(sum(yvec))
	
	return(list(Expected.freqs = Expect,lam.hat=lam, kappa.hat = kappa, P.hat=P, logLhat = logL.hat, Gsq=Gsq,pval=pval,bic=bic))
}





# 13 Wrapper function to fit all the models at once
abund.fit <- function(y1, names.animals="ungulates", name.counts="Patroled",pois.method="BFGS", plot.it=TRUE){
	
	method <- name.counts
	k <- length(y1)
	km1 <- k-1
	xvec <- 0:km1
	tot <- sum(y1)
	lam0 <- sample.mean(y1)
	if(sum(y1)==0){log.lam0<-log(0.0000001)}
	log.lam0 <- log(lam0)	
	if(is.na(log.lam0)){log.lam0<- 10}
	Poisson.Stats <- Poisson.ml(log.lam0,y1,my.method=pois.method)	
	NegBin.Stats <- Negbin.ml(log(c(1.5,4)),y1)
	ZIPoiss.Stats <- zipois.ml(c(log(0.01)-log(1-0.01),log.lam0),y1, my.method=pois.method)
	ZINegBi.Stats <- ZINegbin.ml(c(log(0.51)-log(1-0.51),log(0.99)-log(1-0.99),log(420)),y1)
    HurdNBi.Stats <- Hurdnegbin.ml(c(log(0.51)-log(1-0.51),log(0.99)-log(1-0.99),log(420)),y1)
	PoisNB.Stats <- Pois.NB.ml(c(log.lam0,log(0.99)-log(1-0.99),log(420)),y1)
	NBPois.Stats <- NB.Pois.ml(c(log.lam0,log(0.99)-log(1-0.99),log(420)),y1)
	OIPoiss.Stats <- oipois.ml(c(log(0.01)-log(1-0.01),log.lam0),y1, my.method=pois.method)
	OINegBi.Stats <- OINegbin.ml(c(log(0.51)-log(1-0.51),log(0.99)-log(1-0.99),log(420)),y1)



	Expected.Pois <- Poisson.Stats$Expected.freqs
	Expected.NBin <- NegBin.Stats$Expected.freqs
	Expected.ZIP <- ZIPoiss.Stats$Expected.freqs
	Expected.ZINB<- ZINegBi.Stats$Expected.freqs
	Expected.PoisNB <- PoisNB.Stats$Expected.freqs
	Expected.NBPois <- NBPois.Stats$Expected.freqs
	Expected.OIP <- OIPoiss.Stats$Expected.freqs
	Expected.OINB<- OINegBi.Stats$Expected.freqs



	BIC.Pois <- Poisson.Stats$bic
	BIC.NBin <- NegBin.Stats$bic
	BIC.ZIP <- ZIPoiss.Stats$bic
	BIC.ZINB<- ZINegBi.Stats$bic
	BIC.PoisNB <- PoisNB.Stats$bic
	BIC.NBPois <- NBPois.Stats$bic	
	BIC.OIP <- OIPoiss.Stats$bic
	BIC.OINB<- OINegBi.Stats$bic



	###### The figure
	
	if(plot.it==TRUE){
		
		my.xlab=paste("Number of ", names.animals,sep="");
		my.ylab=paste(method, " Frequency", sep="");
		
		ymaxplot <- max(c(y1,Expected.Pois,Expected.NBin,Expected.ZIP, Expected.ZINB))
		ymaxplot <- ymaxplot+0.20*ymaxplot;
		
		par(mfrow=c(2,2),oma=c(3,2,1,0.1),mar=c(3,4,3,2), cex.axis=1.2,cex.lab=1.2)
		plot(xvec,y1,pch=1,col.axis="white",xlab=my.xlab,ylab=my.ylab, main="Poisson",ylim=c(0,ymaxplot))
		legend("topright",legend=c("Observed","Expected"),pch=c(1,16),bty="n")
		points(xvec,Expected.Pois,pch=16, col.axis="white")
		axis(side=2)
		axis(side=1,labels=c(as.character(xvec[1:km1]),paste(">= ",as.character(xvec[k]))),at=xvec)


		plot(xvec,y1,pch=1,col.axis="white",xlab=my.xlab,ylab=my.ylab, main="Neg. Binomial",ylim=c(0,ymaxplot))
		legend("topright",legend=c("Observed","Expected"),pch=c(1,16),bty="n")
		points(xvec,Expected.NBin,pch=16,col.axis="white")
		axis(side=2)
		axis(side=1,labels=c(as.character(xvec[1:km1]),paste(">= ",as.character(xvec[k]))),at=xvec)

		plot(xvec,y1,pch=1,col.axis="white",xlab=my.xlab,ylab=my.ylab, main= "Zero Inflated Poisson",ylim=c(0,ymaxplot))
		legend("topright",legend=c("Observed","Expected"),pch=c(1,16),bty="n")
		points(xvec,Expected.ZIP,pch=16, col.axis="white")
		axis(side=2)
		axis(side=1,labels=c(as.character(xvec[1:km1]),paste(">= ",as.character(xvec[k]))),at=xvec)

		plot(xvec,y1,pch=1,col.axis="white",xlab=my.xlab,ylab=my.ylab, main="Zero Inflated Neg-Binomial",ylim=c(0,ymaxplot))
		legend("topright",legend=c("Observed","Expected"),pch=c(1,16),bty="n")
		points(xvec,Expected.ZINB,pch=16, col.axis="white")
		axis(side=2)
		axis(side=1,labels=c(as.character(xvec[1:km1]),paste(">= ",as.character(xvec[k]))),at=xvec)
	}
	results <- list(Poisson.Stats=Poisson.Stats,NegBin.Stats=NegBin.Stats,ZIPoiss.Stats=ZIPoiss.Stats,ZINegBi.Stats=ZINegBi.Stats, HurdNBi.Stats=HurdNBi.Stats, PoisNB.Stats=PoisNB.Stats, NBPois.Stats=NBPois.Stats, OIPoiss.Stats=OIPoiss.Stats, OINegBi.Stats=OINegBi.Stats)
	
	
	return(results)
}


########## Other things
# I tested if the reduced expression involving "expected" was identical to the straightforward likelihood 
# expression
Poisson.ml2 <- function(guess,yvec){
	
	max.out  <- optim(par=guess, fn=poiss.negloglike, method="BFGS",y=yvec)
	proplogL.hat <- -max.out$value
	lam.hat  <- exp(max.out$par[1])
	n        <- sum(yvec)
	logL.hat <- proplogL.hat +  lfactorial(n) - sum(lfactorial(yvec))
	k   <- length(yvec)
	km1 <- k-1
	x <- 0:km1
	p <- rep(0,k)
	p[1:km1] <- dpois(x[1:km1],lambda=lam.hat)
	p[k] <- 1-sum(p[1:km1])	
	Expect <- n*p
	phats.full <- yvec/n
	Full.logL <- dmultinom(x=yvec,size=n,prob=phats.full,log=TRUE)
	Gsq       <- -2*(logL.hat-Full.logL)
	pval <- 1-pchisq(Gsq,k-1-1)
		
	return(list(Expected.freqs = Expect,lambda.hat = lam.hat, logLhat = logL.hat, Gsq=Gsq,pval=pval))
		
}

