library(MASS)
library(cubature)


### Functions used:

### A small function to draw 'n' colours:
mycols.ftn <- function(n){
	d <- 360/n
	h <- cumsum(c(15,rep(d,(n-1))))
	return(hcl(h=h, c=100,l=65))
}


# Multivariate Normal simulator:
my.rmvn <- function(n,mu.vec, cov.mat){

	# mu.vec  = vector of means, dimension will be labeled as 'p'
	# cov.mat = (positive definite) variance-covariance matrix (a pxp matrix)
	# n       = desired number of multivariate samples of dimension 'p'  	
	
	p <- length(mu.vec);
	Tau <- chol(cov.mat);
	Zmat <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n);
	out <- matrix(0,nrow=p,ncol=n);
	for(i in 1:n){
		
		Z <- Zmat[,i];
		out[,i] <- t(Tau)%*%Z + mu.vec
		
		}
	
	return(out)
	
}

# # Example (uncomment to run)
#muvec <- c(0,0) # dimension = 2
#Sigma <- matrix(c(0.3,0.1,0.1,0.4), nrow=2,ncol=2,byrow=TRUE)
# Now simulating many many realizations and computing the empirical mean
# and empirical variance-covariance matrix of these realizations to see if we 
# recover the means in muvec and the variance-covariance matrix Sigma:
# each column of the output of the 'my.rmvn' ftn. is a simulated realization of the MVN vector
#nreps <- 100000
#long.trial <- my.rmvn(n=nreps,mu.vec=muvec, cov.mat=Sigma)  
#emp.means <- apply(long.trial,1,mean)
#emp.varcov <- var(t(long.trial)) # compute the variance of the transpose of long.trial
# # Now we compare (not bad, ha?)
# print(emp.means)
# print(muvec)
# print(emp.varcov)
# print(Sigma)

# Multivariate Laplace simulator
rbivlaplace <- function(n,mu.vec, cov.mat){

	# mu.vec  = vector of means, 'p', but p will always be 2
	#           note: if mu.vec is a vector of 0's, this function returns
	#                 samples from the symmetric Laplace distribution;
	#                 otherwise, the samples are bivariate, asymmetric laplace 
	#                 distributed.
	# cov.mat = (positive definite) variance-covariance matrix (a pxp matrix)
	# n       = desired number of multivariate samples of dimension 'p'  	
	
	p <- length(mu.vec);if(p!=2){print("Error, doesn't work with dimension p > 2");break}
	Tau <- chol(cov.mat);
	Zmat <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n);
	W.vec<- rexp(n=n,rate=1);
	out <- matrix(0,nrow=p,ncol=n);
	for(i in 1:n){
		
		Z <- Zmat[,i];
		W <- W.vec[i];
		Y <- t(Tau)%*%Z # mvrnorm, with mean 0
		
		out[,i] <- sqrt(W)*Y + W*mu.vec
		
		}
		
		
	
	return(out)
	
}

dmvnorm <- function(xvec,mu.vec,cov.mat,log=FALSE){
	
	p      <- length(mu.vec)
	Sigma  <- cov.mat
	Sigma.inv <- ginv(Sigma)
	if(log==FALSE){

		pi.det <- ((2*pi)^(p/2))*det(Sigma)^(1/2)
		exp.bit<- exp(- t(xvec-mu.vec)%*%Sigma.inv%*%(xvec-mu.vec)*0.5)	
		f.x    <- (1/pi.det)*exp.bit
		return(f.x)
	}else{
		quad.form<- - t(xvec-mu.vec)%*%Sigma.inv%*%(xvec-mu.vec)*0.5	
		lf.x    <- -(p/2)*log(2*pi) -(1/2)*det(Sigma) + quad.form
		return(lf.x)
	}
}

my.outerbivnorm <- function(x1vals,x2vals, mu.vec,cov.mat){
	
	len1 <- length(x1vals)
	len2 <- length(x2vals)
	out  <- matrix(0,nrow=len1,ncol=len2)
	
	for(i in 1:len1){
		
		x1 <- x1vals[i]
		
		for(j in 1:len2){
			
			x2 <- x2vals[j]
			out[i,j] <- 	dmvnorm(xvec=c(x1,x2),mu.vec=mu.vec,cov.mat=cov.mat)	
		}
		
	}
	
	return(out)
}


dbivlaplace <- function(xvec,mu.vec,cov.mat,log=FALSE){
	
	p      <- length(mu.vec)
	Sigma  <- cov.mat
	Sigma.inv <- ginv(Sigma)
	nu        <- (2-p)/2
	pi.det <- ((2*pi)^(p/2))*sqrt(det(Sigma))
	quadform.x <- t(xvec)%*%Sigma.inv%*%(xvec)
	quadform.mu<- t(mu.vec)%*%Sigma.inv%*%(mu.vec)	
	exp.bit.cross <- exp(t(xvec)%*%Sigma.inv%*%(mu.vec))
	bessel.x  <- sqrt((2+quadform.mu)*(quadform.x))
	part1    <- ((2*exp.bit.cross)/pi.det)
	part2    <- (quadform.x/(2+quadform.mu))^(nu/2)
	part3    <- besselK(x=bessel.x, nu=nu)

	
	if(log==TRUE){
		lf.x      <- log(part1) + log(part2) + log(part3)
		return(lf.x)
	}else{
		f.x      <- part1*part2*part3
		return(f.x)

		}	
}

my.outerbivlapl <- function(x1vals,x2vals, mu.vec,cov.mat, islog=FALSE){
	
	len1 <- length(x1vals)
	len2 <- length(x2vals)
	out  <- matrix(0,nrow=len1,ncol=len2)
	
	for(i in 1:len1){
		
		x1 <- x1vals[i]
		
		for(j in 1:len2){
			
			x2 <- x2vals[j]
			out[i,j] <- dbivlaplace(xvec=c(x1,x2),mu.vec=mu.vec,cov.mat=cov.mat, log=islog)	
		}
		
	}
	
	return(out)
}



GInt <- function(x,mus1,mus2,musg,Sigma1,Sigma2,Sigmag,pwr=1){
	
	 g.x <- dbivlaplace(xvec=x,mu.vec=musg,cov.mat=Sigmag,log=FALSE)
	 lf1.x<- dmvnorm(xvec=x,mu.vec=mus1,cov.mat=Sigma1,log=TRUE)
	 lf2.x<- dmvnorm(xvec=x,mu.vec=mus2,cov.mat=Sigma2, log=TRUE)	 
	 
	 return(g.x*(lf1.x-lf2.x)^pwr)
	 
}

F1Int <- function(x,mus1,mus2,Sigma1,Sigma2,pwr=1){
	
	 f1.x <- dmvnorm(xvec=x,mu.vec=mus1,cov.mat=Sigma1,log=FALSE)
	 lf1.x<- dmvnorm(xvec=x,mu.vec=mus1,cov.mat=Sigma1,log=TRUE)
	 lf2.x<- dmvnorm(xvec=x,mu.vec=mus2,cov.mat=Sigma2, log=TRUE)	 
	 
	 return(f1.x*(lf1.x-lf2.x)^pwr)
	 
}


F2Int <- function(x,mus1,mus2,Sigma1,Sigma2,pwr=1){
	
	 f2.x <- dmvnorm(xvec=x,mu.vec=mus2,cov.mat=Sigma2,log=FALSE)
	 lf1.x<- dmvnorm(xvec=x,mu.vec=mus1,cov.mat=Sigma1,log=TRUE)
	 lf2.x<- dmvnorm(xvec=x,mu.vec=mus2,cov.mat=Sigma2, log=TRUE)	 
	 
	 return(f2.x*(lf2.x-lf1.x)^pwr)
	 
}


NSigmaSets <- function(musg, mus1.0, mus2,Sigma1,Sigma2,Sigmag,angles=c(0, pi/8, (5*pi)/4, (15*pi)/8)){
	
	low.x1lim <- round(musg[1] - 3*sqrt(Sigmag[1,1]),digits=2)
	up.x1lim  <- round(musg[1] + 3*sqrt(Sigmag[1,1]),digits=2)+.01

	low.x2lim <- round(musg[2] - 3*sqrt(Sigmag[2,2]),digits=2)
	up.x2lim  <- round(musg[2] + 3*sqrt(Sigmag[2,2]),digits=2)+.01

	nangles <- length(angles)
	names.cases <- paste0(rep("case ",nangles), 1:nangles)

	n.f1means <- matrix(0,nrow=nangles, ncol=2)
	row.names(n.f1means) <- names.cases
	colnames(n.f1means) <- c("mu1", "mu2")

	n.f2means <- matrix(0,nrow=nangles, ncol=2)
	row.names(n.f2means) <- names.cases
	colnames(n.f2means) <- c("mu1", "mu2")

	hypoth <- sqrt(sum((mus1.0-musg)^2))
	settings.mat1 <- matrix(0,nrow=nangles,ncol=4)
	row.names(settings.mat1) <- names.cases
	colnames(settings.mat1) <- c("K12", "Delta K", "Sigsq1", "Sigsqg")	

	settings.mat2 <- matrix(0,nrow=nangles,ncol=4)
	row.names(settings.mat1) <- names.cases
	colnames(settings.mat2) <- c("K12", "Delta K", "Sigsq1", "Sigsqg")		

	# First we have to determine a new position for 'g':
	# this is because the position for model 2 will be modified according to the case
	musgp <- c(mus2[1]+ hypoth*cos(pi/4), mus2[2] + hypoth*sin(pi/4))


	for(i in 1:nangles){
		
		
		# subcase 1: f1 is closest to truth and f1 is the one being offset by an angle theta	
		if(i==1){mus11<-mus1.0[1];mus12 <- mus1.0[2];}else{
			theta <- angles[i]
			mus11 <- hypoth*cos(theta)
			mus12 <- hypoth*sin(theta)
		}	
		n.f1means[i,1] <- mus11
		n.f1means[i,2] <- mus12
		mus1 <- c(mus11,mus12)
		
		#just checking...
		#plot(rbind(mus1,mus2,musg)[,1], rbind(mus1,mus2,musg)[,2],pch=16)
		
		DeltaK <- hcubature(f=GInt,lowerLimit=c(low.x1lim,low.x2lim), upperLimit=c(up.x1lim,up.x2lim), 
		mus1=mus1, mus2=mus2, musg=musg, Sigma1=Sigma1, Sigma2=Sigma2, Sigmag=Sigmag, pwr=1)$integral;
		
		GIntsq  <- hcubature(f=GInt,lowerLimit=c(low.x1lim,low.x2lim), upperLimit=c(up.x1lim,up.x2lim), 
		mus1=mus1, mus2=mus2, musg=musg, Sigma1=Sigma1, Sigma2=Sigma2, Sigmag=Sigmag, pwr=2)$integral;
		
		Sigsqg <- GIntsq - DeltaK^2 

		K12 <- hcubature(f=F1Int,lowerLimit=c(low.x1lim,low.x2lim), upperLimit=c(up.x1lim,up.x2lim), 
		mus1=mus1, mus2=mus2, Sigma1=Sigma1, Sigma2=Sigma2,pwr=1)$integral; 

		F1Intsq <- hcubature(f=F1Int,lowerLimit=c(low.x1lim,low.x2lim), upperLimit=c(up.x1lim,up.x2lim), 
		mus1=mus1, mus2=mus2, Sigma1=Sigma1, Sigma2=Sigma2, pwr=2)$integral;

		Sigsq1 <- F1Intsq - K12^2 

		settings.mat1[i,1] <- K12;
		settings.mat1[i,2] <- DeltaK;
		settings.mat1[i,3] <- Sigsq1;
		settings.mat1[i,4] <- Sigsqg;

		# subcase 2: f2 is closest to truth and f2 is the one being offset by an angle theta	
		
		mus1 <- mus1.0
		if(i==1){mus21<-mus2[1];mus22 <- mus2[2];}else{
			theta <- angles[i]
			mus21 <- musgp[1]-hypoth*sin(theta)
			mus22 <- musgp[2]-hypoth*cos(theta)
		}	
		n.f2means[i,1] <- mus21
		n.f2means[i,2] <- mus22
		mus2 <- c(mus21,mus22)
		
 		#just checking...
		#plot(rbind(mus1,mus2,musgp)[,1], rbind(mus1,mus2,musgp)[,2],pch=16)

		DeltaK <- hcubature(f=GInt,lowerLimit=c(low.x1lim,low.x2lim), upperLimit=c(up.x1lim,up.x2lim), 
		mus1=mus1, mus2=mus2, musg=musgp, Sigma1=Sigma1, Sigma2=Sigma2, Sigmag=Sigmag, pwr=1)$integral;
		
		GIntsq  <- hcubature(f=GInt,lowerLimit=c(low.x1lim,low.x2lim), upperLimit=c(up.x1lim,up.x2lim), 
		mus1=mus1, mus2=mus2, musg=musgp, Sigma1=Sigma1, Sigma2=Sigma2, Sigmag=Sigmag, pwr=2)$integral;
		
		Sigsqgp <- GIntsq - DeltaK^2 

		K12 <- hcubature(f=F1Int,lowerLimit=c(low.x1lim,low.x2lim), upperLimit=c(up.x1lim,up.x2lim), 
		mus1=mus1, mus2=mus2, Sigma1=Sigma1, Sigma2=Sigma2,pwr=1)$integral; 

		F1Intsq <- hcubature(f=F1Int,lowerLimit=c(low.x1lim,low.x2lim), upperLimit=c(up.x1lim,up.x2lim), 
		mus1=mus1, mus2=mus2, Sigma1=Sigma1, Sigma2=Sigma2,pwr=2)$integral; 

		Sigsq1 <- F1Intsq - K12^2 

		settings.mat2[i,1] <- K12;
		settings.mat2[i,2] <- DeltaK;
		settings.mat2[i,3] <- Sigsq1;
		settings.mat2[i,4] <- Sigsqgp;


	} 
	
	return(list(n.f1means=n.f1means, settings.mat1=settings.mat1,n.f2means=n.f2means, settings.mat2=settings.mat2, musgp=musgp))	

}





Tsq.hotelling <- function(Xn.mat, mu.null,alpha){
	
	n <- nrow(Xn.mat)
	p <- ncol(Xn.mat)
	Xbar <- apply(Xn.mat,2,mean)
	S    <- var(Xn.mat)
	S.inv <- ginv(S)
	Tsq.obs <- n*t(Xbar - mu.null)%*%S.inv%*%(Xbar - mu.null)
	Tsq.theo <- (((n-1)*p)/(n-p))*qf(p=1-alpha,p,(n-p))
	if(Tsq.obs>Tsq.theo){reject<-TRUE}else{reject<-FALSE}
	return(list(Tsq.obs=Tsq.obs, Tsq.theo=Tsq.theo,reject=reject))
}

#source("sweatdata.R")
#Tsq.hotelling(Xn.mat=sweat.data,mu.null=c(4,50,10),alpha=0.10)


alpha.p <- function(nvec,sig.1,sig.g,K12,DeltaK,z.alpha){

	# z.alpha has to be a scalar	
	qtle.vec <- (sqrt(nvec)/sig.g)*(K12-DeltaK) - (sig.1/sig.g)*z.alpha
	return(pnorm(q=qtle.vec))
	
}

beta.p <- function(nvec,sig.1,sig.g,K12,DeltaK,z.alpha){
	
	# z.alpha has to be a scalar	
	qtle.vec <- (sqrt(nvec)/sig.g)*(K12-DeltaK) - (sig.1/sig.g)*z.alpha
	return(pnorm(q=(1-qtle.vec)))
	
}

Mi.p <- function(nvec,k,sig.g,DeltaK){
	
	qtle.vec <- -(sqrt(nvec)/sig.g)*((1/nvec)*log(k) + abs(DeltaK))
	return(pnorm(q=(qtle.vec)))
	
}


Wi.p  <- function(nvec,k,sig.g,DeltaK){

	qtle1.vec <- (sqrt(nvec)/sig.g)*((1/nvec)*log(k) - abs(DeltaK))	
	qtle2.vec <- -(sqrt(nvec)/sig.g)*((1/nvec)*log(k) + abs(DeltaK))
	return(pnorm(q=qtle1.vec) - pnorm(q=qtle2.vec))
	
}
