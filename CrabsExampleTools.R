source("AbundanceToolkit2.0.R")
source("NegbinSim.R")


###  Functions for the crabs example


# Function to compute the expected value of the average number of females per squared meter
# That expected value is equal to alpha/beta.
# The output is then the fraction of 'alpha/beta' needed to achieve a targeted density,  
# advertising frequency and a mean osr.  This is an auxiliary function to the function 
# beta.calc below. 
# Explanation from Appendix 1 in Brockmann et al 2018:
# 1. Targeted density = Number of animals in 1000 squared meters of beach (1000 mts x 1 mt)
# 2. advertising frequency 'd': define as theaverage fraction of advertising females with  
#    respect to average non-advertising females.  Specifically, if we define 
# 						Nf1 = number of non-advertising females,
#						Nf2 = number of advertising females, then
# 						Nf = Nf1 + Nf2 = total number of females at the beach.
#
#
#                       E[Nf1] = alpha/beta
#                       E[Nf2] = (alpha/beta)*d, where d= advertising frequency 
# 						Example: d=1/2 means that the number of advertising females is
#                       half the number of non-advertising females on average.	
# 3. mean osr, mu_s = average operational sex ratio. Defined as 
#   E[(Number of males)/(Number of females)] = E[Nm/Nf] = mu_s.
#   The simulations assume that E[Nm] = E[Nf]*(1+o), where o is a multiplier. Then, from 
#   a T.S. expansion we get that to a first order approx  E[Nm/Nf] = mu_s ~ (1+o).
#   Examples:  o = -1/2, mu_s=1+o = 1/2, so males are half as abundant as females
#              o = 1/2, mu_s=1+o = 3/2, and males are 1.5 times more abundant than females
#              o = 1, mu_s = 1+o = 2, hence there are 2 males per female, on average
#              o = 2, mu_s = 1+o = 3, hence there are 3 males per female, on average  

alpha.beta.calc <- function(target.dens,ad.freq,osr){
	
	return(target.dens/((1+osr)*(1+ad.freq)) )
	
}

# Function to compute the value of the parameter 'beta' (see description of 'alpha.beta.calc()' above)
# given a value of 'alpha' and a total beach area A
# Because (alpha/beta) = target.dens/((1+osr)*(1+ad.freq)),
# one can solve for the beta value given a value of alpha, a required target density, osr and ad.freq

beta.calc <- function(A,alpha.val,single.targetdens,single.adfreq,single.osr){

	albetA <- alpha.beta.calc(target.dens=single.targetdens,ad.freq=single.adfreq, osr=single.osr)
	single.beta <- alpha.val/(albetA/A)
	return(single.beta)
}
# Sample run
#beta.calc(A=1000,alpha.val=0.32,single.targetdens=3000,single.adfreq=0.5,single.osr=2)

#### Short Simulator:
short.sim <- function(adv.freq,osr,beta,a,nsims=10,alpha=0.32,area=1000,hetero.abund=FALSE,male.hetero=FALSE,
					  decision.mode="rand.unif", tracking=FALSE,nonunif=FALSE, ncateg=8, 
					  my.models = c("Poisson", "NegBin", "ZIPoiss", "ZINegBi", "HurdNBi","PoisNB", "NBPois", "OIPoiss", "OINegBi")){
	
	
	
	####### Settings:
	# a = width of male crab vision in meters
	#area <- 1000 # Total study area in squared meters
	#alpha <- 0.32 # Proportional to mean abundance per squared meter
	#####  Set the next two options to TRUE if you want 'Pairs clustered vertical and horizontal' as in Tables 1&2
	#hetero.abund <- FALSE # Are FEMALE in pairs clustered in space? T/F
	#nonunif <- FALSE # Do females land heterogeneously in space? T/F
	
	#####  Are the unpaired males clustered in space?
	#male.hetero <- FALSE #long.mal.het[ell] # Do MALE crabs arrive in clusters at the beach? T/F
	
	##### Unpaired male joining rule
	#decision.mode <- "rand.unif"#"Allee"#"weighted"# Decision rule to pick a female within vision quadrat: 
							   # can be "rand.unif", "nearest", "weighted" or "Allee" (4 choices)
						       # The "weighted" option weights according to number of males already present  
	
	#####  Do males track females?
	#tracking <- FALSE # Do males land more often in areas with more adv. females? T/F


	quad.width <- a # width of male crab vision in meters
	nmales.pfem <- osr-1 #*** number  of males per pair (previously said 'per 'female')
	do.plot <- FALSE # Plot simulation?
	probs.perfemt <- c(0.95,0.75)#c(0.75,0.95) #c(0.75,0.99) #c(0.95,0.75)  
						# **** Probability of joining if 1 or less males 
						# already present (first element), and if 2 or more are present...
	pfail <- 0.001 # *** Prob(failing to join even if "available")
	shape1  <- 2 # beta distrib. param. 1 to set up females in space if non-homogenous space assumed
	shape2 <- 5 # beta distrib. param. 2 to set up females in space if non-homogenous space assumed


	
	nmodels <- length(my.models)
	bics.mat <- matrix(0,nrow=nsims,ncol=nmodels)
	colnames(bics.mat) <- my.models
	counts <- list()
	fits.list <- list()
	ncols.sims <- rep(0,nsims)
	
	for(i in 1:nsims){

		one.sim <- Females.sim(A=area,a=quad.width,ad.freq= adv.freq ,alpha=alpha,beta=beta,nmales=nmales.pfem, 
		hetero=hetero.abund,male.hetero=male.hetero, plot.it=do.plot, decision.rule=decision.mode,
		weights=probs.perfemt,pfail=pfail,follow.ad=tracking,non.unif.space=nonunif, shape1=shape1,shape2=shape2)
	
		stats.per.fem <- one.sim$nmales.perfemme
		# Compute frequencies of counts!!
		tabled.counts <- table(stats.per.fem)
		sum.counts <- sum(tabled.counts)
		len.counts <- length(tabled.counts)

		if((sum.counts==0)|(len.counts==1)){
			one.sim <- Females.sim(A=area,a=quad.width,ad.freq= adv.freq ,alpha=alpha,beta=beta,nmales=nmales.pfem, 
			hetero=hetero.abund,male.hetero=male.hetero, plot.it=do.plot, decision.rule=decision.mode,
			weights=probs.perfemt,pfail=pfail,follow.ad=TRUE,non.unif.space=nonunif, shape1=shape1,shape2=shape2)
			stats.per.fem <- one.sim$nmales.perfemme
		}
	
		# Compute frequencies of counts!!
		tabled.counts <- table(stats.per.fem)
		#print(tabled.counts)
		x.vec <- as.numeric(names(tabled.counts))
		old.len <- length(tabled.counts)
		contig.xs <- 0:(length(x.vec)-1)
		mismatch <- x.vec!=contig.xs
		first0 <- which(mismatch==TRUE,arr.ind=TRUE)
		if(length(first0)>0){ 
			tail.sum <- sum(tabled.counts[first0[1]:old.len])
			tabled.counts <- c(tabled.counts[1:(first0[1]-1)],tail.sum)
			new.len <- length(tabled.counts)
			x.vec <- 0:(new.len-1)
		}
		new.dat <- rbind(x.vec,tabled.counts)
		colnames(new.dat) <- as.character(x.vec)
		ncols.sims[i] <- ncol(new.dat)
		
		
		#print(new.dat)
		#print(paste("Tot.dens = ",one.sim$Tot.dens))
		#print(paste("Realized OSR = ", one.sim$realized.osr))
		print(i)
		counts[[i]] <- new.dat
	
	}
	
	min.k.sims <- min(ncols.sims)
	k <- min(ncateg, min.k.sims)
	km1 <- k-1
	simdat <- matrix(0,nrow=nsims,ncol=k)
	colnames(simdat) <- as.character(0:km1)
	code.string <- paste0(rep("ith.fit$",nmodels),my.models,rep(".Stats$bic",nmodels))
	
	
	for(i in 1:nsims){
		
		new.dat <- counts[[i]]
		len.dat <- ncol(new.dat)
		tail.sum <- sum(new.dat[2,k:len.dat])  
		short.counts <- c(new.dat[2,1:km1],tail.sum)
		simdat[i,] <- short.counts

		ith.fit <- abund.fit(y1=simdat[i,], names.animals="Horse-shoe crabs", 
		name.counts="Satellite males", pois.method="BFGS", plot.it=FALSE)
		fits.list[[i]] <- ith.fit
	
		bics.ithrow <- rep(0,nmodels)
		for(j in 1:nmodels){bics.ithrow[j] <- eval(parse(text=code.string[j]))} 	
		bics.mat[i,] <- bics.ithrow
		#c(ith.fit$Poisson.Stats$bic, ith.fit$NegBin.Stats$bic,ith.fit$ZIPoiss.Stats$bic, ith.fit$ZINegBi.Stats$bic, ith.fit$HurdNBi.Stats$bic)
	}
	
	#print(bics.mat)
	return(list(bics.mat=bics.mat, fits.list=fits.list, counts=counts, simdat=simdat))
}

fits.list.fn <- function(simdat,fitted.models=c("Poisson", "NegBin", "ZIPoiss", "ZINegBi", "HurdNBi")){
	
	nreps <- dim(simdat)[1]	
	nmodels <- length(fitted.models)
	bics.mat <- matrix(0,nrow=nreps,ncol=nmodels)
	colnames(bics.mat) <- my.models

	allfits.list <- list()
	
	code.string <- paste0(rep("ith.fit$",nmodels),fitted.models,rep(".Stats$bic",nmodels))
	
	
	for(i in 1:nreps){
		
		ith.fit <- abund.fit(y1=simdat[i,], names.animals="Horse-shoe crabs", 
		name.counts="Satellite males", pois.method="BFGS", plot.it=FALSE)
		allfits.list[[i]] <- ith.fit
		bics.mat[i,] <- eval(parse(text=code.string))
	}

	return(list(bics.mat=bics.mat, fits.list=allfits.list, simdat=simdat))
}




pred.pis <- function(fits.list, simdat){
	
	nvec <- apply(simdat,1,sum)
	ndat <- nrow(simdat)
	ncateg <- ncol(simdat)
	
	Pois.pis <- matrix(0,nrow=ndat,ncol=ncateg)
	NB.pis <- matrix(0,nrow=ndat,ncol=ncateg)
	ZIP.pis <- matrix(0,nrow=ndat,ncol=ncateg)	
	ZNB.pis <- matrix(0,nrow=ndat,ncol=ncateg)
	HNB.pis <- matrix(0,nrow=ndat,ncol=ncateg)		
	PoisNB.pis <- matrix(0,nrow=ndat, ncol=ncateg)
	NBPois.pis <- matrix(0,nrow=ndat,ncol=ncateg)
	OIPoiss.pis <- matrix(0,nrow=ndat,ncol=ncateg)
	OINegBi.pis <- matrix(0,nrow=ndat,ncol=ncateg) 
	
	for(i in 1:ndat){
		
		Pois.pis[i,] <- fits.list[[i]]$Poisson.Stats$Expected.freqs/nvec[i]
		NB.pis[i,] <- fits.list[[i]]$NegBin.Stats$Expected.freqs/nvec[i]
		ZIP.pis[i,] <- fits.list[[i]]$ZIPoiss.Stats$Expected.freqs/nvec[i]
		ZNB.pis[i,] <- fits.list[[i]]$ZINegBi.Stats$Expected.freqs/nvec[i]
		HNB.pis[i,] <- fits.list[[i]]$HurdNBi.Stats$Expected.freqs/nvec[i]
		PoisNB.pis[i,] <- fits.list[[i]]$PoisNB.Stats$Expected.freqs/nvec[i]
		NBPois.pis[i,] <- fits.list[[i]]$NBPois.Stats$Expected.freqs/nvec[i]
		OIPoiss.pis[i,] <- fits.list[[i]]$OIPoiss.Stats$Expected.freqs/nvec[i]
		OINegBi.pis[i,] <- fits.list[[i]]$OINegBi.Stats$Expected.freqs/nvec[i] 

		
			}
	
	return(list(Pois.pis = Pois.pis, NB.pis=NB.pis, ZIP.pis=ZIP.pis, ZNB.pis=ZNB.pis, HNB.pis=HNB.pis,
		   PoisNB.pis = PoisNB.pis, NBPois.pis = NBPois.pis, OIPoiss.pis = OIPoiss.pis, OINegBi.pis = OINegBi.pis))
}

#trial.list1 <- list(size=nvec[1], pis.g=Pois.pis[1,], type="Sgg", pis.f=NULL)
#trial.list2 <- list(size=nvec[1], pis.g=Pois.pis[1,], type="Sgf", pis.f=NB.pis[1,])
#H.multinom.loop(par.lst=trial.list1)

new.pis2add <- function(newfits.list,simdat,newdists.names){
	
	nvec <- apply(simdat,1,sum)
	ndat <- nrow(simdat)
	ncateg <- ncol(simdat)

	nnew.dists <- length(newdists.names)
	
	#creating as many empty matrices with pis as new models there are
	for(i in 1:nnew.dists){
		code.string <- paste0(newdists.names[i],".pis <- matrix(0,nrow=ndat,ncol=ncateg)")
		eval(parse(text=code.string))
	}
	
	for(j in 1:ndat){
		for(i in 1:nnew.dists){
			code.string <- paste0(newdists.names[i],".pis[j,] <- newfits.list[[j]]$",newdists.names[i],".Stats$Expected.freqs/nvec[j]")
			eval(parse(text=code.string))
		}
	}
	
	code.string2 <- paste0(newdists.names,rep(".pis = ",nnew.dists), newdists.names,rep(".pis",nnew.dists))
	
	outlist <- list()
	names.out <- rep("NA",nnew.dists)
	
	for(i in 1:nnew.dists){
		
		outlist[[i]] <- eval(parse(text=code.string2[i]))
		names.out[i] <- paste0(newdists.names[i],".pis")
	}
	names(outlist) <- names.out
	
	return(outlist)		
	
}

aics.calc <- function(fits.list, my.models=c("Poisson", "NegBin", "ZIPoiss", "ZINegBi", "HurdNBi", "PoisNB", "NBPois", "OIPoiss", "OINegBi"),nparms=c(1,2,2,3,3,3,3,2,3)){
	
	N <- length(fits.list)
	nmods <- length(fits.list[[1]])
	
	aics.mat <- matrix(0,nrow=N,ncol=nmods)
	colnames(aics.mat) <- my.models
	
	for(i in 1:N){
		
		stats.minilist <- fits.list[[i]]
		for(j in 1:nmods){
			
			loglike <- stats.minilist[[j]]$logLhat[1]
			npars   <- nparms[j]
			jth.aic <- -2*loglike + 2*npars
			aics.mat[i,j] <- jth.aic
		}	
	}
	
	return(aics.mat)
	
}

