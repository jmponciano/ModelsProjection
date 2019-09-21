#Function to simulate the Negative Binomial (Poisson-Gamma) process

Females.sim <- function(simple.draw=FALSE,A,a,ad.freq,alpha,beta,nmales,hetero=TRUE,male.hetero=FALSE,plot.it=TRUE,decision.rule="rand.unif",weights=c(0.80,0.20),
pfail=0.20,follow.ad=TRUE,non.unif.space=TRUE, shape1=2,shape2=7){
	
	# nmales = average number of arriving males per female per unit area

	# 3. For each male:
 	#    If there are 1 or more signaling females in its entourage, choose one of three rules to decide 
 	#    to which female the male is attaching.
 	#    Rule 1: choose a signaling female at random 
 	#    Rule 2: choose the closest female
 	#    Rule 3: choose female at random, but the probability of joining every female is weighted by the 
 	#            number of dudes attached to her. 1 or 0 have same (high) chances.  
 	#            2 or more should have a lower prob.
	pjoin <- 1-pfail;
	k<- 2 # number of types of females, it's 2 always
	inv.beta <- beta^(-1);
  	# Setting the average number of females per unit area
  	lams.const <- c(alpha/beta, (alpha/beta)*ad.freq)*A;
  	lams.rand  <- rgamma(n=k,shape=c(alpha, alpha*ad.freq)*A,scale=rep(inv.beta,2));
  	if(hetero==FALSE){
  		lambdas <- lams.const}else{
  		  	lambdas <- lams.rand;
	  }
	
	# Drawing at random the total number of males using an average abundance (nmales)
	
	if(male.hetero==FALSE){
		
		ave.males <- sum(lams.const)*nmales;
		}else if(male.hetero==TRUE){
		ave.males <- sum(lams.rand)*nmales;
	}
	qu <-rpois(n=1,lambda=ave.males);		
	
	# Drawing the number of females at random 
	pois.draw <- rep(0,k);
	for(i in 1:k){pois.draw[i]<- rpois(n=1,lambda=lambdas[i]);}
	nis <- pois.draw; # totals per species drawn randomly
 	Ntot <- sum(nis);     # total number of individuals in A

	Tot.males <- Ntot+qu;
	Tot.females <- Ntot
	Tot.dens <- Tot.males+Tot.females
	realized.osr <- Tot.males/Tot.females

	
	if(simple.draw==TRUE){
		
		return(list(Tot.males=Tot.males,Tot.females=Tot.females,Tot.dens=Tot.dens,realized.osr=realized.osr))
	}else if(simple.draw==FALSE){

		# Placing the females at uniformly random
	 	if(non.unif.space==TRUE){
 			xs <- rbeta(n=Ntot,shape1=shape1,shape2=shape2)*A 
 		}else{
 		xs<- runif(n=Ntot,min=0,max=A);
 		}
 		ys <- runif(n=Ntot,min=0,max=1);
 		cls <- rainbow(k);
 		marks <- rep(cls[1:k],nis);
 		spp.numtag <- rep(1:k,nis);
 	 
		# Plotting the females in space 
 		if (plot.it == T)
 		{
   		par(mai=c(0.85,0.85,0.35,0.70),oma=c(1.45,1.25,0.5,0.75));
   		plot(xs,ys,type="n", main="Study Area",xlab="West-East distance (kms)", 
   		ylab="North-South distance (kms)", cex=1.5,cex.lab=1.5, cex.axis=1.25,cex.main=1.5);
   		for(i in 1:Ntot)
   		{
			points(xs[i],ys[i],pch=20,col=marks[i]);
   		}
 	}

		fem.supclass <-list()
	   	for(i in 1:Ntot)
   			{
   			ith.fem <- list();
   			coords <- c(xs[i],ys[i]);
   			is.she.adv <- spp.numtag[i];
   			dudes.attached <- 0;# list of the male identifiers attached to this female. 
   								# length(dudes.attached)=num.males attached
   			which.qu.ishe <- 0;
   			ith.fem[[1]] <- coords;
   			ith.fem[[2]] <- is.she.adv;
   			ith.fem[[3]] <- dudes.attached;
   			ith.fem[[4]] <- which.qu.ishe
   			names(ith.fem) <- c("coords","is.she.adv","dudes.attached","which.qu.is.she");
   			#print(ith.fem)
   			fem.supclass[[i]] <- ith.fem;
   			
   		}

		mal.supclass <- list()
		for(i in 1:qu){

   			ith.male <- list();
   			coords.square <- matrix(c(0,0,0,0),nrow=1,ncol=4);
   			is.dude.attached <- 0; # 0=no, 1=first female, 2=second female in quadrat, and so on
   			coords.girls <- matrix(c(0,0),nrow=1,ncol=2);
   			id.girls <- 0;
   			ith.male[[1]] <- coords.square;
   			ith.male[[2]] <- is.dude.attached;
   			ith.male[[3]] <- coords.girls;
			ith.male[[4]] <- id.girls
			names(ith.male) <- c("coords.square","is.dude.attached","coords.girls","id.girls")
   			mal.supclass[[i]] <- ith.male
		
	}
 	
	 	# Locating males arriving on the beach at random along the horizontal axis
		# 	Done by locating at random the lower left corners of the quadrats
 	
	 	if(follow.ad==TRUE){
 			xs.hist <- hist(xs[(nis[2]+1):Ntot],plot=FALSE)
 			rel.freqs <- xs.hist$counts/sum(xs.hist$counts)
 			mid.xs <- xs.hist$mids
 			band.corner <- runif(n=qu,min=-50,max=49)#runif(n=qu,min=-50,max=49)
 			band.location <- sample(mid.xs,size=qu,replace=TRUE, prob=rel.freqs)
 			low.left.corners <- band.corner+band.location
 		}else{low.left.corners <- runif(n=qu,min=1,max=(A-1))}
 	
		n.sampled <- matrix(rep(0,qu*k),nrow=qu,ncol=k);
 		quadrats.xs <- matrix(rep(0,4*qu),nrow=qu,ncol=4);
 		quadrats.ys <- quadrats.xs;

		#Attributes females: 1. coords, 2. is.she.adv, 3. dudes.attached, 4.which.qu.ishe 
		#Attributes males: 1, coords.square, 2. is.dude.attached, 3. coords.girl, 4. id.girls

			for(i in 1:qu){

				llc <- low.left.corners[i];
				x.left<- llc;
				y.bott<- 0;
				x.right<- llc+a;
				y.top  <- 1;
				x.mid <- llc+(a/2)
				mal.supclass[[i]][[1]][1,] <- c(x.left, y.bott, x.right, y.top)
			
				arefem.adv   <- rep(0,Ntot)
				femmes.inbox <- rep(0,Ntot)
				dist2fems    <- rep(0,Ntot)
				
				for(j in 1:Ntot){

					# For every female, check if it falls within one of the quadrats
 					# Identifying the females within a 2x1 square around each male and keep track of 
 					# their coordinates
			
					is.inbox <- ((xs[j]>=x.left)&&(xs[j]<=x.right))&&((ys[j]>=y.bott)&&(ys[j]<=y.top));
					femmes.inbox[j] <- is.inbox
		       		arefem.adv[j] <- spp.numtag[j]==2; # if TRUE, female is advertising
		       		dist2fems[j] <- abs(xs[j]-x.mid)
		       		
				}# end of first j-for loop over all females for each male's quadrat
				
				dist2fems.ad <- dist2fems*(spp.numtag==2)
				dist2fems.ad[dist2fems.ad==0]<-A
					
				num.closest.fem <- which(dist2fems.ad==min(dist2fems.ad),arr.ind=TRUE)
				#tot.inbox<- sum(femmes.inbox)
				tot.advinbox <- sum(arefem.adv*femmes.inbox)
				
				# If there's no adv. in quadrat, go to closest female and set quadrat there
				if((tot.advinbox==0)&(follow.ad==TRUE)){
					
					x.left  <- xs[num.closest.fem]-(a/2)
					x.right <- xs[num.closest.fem]+(a/2)
					mal.supclass[[i]][[1]][1,] <- c(x.left, y.bott, x.right, y.top)

					}

	      		if(plot.it==T)
					{  
	  					rect(xleft=x.left,ybottom=y.bott,xright =x.right,ytop=y.top,
	  					density=NA, col=rgb(0.3,0.2,0.1,alpha=0.1), border="blue");
					}


				# Now loop over the females again and save characteristics in the male's class attributes
				for(j in 1:Ntot){
	   			
	   				is.inbox2 <- ((xs[j]>=x.left)&&(xs[j]<=x.right))&&((ys[j]>=y.bott)&&(ys[j]<=y.top));
	   				
	   					if(is.inbox2==T){	
				       		fem.coords <- c(xs[j],ys[j]);
		    			   		fem.supclass[[j]][[4]] <- c(fem.supclass[[j]]$which.qu.is.she,i) 
							mal.supclass[[i]]$coords.girls <- rbind(mal.supclass[[i]]$coords.girls,fem.coords);
							mal.supclass[[i]]$id.girls <- c(mal.supclass[[i]]$id.girls,j);
							# Just counting females in quadrats, and of which type...
							which.spp  <- spp.numtag[j];
							n.sampled[i,which.spp] <- n.sampled[i,which.spp] + 1;

	       			
    	   			}# end if it's in box	

				}# end of second j-for loop over all females for each male's quadrat
			
				
			}
			
			# Now that we've given a chance to males to re-locate according to whether they 
			# landed close to an advertising female or not, we do the attachment decision tree			
			
			for(i in 1:qu){
				# Build decision tree for attachment depending on number of advertising females present
				# and on number of males attached to each female.  Consider the decision rule too...
				# Case 1: no adv females, Case 2: exactly 1 adv. female, Case 3, more than 1 adv. female
				# Case 1

				# First, we count num. of adv. females in quadrat
				nadv.fem <- n.sampled[i,2];
				#	- from the list of females in quadrat, gather ids
				local.femmes.ids <- mal.supclass[[i]]$id.girls[-1];
				are.local.adv <- spp.numtag[local.femmes.ids]
				index.ofadv <- which(are.local.adv==2,arr.ind=TRUE);
				ids.ofadv <- local.femmes.ids[index.ofadv];
				
				
				U2 <- runif(n=1)
				if(nadv.fem==0){
				
					mal.supclass[[i]]$is.dude.attached <- 0;
			
				}else if(nadv.fem==1){
				# Case 2	
				
				
					# which female is it?
					#   - from these ids, get the only female that is advertising
					the.one <- ids.ofadv[1];
					# count the number of males already attached to this female
					ndudes.there <- length(fem.supclass[[the.one]]$dudes.attached[-1])
					if(ndudes.there<=1){weight <- weights[1]}else{weight<- weights[2]}
					# Decide according to rule whether to join this female	
					if(decision.rule=="rand.unif"){
	
							#Since only one, just join with prob pjoin!
							if(U2<=pjoin){
								mal.supclass[[i]]$is.dude.attached <- 1;
								fem.supclass[[the.one]]$dudes.attached <- c(fem.supclass[[the.one]]$dudes.attached,i) 
							}			
						}else if(decision.rule=="weighted"){
					
						u <- runif(n=1);
						#print(u)
						if(u<=weight){
							
							if(U2<=pjoin){
								#attach it!
								mal.supclass[[i]]$is.dude.attached <- 1;
								fem.supclass[[the.one]]$dudes.attached <- c(fem.supclass[[the.one]]$dudes.attached,i) 
								}
							}
						# end if "weighted" ############# ADD new decision here
						}else if(decision.rule=="Allee"){
							
						nsatels <- ndudes.there
						A.weight <- weight+((1-weight)/(nsatels+1))*nsatels
						u <- runif(n=1);
						#print(u)
						if(u<=A.weight){
							
							if(U2<=pjoin){
								#attach it!
								mal.supclass[[i]]$is.dude.attached <- 1;
								fem.supclass[[the.one]]$dudes.attached <- c(fem.supclass[[the.one]]$dudes.attached,i) 
								}
							}
							
						}else{
							
							if(U2<=pjoin){
								#Nearest: since only one, just join!
								mal.supclass[[i]]$is.dude.attached <- 1;
								fem.supclass[[the.one]]$dudes.attached <- c(fem.supclass[[the.one]]$dudes.attached,i) 
							}
					}# End attaching decision tree for Case 2
					
					
				}else if(nadv.fem>1){
				
				# Case 3	
					if(decision.rule=="nearest"){	
						coords.advert <- matrix(0,nrow= nadv.fem,ncol=2)
						dist.2male <- rep(0,nadv.fem)
						xm <- x.left + (1/2)*a
						ym <- 0
					
						for(h in 1:nadv.fem){
							h.coords <- fem.supclass[[ids.ofadv[h]]]$coords
							coords.advert[h,] <- h.coords;
							dist.2male[h] <- sqrt(sum((c(xm,ym)-h.coords)^2))
						}
				
						loc.smallest <- which(dist.2male==min(dist.2male),arr.ind=TRUE)
						the.one <- ids.ofadv[loc.smallest]	
						
						if(U2<=pjoin){
							#attach it!
							mal.supclass[[i]]$is.dude.attached <- 1;
							fem.supclass[[the.one]]$dudes.attached <- c(fem.supclass[[the.one]]$dudes.attached,i) 
						}
					}else if(decision.rule=="rand.unif"){
					
						rand.pick <- sample(x=1:nadv.fem, size=1,prob=rep((1/nadv.fem),nadv.fem))
						the.one <- ids.ofadv[rand.pick];
						
						if(U2<=pjoin){	
							#attach it!
							mal.supclass[[i]]$is.dude.attached <- 1;
							fem.supclass[[the.one]]$dudes.attached <- c(fem.supclass[[the.one]]$dudes.attached,i) 
						}
						
					}else if(decision.rule=="weighted"){
						
						weights0.vec <- rep(0,nadv.fem)
						for(h in 1:nadv.fem){
							
							this.one <- ids.ofadv[h]
							ndudes.there <- length(fem.supclass[[this.one]]$dudes.attached[-1])
							if(ndudes.there<=1){weights0.vec[h] <- weights[1]}else{weights0.vec[h]<- weights[2]}
						}
						weights.vec <- weights0.vec/sum(weights0.vec)					
						
						rand.pick <- sample(x=1:nadv.fem, size=1,prob=weights.vec)
						the.one <- ids.ofadv[rand.pick];
						if(U2<=pjoin){
							mal.supclass[[i]]$is.dude.attached <- 1;
							fem.supclass[[the.one]]$dudes.attached <- c(fem.supclass[[the.one]]$dudes.attached,i) 
						}
					}else if(decision.rule=="Allee"){
						
						weights0.vec <- rep(0,nadv.fem)
						for(h in 1:nadv.fem){
							
							this.one <- ids.ofadv[h]
							ndudes.there <- length(fem.supclass[[this.one]]$dudes.attached[-1])
							if(ndudes.there<=1){weights0.vec[h] <- weights[1]}else{
								
								A.weight <- weights[2]+((1-weights[2])/(ndudes.there+1))*ndudes.there
								weights0.vec[h]<- A.weight}
						}
						weights.vec <- weights0.vec/sum(weights0.vec)					
						
						rand.pick <- sample(x=1:nadv.fem, size=1,prob=weights.vec)
						the.one <- ids.ofadv[rand.pick];
						if(U2<=pjoin){
							mal.supclass[[i]]$is.dude.attached <- 1;
							fem.supclass[[the.one]]$dudes.attached <- c(fem.supclass[[the.one]]$dudes.attached,i) 
						}
						
						
					}
				
				
				}# End "else if" section for Case 3 (more than 1 adv.female in quadrat)
			
	
  		}# end of for loop over all males' quadrats
 	
	 	#   return the number of males attached to each female.
		nmales.perfemme <- rep(0,Ntot)
		for(ell in 1:Ntot){
		
		dudes.there <- length(fem.supclass[[ell]]$dudes.attached[-1])
		nmales.perfemme[ell] <- dudes.there		
	}	

	
	
		#Calculate Nq's = number of adv. females per quadrat
		Nadfem.perq <- n.sampled[,2];
	
		result <- list(n.sampled=n.sampled, Nadfem.perq=Nadfem.perq, lambdas=lambdas, mal.supclass=mal.supclass, fem.supclass=fem.supclass,nmales.perfemme=nmales.perfemme,Tot.males=Tot.males,Tot.females=Tot.females,Tot.dens=Tot.dens,realized.osr=realized.osr);
		return(result);
	}


} #end of the function



#-----------

# Simulating an observation
# Decision rule can be "rand.unif" or "nearest" or "weighted"
#sim.sa <- Females.sim(A=1000,a=2,ad.freq= 2,alpha=0.15,beta=(1/2),nmales=1/3, hetero=FALSE,plot.it=TRUE, decision.rule="nearest",weights=c(0.90,0.10));


