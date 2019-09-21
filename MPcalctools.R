library(ape)
library(phangorn)
library(combinat)
library(progress)
library(smacof)
library(dfoptim)
library(plyr)
library(pracma)
library(seqinr)
library(scatterplot3d)
library(plot3D)
library(mvtnorm)


# Model Projection Tools
# comment for later: use doFuture and foreach packages for parallel processing

# Analytical Neg entropy for a Multivariate Normal (MVN) distribution
H.MVnorm <- function(par.lst=list(mu=c(0,0),sigma=matrix(c(1,0,0,1),nrow=2))) {
  sigma=par.lst$sigma
  d=dim(sigma)[1]
  out <- log(((2*pi*exp(1))^d)*det(sigma))/2
  return(out)
}

# Auxiliary function for the function Hse.wKL which computes Berrett et al's estimator:
Find.wkdN <- function(X,k=0){

  X=as.matrix(X)
  if (dim(X)[1]<dim(X)[2]) { #checking that observations are rows
    X <- t(X)
  }
  numobs <- dim(X)[1]
  d <- dim(X)[2]
  if(k==0){k <- floor(numobs^(1/3))}
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
  if (d<=7){w<-w;k=k}else if(d>7){
  	
  	# then optimize the weights according to Berrett's email
  	# setting k = n^(1/3)
  	
	ell.vec <- 1:floor(d/4)
  	ws <- w
	constr.f <- function(ws,ell.vec,d,k){
		nells <- length(ell.vec)
		cond1 <- all.equal(target=sum(ws), current=1)
		Gamma.test <- rep(0,nells)
		for(ell in ell.vec){
			jsum.elems <- rep(0,k)
			
			for(j in 1:k){jsum.elems[j] <- ws[j]*gamma((j+(2*ell)/d))/gamma(j)}
				Gamma.test[ell] <- (sum(jsum.elems)-0)^2
		}
		tot.gtest <- sum(Gamma.test)		
		if(cond1==1){
			eucl.norm <- sqrt(sum(ws^2))
			ofn <- eucl.norm+tot.gtest
			return(ofn)
		}else{return(.Machine$double.xmax)}
	}
	
	constr.f(ws=ws,ell.vec=ell.vec,d=d,k=k)
	best.ws <- optim(par=ws, fn=constr.f,method="BFGS", ell.vec=ell.vec,d=d,k=k)
	best.ws <- optim(par=best.ws$par, fn=constr.f,method="Nelder-Mead", ell.vec=ell.vec,d=d,k=k)	  		
	best.ws <- optim(par=best.ws$par, fn=constr.f,method="BFGS", ell.vec=ell.vec,d=d,k=k)	    	
  	w <- best.ws$par
  	}

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
  # X is a matrix of multivariate observations with rows being observations
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



### A small function to draw 'n' colours:
mycols.ftn <- function(n){
	d <- 360/n
	h <- cumsum(c(15,rep(d,(n-1))))
	return(hcl(h=h, c=100,l=65))
}


find_mg_xyz=function(mgKLvec,points){
  # Points is a matrix of coordinates of the models from the (re-scaled) NMDS
  # this function assumes that the points matrix has colnames D1..Ddim(NMDS)
  # projects generating model onto nmds space 
  # mgKLvec is a vector of KL distances to generating model length dim(points)[1]
  # 'z' is 'h'
  # What you want is the element xyz of the list 'out':
  # these are the coordinates of 'm' in the model projection space.
  
  npts=dim(points)[[1]]
  dmnz=dim(points)[[2]]
  xyz=rep(0,(dmnz+1))
  xyz[1:dmnz]=aaply(.data = points,.margins = 2,.fun = mean)
  names(xyz) <- c(colnames(points),"z")
  notop <- list(points=points,mgKLvec=mgKLvec,scale.vec=rep(1,length(xyz)))
  lossxyz=function(xyz,notop){
    
    len <- length(xyz)
    points=notop$points
    mgKLvec=notop$mgKLvec
    xyz <- xyz/notop$scale.vec
    z=xyz[len]
    #    d=pdist2(X = points,Y = xyz[1:2])@dist
    d=pdist2(X = points,Y = xyz[1:(len-1)])  #using pdist2 from pracma
    out=mean((mgKLvec - d^2 - z^2)^2) # removing ^2 from mgKLvec
    return(out)
  }
  bst.xyz=nmk(par = xyz,fn = lossxyz,notop=notop)
  scale.vec <- 10/bst.xyz$par
  notop <- list(points=points,mgKLvec=mgKLvec,scale.vec=scale.vec)
  bst.xyz=nmk(par = bst.xyz$par,fn = lossxyz,notop=notop)  
  scale.vec <- 10/bst.xyz$par
  notop <- list(points=points,mgKLvec=mgKLvec,scale.vec=scale.vec)
  bst.xyz=nmk(par = bst.xyz$par,fn = lossxyz,notop=notop)  
  scale.vec <- 10/bst.xyz$par
  notop <- list(points=points,mgKLvec=mgKLvec,scale.vec=scale.vec)
  bst.xyz=nmk(par = bst.xyz$par,fn = lossxyz,notop=notop)  
  names(bst.xyz$par) <- c(colnames(points),"z")
  xyz <- bst.xyz$par/scale.vec
  out <- list(xyz=xyz,bst.xyz=bst.xyz,notop=notop)
  
  return(out)
}

####  Function to compute the multinomial self entropy
H.multinom.loop <- function(par.lst=list(size=15,pis.g=c(.1,.3,.6),type="Sgg",pis.f=NULL)) {
  # https://stats.stackexchange.com/questions/207893/entropy-of-the-multinomial-distribution/208075#208075
  # function will tolerate 0s in prob vector
  # Returns -Sgg or -Sgf
  
  n <- par.lst$size
  p1 <- par.lst$pis.g
  p2 <- p1
  if(par.lst$type=="Sgf"){p2 <- par.lst$pis.f}
  k <- length(p1)
  piece1 <- -lfactorial(n)
  piece2 <- -n*sum(p1*log(p2),na.rm = T) #0*log(0) is NaN should be 0 so na.rm
  piece3 <-0
  for (j in 1: 1:k){
    for (xi in 0:n) {
      piece3 <- piece3 + dbinom(x = xi,size = n,prob = p1[j])*lfactorial(xi)
    }
  }
  H <- piece1+piece2+piece3
  out <- list(H=H)
  return(out)  
} 


g2mods.lines <- function(allcoords=ThreeDcoords2[-c(11,13),],refpoint.name="g"){
		
	g.row <- which(rownames(allcoords)==refpoint.name, arr.ind=TRUE)

	g.pt <- allcoords[g.row,]		
	all.otherpts <- allcoords[-g.row,]
		
	nmods <- nrow(all.otherpts)
	xo.vec <- rep(g.pt[1],nmods)
	yo.vec <- rep(g.pt[2],nmods)	
	zo.vec <- rep(g.pt[3],nmods)	
		
	x1.vec <- all.otherpts[,1]
	y1.vec <- all.otherpts[,2]
	z1.vec <- all.otherpts[,3]
		
	return(list(xo.vec=xo.vec,yo.vec=yo.vec,zo.vec=zo.vec,x1.vec=x1.vec,y1.vec=y1.vec,z1.vec=z1.vec))
		
}


MP.coords <- function(allKLs=FALSE,Entropies, N=1,nmds.dim=3){
	
	# If allKLs=FALSE, input 'Entropies' is a list with a matrix of model cross-entropies, an estimate of Sgg
	# (scaled by sample size), a vector of estimated cross-entropies between g and all models, Sgfis, 
	# and its counterpart, the vector of cross-entropies from the models to g, Sfisg
	# Also, the function assumes the estimate of Sgg is already scaled by N.
	# Alternatively,if allKLs=TRUE,input 'Entropies' is a matrix of KL divergences 
	# between all models and also between all models and g.  
	# (or the true Sgg if testing a simulated data set, for example)
	# N = sample size. If allKLs==FALSE and do not want to scale KLs by sample size, use N=1

	if(allKLs==FALSE){
	
		Sfifjs.hat <- Entropies[[1]];
		Sgg.hat <- Entropies[[2]];
		Sgfis.hat <- Entropies[[3]]; 
		Sfisg.hat <- Entropies[[4]];
		
		nmods <- dim(Sfifjs.hat)[1]
		nmodsp1 <- nmods +1
		
		Self.entropies <- diag(Sfifjs.hat)/N
		Self.entropies.mat0 <- matrix(rep(Self.entropies,each=nmods),nrow=nmods,ncol=nmods,byrow=TRUE)
		Self.entropies.mat <- rbind(cbind(Self.entropies.mat0,Self.entropies.mat0[,1]),rep(Sgg.hat,nmodsp1))
		Sfifjs.full <- rbind(cbind(Sfifjs.hat/N,Sfisg.hat/N),c(Sgfis.hat/N,Sgg.hat))
		row.names(Self.entropies.mat) <- c(rownames(Sfifjs.hat),"g")
		colnames(Self.entropies.mat) <- c(colnames(Sfifjs.hat),"g")
		row.names(Sfifjs.full) <- c(rownames(Sfifjs.hat),"g")
		colnames(Sfifjs.full) <- c(colnames(Sfifjs.hat),"g")
		
		All.KLs <- 2*(Self.entropies.mat - Sfifjs.full)
		Sym.AllKLs <- (All.KLs + t(All.KLs))/2
		Sym.AllKLs[Sym.AllKLs<0] <- .Machine$double.xmin

	}else if(allKLs==TRUE){
		
		All.KLs <- Entropies
		nmodsp1  <- dim(All.KLs)[1]
		nmods <- nmodsp1-1
		Sym.AllKLs <- (All.KLs + t(All.KLs))/2
		Sym.AllKLs[Sym.AllKLs<0] <- .Machine$double.xmin
		
	}

	# KL divergences between g and all models
	Estim.KLs.gfjs <- Sym.AllKLs[nmodsp1,1:nmods]

	# Estimated KL divergences between all models
	KLfifjs.4mnds <- Sym.AllKLs[1:nmods,1:nmods]
	diag(KLfifjs.4mnds) <- rep(0,nmods)	
	Data.KL.dists <- as.dist(KLfifjs.4mnds)

	# Data-based NMDS coordinates of all models
	Data.nmds <- smacofSym(delta=Data.KL.dists, ndim=nmds.dim, type="ratio") # type="ordinal" was before
	Data.fis.coords <- Data.nmds$conf
	
	# Euclidean distance between points in NMDS space
	Data.nmds.meandists <- mean(dist(Data.fis.coords,Data.fis.coords,method="euclidean"))
	Data.resc.nmdscoords <- Data.fis.coords*mean(sqrt(Data.KL.dists))/Data.nmds.meandists
	
	# Solving for h.hat and m.hat 
	Estim.m.loc <- find_mg_xyz(mgKLvec=Estim.KLs.gfjs,points=Data.resc.nmdscoords)
	Estim.m.hat.coords <- Estim.m.loc$xyz[1:nmds.dim]

	# Putting the 2-dimnl' coordinates of all the models, m.hat, and g.hat in a matrix

	TwoDcoords <- rbind(Data.resc.nmdscoords[,1:2],Estim.m.hat.coords[1:2])
	ThreeDcoords1 <- cbind(TwoDcoords,rep(0,nmodsp1)) 
	ThreeDcoords2 <- rbind(ThreeDcoords1, c(Estim.m.hat.coords[1:2],abs(Estim.m.loc$xyz[(nmds.dim+1)])))
	row.names(ThreeDcoords2)[nmodsp1:(nmods+2)] <- c( "M","g")
	colnames(ThreeDcoords2) <- c("MPxs", "MPys","MPzs")
	
	out <- list(Data.KL.dists = Data.KL.dists,Data.nmds=Data.nmds,Estim.m.loc=Estim.m.loc, XYs.mat=ThreeDcoords2)
	
}


plot.MP <- function(XYs.mat, phi=35,theta=-295,r=60,axis.extent=0.10,my.main="my.main", true.compare=FALSE,XYs.true=NULL){

	hyps <- g2mods.lines(allcoords=XYs.mat,refpoint.name="g")
	mods2m <- g2mods.lines(allcoords=XYs.mat,refpoint.name="M")
	
	xyz.range <- apply(XYs.mat,2,FUN=function(x){range(x)})
	xyz.plot.range <- xyz.range
	xs.abs.range <- abs(xyz.plot.range[2,1]-xyz.plot.range[1,1])
	ys.abs.range <- abs(xyz.plot.range[2,2]-xyz.plot.range[1,2])	

	xyz.plot.range[,1] <- xyz.plot.range[,1] +axis.extent*c(-xs.abs.range,xs.abs.range)
	xyz.plot.range[,2] <- xyz.plot.range[,2] +axis.extent*c(-ys.abs.range,ys.abs.range)
	xyz.plot.range[,3] <- xyz.plot.range[,3] +axis.extent*c(0,xyz.plot.range[,3][2])
		
	plot.cols <- mycols.ftn(n=nrow(XYs.mat))
	plot.mat <- cbind(matrix(rep(1,9), nrow=3,ncol=3),rep(2,3))
	if(true.compare==TRUE){plot.mat <- cbind(plot.mat, plot.mat+2)}
	layout(mat=plot.mat, widths=c(1,1,0.75,1,1,1,0.75,1))
	scatter3D(x=XYs.mat[,1], y=XYs.mat[,2], z=XYs.mat[,3], pch=19,col= 	
		plot.cols, xlab="",ylab="",zlab="",phi=phi,theta=theta,r=r,colkey=F,colvar=c(1:nrow(XYs.mat)), 
		cex=2,ticktype="detailed",cex.axis=1.3,bty="u", col.axis="gray", col.panel="white", col.grid="gray",
		lwd.axis=0.5,lwd.grid=0.5,adj=0, xlim=xyz.plot.range[,1],ylim=xyz.plot.range[,2],zlim=xyz.plot.range[,3])
	segments3D(x0=hyps$xo.vec, y0=hyps$yo.vec,z0=hyps$zo.vec, x1=hyps$x1.vec,
				 y1=hyps$y1.vec,z1=hyps$z1.vec,col="gray",lwd=3,lty=1,add=T,plot=T)
	segments3D(x0=mods2m$xo.vec, y0=mods2m$yo.vec,z0=mods2m$zo.vec, x1=mods2m$x1.vec,
				 y1=mods2m$y1.vec,z1=mods2m$z1.vec,col="black",lwd=1,lty=2,add=T,plot=T)
	mtext(text=my.main[1], side=3,cex=2, adj=0.85)
	mtext(text="A    ", side=3,cex=1.75, adj=0)
	plot(0,0, axes=FALSE, type="n", xlab="", ylab="")
	legend(-1.25,0.5, col=plot.cols, legend=row.names(XYs.mat), bty="n",cex=2,pch=19)	
	
	if(true.compare==TRUE){
		XYs.mat <- XYs.true
		hyps <- g2mods.lines(allcoords=XYs.mat,refpoint.name="g")
		mods2m <- g2mods.lines(allcoords=XYs.mat,refpoint.name="M")
	
		xyz.range <- apply(XYs.mat,2,FUN=function(x){range(x)})
		xyz.plot.range <- xyz.range
		xs.abs.range <- abs(xyz.plot.range[2,1]-xyz.plot.range[1,1])
		ys.abs.range <- abs(xyz.plot.range[2,2]-xyz.plot.range[1,2])	

		xyz.plot.range[,1] <- xyz.plot.range[,1] +axis.extent*c(-xs.abs.range,xs.abs.range)
		xyz.plot.range[,2] <- xyz.plot.range[,2] +axis.extent*c(-ys.abs.range,ys.abs.range)
		xyz.plot.range[,3] <- xyz.plot.range[,3] +axis.extent*c(0,xyz.plot.range[,3][2])

		scatter3D(x=XYs.mat[,1], y=XYs.mat[,2], z=XYs.mat[,3], pch=19,col= 	
		plot.cols, xlab="",ylab="",zlab="",phi=phi,theta=theta,r=r,colkey=F,colvar=c(1:nrow(XYs.mat)), 
		cex=2,ticktype="detailed",cex.axis=1.3,bty="u", col.axis="gray", col.panel="white", col.grid="gray",
		lwd.axis=0.5,lwd.grid=0.5,adj=0, xlim=xyz.plot.range[,1],ylim=xyz.plot.range[,2],zlim=xyz.plot.range[,3])
		segments3D(x0=hyps$xo.vec, y0=hyps$yo.vec,z0=hyps$zo.vec, x1=hyps$x1.vec,
				 y1=hyps$y1.vec,z1=hyps$z1.vec,col="gray",lwd=3,lty=1,add=T,plot=T)
		segments3D(x0=mods2m$xo.vec, y0=mods2m$yo.vec,z0=mods2m$zo.vec, x1=mods2m$x1.vec,
				 y1=mods2m$y1.vec,z1=mods2m$z1.vec,col="black",lwd=1,lty=2,add=T,plot=T)
		mtext(text=my.main[2], side=3,cex=2, adj=0.85)
		mtext(text="B  ", side=3,cex=1.75, adj=0)
		plot(0,0, axes=FALSE, type="n", xlab="", ylab="")
		legend(-1.25,0.5, col=plot.cols, legend=row.names(XYs.mat), bty="n",cex=2,pch=19)	
	}
	
}


plot.MP4eps <- function(XYs.mat, phi=35,theta=-295,r=60,axis.extent=0.10,my.main="my.main", true.compare=FALSE,XYs.true=NULL){

	hyps <- g2mods.lines(allcoords=XYs.mat,refpoint.name="g")
	mods2m <- g2mods.lines(allcoords=XYs.mat,refpoint.name="M")
	
	xyz.range <- apply(XYs.mat,2,FUN=function(x){range(x)})
	xyz.plot.range <- xyz.range
	xs.abs.range <- abs(xyz.plot.range[2,1]-xyz.plot.range[1,1])
	ys.abs.range <- abs(xyz.plot.range[2,2]-xyz.plot.range[1,2])	

	xyz.plot.range[,1] <- xyz.plot.range[,1] +axis.extent*c(-xs.abs.range,xs.abs.range)
	xyz.plot.range[,2] <- xyz.plot.range[,2] +axis.extent*c(-ys.abs.range,ys.abs.range)
	xyz.plot.range[,3] <- xyz.plot.range[,3] +axis.extent*c(0,xyz.plot.range[,3][2])
		
	plot.cols <- mycols.ftn(n=nrow(XYs.mat))
	plot.mat <- cbind(matrix(rep(1,9), nrow=3,ncol=3),rep(2,3))
	if(true.compare==TRUE){plot.mat <- cbind(plot.mat, plot.mat+2)}
	layout(mat=plot.mat, widths=c(1,1,0.75,1.5,1,1,0.75,1.5))
	scatter3D(x=XYs.mat[,1], y=XYs.mat[,2], z=XYs.mat[,3], pch=19,col= 	
		plot.cols, xlab="",ylab="",zlab="",phi=phi,theta=theta,r=r,colkey=F,colvar=c(1:nrow(XYs.mat)), 
		cex=1.5,ticktype="detailed",cex.axis=1.3,bty="u", col.axis="gray", col.panel="white", col.grid="gray",
		lwd.axis=0.5,lwd.grid=0.5,adj=0, xlim=xyz.plot.range[,1],ylim=xyz.plot.range[,2],zlim=xyz.plot.range[,3])
	segments3D(x0=hyps$xo.vec, y0=hyps$yo.vec,z0=hyps$zo.vec, x1=hyps$x1.vec,
				 y1=hyps$y1.vec,z1=hyps$z1.vec,col="gray",lwd=3,lty=1,add=T,plot=T)
	segments3D(x0=mods2m$xo.vec, y0=mods2m$yo.vec,z0=mods2m$zo.vec, x1=mods2m$x1.vec,
				 y1=mods2m$y1.vec,z1=mods2m$z1.vec,col="black",lwd=1,lty=2,add=T,plot=T)
	mtext(text=my.main[1], side=3,cex=1.75, adj=0.85)
	mtext(text="A    ", side=3,cex=1.75, adj=-0.20)
	plot(0,0, axes=FALSE, type="n", xlab="", ylab="")
	legend(-1.25,0.5, col=plot.cols, legend=row.names(XYs.mat), bty="n",cex=1.5,pch=19)	
	
	if(true.compare==TRUE){
		XYs.mat <- XYs.true
		hyps <- g2mods.lines(allcoords=XYs.mat,refpoint.name="g")
		mods2m <- g2mods.lines(allcoords=XYs.mat,refpoint.name="M")
	
		xyz.range <- apply(XYs.mat,2,FUN=function(x){range(x)})
		xyz.plot.range <- xyz.range
		xs.abs.range <- abs(xyz.plot.range[2,1]-xyz.plot.range[1,1])
		ys.abs.range <- abs(xyz.plot.range[2,2]-xyz.plot.range[1,2])	

		xyz.plot.range[,1] <- xyz.plot.range[,1] +axis.extent*c(-xs.abs.range,xs.abs.range)
		xyz.plot.range[,2] <- xyz.plot.range[,2] +axis.extent*c(-ys.abs.range,ys.abs.range)
		xyz.plot.range[,3] <- xyz.plot.range[,3] +axis.extent*c(0,xyz.plot.range[,3][2])

		scatter3D(x=XYs.mat[,1], y=XYs.mat[,2], z=XYs.mat[,3], pch=19,col= 	
		plot.cols, xlab="",ylab="",zlab="",phi=phi,theta=theta,r=r,colkey=F,colvar=c(1:nrow(XYs.mat)), 
		cex=1.5,ticktype="detailed",cex.axis=1.3,bty="u", col.axis="gray", col.panel="white", col.grid="gray",
		lwd.axis=0.5,lwd.grid=0.5,adj=0, xlim=xyz.plot.range[,1],ylim=xyz.plot.range[,2],zlim=xyz.plot.range[,3])
		segments3D(x0=hyps$xo.vec, y0=hyps$yo.vec,z0=hyps$zo.vec, x1=hyps$x1.vec,
				 y1=hyps$y1.vec,z1=hyps$z1.vec,col="gray",lwd=3,lty=1,add=T,plot=T)
		segments3D(x0=mods2m$xo.vec, y0=mods2m$yo.vec,z0=mods2m$zo.vec, x1=mods2m$x1.vec,
				 y1=mods2m$y1.vec,z1=mods2m$z1.vec,col="black",lwd=1,lty=2,add=T,plot=T)
		mtext(text=my.main[2], side=3,cex=1.75, adj=0.85)
		mtext(text="B  ", side=3,cex=1.75, adj=0)
		plot(0,0, axes=FALSE, type="n", xlab="", ylab="")
		legend(-1.25,0.5, col=plot.cols, legend=row.names(XYs.mat), bty="n",cex=1.5,pch=19)	
	}
	
}



#trying the functions MP.coords and plot.MP
# all.crossentropies <- as.matrix(read.table(file="SimAllSfifjs.txt", header=TRUE))
# row.names(all.crossentropies) <- colnames(all.crossentropies)
# mod.crossentropies <- all.crossentropies[1:9,1:9]
# true.selfentropy <- all.crossentropies[10,10]/1000
# truth2mods <- all.crossentropies[10,1:9]
# mods2truth <- all.crossentropies[1:9,10]
# Entropy1 <- list(Sfifjs.hat=mod.crossentropies,Sgg.hat=true.selfentropy,Sgfis.hat=truth2mods,Sfisg.hat = mods2truth)
# FullKLs.examp <- as.matrix(read.table(file="SimSym.AllKLs.txt", header=TRUE))
# row.names(FullKLs.examp) <- colnames(FullKLs.examp)
# Entropy2 <- FullKLs.examp

# # Two versions: using all the cross and self entropies to come up internally with the KLs or
# # inputing directly the KLs matrix that includes all the models and the truth nonparameric
# # estimates of KLgfis (symmetrized)
# trial.XYs1 <- MP.coords(allKLs=FALSE,Entropies=Entropy1, N=1000,nmds.dim=3)
# trial.XYs2 <- MP.coords(allKLs=TRUE,Entropies=Entropy2, N=1000,nmds.dim=3)
# trial.XYs1$XYs.mat
# trial.XYs2$XYs.mat

# plot.MP(XYs.mat=trial.XYs1$XYs.mat)



entropies.matcalc <- function(X,models.mat,modspis.persamp, pis.true=NULL){
	
	# X= the matrix of all multinomial data samples
	#    with one sample per row. nrow= N,ncol=ncateg, the number of MNOM categories
	#
	# models.mat= a matrix with as many columns as models I want to test
	#             and as many rows as multinomial samples I have in my entire data set. 
	#             So the number of rows in models.mat should match the number of rows
	#             in X.  For sample 'i' and model 'j', models.mat[i,j] contains the
	#             probability model name assigned to that sample by model 'j' 
	#
	# modspis.persamp = a list with the predicted pis per sample, per probability model.  
	#                   Each element of the list is a matrix of predicted pis 
	#                   per sample according to the probability model fit to each sample.
	#                   The dimension of this list is equal to the number of probability models, 
	#                   which is not necessarily equal to the number of columns of 
	#                   'models.mat' above.  Each 'model' (column) in 'models.mat' can be composed of a single
	#                   or multiple probability models
	#             
	#                   
	# pis.true = if available from simulation, the vector of true pis from the generating model
	#            This is used to test the method via simulations
	#
	# Note that the dimension of the list modspis.persamp HAS to match the maximum entry number in
	# the matrix 'models.mat'
	
	k <- ncol(X)
	N <- nrow(X)
	tots.persamp <- apply(X,1,sum)
	nmods <- ncol(models.mat)
	nmodsp1 <- nmods+1
	
	
	if(is.null(pis.true)==TRUE){
		
		truepis.mat <- t(apply(X,1,FUN=function(x){x/sum(x)}))
		
		}else if(is.null(pis.true)==FALSE){
			
			truepis.mat <- matrix(rep(pis.true,N), nrow=N,ncol=k,byrow=TRUE)
	}
	
	pis.list <- list()# a list of dimenson nmodsp1 with the pis per sample, per probability
				      # model combination
	
	for(j in 1:nmods){
		
		pismod.mat <- matrix(0,nrow=N,ncol=k)
		for(i in 1:N){
			
			jth.list2pick <- models.mat[i,j]
			pismod.mat[i,] <- modspis.persamp[[jth.list2pick]][i,]
			
		}
		pis.list[[j]] <- pismod.mat
		
	}
	names(pis.list) <- colnames(models.mat)
	pis.list[[nmodsp1]] <- truepis.mat
	names(pis.list) <- c(colnames(models.mat),"g")	
	
	
	Sgfs.mat <- matrix(0,nmodsp1,nmodsp1)
	for(i in 1:nmodsp1){
		for(j in 1:nmodsp1){
			
			if(i==j){
				Sggvec <- rep(0,N)
				for(n in 1:N){
					mnom.size <- tots.persamp[n]
					pis.g <- pis.list[[i]][n,]
					parslist <- list(size=mnom.size, pis.g=pis.g, type="Sgg")
					Sggvec[n] <- -H.multinom.loop(par.lst=parslist)$H	
				}
				
				Sgfs.mat[i,j] <- sum(Sggvec)
			
			}else{
				
				Sgfsvec <- rep(0,N)
				for(n in 1:N){
					
					mnom.size <- tots.persamp[n]
					pis.g <- pis.list[[i]][n,]
					pis.f <- pis.list[[j]][n,]
					parslist <- list(size=mnom.size,pis.g=pis.g, type="Sgf", pis.f=pis.f)
					Sgfsvec[n] <- -H.multinom.loop(par.lst=parslist)$H
				}		
				Sgfs.mat[i,j] <- sum(Sgfsvec)
				
			}#End if-else
			
		}# End of internal models' loop
	}# End of external models' loop
	
	colnames(Sgfs.mat) <- names(pis.list)
	row.names(Sgfs.mat) <- names(pis.list)

	entropies.list <- list(Sfifjs.hat=Sgfs.mat[1:nmods,1:nmods], Sgg.hat = Sgfs.mat[nmodsp1,nmodsp1]/N, 
	Sgfis.hat=Sgfs.mat[nmodsp1,1:nmods], Sfisg.hat = Sgfs.mat[1:nmods,nmodsp1])
	return(entropies.list)
}



