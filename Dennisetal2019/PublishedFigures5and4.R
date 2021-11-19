# Figures 5 and 4 (in that order) of Dennis et al 2019
# 
source("LaplaceSimsTools2.0.R")



mult.factor <- 10 #how many times bigger is the variance of g (both dimensions)
musg  <- c(0,0) 
mus1.f1 <- c(1.5,1.5)
mus.f2 <- c(3.5,3.5)


vars <- 20 # The size of the variances of f1 and f2.  Assume same variance in both dimensions
vars.g <- vars*mult.factor 

# covariance matrix of g (Sigma0) and of f1 and f2 (f1 and f2 have the same cov. matrix)
Sigma0 <- matrix(c(vars.g,0,0,vars.g), nrow=2,ncol=2,byrow=TRUE)
Sigma <- matrix(c(vars,0,0,vars), nrow=2,ncol=2,byrow=TRUE) 

# Compute four different sets of means and four different sets of sigma1.sq, sigmag.sq, K12s and DeltaK's
# each one of these four sets corresponds to a different geometric configuration (f1 closer to g):
# from all models aligned and mis-aligned 
# To compute these settings I simply started with the models aligned g, f1, f2 and rotated f1 or f2 clockwise


twoangles.params <- NSigmaSets(musg=musg,mus1.0=mus1.f1,mus2=mus.f2,Sigma1=Sigma, Sigma2=Sigma, Sigmag=Sigma0,angles=c(0, (15*pi)/8))
print(twoangles.params)

settings.mat1 <- twoangles.params$settings.mat1
settings.mat2 <- twoangles.params$settings.mat2


### Now do the plot:
plot.fname <- "Figure5.eps"
samp.sizes <- 2:50
alpha <- 0.05
k <- 2

ncases <- nrow(settings.mat1)
z.alpha <- qnorm(p=1-alpha)
#tiff(plot.fname,width=8,height=8,units="in",res=600,compression="lzw",type="cairo")
postscript(plot.fname, horizontal=FALSE,width=8,height=8, paper="letter") 

one <-matrix(1,nrow=4,ncol=4)
two <-matrix(2,nrow=4,ncol=4)
three <-matrix(3,nrow=4,ncol=4)
four <-matrix(4,nrow=4,ncol=4)
mat <- rbind(cbind(one,two), cbind(three,four),rep(5,8))
par(mar=c(4,4,2,3), oma=c(3,3,2,4))
layout(mat)
cases1 <- c("A","B") 
cases2 <- c("C","D")

for(i in 1:ncases){
  
  K12 <- settings.mat1[i,1]
  DeltaK <- settings.mat1[i,2]
  Sig.1 <- sqrt(settings.mat1[i,3])
  Sig.g <- sqrt(settings.mat1[i,4])	
  
  a.prime <- alpha.p(nvec=samp.sizes, sig.1=Sig.1,sig.g=Sig.g,K12 = K12,DeltaK=DeltaK,z.alpha=z.alpha)
  #b.prime <- beta.p(nvec=samp.sizes, sig.1=Sig.1,sig.g=Sig.g,K12 = K12,DeltaK=DeltaK,z.alpha=z.alpha)	
  Mi.prime<- Mi.p(nvec=samp.sizes,k=k,sig.g=Sig.g,DeltaK=DeltaK)
  Wi.prime<- Wi.p(nvec=samp.sizes,k=k,sig.g=Sig.g,DeltaK=DeltaK)
  
  plot(samp.sizes,a.prime, type="l", col="red", lwd=2, cex.lab=1.5,ylim=c(0,1),ylab="",axes=F,xlab="n")
  axis(1, cex.axis=1.5, outer=FALSE) 
  axis(2, cex.axis=1.5, outer=FALSE) 
  box() #box(bty="l)
  
  #points(samp.sizes,b.prime, type="l", col="blue", lwd=2, cex.lab=1.5,cex.axis=1.5)
  points(samp.sizes,Mi.prime, type="l", col="lightgreen", lwd=3, cex.lab=1.5,cex.axis=1.5,lty=2)
  points(samp.sizes,Wi.prime, type="l", col="gray", lwd=3, cex.lab=1.5,cex.axis=1.5, lty=3)
  
  text(x=45,y=0.8,cases1[i],cex=2)			
  
}


for(i in 1:ncases){
  
  K12 <- settings.mat2[i,1]
  DeltaK <- settings.mat2[i,2]
  Sig.1 <- sqrt(settings.mat2[i,3])
  Sig.g <- sqrt(settings.mat2[i,4])	
  
  b.prime <- beta.p(nvec=samp.sizes, sig.1=Sig.1,sig.g=Sig.g,K12 = K12,DeltaK=DeltaK,z.alpha=z.alpha)	
  Mi.prime<- Mi.p(nvec=samp.sizes,k=k,sig.g=Sig.g,DeltaK=DeltaK)
  Wi.prime<- Wi.p(nvec=samp.sizes,k=k,sig.g=Sig.g,DeltaK=DeltaK)
  
  plot(samp.sizes,b.prime, type="l", col="blue", lwd=2, cex.lab=1.5,ylim=c(0,1),ylab="",axes=F,xlab="n")
  axis(1, cex.axis=1.5, outer=FALSE) 
  axis(2, cex.axis=1.5, outer=FALSE) 
  box() #box(bty="l)
  
  points(samp.sizes,Mi.prime, type="l", col="lightgreen", lwd=3, cex.lab=1.5,cex.axis=1.5,lty=2)
  points(samp.sizes,Wi.prime, type="l", col="gray", lwd=3, cex.lab=1.5,cex.axis=1.5,lty=3)
  
  text(x=45,y=0.8,cases2[i],cex=2)			
  
}



mtext(text=expression(paste(alpha,"', ",beta,"',",M[i],"',",W[i],"'")), side=2, cex=1.5, outer=TRUE, at=0.6)
plot(0,0,col="white",pch=".", axes=F,xlab="",ylab="",mar=c(1,0,1,0))
legend(-0.95,-0.3, legend=c(expression(paste(alpha,"'")), expression(paste(beta,"'")), expression(paste(M[i],"'")), 
                            expression(paste(W[i],"'"))), col=c("red", "blue", "lightgreen", "gray"), lwd=c(3,3,3,3), lty=c(1,1,2,3), bty="n", horiz=T, xpd=T, inset=c(0,0), cex=2.5)

dev.off()


######## FIGURE 4 ################
## Now the contours plot:
postscript("Figure4.eps", horizontal=FALSE,width=8,height=8, paper="letter") 
#tiff("PrimesPlot2.0Geometry.tiff",width=8,height=8,units="in",res=600,compression="lzw",type="cairo")
par(mfrow=c(2,2), mar=c(5,5,2,2), oma=c(2,2,2,3))

###  Upper left (case a) 
mult.factor <- 1.5 #how many times bigger is the variance of g (both dimensions)
musg  <- c(0,0) 
mus1.f1 <- twoangles.params$n.f1means[1,]
mus.f2 <- twoangles.params$n.f2means[1,]
musgp <- twoangles.params$musgp

## Use small variances just to create a schematic, for the sake of simplicity:
vars <- 2.5#0.4#20 # The size of the variances of f1 and f2.  Assume same variance in both dimensions
vars.g <- vars*mult.factor 

# covariance matrix of g (Sigma0) and of f1 and f2 (f1 and f2 have the same cov. matrix)
Sigma0 <- matrix(c(vars.g,0,0,vars.g), nrow=2,ncol=2,byrow=TRUE)
Sigma <- matrix(c(vars,0,0,vars), nrow=2,ncol=2,byrow=TRUE) 

low.x1lim <- round(musg[1] - 3*sqrt(vars.g),digits=2)
up.x1lim  <- round(musgp[1] + 3*sqrt(vars.g),digits=2)+.01

low.x2lim <- round(musg[2] - 3*sqrt(vars.g),digits=2)
up.x2lim  <- round(musgp[2] + 3*sqrt(vars.g),digits=2)+.01

x1.vals <- seq(low.x1lim,up.x1lim,by=0.1) + musg[1] # set lower imit at -4.5 for smaller variances
x2.vals <- seq(low.x2lim,up.x2lim,by=0.1) + musg[2]


mvnorm1.vals <- my.outerbivnorm(x1vals=x1.vals,x2vals=x2.vals,mu.vec=mus1.f1, cov.mat = Sigma )
mvnorm2.vals <- my.outerbivnorm(x1vals=x1.vals,x2vals=x2.vals,mu.vec=mus.f2, cov.mat = Sigma )
laplac.vals <- my.outerbivlapl(x1vals=x1.vals,x2vals=x2.vals,mu.vec=musg, cov.mat = Sigma0, islog=TRUE)

nlevs <- 8
mvnormlevs <- pretty(range(mvnorm1.vals,finite=TRUE), nlevs)

contour(x=x1.vals,y=x2.vals,z=mvnorm1.vals, nlevels=nlevs,levels=mvnormlevs, col="blue", xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=1.5,cex.axis=1.5,drawlabels=FALSE)
par(new=TRUE)
contour(x=x1.vals,y=x2.vals,z=mvnorm2.vals,add=T, nlevels=nlevs,levels=mvnormlevs, col="red",drawlabels=FALSE)
par(new=TRUE)
contour(x=x1.vals,y=x2.vals,z=laplac.vals,add=T, nlevels=10, col="black",drawlabels=FALSE)
points(mus1.f1[1],mus1.f1[2], pch=16, col="blue")
points(mus.f2[1],mus.f2[2], pch=16, col="red")
points(musg[1],musg[2], pch=16, col="black")
text(9.75,9.75, "A", cex=2)

###  Upper right (case b) 
musg  <- c(0,0) 
mus1.f1 <- twoangles.params$n.f1means[2,]
mus.f2 <- twoangles.params$n.f2means[1,]


mvnorm1.vals <- my.outerbivnorm(x1vals=x1.vals,x2vals=x2.vals,mu.vec=mus1.f1, cov.mat = Sigma )
mvnorm2.vals <- my.outerbivnorm(x1vals=x1.vals,x2vals=x2.vals,mu.vec=mus.f2, cov.mat = Sigma )
laplac.vals <- my.outerbivlapl(x1vals=x1.vals,x2vals=x2.vals,mu.vec=musg, cov.mat = Sigma0,islog=TRUE )

nlevs <- 8
mvnormlevs <- pretty(range(mvnorm1.vals,finite=TRUE), nlevs)

contour(x=x1.vals,y=x2.vals,z=mvnorm1.vals, nlevels=nlevs,levels=mvnormlevs, col="blue", xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=1.5,cex.axis=1.5,drawlabels=FALSE)
par(new=TRUE)
contour(x=x1.vals,y=x2.vals,z=mvnorm2.vals,add=T, nlevels=nlevs,levels=mvnormlevs, col="red",drawlabels=FALSE)
par(new=TRUE)
contour(x=x1.vals,y=x2.vals,z=laplac.vals,add=T, nlevels=10, col="black",drawlabels=FALSE)
points(mus1.f1[1],mus1.f1[2], pch=16, col="blue")
points(mus.f2[1],mus.f2[2], pch=16, col="red")
points(musg[1],musg[2], pch=16, col="black")
text(9.75,9.75, "B", cex=2)

###  lower left  (case c) 
musg  <- musgp 
mus1.f1 <- twoangles.params$n.f1means[1,]
mus.f2 <- twoangles.params$n.f2means[1,]

mvnorm1.vals <- my.outerbivnorm(x1vals=x1.vals,x2vals=x2.vals,mu.vec=mus1.f1, cov.mat = Sigma )
mvnorm2.vals <- my.outerbivnorm(x1vals=x1.vals,x2vals=x2.vals,mu.vec=mus.f2, cov.mat = Sigma )
laplac.vals <- my.outerbivlapl(x1vals=x1.vals,x2vals=x2.vals,mu.vec=musg, cov.mat = Sigma0,islog=TRUE)

#contour(x=x1.vals,y=x2.vals,z=laplac.vals, nlevels=20, col="red", add=T)
#lapl.sims <- rbivlaplace(n=20000,mu.vec=musg, cov.mat=Sigma0)
#dim(lapl.sims)
#plot(lapl.sims[1,],lapl.sims[2,], pch=16)
#contour(x=x1.vals,y=x2.vals,z=laplac.vals, nlevels=20, col="red", add=T)
#mean(lapl.sims[1,])
#mean(lapl.sims[2,])


nlevs <- 8
mvnormlevs <- pretty(range(mvnorm1.vals,finite=TRUE), nlevs)

contour(x=x1.vals,y=x2.vals,z=mvnorm1.vals, nlevels=nlevs,levels=mvnormlevs, col="blue", xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=1.5,cex.axis=1.5,drawlabels=FALSE)
par(new=TRUE)
contour(x=x1.vals,y=x2.vals,z=mvnorm2.vals,add=T, nlevels=nlevs,levels=mvnormlevs, col="red",drawlabels=FALSE)
par(new=TRUE)
contour(x=x1.vals,y=x2.vals,z=laplac.vals,add=T, nlevels=15, col="black",drawlabels=FALSE)
points(mus1.f1[1],mus1.f1[2], pch=16, col="blue")
points(mus.f2[1],mus.f2[2], pch=16, col="red")
points(musg[1],musg[2], pch=16, col="black")
text(9.75,9.75, "C", cex=2)



###  lower right  (case d) 
musg  <- musgp 
mus1.f1 <- twoangles.params$n.f1means[1,]
mus.f2 <- twoangles.params$n.f2means[2,]

mvnorm1.vals <- my.outerbivnorm(x1vals=x1.vals,x2vals=x2.vals,mu.vec=mus1.f1, cov.mat = Sigma )
mvnorm2.vals <- my.outerbivnorm(x1vals=x1.vals,x2vals=x2.vals,mu.vec=mus.f2, cov.mat = Sigma )
laplac.vals <- my.outerbivlapl(x1vals=x1.vals,x2vals=x2.vals,mu.vec=musg, cov.mat = Sigma0, islog=TRUE)

nlevs <- 8
mvnormlevs <- pretty(range(mvnorm1.vals,finite=TRUE), nlevs)

contour(x=x1.vals,y=x2.vals,z=mvnorm1.vals, nlevels=nlevs,levels=mvnormlevs, col="blue", xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=1.5,cex.axis=1.5,drawlabels=FALSE)
par(new=TRUE)
contour(x=x1.vals,y=x2.vals,z=mvnorm2.vals,add=T, nlevels=nlevs,levels=mvnormlevs, col="red",drawlabels=FALSE)
par(new=TRUE)
contour(x=x1.vals,y=x2.vals,z=laplac.vals,add=T, nlevels=15, col="black",drawlabels=FALSE)
points(mus1.f1[1],mus1.f1[2], pch=16, col="blue")
points(mus.f2[1],mus.f2[2], pch=16, col="red")
points(musg[1],musg[2], pch=16, col="black")
text(9.75,9.75, "D", cex=2)

dev.off()

# Now a Monte Carlo approximation of the same...
