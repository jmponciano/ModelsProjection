#### Example of the MPMS calculations with the horseshoe crabs simulator:
source("CrabsExampleTools.R")
source("MPcalctools.R")



#####  First,  setting the simulator parameters according to one of the most 
#####  plausible scenarios in Brockmann et al's 2018
##### Setting up beta, for a target density 'high'

Targetdens <- 1500
alpha.val <- 0.32
A <- 1000
osr <- 2.5
adv.freq <- 10
beta.val <- beta.calc(A=A,alpha.val=alpha.val,single.targetdens=Targetdens,single.adfreq=adv.freq,single.osr=osr)

### Setting up the width of male search
width.search <- 34 # per table 1 in paper

### Setting up the models to be fit:
my.models <- c("Poisson", "NegBin", "ZIPoiss", "ZINegBi", "HurdNBi", "PoisNB", "NBPois", "OIPoiss", "OINegBi")
my.modsf <- 1:9 # 1= Poi, 2=NB, 3=ZIP,  4= ZNB, 5= HNB, 6=PoisNB, 7=NBPois, 8=OIPoiss,9=OINegBi


###### Cranking up a simulation with data for 300 tides:
###### One simulation with 300 tides

ntides <- 300
ten.tides <- short.sim(adv.freq = adv.freq, osr=osr, beta=beta.val, a=width.search, nsims=ntides, alpha=0.32, area=1000, hetero.abund=FALSE, male.hetero=FALSE, decision.mode="rand.unif", tracking=FALSE,nonunif=FALSE,ncateg=8,my.models=my.models)

#### Simulation results:  for each tide the output is:
#### object$simdat" The simulated counts, one row per tide. Data under column 'j' and row 'i' is the number of
#### pairs with 'j' satellites observed during tide 'i'
#### object$bics.mat:  For tide 'i' (rows) and model j'(cols), the value of the bic score 
#### object$fits.list: the list of the fit stats for every model and every tide.
#### if data for 'n' tides is simulated, fits.list is of length 'n' and
#### every element 'i' of that list is itself a list with the fit statistics for every model, for tide 'i'


###  Now using my function to take the list of fitted stats and return the 'pis' for all probability models 
###  and all the reps of vectors of counts
###  function name:  "pred.pis"

fits.list <- ten.tides$fits.list
simdat <-  ten.tides$simdat


test.predpis <- pred.pis(fits.list= fits.list,simdat=simdat)


### Creating a matrix of models:
nmods <- length(my.models)
tested.models1 <- matrix(0, nrow = ntides,ncol=nmods)
for(i in 1:nmods){tested.models1[,i] <- rep(i,ntides)}
colnames(tested.models1) <- my.models

##### Computation of true pis for this setting, keeping 6 categories as the simulations:
Ntides <- 20000
Inf.tides <- short.sim(adv.freq = adv.freq, osr=osr, beta=beta.val, a=width.search, nsims=Ntides, alpha=0.32, area=1000, hetero.abund=FALSE, male.hetero=FALSE, decision.mode="rand.unif", tracking=FALSE,nonunif=FALSE,ncateg=8)
N.tidessims <- Inf.tides$simdat
dim(N.tidessims)
save.image("twenty.th.tides.RData")


N.tides.counts <- Inf.tides$counts
N.tides.counts[[1]]

k <- ncol(ten.tides$simdat)
km1 <- k-1
truelongdat <- matrix(0,nrow=Ntides,ncol=k)
colnames(truelongdat) <- as.character(0:km1)
	
for(i in 1:Ntides){
		
	new.dat <- Inf.tides$counts[[i]]
	len.dat <- ncol(new.dat)
	if(len.dat>k){
		tail.sum <- sum(new.dat[2,k:len.dat])  
		short.counts <- c(new.dat[2,1:km1],tail.sum)
		truelongdat[i,] <- short.counts
	}else if(len.dat ==k){
		truelongdat[i,] <- new.dat[2,]
	}else if(len.dat<k){
		nzeros <- k-len.dat
		truelongdat[i,] <- c(new.dat[2,],rep(0,nzeros))
	}	
}

dim(truelongdat)

TRUE.PIS <- apply(truelongdat,2,sum)/sum(apply(truelongdat,2,sum))
sum(TRUE.PIS)

save.image("July9th.RData")

#####  Now, finally, gathering the input for the model projection calculation 
#####  First, knowing the true model probabilities
out4MPcoords1 <- entropies.matcalc(X=simdat, models.mat=tested.models1, modspis.persamp=test.predpis, pis.true=TRUE.PIS)

#####  Second, not knowing the true model probabilities, but estimating Sgg from data
out4MPcoords2 <- entropies.matcalc(X=simdat, models.mat=tested.models1, modspis.persamp=test.predpis, pis.true=NULL)

true.MP <- MP.coords(allKLs=FALSE, Entropies=out4MPcoords1, N=ntides, nmds.dim=3)
estim.MP <- MP.coords(allKLs=FALSE, Entropies=out4MPcoords2, N=ntides, nmds.dim=3)

round(dist(estim.MP$XYs.mat),digits=5)
round(dist(true.MP$XYs.mat),digits=5)

estim.true.Mg <- rbind(estim.MP$XYs.mat[10:11,],true.MP$XYs.mat[10:11,])
row.names(estim.true.Mg) <- c("hat.m","hat.g","true.m", "true.g")
dist(estim.true.Mg)

#> dist(estim.true.Mg)
#              hat.m        hat.g       true.m
#hat.g  0.0003226238                          
#true.m 0.3830739458 0.3830740817             
#true.g 0.3830741268 0.3830739490 0.0003723332


plot.fname <- "Figure5-v3-5.tiff"
tiff(plot.fname,width=16,height=8,units="in",compression="lzw",type="cairo",res=600,family="times")
par(mar=c(3,3,3,1), oma=c(3,3,3,1))
plot.MP(XYs.mat=estim.MP$XYs.mat, my.main=c("Estimated Projection","True Projection"), true.compare=TRUE, XYs.true=true.MP$XYs.mat);
dev.off()

#####  This is the figure 4 for the final manuscript
postscript("Figure5-v3-5.eps",width=16,height=8, paper="a4")
par(mar=c(3,3,3,1), oma=c(3,3,3,1))
plot.MP4eps(XYs.mat=estim.MP$XYs.mat, my.main=c("Estimated Projection","True Projection"), true.compare=TRUE, XYs.true=true.MP$XYs.mat);
dev.off()

save.image("July10th.RData")




########  Doing Figure 5: including weighted AICs in here


code.string <- paste0(rep("ith.fit$",nmods),my.models,rep(".Stats$logLhat",nmods))

logLhats.mat <- matrix(0,nrow=ntides,ncol=nmods)
for(i in 1:ntides){

	ith.fit <- fits.list[[i]]
	logLhats.ithrow <- rep(0,nmods)
	for(j in 1:nmods){logLhats.ithrow[j] <- eval(parse(text=code.string[j]))} 	
	logLhats.mat[i,] <- logLhats.ithrow
	
}
nmodpars <- c(1,2,2,3,3,3,3,2,3)
logLtots <- apply(logLhats.mat,2,sum)
AICs <- -2*logLtots + 2*nmodpars
samp.sizes <- apply(simdat,2,sum) #don't really need this

delta.is <- AICs - min(AICs)
exp.deltais <- exp(-delta.is/2)
w.is <- exp.deltais/sum(exp.deltais) 

mods.coords <- estim.MP$XYs.mat[1:nmods,1:2]
wAIC.coords <- c(sum(mods.coords[,1]*w.is),sum(mods.coords[,2]*w.is))
little.m <- estim.MP$XYs.mat[(nmods+1),1:2]
big.m <- true.MP$XYs.mat[(nmods+1),1:2]

dist(rbind(wAIC.coords,little.m,big.m))

#dist(rbind(wAIC.coords,little.m,big.m))
#         wAIC.coords  little.m
#little.m   1.3783360          
#big.m      1.4760623 0.3739722

#> dist(rbind(wAIC.coords,little.m,big.m))
#         wAIC.coords  little.m
#little.m   1.7614023          
#big.m      1.7845545 0.3830739


fig5cols <- mycols.ftn(n=nmods+3) # nmods + little m + big m + wAIC

postscript("Figure6-v3-5.eps",width=10,height=8, paper="letter")
plot.mat <- cbind(matrix(rep(1,9), nrow=3,ncol=3),rep(2,3))
par(mar=c(5,5,3,1), oma=c(2,2,3,1))
layout(mat=plot.mat, widths=c(1,1,0.75,1.5))
plot(c(estim.MP$XYs.mat[1:nmods,1], little.m[1], big.m[1], wAIC.coords[1]), 
c(estim.MP$XYs.mat[1:nmods,2], little.m[2], big.m[2], wAIC.coords[2]), pch=19, col=fig5cols, cex=2,
xlab="axis 1", ylab="axis 2", main="Models Hyperplane", cex.main=2,cex.lab=2,cex.axis=1.5,xlim=c(-2.5,3.5))
text(little.m[1],little.m[2],"m", cex=1.5,pos=3)
text(big.m[1],big.m[2],"M", cex=1.5,pos=3)
text(wAIC.coords[1],wAIC.coords[2],"wAIC", cex=1.5,pos=3)

plot(0,0, axes=FALSE, type="n", xlab="", ylab="")
legend(-1.,0.5, col=fig5cols, legend=c(my.models,"m","M","wAIC"), bty="n",cex=2,pch=19)	
dev.off()

tiff("Figure6-v3-5.tiff",width=10,height=8, units="in",compression="lzw",type="cairo",res=600,family="times")
plot.mat <- cbind(matrix(rep(1,9), nrow=3,ncol=3),rep(2,3))
par(mar=c(5,5,3,1), oma=c(2,2,3,1))
layout(mat=plot.mat, widths=c(1,1,0.75,1.5))
plot(c(estim.MP$XYs.mat[1:nmods,1], little.m[1], big.m[1], wAIC.coords[1]), 
c(estim.MP$XYs.mat[1:nmods,2], little.m[2], big.m[2], wAIC.coords[2]), pch=19, col=fig5cols, cex=2,
xlab="axis 1", ylab="axis 2", main="Models Hyperplane", cex.main=2,cex.lab=2,cex.axis=1.5,xlim=c(-2.5,3.5))
text(little.m[1],little.m[2],"m", cex=1.5,pos=3)
text(big.m[1],big.m[2],"M", cex=1.5,pos=3)
text(wAIC.coords[1],wAIC.coords[2],"wAIC", cex=1.5,pos=3)

plot(0,0, axes=FALSE, type="n", xlab="", ylab="")
legend(-1.,0.5, col=fig5cols, legend=c(my.models,"m","M","wAIC"), bty="n",cex=2,pch=19)	
dev.off()

save.image("Sept6th.RData")


## SMACOF object: trying to plot the models with a confidence ellipse

## Re-scaling of axis to match previous figure:
xs.all <- estim.MP$XYs.mat[1:nmods,1]
ys.all <- estim.MP$XYs.mat[1:nmods,2]
smacof.estim <- estim.MP$Data.nmds
conf.estim <- confEllipse(smacof.estim)

xs.smacof <- conf.estim$X[,1]
ys.smacof <- conf.estim$X[,2]

lm.xs <- lm(xs.all~xs.smacof)
xs.4axis <- c(-0.1,0,0.1,0.2)
labs.4axis <- round(lm.xs$coef[1] + lm.xs$coef[2]*xs.4axis,digits=2)

lm.ys <- lm(ys.all~ys.smacof)
ys.4axis <- c(-.15,-.10,-.05,0,0.05,0.10,0.15)
labs.4axis2 <- round(lm.ys$coef[1] + lm.ys$coef[2]*ys.4axis,digits=2)

postscript("Figure4p5-v3p2.eps",width=8,height=8, paper="letter")
par(mar=c(5,5,3,1), oma=c(2,2,3,1))
plot(conf.estim, plot.dim=c(1,2), label.conf=list(label=TRUE, pos=5,cex=0.95),ell=list(lty=1:9, lwd=2, col="grey",adj=c(0.8,1)),cex.lab=1.5,cex.axis=1.25, col="black", cex.main=1.5,pch=19,axes=FALSE)
axis(side=1, at=xs.4axis,labels=labs.4axis)
axis(side=2,at=ys.4axis,labs.4axis2)
dev.off()


tiff("Figure4p5-v3p2.tiff",width=8,height=8, units="in",compression="lzw",type="cairo",res=600,family="times")
par(mar=c(5,5,3,1), oma=c(2,2,3,1))
plot(conf.estim, plot.dim=c(1,2), label.conf=list(label=TRUE, pos=5,cex=0.95),ell=list(lty=1:9, lwd=2, col="grey",adj=c(0.8,1)),cex.lab=1.5,cex.axis=1.25, col="black", cex.main=1.5,pch=19,axes=FALSE)
axis(side=1, at=xs.4axis,labels=labs.4axis)
axis(side=2,at=ys.4axis,labs.4axis2)
dev.off()




save.image("Sept19th.RData")

###### OUTPUT TO PUT ON PAPER:
> dist(estim.MP$XYs.mat)
             Poisson       NegBin      ZIPoiss      ZINegBi      HurdNBi       PoisNB       NBPois      OIPoiss      OINegBi            M
NegBin  4.4671083507                                                                                                                     
ZIPoiss 5.0935841333 1.5161583938                                                                                                        
ZINegBi 5.2142191416 1.2606895949 0.4278362764                                                                                           
HurdNBi 5.2244412255 1.1969917328 0.5234683040 0.0977307769                                                                              
PoisNB  5.1085779359 0.8956082347 0.8130911946 0.4254229023 0.3369049563                                                                 
NBPois  5.1988332063 1.2434783761 0.4320785679 0.0179829239 0.0913943141 0.4122826867                                                    
OIPoiss 0.0018115021 4.4685867505 5.0953003820 5.2158838201 5.2260923142 5.1101803220 5.2004966052                                       
OINegBi 1.6949971467 2.8533169874 3.7387517875 3.7681794485 3.7592555517 3.5874556275 3.7512512834 1.6961761617                          
M       4.5123514073 2.3328017336 1.2656600692 1.6745139126 1.7614023164 1.9744392401 1.6738741433 4.5141554613 3.5030578142             
g       4.5123514189 2.3328017559 1.2656601103 1.6745139437 1.7614023460 1.9744392665 1.6738741744 4.5141554729 3.5030578291 0.0003226238
> dist(true.MP$XYs.mat)
             Poisson       NegBin      ZIPoiss      ZINegBi      HurdNBi       PoisNB       NBPois      OIPoiss      OINegBi            M
NegBin  4.4671083507                                                                                                                     
ZIPoiss 5.0935841333 1.5161583938                                                                                                        
ZINegBi 5.2142191416 1.2606895949 0.4278362764                                                                                           
HurdNBi 5.2244412255 1.1969917328 0.5234683040 0.0977307769                                                                              
PoisNB  5.1085779359 0.8956082347 0.8130911946 0.4254229023 0.3369049563                                                                 
NBPois  5.1988332063 1.2434783761 0.4320785679 0.0179829239 0.0913943141 0.4122826867                                                    
OIPoiss 0.0018115021 4.4685867505 5.0953003820 5.2158838201 5.2260923142 5.1101803220 5.2004966052                                       
OINegBi 1.6949971467 2.8533169874 3.7387517875 3.7681794485 3.7592555517 3.5874556275 3.7512512834 1.6961761617                          
M       4.1468776933 2.1494228201 1.3458667302 1.7095859193 1.7845545353 1.9397031373 1.7049315000 4.1486767599 3.1205851558             
g       4.1468777101 2.1494228524 1.3458667817 1.7095859598 1.7845545741 1.9397031730 1.7049315407 4.1486767766 3.1205851780 0.0003723332
> 0.0003226238/0.0003723332
[1] 0.8664922
> 1.11/0.77
[1] 1.441558
> estim.true.Mg <- rbind(estim.MP$XYs.mat[10:11,],true.MP$XYs.mat[10:11,])
> row.names(estim.true.Mg) <- c("hat.m","hat.g","true.m", "true.g")
> dist(estim.true.Mg)
              hat.m        hat.g       true.m
hat.g  0.0003226238                          
true.m 0.3830739458 0.3830740817             
true.g 0.3830741268 0.3830739490 0.0003723332


> dist(rbind(wAIC.coords,little.m,big.m))
         wAIC.coords  little.m
little.m   1.7614023          
big.m      1.7845545 0.3830739
