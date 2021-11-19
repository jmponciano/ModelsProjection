
# Figure 8


figpath1 <- getwd()
figname  <- "Figure7.eps" #OLD NAME: "Figure5.tiff"
figpath <- paste0(figpath1,"/",figname)

#tiff(figpath, width=6,height=7, units="in", res=600, compression="lzw", type="cairo")

postscript(figpath, horizontal=TRUE, paper="letter", width=11,height=11)


xb=seq(.01,10,.01)   # axis
adj1=.4                    # adjustment for -2 label
adj2=.595                    # adjustment for 0 label
adj3=.78                   # adjustment for +2 label

par(lwd=1.5, oma=c(2,1,1,1))
n=6
lambda=n*.25
df=3
y1=dchisq(xb,df)
y2=dchisq(xb,df,lambda)
ht=max(y1)
y1=y1+ht*1.1
xb1=xb-2*df
xb2=xb-2*df

y1=c(ht*1.1,y1)     # Tie down the pdfs to the axes.
y2=c(0,y2)
xb1=c(-2*df,xb1)
xb2=c(-2*df,xb2)

max(y2)

y=cbind(y1,y2)
x=cbind(xb1,xb2)
matplot(x,y,type="l",ylab="Density",xlab="AIC difference",lty="solid",
        col="black",bty="n",xaxt="n",yaxt="n", cex.lab=1.5,cex.axis=1.25)
abline(v=c(-2,0,2))
abline(h=ht*1.1)
abline(h=0)

mtext(text="-2",side=1,adj=adj1,cex=1.15)
mtext(text="0",side=1,adj=adj2,cex=1.15)
mtext(text="+2",side=1,adj=adj3,cex=1.15)

n=12
lambda=n*.25
df=3
y2=dchisq(xb,df,lambda)
y2=c(0,y2)

max(y2)

points(xb2,y2,type="l",lty=2)

n=24
lambda=n*.25
df=3
y2=dchisq(xb,df,lambda)
y2=c(0,y2)

max(y2)

points(xb2,y2,type="l",lty=3)

n=48
lambda=n*.25
df=3
y2=dchisq(xb,df,lambda)
y2=c(0,y2)

max(y2)

points(xb2,y2,type="l",lty=4)

legend(x=2.3,y=.23,c("n=6","n=12","n=24","n=48"),lty=c(1,2,3,4), bty="n",cex=1.25)
text(c(-6,-6),c(.23+ht*1.16,.2),c("A","B"), cex=2)

dev.off()