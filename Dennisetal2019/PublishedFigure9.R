# Figure 9


figpath1 <- getwd()
figname  <- "Figure9.eps" #OLD NAME: "Figure6.tiff"
figpath <- paste0(figpath1,"/",figname)

#tiff(figpath, width=6,height=7, units="in", res=600, compression="lzw", type="cairo")

postscript(figpath, horizontal=TRUE, paper="letter", width=11,height=11)


xb=seq(.01,10,.01)   # axis
adj1=.58                    # adjustment for -2 label
adj2=.693                    # adjustment for 0 label
adj3=.80                   # adjustment for +2 label

par(lwd=1.5, oma=c(2,2,1,1))

n=6
lambda=n*.25
df=3
y1=dchisq(xb,df)
y2=dchisq(xb,df,lambda)
ht=max(y1)
y1=y1+ht*1.1
xb1=xb-df*log(n)
xb2=xb-df*log(n)

y1=c(ht*1.1,y1)     # Tie down the pdfs to the axes.
y2=c(0,y2)
xb1=c(-df*log(n),xb1)
xb2=c(-df*log(n),xb2)

max(y2)

y=cbind(y1,y2)
x=cbind(xb1,xb2)
matplot(x,y,type="l",xlim=c(-12,5),ylab="Density",xlab="SIC difference",
        lty="solid", col="black", bty="n",xaxt="n",yaxt="n", cex.lab=1.5, cex.axis=1.25)
abline(v=c(-2,0,2))
abline(h=ht*1.1)
abline(h=0)

mtext(text="-2",side=1,adj=adj1,cex=1.15)
mtext(text="0",side=1,adj=adj2,cex=1.15)
mtext(text="+2",side=1,adj=adj3,cex=1.15)

n=12
xb=seq(.01,10+3*(log(n)-log(6)),.01)
lambda=n*.25
df=3
y1=dchisq(xb,df)
y2=dchisq(xb,df,lambda)
y1=y1+ht*1.1

y1=c(ht*1.1,y1)
y2=c(0,y2)

xb1=xb-df*log(n)
xb2=xb-df*log(n)

xb1=c(-df*log(n),xb1)
xb2=c(-df*log(n),xb2)

points(xb1,y1,type="l",lty=2)
points(xb2,y2,type="l",lty=2)

n=24
xb=seq(.01,10+3*(log(n)-log(6)),.01)
lambda=n*.25
df=3
y1=dchisq(xb,df)
y2=dchisq(xb,df,lambda)
y1=y1+ht*1.1

y1=c(ht*1.1,y1)
y2=c(0,y2)

xb1=xb-df*log(n)
xb2=xb-df*log(n)

xb1=c(-df*log(n),xb1)
xb2=c(-df*log(n),xb2)

points(xb1,y1,type="l",lty=3)
points(xb2,y2,type="l",lty=3)

n=48
xb=seq(.01,10+3*(log(n)-log(6)),.01)
lambda=n*.25
df=3
y1=dchisq(xb,df)
y2=dchisq(xb,df,lambda)
y1=y1+ht*1.1

y1=c(ht*1.1,y1)
y2=c(0,y2)

xb1=xb-df*log(n)
xb2=xb-df*log(n)

xb1=c(-df*log(n),xb1)
xb2=c(-df*log(n),xb2)

points(xb1,y1,type="l",lty=4)
points(xb2,y2,type="l",lty=4)


legend(2.1,.21,c("n=6","n=12","n=24","n=48"),lty=c(1,2,3,4), bty="n", cex=1.25)
legend(2.1,.21+ht*1.1,c("n=6","n=12","n=24","n=48"),lty=c(1,2,3,4), bty="n", cex=1.25)
text(c(-12,-12),c(.23+ht*1.1,.23),c("A","B"),cex=2)

dev.off()

