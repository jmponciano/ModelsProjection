library(plot3D)
library(pracma)
library(plyr)
Drange=function(dis,incld0=F,lmts=T,lngth=T){
  v=as.vector(dis)
  if (incld0) idx=1:length(dis)
  else idx=v != 0
  v=v[idx]
  limits=range(v,na.rm=T,finite=T)
  length=limits[2]-limits[1]
  out=list()
  if (lmts) out=c(out,list(limits=limits))
  if (lngth) out=c(out,list(length=length))
  return(out)
}


Range=function(mat,margins=c(2),lmts=T, lngth=T){
  out=alply(.data = mat,.margins = margins,.fun = Drange,lmts=lmts,lngth=lngth,incld0=T)
  return(out)
}


figf=matrix(c(11.5,14.6,0,11.4,10.5,0,12.1,7.5,0,17.5,10.1,0,15.8,2,0),ncol=3,byrow=T)
figf[,1:2]=scale(figf[,1:2],center = F,scale=unlist(Range(figf[,1:2],lmts=F)))
figf[,1]=figf[,1]-min(figf[,1])
figf[,2]=figf[,2]-min(figf[,2])

h=.4
p=c(.3,.5,0)
g=c(.3,.5,h)


figfpg=rbind(figf,p,g)
x=figfpg[,1]
y=figfpg[,2]
z=figfpg[,3]
f1=figf[1,]
f2=figf[2,]
f3=figf[3,]
f4=figf[4,]
f5=figf[5,]
mdnms=c("f1","f2","f3","f4","f5","m","g")
rownames(figfpg)=mdnms

#pdf(file = "modelProjectionFig.pdf",)
#setEPS(); postscript(file = "modelProjectionFig.eps")
#tiff("Figure4.tiff",width=22,height=22,units="cm",res=600,compression="lzw",type="cairo") # w/l 22/22
postscript("Figure4.eps", width=22,height=22)


par(mfrow = c(1, 1))
fig=scatter3D(x,y,z,type="p",bty=c("g"),xlab="Hyperplane axis 1 ",ylab="Hyperplane axis 2",zlab="Z",
              xlim=c(-0.1,1),ylim=c(-0.1,1),zlim=c(-0.1,1), pch=rep(16,7),colvar=c(1,1,1,1,1,2,2), cex=4,
              col=c("black","darkgray"),colkey=F,labels=mdnms,phi=20,theta=-20,r=2,cex.axis=2.5) 

segments3D(figf[1, 1],figf[1, 2],figf[1, 3],g[1],g[2],g[3],col="black",lwd=4,add=T,plot=T)
segments3D(f2[1],f2[2],f2[3],g[1],g[2],g[3],col="black",lwd=4,add=T,plot=T)
segments3D(f3[1],f3[2],f3[3],g[1],g[2],g[3],col="black",lwd=4,add=T,plot=T)
segments3D(f4[1],f4[2],f4[3],g[1],g[2],g[3],col="black",lwd=4,add=T,plot=T)
segments3D(f5[1],f5[2],f5[3],g[1],g[2],g[3],col="black",lwd=4,add=T,plot=T)

segments3D(f1[1], f1[2], f1[3],f2[1], f2[2], f2[3],col="black",lwd=1, lty=2,add=T,plot=T)
segments3D(f1[1], f1[2], f1[3],f3[1], f3[2], f3[3],col="black",lwd=1, lty=2,add=T,plot=T)
segments3D(f1[1], f1[2], f1[3],f4[1], f4[2], f4[3],col="black",lwd=1, lty=2,add=T,plot=T)
segments3D(f1[1], f1[2], f1[3],f5[1], f5[2], f5[3],col="black",lwd=1, lty=2,add=T,plot=T)
segments3D(f2[1], f2[2], f2[3],f3[1], f3[2], f3[3],col="black",lwd=1, lty=2,add=T,plot=T)
segments3D(f2[1], f2[2], f2[3],f4[1], f4[2], f4[3],col="black",lwd=1, lty=2,add=T,plot=T)
segments3D(f2[1], f2[2], f2[3],f5[1], f5[2], f5[3],col="black",lwd=1, lty=2,add=T,plot=T)
segments3D(f3[1], f3[2], f3[3],f4[1], f4[2], f4[3],col="black",lwd=1, lty=2,add=T,plot=T)
segments3D(f3[1], f3[2], f3[3],f5[1], f5[2], f5[3],col="black",lwd=1, lty=2,add=T,plot=T)
segments3D(f4[1], f4[2], f4[3],f5[1], f5[2], f5[3],col="black",lwd=1, lty=2,add=T,plot=T)

segments3D(p[1],p[2],p[3],g[1],g[2],g[3],col="gray",lwd=4,add=T,plot=T)

d2r=function(d){d*pi/180}
r2d=function(r){r*180/pi}
norm_vec <- function(x){sqrt(crossprod(x))}
EYE=function(theta,phi,r,xlim,ylim,zlim,degree=T,bndP=c(0,0,0)){
  if (degree) {
    t=d2r(theta)
    p=d2r(phi)
  } else {
    t=theta
    p = phi
  }
  dx=(xlim[2]-xlim[1])*bndP[1]
  dy=(ylim[2]-ylim[1])*bndP[2]
  dz=(zlim[2]-zlim[1])*bndP[3]
  cntr=(c(sum(xlim),sum(ylim),sum(zlim))/2)-bndP/2
  vpxyz=sph2cart(tpr=c(t-(pi/2) ,p,r))+cntr
  out=list(vpxyz=vpxyz,cntr=cntr)
  return(out)
}
xlim=c(-0.1,1)
ylim=c(-0.1,1)
zlim=c(-0.1,1)


bp=0
ey=EYE(theta=-20,phi=20,r=2,xlim,ylim,zlim,bndP=c(bp,bp,bp))$vp
cntr=EYE(theta=-20,phi=20,r=2,xlim,ylim,zlim,bndP=c(bp,bp,bp))$cntr




scl=.4
names(x) <- NULL
names(y) <- NULL
names(z) <- NULL
sft=scl*((matrix(ey,ncol=3,nrow=length(x),byrow=T)-cbind(x,y,z))/aaply(.data = (matrix(ey,ncol=3,nrow=length(x),byrow=T)-cbind(x,y,z)),.margins = 1,.fun = norm_vec))
colnames(sft)=c("x","y","z")
text3D(x+sft[ ,1],y+sft[, 2],z+sft[ ,3],labels=mdnms,colvar=c(rep(1,5),2,2),col=c("white","white"),add=T,plot=T,colkey=F,adj=0.5,font=2)


dev.off()
