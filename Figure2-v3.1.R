## FigurE 2
library(shape)


######################## New Figure 2, August 22nd 2019 #####################
#postscript("Figure2-v3.1.eps",width=37,height=24)#30/25
tiff("Figure2.tiff",width=37,height=24,units="cm",res=600,compression="lzw",type="cairo")

g   <- c(5,9.3)
f2  <- c(2.4,9) # ZIP
f4  <- c(2,6.75) #ZINEGBIN
f5  <- c(1.9,4.5) # HURDNEGBIN 
f6  <- c(3,2.75) # POINB
f3  <- c(6,3.5) # NEGBIN
f1  <- c(8,4.5) # POISSON

#mod.labs <- c("POIS", "ZIP", "NBIN", "ZINBI", "HUNBI", "POINB") 


# Panel A
par(mar=c(0,0,0,0),mfrow=c(1,2))
plot(0,0,pch=".",col="white",axes=F,xlim=c(1,10),ylim=c(1,11))
plotcircle(mid=g,r=1,col="grey90",lwd=2)
text(g[1],g[2],"g",font=3,cex=2)


center <- f2 
plotcircle(mid=center,r=1,lwd=4,col="grey90")
text(center[1],center[2],expression("f"[2]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.98 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/2)*0.95,center[2]+(v1[2]/2)*3.3,expression("d"[2]),cex=2) # manually change model num

center <- f4 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[4]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.95 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/3)*0.8,center[2]+(v1[2]/3)*1.3,expression("d"[4]),cex=2) # manually change model num


center <- f5 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[5]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/2)*0.84,center[2]+(v1[2]/2)*1.1,expression("d"[5]),cex=2) # manually change model num


center <- f6 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[6]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.85 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/2)*0.65,center[2]+v1[2]/2,expression("d"[6]),cex=2) # manually change model num

center <- f3 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[3]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.85 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/2)*1.7,center[2]+v1[2]/2,expression("d"[3]),cex=2) # manually change model num

center <- f1 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[1]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.88 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/2)*1.2,center[2]+v1[2]/2,expression("d"[1]),cex=2) # manually change model num

text(5.5,10.5,expression("K-L Information"),font=2,cex=2.5)
text(1,10.5,"A", font=1, cex=3)


#######  PANEL B

plot(0,0,pch=".",col="white",axes=F,xlim=c(1,10),ylim=c(1,11))
plotcircle(mid=g,r=1,col="grey90",lwd=4)
text(g[1],g[2],expression("f"[2]),font=3,cex=2)


center <- f4 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[4]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.95 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/3)*1.2,center[2]+(v1[2]/3)*1.6,expression("d"[4]-"d"[2]),cex=2,srt=48) # manually change model num


center <- f5 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[5]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/2)*1.05,center[2]+(v1[2]/2)*1.25,expression("d"[5]-"d"[2]),cex=2,srt=62) # manually change model num



center <- f6 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[6]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.85 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/2)*0.65,center[2]+v1[2]/2,expression("d"[6]-"d"[2]),cex=2,srt=76) # manually change model num


center <- f3 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[3]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.85 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/2)*1.7,center[2]+v1[2]/2,expression("d"[3]-"d"[2]),cex=2,srt=97) # manually change model num


center <- f1 
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[1]),font=3,cex=2) # manually change model num
v1	<-	g-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.88 # manually change fraction multiplier
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=g[1]-u[1],y1=g[2]-u[2],lwd=2,length=0.1)
text(center[1]+(v1[1]/2)*1.2,center[2]+v1[2]/2,expression("d"[1]-"d"[2]),cex=2,srt=118) # manually change model num

text(4.8,10.5,expression(paste(Delta," AIC")),font=2,cex=2.5)
text(2,10.5,"B", font=1, cex=3)


dev.off()


########################  OLD FIGURE 2 BELOW, FROM MANUSCRIPT VERSION Mod.Proj.Frontiers.3.0...########
#tiff("Figure2.tiff",width=37,height=24,units="cm",res=600,compression="lzw",type="cairo")
postscript("Figure2.eps",width=37,height=24)#30/25
# Panel A
par(mar=c(0,0,0,0),mfrow=c(1,2))
plot(0,0,pch=".",col="white",axes=F,xlim=c(1,10),ylim=c(1,11))
plotcircle(mid=c(5,9),r=1,col="grey90",lwd=2)
text(5,9,"g",font=3,cex=2)

plotcircle(mid=c(2,7),r=1,lwd=4,col="grey90")
text(2,7,expression("f"[2]),font=3,cex=2)
v1	<-	c(5,9)-c(2,7)
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=2+u[1],y0=7+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(2+u[1],7+u[1],expression("d"[2]),adj=c(0,-0.75),cex=2)


center	<-	c(3.5,5.5)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[5]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[5]),adj=c(0.75,-3),cex=2)

center	<-	c(5,3)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[4]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[4]),adj=c(1.25,-9),cex=2)

center	<-	c(7,4.5)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[1]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[1]),adj=c(2.5,-6),cex=2)

center	<-	c(8,6.5)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[3]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[3]),adj=c(2,-3),cex=2)
text(5.5,10.5,expression("K-L Information"),font=2,cex=3)
text(1,10.5,"A", font=1, cex=3)

### Panel B

plot(0,0,pch=".",col="white",axes=F,xlim=c(1,9),ylim=c(1,11))
plotcircle(mid=c(5,9),r=1,col="grey90",lwd=4)
text(5,9,expression("f"[2]),font=3,cex=2)

plotcircle(mid=c(2,7),r=1,lwd=2,col="grey90")
text(2,7,expression("f"[5]),font=3,cex=2)
v1	<-	c(5,9)-c(2,7)
u	<-	(v1/(sqrt(sum(v1^2))))
arrows(x0=2+u[1],y0=7+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(2+u[1],7+u[1],expression("d"[5]-"d"[2]),adj=c(0.15,0.15),cex=2,srt=35)


center	<-	c(3.5,5.5)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[4]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[4]-"d"[2]),adj=c(-0.6,-0.8),cex=2,srt=67)

center	<-	c(5,3)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[1]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[1]-"d"[2]),adj=c(-2,1.1),cex=2,srt=90)

center	<-	c(7,4.5)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[3]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[3]-"d"[2]),adj=c(-1.8,2.8),cex=2,srt=114)

text(4.8,10.5,expression(paste(Delta," AIC")),font=2,cex=3)
text(2,10.5,"B", font=1, cex=3)


dev.off() 

