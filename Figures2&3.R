#### Figuras para Jose Miguel 

## Figure 1: THIS IS THE OLD VERSION OF FIGURE 1, GO DIRECTLY TO "FIGURE 2" BELOW, LINE 66
postscript("~/Desktop/figure1.eps",width=25,height=17)
#tiff("~/Desktop/figure1.tiff",width=25,height=17,units="cm",res=300,compression="lzw",type="cairo")
par(mar=c(0,0,0,0))
plot(0,0,pch=".",col="white",xlim=c(0.9,14.1),ylim=c(0.9,14.1),axes=F)

poly.mat	<-	matrix(NA,4,2)

poly.mat[1,]	<-	c(1,1)
poly.mat[2,]	<-	c((1+3),(sqrt(8^2-3^2)+1))
poly.mat[3,]	<-	c(poly.mat[2,1]+10,poly.mat[2,2])
poly.mat[4,]	<-	c(11,1)

polygon(poly.mat)

para.center	<-	apply(poly.mat,2,sum)/4

triangle1.mat	<-	matrix(NA,3,2)
triangle1.mat[1,]	<-	c(para.center[1]-1.5,para.center[2])
triangle1.mat[2,]	<-	c(triangle1.mat[1,1],triangle1.mat[1,2]+6)
triangle1.mat[3,]	<-	c(para.center[1]+1.5,para.center[2])

polygon(triangle1.mat,col="white")

triangle2.mat	<-	matrix(NA,3,2)
triangle2.mat[1,]	<-	c(triangle1.mat[1,1],triangle1.mat[1,2]+6)
triangle2.mat[2,]	<-	c(para.center[1]+1.5,triangle2.mat[1,2])
triangle2.mat[3,]	<-	c(para.center[1]+1.5,triangle1.mat[1,2])

polygon(triangle2.mat,col="white",lty=2)

lines(c(para.center[1]-1.5,para.center[1]+1.5),c(para.center[2],para.center[2]+6),lty=2)
lines(c(para.center[1]-1.5,para.center[1]+1.5),c(para.center[2]+6,para.center[2]))

square.mat	<-	matrix(NA,4,2)

square.mat[1,]	<-	c(para.center[1]-1.5,para.center[2])
square.mat[2,]	<-	c(square.mat[1,1],square.mat[1,2]+6)
square.mat[3,]	<-	c(para.center[1]+1.5,square.mat[2,2])
square.mat[4,]	<-	c(square.mat[3,1],para.center[2])

points(square.mat,pch=19)

text(para.center[1],para.center[2],"a",adj=c(0,1),cex=2.5)
text(para.center[1],para.center[2]+6,"c",adj=c(0,-0.5),cex=2.5)
text(para.center[1]-1.5,para.center[2],expression(theta["0k"]),adj=c(1,1),cex=2.5)
text(para.center[1]+1.5,para.center[2],expression(hat(theta)["k"]),adj=c(-0.25,0.75),cex=2.5)
text(para.center[1]-1.5,para.center[2]+6,expression(paste("truth\n",theta["0"])),adj=c(1,-0.5),cex=2.5)
text(para.center[1]+1.5,para.center[2]+6,expression(hat(theta)["0"]),adj=c(-0.25,-0.25),cex=2.5)
text(para.center[1],para.center[2]+6/2,"b",adj=c(-2,2.25),cex=2.5)
text(para.center[1],para.center[2]+6/2,"d",adj=c(-0.75,-2.25),cex=2.5)
text(para.center[1]-1.5,para.center[2]+6/2,"h",adj=c(1,0),cex=2.5)
text(para.center[1]+1.5,para.center[2]+6/2,"e",adj=c(-0.5,0),cex=2.5)


text(1,1,expression(paste("Set W(",hat(theta)[k],",",theta[0],") = b")),cex=2.5,adj=c(-0.15,-0.15))
text(11,1,expression(paste("plane ",Theta[k])),cex=2.5,adj=c(1,-0.15))
text(11,6,"MLE with\n approximating\n model",cex=2)
dev.off()

## FigurE 2
library(shape)

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


####  SKIP THIS FIGURE 2A, I DIDN'T PUT IT IN THE PAPER, GO DIRECTLY TO FIGURE 3, LINE 252
####  This is an alternative to figure 2, didn't put it in... 
### Figure 2A
#tiff("~/Desktop/figure2A.tiff",width=20,height=22,units="cm",res=300,compression="lzw",type="cairo")
postscript("~/Desktop/figure2A.eps",width=20,height=22,horizontal=F)
# Panel A
par(mar=c(0,0,0,0))
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
text(center[1]+u[1],center[2]+u[1],expression("d"[4]),adj=c(1.25,-10),cex=2)

center	<-	c(7,4.5)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[1]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[1]),adj=c(2.5,-7),cex=2)


center	<-	c(8,6.5)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[3]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[3]),adj=c(2,-4.5),cex=2)
text(5,10.5,expression("K-L Information"),font=2,cex=3)

dev.off()

### Figure 2B
#tiff("~/Desktop/figure2B.tiff",width=20,height=22,units="cm",res=300,compression="lzw",type="cairo")
postscript("~/Desktop/figure2B.eps",width=20,height=22,horizontal=F)
# Panel A
par(mar=c(0,0,0,0))
plot(0,0,pch=".",col="white",axes=F,xlim=c(1,10),ylim=c(1,11))
plotcircle(mid=c(5,9),r=1,col="grey90",lwd=2)
text(5,9,"g",font=3,cex=2)

plotcircle(mid=c(2,7),r=1,lwd=4,col="grey90")
text(2,7,expression("f"[2]),font=3,cex=2)
v1	<-	c(5,9)-c(2,7)
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=2+u[1],y0=7+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(2+u[1],7+u[1],expression("d"[2]),adj=c(0,-0.75),cex=2)


center	<-	c(5,3)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[4]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[4]),adj=c(1.25,-9),cex=2)


center	<-	c(8,6.5)
plotcircle(mid=center,r=1,lwd=2,col="grey90")
text(center[1],center[2],expression("f"[3]),font=3,cex=2)
v1	<-	c(5,9)-center
u	<-	(v1/(sqrt(sum(v1^2))))*0.9
arrows(x0=center[1]+u[1],y0=center[2]+u[2],x1=5-u[1],y1=9-u[2],lwd=2,length=0.1)
text(center[1]+u[1],center[2]+u[1],expression("d"[3]),adj=c(2,-4),cex=2)
text(5,10.5,expression("K-L Information"),font=2,cex=3)

dev.off()

####### Figure 3
library(shape)

CurlyBraces <- function(x0, x1, y0, y1, pos = 1, direction = 1, depth = 1,lwd=1,lty=1) {

    a=c(1,2,3,48,50)    # set flexion point for spline
    b=c(0,.2,.28,.7,.8) # set depth for spline flexion point

    curve = spline(a, b, n = 50, method = "natural")$y * depth

    curve = c(curve,rev(curve))

    if (pos == 1){
        a_sequence = seq(x0,x1,length=100)
        b_sequence = seq(y0,y1,length=100)  
    }
    if (pos == 2){
        b_sequence = seq(x0,x1,length=100)
        a_sequence = seq(y0,y1,length=100)      
    }

    # direction
    if(direction==1)
        a_sequence = a_sequence+curve
    if(direction==2)
        a_sequence = a_sequence-curve

    # pos
    if(pos==1)
        lines(a_sequence,b_sequence,   xpd=NA,lwd=lwd,lty=lty) # vertical
    if(pos==2)
        lines(b_sequence,a_sequence, xpd=NA,lwd=lwd,lty=lty) # horizontal

}

#tiff("Figure3.tiff",width=31,height=15,units="cm",res=600,compression="lzw",type="cairo") # w/l 22/22
postscript("Figure3.eps",width=31,height=15) #62/30

par(mfrow=c(1,2), mar=c(0,0,2,1), oma=c(0,2,0,4))
#Figure 3a
plot(0,0,pch=".",col="white",xlim=c(-1,17),ylim=c(0,10),axes=F)

# Plane where models are located
poly.mat	<-	matrix(NA,4,2)
poly.mat[1,]	<-	c(-4,1)
poly.mat[2,]	<-	c((1+3),(sqrt(8^2-3^2)+1))
poly.mat[3,]	<-	c(poly.mat[2,1]+19,poly.mat[2,2])
poly.mat[4,]	<-	c(15,1)
polygon(poly.mat,border="grey", lwd=2)

m     <- c(10,5) # 'm' is located at (10,5)
g     <- c(10,8) # location of 'g'
f2    <- c(3.5,5)
f3    <- c(14.2,5)

# the 'h' edge, dotted line
lines(c(g[1],m[1]),c(g[2],m[2]),lty=2)

# horizontal line aligning 'm', 'f2' and 'f3'
# 0.37 is the left border of the plane at y=5 (m's 'y' coordinate)
lines(c(0.37,f3[1]),c(m[2],m[2]),lty=2)

# dotted line marking the 'y1*' ordinate of 'm'
lines(c(m[1],(m[1]-4)),c(m[2],1),lty=2)

# Marking the y1 and y2 star coordinates
text(6,1,expression(paste("y"[1],"*")),adj=c(0.5,1.5),cex=1.5)
text(0.3,5,expression(paste("y"[2],"*")),adj=c(1,0.5),cex=1.5)

# Marking the sub-plot letter
mtext("A",side=3,adj=0.1,cex=2)


# Labeling 'm' and 'h'
text(m[1],m[2],expression(paste("m(","y"[1],"*",",","y"[2],"*",")")),adj=c(-0.03,-0.5),cex=1.2)
text(m[1],(m[2]+1),"h",adj=c(-0.5,-0.5),cex=1.5)

# plotting line 'g-f2' and 'f2'
arrows(x0=f2[1], y0=f2[2], x1=g[1],y1=g[2], lwd=2,length=0)
plotcircle(f2,r=0.6,col="grey90",lwd=2)
text(f2[1],f2[2],expression("f"[2]),font=3,cex=1.5)

# plotting line 'g-f3' and 'f3'
arrows(x0=(f3[1]), y0=f3[2], x1=g[1],y1=g[2], lwd=2,length=0)
plotcircle(f3,r=0.6,col="grey90",lwd=2)
text(f3[1],f3[2],expression("f"[3]),font=3,cex=1.5)

# Plotting 'g'
main <- g
plotcircle(main,r=0.6,col="grey90",lwd=2)
text(main[1],main[2],"g",cex=1.5,font=3)
text(main[1]+4,main[2]+0.2,"Generating model",cex=1.2)
#text(main[1]-10,main[2]+1.2,"A",cex=3)


# Labeling the KL divergences
text((f3[1]-3.1),7.25,expression(paste("KL(g,","f"[3],")","=","S"[gg]-"S"[gf[3]])),adj=c(0,0),cex=1.2)
text((f2[1]+3.2),6.5,expression(paste("KL(g,","f"[2],")","=","S"[gg]-"S"[gf[2]])),adj=c(1,0),cex=1.2)

# Labeling 'm-f3'
CurlyBraces(m[1],(f3[1]-0.6),5,5,pos=2,direction=2,depth=0.5,lwd=1.5)
text(11.5,4.3,expression(paste("d(","f"[3],",m)")),cex=1.2,adj=c(0.5,1))
# Labeling 'm-f2'
CurlyBraces((f2[1]+0.6),m[1],5,5,pos=2,direction=2,depth=0.5,lwd=1.5)
text(7,4.3,expression(paste("d(","f"[2],",m)")),cex=1.2,adj=c(0.5,1))

# Labeling 'f2-f3'
CurlyBraces(f2[1],f3[1],4,4,pos=2,direction=2,depth=0.5,lwd=1.5)
text(8.7,3.6,expression(paste("d(","f"[2],",f"[3],")")),cex=1.2,adj=c(0,1.1))


######-------- Figure 3b
#par(oma=c(0,5,0,6))
plot(0,0,pch=".",col="white",xlim=c(-1,17),ylim=c(0,10),axes=F,oma=c(0,5,0,6))

# Plane where models are located
poly.mat	<-	matrix(NA,4,2)
poly.mat[1,]	<-	c(-4,1)
poly.mat[2,]	<-	c((1+3),(sqrt(8^2-3^2)+1))
poly.mat[3,]	<-	c(poly.mat[2,1]+19,poly.mat[2,2])
poly.mat[4,]	<-	c(15,1)
polygon(poly.mat,border="grey", lwd=2)

m     <- c(4,5) # 'm' is located at (4,5)
g     <- c(4,7.5) # location of 'g'
f2    <- c(11,5)
f3    <- c(15,5)

# the 'h' edge, dotted line
lines(c(g[1],m[1]),c(g[2],m[2]),lty=2)

# horizontal line aligning 'm', 'f2' and 'f3'
# 0.37 is the left border of the plane at y=5 (m's 'y' coordinate)
lines(c(0.37,f3[1]),c(m[2],m[2]),lty=2)

# dotted line marking the 'y1*' ordinate of 'm'
lines(c(m[1],(m[1]-4)),c(m[2],1),lty=2)

# Marking the y1 and y2 star coordinates
text((m[1]-4),1,expression(paste("y"[1],"*")),adj=c(0.5,1.5),cex=1.5)
text(0.3,5,expression(paste("y"[2],"*")),adj=c(1,0.5),cex=1.5)

# Marking the sub-plot letter
mtext("B",side=3,adj=0.1,cex=2)

# Labeling 'm' and 'h'
text(m[1],m[2],expression(paste("m(","y"[1],"*",",","y"[2],"*",")")),adj=c(-0.1,-0.5),cex=1.2)
text(m[1],(m[2]+0.5),"h",adj=c(1.2,-0.5),cex=1.5)

# plotting line 'g-f2' and 'f2'
arrows(x0=f2[1], y0=f2[2], x1=g[1],y1=g[2], lwd=2,length=0)
plotcircle(f2,r=0.6,col="grey90",lwd=2)
text(f2[1],f2[2],expression("f"[2]),font=3,cex=1.5)

# plotting line 'g-f3' and 'f3'
arrows(x0=(f3[1]), y0=f3[2], x1=g[1],y1=g[2], lwd=2,length=0)
plotcircle(f3,r=0.6,col="grey90",lwd=2)
text(f3[1],f3[2],expression("f"[3]),font=3,cex=1.5)

# Plotting 'g'
main <- g
plotcircle(main,r=0.6,col="grey90",lwd=2)
text(main[1],main[2],"g",cex=1.5,font=3)
text(main[1]+4,main[2]+0.2,"Generating model",cex=1.2)#,adj=c(0.5,-2)

# Labeling the KL divergences
text((f3[1]-5.5),6.5,expression(paste("KL(g,","f"[3],")","=","S"[gg]-"S"[gf[3]])),adj=c(0,0),cex=1.2)
text((f2[1]-5.25),6.5,expression(paste("KL(g,","f"[2],")","=","S"[gg]-"S"[gf[2]])),adj=c(1,0),cex=1.2)

# Labeling 'm-f2'
CurlyBraces((f2[1]-0.6),m[1],5,5,pos=2,direction=2,depth=0.5,lwd=1.5)
text(7,4.3,expression(paste("d(","f"[2],",m)")),cex=1.2,adj=c(0.5,1))

# Labeling 'm-f3'
CurlyBraces(m[1],(f3[1]-0.6),4,4,pos=2,direction=2,depth=0.5,lwd=1.5)
text(8.8,3.4,expression(paste("d(","f"[3],",m)")),cex=1.2,adj=c(0.5,1))

# Labeling 'f2-f3'
CurlyBraces((f2[1]+0.6),(f3[1]-0.6),5,5,pos=2,direction=2,depth=0.5,lwd=1.5)
text(12.2,4.3,expression(paste("d(","f"[2],",f"[3],")")),cex=1.2,adj=c(0,1.1))

dev.off()