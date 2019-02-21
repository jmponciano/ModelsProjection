# Akaike's Geometry

fig.fname <- "Figure1.tiff"

postscript("Figure1.eps",width=28,height=19)# 50/34
#tiff(fig.fname,width=28,height=19,units="cm",res=600,compression="lzw",type="cairo")
par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(1,1,1,3))

# First drawing
plot(0,0,pch=".",col="white",xlim=c(0.9,14.1),ylim=c(0.9,14.1),axes=F)


# Drawing the plane THETA.K
poly.mat	<-	matrix(NA,4,2)

poly.mat[1,]	<-	c(1,1)
poly.mat[2,]	<-	c((1+3),(sqrt(8^2-3^2)+1))
poly.mat[3,]	<-	c(poly.mat[2,1]+10,poly.mat[2,2])
poly.mat[4,]	<-	c(11,1)
polygon(poly.mat)
text(11,1,expression(paste("plane ",Theta[k])),cex=1.5,adj=c(1,-0.15))

# Locate the center of the plot 
para.center	<-	apply(poly.mat,2,sum)/4

# Plotting the triangle with b=hypothenuse and sides 'a' and 'h'
# Coordinates and drawing of the first triangle, using dotted lines
triangle1.mat	<-	matrix(NA,3,2)
triangle1.mat[1,]	<-	c(para.center[1]-1.5,para.center[2])
triangle1.mat[2,]	<-	c(triangle1.mat[1,1],triangle1.mat[1,2]+6)
triangle1.mat[3,]	<-	c(para.center[1]+1.5,para.center[2])

polygon(triangle1.mat,col="white", lty=2)

# Drawing 'h' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]-1.5), c(para.center[2],para.center[2]+6),lty=1,lwd=2)


# Coordinates of the 4 corners of the polygon:
square.mat	<-	matrix(NA,5,2)

# Theta_0_k
square.mat[1,]	<-	c(para.center[1]-1.5,para.center[2])
# Theta_0
square.mat[2,]	<-	c(square.mat[1,1],square.mat[1,2]+6)
# Theta_hat_k
square.mat[3,]	<-	c(para.center[1]+1.5,para.center[2])
# Theta_0_hat_prime  
square.mat[4,]	    <-	c(para.center[1]+1.5,para.center[2]+6+0.45)
# Theta_0_hat
square.mat[5,]	    <-	c(para.center[1]+1.5,para.center[2]+6)

# Plotting the 4 corners and labelling them:
points(square.mat[1:3,],pch=19)
text(para.center[1]-1.5,para.center[2],expression(theta["0k"]),adj=c(1,1),cex=1.5)
text(para.center[1]+1.5,para.center[2],expression(hat(theta)["k"]),adj=c(-0.25,0.75),cex=1.5)
text(para.center[1]-1.5,para.center[2]+6,expression(paste("Generating model\n",theta[0])),adj=c(1,-0.5),cex=1.5)
#text(para.center[1]+1.5,para.center[2]+6+1,expression(hat(theta)["0"]),adj=c(-0.25,0.6),cex=2.5)

text(para.center[1],para.center[2],"a",adj=c(0,1),cex=1.5)
text(para.center[1],para.center[2]+6/2,"b",adj=c(-2,2.25),cex=1.5)
text(para.center[1]-1.55,para.center[2]+6/2,"h",adj=c(1,0),cex=1.5)
text(11,5.1,"MLE with\n approximating\n model",cex=1.2)
text(1,1,expression(paste("W(",hat(theta)[k],",",theta[0],") = ",b^2)),cex=1.5,adj=c(-0.15,-0.15))

mtext("A",adj=0.1,cex=2)


############################# Second Drawing: #############################
plot(0,0,pch=".",col="white",xlim=c(0.9,14.1),ylim=c(0.9,14.1),axes=F)


# Drawing the plane THETA.K
poly.mat	<-	matrix(NA,4,2)

poly.mat[1,]	<-	c(1,1)
poly.mat[2,]	<-	c((1+3),(sqrt(8^2-3^2)+1))
poly.mat[3,]	<-	c(poly.mat[2,1]+10,poly.mat[2,2])
poly.mat[4,]	<-	c(11,1)
polygon(poly.mat)
text(11,1,expression(paste("plane ",Theta[k])),cex=1.5,adj=c(1,-0.15))

# Locate the center of the plot 
para.center	<-	apply(poly.mat,2,sum)/4

# Plotting the polygon
# Coordinates of the 4 corners of the polygon:
square.mat	<-	matrix(NA,5,2)


# Theta_0_k
square.mat[1,]	<-	c(para.center[1]-1.5,para.center[2])
# Theta_0
square.mat[2,]	<-	c(square.mat[1,1],square.mat[1,2]+6)
# Theta_0_hat_prime  
square.mat[3,]	    <-	c(para.center[1]+1.5,para.center[2]+6+0.45)
# Theta_hat_k
square.mat[4,]	<-	c(para.center[1]+1.5,para.center[2])
# Theta_0_hat
square.mat[5,]	    <-	c(para.center[1]+1.5,para.center[2]+6)

#Plotting the polygon
polygon(square.mat[1:4,],col="white", lty=2)


# Drawing 'h' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]-1.5), c(para.center[2],para.center[2]+6),lty=1,lwd=2)


# Plotting the 4 corners and labelling them:
points(square.mat[1:4,],pch=19)
text(para.center[1]-1.5,para.center[2],expression(theta["0k"]),adj=c(1,1),cex=1.5)
text(para.center[1]+1.5,para.center[2],expression(hat(theta)["k"]),adj=c(-0.25,0.75),cex=1.5)
text(para.center[1]-1.5,para.center[2]+6,expression(paste("Generating model\n",theta[0])),adj=c(1,-0.5),cex=1.5)
text(para.center[1]+1.5,para.center[2]+6+1,expression(hat(theta)["0"]),adj=c(-0.25,0.6),cex=1.5)
text(11.2,para.center[2]+7.2,"MLE with\n generating\n model",cex=1.2)

text(para.center[1],para.center[2],"a",adj=c(0,1),cex=1.5)
text(para.center[1],para.center[2]+6+0.25,"c",adj=c(0,-0.5),cex=1.5)
text(para.center[1],para.center[2]+6/2,"b",adj=c(-2,2.25),cex=1.5)
text(para.center[1],para.center[2]+3.2,"d",adj=c(-0.75,-2.25),cex=1.5)
text(para.center[1]-1.55,para.center[2]+6/2,"h",adj=c(1,0),cex=1.5)
text(para.center[1]+1.5,para.center[2]+6/2,"e",adj=c(-0.5,0),cex=1.5)
text(11,5.1,"MLE with\n approximating\n model",cex=1.2)
text(1,1,expression(paste("W(",hat(theta)[k],",",theta[0],") = ",b^2)),cex=1.5,adj=c(-0.15,-0.15))

# Marking the angle phi:
circ.eq <- function(x,a,b,r){
	
	y <- sqrt(r*r-(x-a)^2) + b
	return(y)
}

a <- para.center[1]-1.5
b <- para.center[2]+6
xs <- seq(a,a+0.55,by=0.01)
ys <- b-(circ.eq(x=xs,a=a,b=b,r=0.55) - b)

points(xs,ys, type="l",col="blue",lwd=2, lty=2)
xs <- seq(a+0.54,a+0.55,by=0.01)
ys <- circ.eq(x=xs,a=a,b=b,r=0.55)
points(xs,ys, type="l",col="blue", lwd=2, lty=2)

text(a+0.62,10.1, expression(phi), cex=1.5)

# Drawing 'd' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]+1.5), c(para.center[2],para.center[2]+6+0.45),lty=2)

# Drawing 'b' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]+1.5), c(para.center[2]+6,para.center[2]),lty=2)

#text(11.25,para.center[2]+2,expression(paste(e^2," ~ LLN likelihood ratio")),cex=1.2)


# Just add to the previous plot the law of cosine:
text(7.6,2.5,expression(paste(d^2," = ",h^2," + ",c^2," - ", 2%.%h%.%c%.%cos(phi))),cex=1.2,adj=c(1,-0.15))

mtext("B",adj=0.1,cex=2)



############################# Third Drawing: #############################
plot(0,0,pch=".",col="white",xlim=c(0.9,14.1),ylim=c(0.9,14.1),axes=F)


# Drawing the plane THETA.K
poly.mat	<-	matrix(NA,4,2)

poly.mat[1,]	<-	c(1,1)
poly.mat[2,]	<-	c((1+3),(sqrt(8^2-3^2)+1))
poly.mat[3,]	<-	c(poly.mat[2,1]+10,poly.mat[2,2])
poly.mat[4,]	<-	c(11,1)
polygon(poly.mat)
text(11,1,expression(paste("plane ",Theta[k])),cex=1.5,adj=c(1,-0.15))

# Locate the center of the plot 
para.center	<-	apply(poly.mat,2,sum)/4

# Plotting the polygon
# Coordinates of the 4 corners of the polygon:
square.mat	<-	matrix(NA,5,2)


# Theta_0_k
square.mat[1,]	<-	c(para.center[1]-1.5,para.center[2])
# Theta_0
square.mat[2,]	<-	c(square.mat[1,1],square.mat[1,2]+6)
# Theta_0_hat_prime  
square.mat[3,]	    <-	c(para.center[1]+1.5,para.center[2]+6+0.45)
# Theta_hat_k
square.mat[4,]	<-	c(para.center[1]+1.5,para.center[2])
# Theta_0_hat
square.mat[5,]	    <-	c(para.center[1]+1.5,para.center[2]+6)

#Plotting the polygon
polygon(square.mat[1:4,],col="white", lty=2)
polygon(rbind(square.mat[1:2,],square.mat[5,],square.mat[4,]),col="white", lty=2)



# Drawing 'h' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]-1.5), c(para.center[2],para.center[2]+6),lty=1,lwd=2)


# Plotting the 4 corners and labelling them:
points(square.mat[1:5,],pch=19)
text(para.center[1]-1.5,para.center[2],expression(theta["0k"]),adj=c(1,1),cex=1.5)
text(para.center[1]+1.5,para.center[2],expression(hat(theta)["k"]),adj=c(-0.25,0.75),cex=1.5)
text(para.center[1]-1.5,para.center[2]+6,expression(paste("Generating model\n",theta[0])),adj=c(1,-0.5),cex=1.5)
text(para.center[1]+1.5,para.center[2]+6+1,expression(hat(theta)["0"]),adj=c(-0.25,0.6),cex=1.5)


text(para.center[1],para.center[2],"a",adj=c(0,1),cex=1.5)
text(para.center[1],para.center[2]+6+0.25,"c",adj=c(0,-0.5),cex=1.5)
text(para.center[1],para.center[2]+6/2,"b",adj=c(-2,2.25),cex=1.5)
text(para.center[1]+0.2,para.center[2]+2.9,"d",adj=c(-0.75,-2.25),cex=1.5)
text(para.center[1]-1.55,para.center[2]+6/2,"h",adj=c(1,0),cex=1.5)
text(para.center[1]+1.5,para.center[2]+6/2,"e",adj=c(-0.5,0),cex=1.5)
text(11,6,expression(paste(e^2," ~ LLN LLR")),cex=1.2) 
text(1,1,expression(paste("W(",hat(theta)[k],",",theta[0],") = ",b^2)),cex=1.5,adj=c(-0.15,-0.15))


# Drawing 'd' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]+1.5), c(para.center[2],para.center[2]+6),lty=2)

# Drawing 'b' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]+1.5), c(para.center[2]+6,para.center[2]),lty=2)


# Marking the angle phi:
circ.eq <- function(x,a,b,r){
	
	y <- sqrt(r*r-(x-a)^2) + b
	return(y)
}

a <- para.center[1]-1.5
b <- para.center[2]+6
xs <- seq(a,a+0.55,by=0.01)
ys <- b-(circ.eq(x=xs,a=a,b=b,r=0.55) - b)

points(xs,ys, type="l",col="blue",lwd=2, lty=1)
#xs <- seq(a+0.54,a+0.55,by=0.01)
#ys <- circ.eq(x=xs,a=a,b=b,r=0.55)
#points(xs,ys, type="l",col="blue", lwd=2, lty=2)

text(a+0.8,9.7, expression(paste(phi," = ", pi/2)), cex=1.15)

text(4.7,2.8,expression(paste(d^2," ","="," ",c^2," + ",h^2)),cex=1.2,adj=c(1,-0.15))

mtext("C",adj=0.1,cex=2)

########################### Fourth Drawing: ##################################
plot(0,0,pch=".",col="white",xlim=c(0.9,14.1),ylim=c(0.9,14.1),axes=F)


# Drawing the plane THETA.K
poly.mat	<-	matrix(NA,4,2)

poly.mat[1,]	<-	c(1,1)
poly.mat[2,]	<-	c((1+3),(sqrt(8^2-3^2)+1))
poly.mat[3,]	<-	c(poly.mat[2,1]+10,poly.mat[2,2])
poly.mat[4,]	<-	c(11,1)
polygon(poly.mat)
text(11,1,expression(paste("plane ",Theta[k])),cex=1.5,adj=c(1,-0.15))

# Locate the center of the plot 
para.center	<-	apply(poly.mat,2,sum)/4

# Plotting the polygon
# Coordinates of the 4 corners of the polygon:
square.mat	<-	matrix(NA,5,2)


# Theta_0_k
square.mat[1,]	<-	c(para.center[1]-1.5,para.center[2])
# Theta_0
square.mat[2,]	<-	c(square.mat[1,1],square.mat[1,2]+6)
# Theta_0_hat_prime  
square.mat[3,]	    <-	c(para.center[1]+1.5,para.center[2]+6+0.45)
# Theta_hat_k
square.mat[4,]	<-	c(para.center[1]+1.5,para.center[2])
# Theta_0_hat
square.mat[5,]	    <-	c(para.center[1]+1.5,para.center[2]+6)

#Plotting the polygon
#polygon(square.mat[1:4,],col="white", lty=2)
polygon(rbind(square.mat[1:2,],square.mat[5,],square.mat[4,]),col="white", lty=2)

# Drawing 'h' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]-1.5), c(para.center[2],para.center[2]+6),lty=1,lwd=2)

# Plotting the 4 corners and labelling them:
points(rbind(square.mat[1:2,],square.mat[5,],square.mat[4,]),pch=19)
text(para.center[1]-1.5,para.center[2],expression(theta["0k"]),adj=c(1,1),cex=1.5)
text(para.center[1]+1.5,para.center[2],expression(hat(theta)["k"]),adj=c(-0.25,0.75),cex=1.5)
text(para.center[1]-1.5,para.center[2]+6,expression(paste("Generating model\n",theta[0])),adj=c(1,-0.5),cex=1.5)
text(para.center[1]+1.5,para.center[2]+6+0.5,expression(hat(theta)["0"]),adj=c(-0.25,0.6),cex=1.5)

text(para.center[1],para.center[2],"a",adj=c(0,1),cex=1.5)
text(para.center[1],para.center[2]+6+0.25,"c",adj=c(0,-0.5),cex=1.5)
text(para.center[1],para.center[2]+6/2,"b",adj=c(-2,2.25),cex=1.5)
text(para.center[1],para.center[2]+3.2,"d",adj=c(-0.75,-2.25),cex=1.5)
text(para.center[1]-1.55,para.center[2]+6/2,"h",adj=c(1,0),cex=1.5)
text(para.center[1]+1.5,para.center[2]+6/2,"e",adj=c(-0.5,0),cex=1.5)
#text(11,4.7,"MLE with\n approximating\n model",cex=1.2)
text(1,1,expression(paste("W(",hat(theta)[k],",",theta[0],") = ",b^2," = ",e^2,"-",2*a^2,"-",c^2)),cex=1.5,adj=c(-0.15,-0.15))

# Drawing 'd' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]+1.5), c(para.center[2],para.center[2]+6),lty=2)

# Drawing 'b' as a line, not dotted
lines(c(para.center[1]-1.5, para.center[1]+1.5), c(para.center[2]+6,para.center[2]),lty=2)


text(11.5,para.center[2]+8,expression(paste(b^2," ","="," ",h^2," + ",a^2,";")),cex=1.2)
text(11.5,para.center[2]+7,expression(paste(e^2," ","="," ",d^2," - ",a^2,";")),cex=1.2)
text(11.5,para.center[2]+6,expression(paste(d^2," ","="," ",c^2," + ",h^2,";")),cex=1.2)
text(11.5,para.center[2]+5,expression(paste(b^2," - ",e^2," ","="," ",2*a^2," - ",c^2,".")),cex=1.2)

mtext("D",adj=0.1,cex=2)

dev.off()



