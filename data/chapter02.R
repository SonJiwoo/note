#### Fig 2.1 
# pdf("fig2_1.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

x<-0:10
plot(x,dpois(x,2.1), type="h",lwd=1, xlab=expression(italic(y)),ylab=
       #    expression(paste(italic("p(y"),"|",theta==2.1,")")) )
       expression(paste(italic("p"),"(",italic("y"),"|",theta==2.1,")")) )

x<-0:100
plot(x,dpois(x,21), type="h",lwd=1, xlab=expression(italic(y)),ylab=
       expression(paste(italic("p"),"(",italic("y"),"|",theta==21,")")) )
dev.off()


#### Fig 2.2
# pdf("fig2_2.pdf",height=3.5,width=7,family="Times")

mu<-10.75
sig<- .8

par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

x<- seq(7.9,13.9,length=500)
plot(x,pnorm(x,mu,sig),type="l",ylab=expression(paste(italic("F"),"(",italic("y"),")")),xlab=
       expression(italic(y)),lwd=1)
abline(h=c(0,.5,1),col="gray")
plot(x,dnorm(x,mu,sig),type="l",ylab=expression(paste(italic("p"),"(",italic("y"),")")),
     xlab=
       expression(italic(y)),lwd=1)
abline(v=mu,col="gray")
dev.off()


#### Fig 2.3
# pdf("fig2_3.pdf",height=3.5,width=7,family="Times")
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
x<-seq(7.75,13.75,length=100)
mu<-10.75 ; sig<-.8

par(mfrow=c(1,2))
plot(x, dnorm(x,mu,sig),type="l",xlab=expression(italic(y)),
     ylab= expression(paste(italic("p"),"(",italic("y"),")")) )
abline(v=mu,lty=1,col=gray(0))
abline(v=mu,lty=2,col=gray(.33))
abline(v=mu,lty=4,col=gray(.66))

x<-seq( 0,300000,length=200)
mu<-10.75 ; sig<-.8
plot(x, dlnorm(x,mu,sig)*1e5,type="l", xlab=expression(italic(y)),
     ylab= expression( 10^5*paste(italic("p"),"(",italic("y"),")")) )
abline(v=24600,col=gray(0))
abline( v=qlnorm(.5,mu,sig),col=gray(.3),lty=2)
abline(v=exp(mu+.5*sig^2) , col=gray(.7),lty=4)
legend(150000,1.0,c("mode","median","mean"),
       lty=c(1,2,4),col=gray(c(0,.33,.66)),
       bty="n",cex=.85)

dev.off()

#### Happiness
x<-c(467,667,113,25)
xc<-c(x[1]+x[2],x[3]+x[4])
xc[1]/sum(xc) 

