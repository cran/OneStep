
library(OneStep)

n<-200
theta<-c(0.5,1.5)

M<-100
M <- 5

f <- function()
{
  o.sample<-rbeta(n,shape1=theta[1],shape2=theta[2])
  benchonestep(o.sample, "beta", methods=c("mle", "mme", "onestep")) 
}  

resall <- replicate(M, f())  
  
tabMLE <- resall[1:2,"mle",]
tabME <- resall[1:2,"mme",]
tabLCE <- resall[1:2,"onestep",]


###Visualisation et Information de Fisher

res<-sqrt(n)*(tabMLE-matrix(rep(theta,M),2,M))
res2<-sqrt(n)*(tabME-matrix(rep(theta,M),2,M))
res3<-sqrt(n)*(tabLCE-matrix(rep(theta,M),2,M))


tabtime<-t(resall[3,,])


I<-matrix(0,2,2)
I[1,1]<-trigamma(theta[1]) -trigamma(theta[1]+theta[2])
I[2,1]<--trigamma(theta[1]+theta[2])
I[1,2]<--trigamma(theta[1]+theta[2])
I[2,2]<-trigamma(theta[2]) -trigamma(theta[1]+theta[2])

covlim<-solve(I)

## Moments variance
J11<-theta[2]/(theta[1]+theta[2])^2
J12<--theta[1]/(theta[1]+theta[2])^2
J21<-(2*theta[1]^2*theta[2]+2*theta[1]*theta[2]^2+2*theta[1]*theta[2]+theta[2]^2+theta[2])/(theta[1]+theta[2])^2/(theta[1]+theta[2]+1)^2
J22<--theta[1]*(theta[1]+1)*(2*theta[1]+2*theta[2]+1)/(theta[1]+theta[2])^2/(theta[1]+theta[2]+1)^2

J<-matrix(c(J11,J12,J21,J22),2,2,byrow=TRUE)

A11<-theta[1]*theta[2]/(theta[1]+theta[2])^2/(theta[1]+theta[2]+1)
A12<-2*theta[2]*theta[1]*(theta[1]+1)/(theta[1]+theta[2])^2/(theta[1]+theta[2]+1)/(theta[1]+theta[2]+2)
#A22<-theta[1]*(theta[1]+1)*(2*theta[1]^3+6*theta[1]^2*theta[2]+4*theta[1]*theta[2]^2+14*theta[1]*theta[2]+4*theta[1]^2+4*theta[1]+6*theta[2]^2+6*theta[2])/(theta[1]+theta[2])^2/(theta[1]+theta[2]+1)^2/(theta[1]+theta[2]+2)/(theta[1]+theta[2]+3)
A22<-theta[1]*(theta[1]+1)*(theta[1]+2)*(theta[1]+3)/(theta[1]+theta[2])/(theta[1]+theta[2]+1)/(theta[1]+theta[2]+2)/(theta[1]+theta[2]+3)-theta[1]^2*(theta[1]+1)^2/(theta[1]+theta[2])^2/(theta[1]+theta[2]+1)^2


A<-matrix(c(A11,A12,A12,A22),2,2,byrow=TRUE)
Jinv<-solve(J)
Sigma<-Jinv%*%A%*%t(Jinv)


layout(matrix(1:6,2,3,byrow=TRUE))
hist(res[1,],freq=FALSE,nclass=40,xlim=c(-3,3),ylim=c(0,0.8),main="MLE theta1",xlab="")
x<-seq(-3,3,length=100)
y<-dnorm(x,mean=0,sd=sqrt(covlim[1,1]))
lines(x,y,col="red")
hist(res2[1,],freq=FALSE,nclass=40,xlim=c(-3,3),ylim=c(0,0.8),main="ME theta1",xlab="")
lines(x,y,col="red")
y2<-dnorm(x,mean=0,sd=sqrt(Sigma[1,1]))
lines(x,y2,col="blue")
hist(res3[1,],freq=FALSE,nclass=40,xlim=c(-3,3),ylim=c(0,0.8),main="LCE theta1",xlab="")
lines(x,y,col="red")

hist(res[2,],freq=FALSE,nclass=40,xlim=c(-8,8),ylim=c(0,0.25),main="MLE theta2",xlab="")
x<-seq(-10,10,length=100)
y<-dnorm(x,mean=0,sd=sqrt(covlim[2,2]))
lines(x,y,col="red")
hist(res2[2,],freq=FALSE,nclass=40,xlim=c(-8,8),ylim=c(0,0.25),main="ME theta2",xlab="")
lines(x,y,col="red")
y2<-dnorm(x,mean=0,sd=sqrt(Sigma[2,2]))
lines(x,y2,col="blue")
hist(res3[2,],freq=FALSE,nclass=40,xlim=c(-8,8),ylim=c(0,0.25),main="LCE theta2",xlab="")
y<-dnorm(x,mean=0,sd=sqrt(covlim[2,2]))
lines(x,y,col="red")


