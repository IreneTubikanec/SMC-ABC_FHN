
#-----------------------------------------------------------------------------------
# Author: Irene Tubikanec
# Date:  2024-05-24
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Visualizaton of SBP SMC-ABC results for the stochastic FHN model
#-----------------------------------------------------------------------------------

#--------------------------------------------------------
#read results
#--------------------------------------------------------

filename<-"ABC_Results" 

#weights
weights_kept<-as.vector(t(read.table(paste(filename,"/weightskept.txt",sep=""),header=F)))
weights_kept<-weights_kept/sum(weights_kept)

#model parameters
eps_kept<-as.vector(t(read.table(paste(filename,"/epskept.txt",sep=""),header=F)))
gamma_kept<-as.vector(t(read.table(paste(filename,"/gammakept.txt",sep=""),header=F)))
beta_kept<-as.vector(t(read.table(paste(filename,"/betakept.txt",sep=""),header=F)))
sig_kept<-as.vector(t(read.table(paste(filename,"/sigkept.txt",sep=""),header=F)))

#prior distribution
Pr<-matrix(0,nrow=4,ncol=2) 
Pr[1,]<-c(0,0.5) #eps
Pr[2,]<-c(0,6) #gamma
Pr[3,]<-c(0,6) #beta
Pr[4,]<-c(0,1) #sig

#--------------------------------------------------------
#plot results
#--------------------------------------------------------

par(mfrow=c(2,2),mai=c(0.35,0.3,0.25,0.15)) 

#-----------------------------
#eps
#-----------------------------

plot(density(eps_kept,weights=weights_kept,from=Pr[1,1],to=Pr[1,2]),col="blue",xlim=c(Pr[1,1],Pr[1,2]),xlab="",main="")
curve(dunif(x,min=Pr[1,1],max=Pr[1,2]),add=TRUE,col="grey",lty=1)
abline(v=0.1,col="green",lty=2)
mtext(expression(epsilon),side=3,line=-1.7,at=Pr[1,2],cex=1.3)

#-----------------------------
#gamma
#-----------------------------

plot(density(gamma_kept,weights=weights_kept,from=Pr[2,1],to=Pr[2,2]),col="blue",xlim=c(Pr[2,1],Pr[2,2]),xlab="",main="")
curve(dunif(x,min=Pr[2,1],max=Pr[2,2]),add=TRUE,col="grey",lty=1)
abline(v=1.5,col="green",lty=2)
mtext(expression(gamma),side=3,line=-1.7,at=Pr[2,2],cex=1.3)

#-----------------------------
#beta
#-----------------------------

plot(density(beta_kept,weights=weights_kept,from=Pr[3,1],to=Pr[3,2]),col="blue",xlim=c(Pr[3,1],Pr[3,2]),xlab="",main="")
curve(dunif(x,min=Pr[3,1],max=Pr[3,2]),add=TRUE,col="grey",lty=1)
abline(v=0.8,col="green",lty=2)
mtext(expression(beta),side=3,line=-1.7,at=Pr[3,2],cex=1.3)

#-----------------------------
#sig
#-----------------------------

plot(density(sig_kept,weights=weights_kept,from=Pr[4,1],to=Pr[4,2]),col="blue",xlim=c(Pr[4,1],Pr[4,2]),xlab="",main="")
curve(dunif(x,min=Pr[4,1],max=Pr[4,2]),add=TRUE,col="grey",lty=1)
abline(v=0.3,col="green",lty=2)
mtext(expression(sigma),side=3,line=-1.7,at=Pr[4,2],cex=1.3)




