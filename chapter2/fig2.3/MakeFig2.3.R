


pdf(file="fig2.3.pdf",width=4,height=1)

layout(matrix(1:8,nrow=2))

load("ScenA.rep1.endmatrices.times")

matplot(repro.gen[,tmp],type="p",pch=".",col= rgb(0,0,0,0.01),cex=0.5,xlab="generation",ylab="repro trait",ylim=c(-0.15,0.25))	  
matplot(eco.gen[,tmp],type="p",pch=".",col=rgb(0,0,0,0.01),cex=0.5,xlab="generation",ylab="eco trait",ylim=c(0.15,0.5))

##INDEX OF WHICH POPULATION EACH INDIVIDUAL IS IN
	pops<-rep(1:2,each=500)

tmp<-c(sample(1:500,50),sample(501:1000,50))
cols<-array(NA,100)
cols[pops[tmp]==1]<-rgb(0,0,0,0.01)
cols[pops[tmp]==2]<-rgb(0.5,0.5,0.5,0.01)

load("ScenB.rep1.endmatrices.times")

matplot(repro.gen[,tmp],type="p",pch=".",col= cols,cex=0.5,xlab="generation",ylab="repro trait",ylim=c(-0.15,0.25))	  
matplot(eco.gen[,tmp],type="p",pch=".",col=cols,cex=0.5,xlab="generation",ylab="eco trait",ylim=c(0.15,0.5))

load("ScenC.rep1.endmatrices.times")

matplot(repro.gen[,tmp],type="p",pch=".",col= cols,cex=0.5,xlab="generation",ylab="repro trait",ylim=c(-0.15,0.25))	  
matplot(eco.gen[,tmp],type="p",pch=".",col=cols,cex=0.5,xlab="generation",ylab="eco trait",ylim=c(0.15,0.5))

load("ScenD.rep1.endmatrices.times")

matplot(repro.gen[,tmp],type="p",pch=".",col= cols,cex=0.5,xlab="generation",ylab="repro trait",ylim=c(-0.15,0.25))	  
matplot(eco.gen[,tmp],type="p",pch=".",col=cols,cex=0.5,xlab="generation",ylab="eco trait",ylim=c(0.15,0.5))

dev.off()

