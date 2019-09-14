###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##Figure 2.2, loss of diversity with random interbreeding

##number of generations to simulate for
ngen<-6

##matrix to store histogram data
h<-matrix(NA,ncol=ngen,nrow=99)

##starting phenotype values: 
##600000 individuals from species 1, mean trait=0; 400000 individuals species 2, mean trait 0.1)
x<-c(rnorm(600000,0,0.01),rnorm(400000,0.1,0.01))

##histogram of trait distribution at start before interbreeding
h[,1]<-hist(x,breaks=seq(-0.05,0.15,length.out=100),plot=F)$counts

##now loop for further generations removing the barrier to interbreeding
for (i in (2:ngen)) {

##random mating, irrespective of individual
father<-sample(x,replace=F)
mother<-sample(x,replace=F)

##additive trait, offspring has average phenotype of parents
  x<-(father+mother)/2

##store histogram of trait variation
h[,i]<-hist(x,breaks=seq(-0.05,0.15,length.out=100),plot=F)$counts

}

##set up plot window
layout(matrix(c(1,2,3,4,5,5,5,5,5,5,5,5),ncol=4,byrow=T))

##plot four panels of loss of bimodal distribution over generations 0, 1, 2 and 5
for (i in c(1,2,3,6)) {
	plot(hist(x,breaks=seq(-0.05,0.15,length.out=100),plot=F)$mids,h[,i],type="l",main=paste("t=",i-1,sep=""),ylab="Abundance",xlab="trait value",ylim=c(0,max(h[,c(1,2,3,6)])))
	lines(c(0,0),c(0,10^8),col="grey")
	lines(c(0.1,0.1),c(0,10^8),col="grey")

	}
	
##plot deterministic loss of rarer species
	curve(0.5^(2^x),from=0,to=5,ylab="Frequency of rarer species",xlab="Generations since contact")
	for (i in seq(0,0.45,0.05)) {
	curve(i^(2^x),from=0,to=5,add=T)}
