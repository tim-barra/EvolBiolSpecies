###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##Figure 2.3, panel D - assortative mating + geographic isolation + div selection, eco and repro traits are uncorrelated

library(fields)
library(mclust)
#library(splits)


##1) SET GENERAL PARAMETERS

##NUMBER OF GENERATIONS TO RUN SIMULATION
	ngen<-5000
##NUMBER OF REPLICATES OF THE SIMULATION
	nreps<-1
##SEXUAL REPRODUCTION?
	sexual<-TRUE
##TOTAL NUMBER OF INDIVIDUALS
	nindiv<-1000
##NUMBER OF GEOGRAPHICAL REGIONS
	nreg<-2
##NUMBER OF INDIVIDUALS PER SUB-POPULATION
	npop<-nindiv/nreg
##INDEX OF WHICH POPULATION EACH INDIVIDUAL IS IN
	pops<-rep(1:nreg,each=npop)
##DISPERSAL RATE PER INDIVIDUAL PER GENERATION
	disp.rate<-0.000
##FECUNDITY PER INDIVIDUAL
	fecundity<-5
#PLOIDY LEVEL OF NUCLEAR LOCI = 2 FOR DIPLOID
	ncopy<-2

##2) ECOLOGICAL TRAIT PARAMETERS

##NUMBER OF LOCI CONTRIBUTING TO TRAIT DETERMINING FITNESS IN HABITATS
	nloc.eco<-10
##MUTATION RATE
	mu.eco<-0.0005
##VARIANCE OF EFFECT SIZES OF MUTATIONS 
	alpha.sq.eco<-0.0001
##STRENGTH OF PURIFYING SELECTION - IF BOTH VARIANCES = 0.1, MEAN S ~ 30%, IF ALPHASQ=0.01 AND OMEGASQ=0.1, MEAN S~5%
	omega.sq.eco<-0.003
##NUMBER OF HABITATS
	nopt<-1
##OPTIMUM PHENOTYPES IN EACH HABITAT - ROW, IN EACH REGION - COLUMNS
	optimum.pheno<-matrix(c(0.2,0.4),nrow=nopt,ncol=nreg) ##BEST TRAIT VALUE FOR EACH HABITAT
	colnames(optimum.pheno)<-paste("reg",1:nreg,sep=".")
	rownames(optimum.pheno)<-paste("opt",1:nopt,sep=".")
##RELATIVE AMOUNT OF EACH TYPE OF HABITAT IN EACH REGION
	habitat.size<-matrix(1,nrow=nopt,ncol=nreg)   
	colnames(habitat.size)<-paste("reg",1:nreg,sep=".")
	rownames(habitat.size)<-paste("opt",1:nopt,sep=".")
##PHENOTYPIC VARIANCE
	pheno.var.eco<-0.0015
	
	
##3) REPRODUCTIVE TRAIT PARAMETERS

##CORRELATED WITH ECOLOGICAL TRAIT -  5 INDEPENDENT LOCI ADD TO
##THE FIRST 5 LOCI OF THE ECOLOGICAL TRAIT TO GIVE REPRODUCTIVE TRAIT

##NUMBER OF LOCI CONTRIBUTING TO TRAIT DETERMINING FITNESS IN HABITATS
	nloc.repro<-5
##MUTATION RATE
	mu.repro<-0.0005  ##0.0005
##VARIANCE OF EFFECT SIZES OF MUTATIONS 
	alpha.sq.repro<-0.0001
##STRENGTH OF ASSORTATIVE MATING
	omega.sq.repro<-0.003  ##0.003
##PHENOTYPIC VARIATION
	pheno.var.repro<-0.0015


##4) SEQUENCE DATA

##NUMBER OF SITES - ASSUME NEUTRAL, E.G. 3rd position sites = COULD HAVE SLOW AND FAST SITES TO BETTER RECOVER PHYLOGENY
	lseq<-100
## SUBSTITUTION PARAMETERS FOR DNA
	mu<-0.002*4/3	 #scaling factor for simple mutation model


##5) FUNCTION FOR SEXUAL REPRODUCTION, RANDOM ASSORTMENT

genet.trait<-function(trait) {                                       ##trait = eco or repro, matrix of allele values
	##extract number of loci
	nloc<-ncol(trait)/ncopy
	#set up temporary arrays to house new genotypes                  
		temp1<-trait[parents[,1],]                                  
		temp2<-trait[parents[,2],]                                 
		temp<-temp1                                                 			
	##shuffle copies in parent 1 to determine which copy is passed on
		x<-matrix(sample(c(TRUE,FALSE),nloc*nindiv*fecundity,replace=TRUE),ncol=nloc)       
		temp1[,1:nloc][x]<-temp1[,(1+nloc):(nloc*ncopy)][x]                     
  	##shuffle copies in parent 2 to determine which copy is passed on
		x<-matrix(sample(c(TRUE,FALSE),nloc*nindiv*fecundity,replace=TRUE),ncol=nloc)
		temp2[,1:nloc][x]<-temp2[,(1+nloc):(nloc*ncopy)][x]
	##give the copies to the offspring
		temp[,1:nloc]<-temp1[,1:nloc]
		temp[,(1+nloc):(nloc*ncopy)]<-temp2[,1:nloc]			
	return(temp)}
	   
##6) FUNCTION USED TO PICK MALES FOR EACH FEMALE

pick.males<-function(x) {sample(which(pops==pops[x]),fecundity,prob=p.repro[x,which(pops==pops[x])],replace=FALSE)}


##8) SIMULATION STARTS HERE


start.time<-Sys.time()

results.reps<-matrix(NA,nrow=nreps,ncol=7)

for (reps in (1:nreps)) {

disp.gen<-array(0,ngen)
eco.gen<-matrix(NA,nrow=ngen,ncol=nindiv)
repro.gen<-matrix(NA,nrow=ngen,ncol=nindiv)
fitness.gen<-matrix(NA,nrow=ngen,ncol=nindiv)
effectsize.gen<-matrix(NA,nrow=ngen,ncol=nloc.eco)
effectsize.var.gen<-effectsize.gen
num.clust<-array(NA,ngen)

##FILL IN STARTING ARRAYS
load(paste("Scen0.rep",reps,".endmatrices",sep=""))
rownames(eco)<-pops
rownames(repro)<-pops

#mtdna<-matrix(sample(c("a","c","g","t"),lseq,replace=TRUE),nrow=npop*nreg,ncol=lseq,byrow=TRUE)
#rownames(mtdna)<-rep(1:nreg,each=npop)


for (gen in (1:ngen)) {
		
	##a) random mating within each geographical region
	##Each individual lays fecundity eggs fathered by random other individual 
	females<-rep(1:nindiv,each=fecundity)
	##pick random mate based on difference in the reproductive trait phenotype
	tmp<-rdist(rowSums(cbind(repro[,1:5],eco[,1:5],repro[,6:10],eco[,11:15]))+rnorm(nindiv,0,sqrt(pheno.var.repro)))	
	##or, can use difference across all loci - trait average same, but each locus different, will be incompatible still
	#tmp<-as.matrix(abs(dist(repro,upper=T,diag=T)))	
	##matrix of probabilities of interbreeding
	p.repro<-exp(-((tmp)^2)/(2*omega.sq.repro))
	##pick dads using pick.males function - replace=TRUE means polygyny, replace=FALSE = strict monogamy
	males<-as.vector(sapply(1:nindiv,pick.males))                    
	parents<-cbind(females,males)
	
	##b) Inheritance of traits and mtDNA
	eco<-genet.trait(eco)
	repro<-genet.trait(repro)
	#mtdna<-mtdna[parents[,1],]
	
	##c) Add mutation
	muts<-runif(nindiv*fecundity*nloc.eco*ncopy)<mu.eco
	eco[muts]<- rnorm(sum(muts),eco[muts],sqrt(alpha.sq.eco))
	muts<-runif(nindiv*fecundity*nloc.repro*ncopy)<mu.repro                   
	repro[muts]<-rnorm(sum(muts), repro[muts],sqrt(alpha.sq.repro))
	#muts<-runif(nindiv*fecundity*lseq)<mu                               
	#mtdna[muts]<-sample(c("a","c","g","t"),sum(muts),replace=TRUE)
	
	##d) Dispersal - disp.rate offspring enter disperser pool, then refil empty spaces left by them at random, i.e. some return to home region
	disps<-which(runif(nindiv*fecundity)<disp.rate)
	##only bother if there's more than one disperser picked
	if (length(disps)>1) {
		new.place<-sample(disps)
		eco[disps,]<-eco[new.place,]
		repro[disps,]<-repro[new.place,]
		#mtdna[disps,]<-mtdna[new.place,]
		disp.gen[gen]<-sum(rownames(eco)[disps]!=rownames(eco)[new.place])}
	
	##e) Survival - selection happens at this step
	##calculate phenotype of ecological trait +/- phenotypic variation
		pheno<-rowSums(eco)+rnorm(nindiv,0,sqrt(pheno.var.eco)) 
	##fitness of each individual in each habitat based on exp(-((pheno-opt.pheno)^2)/(2*omega.sq.eco))
		fitness<-exp(-(sweep(matrix(optimum.pheno[,as.integer(rownames(eco))],nrow=1),2,pheno,"-")^2)/(2*omega.sq.eco))
	##fitness of each individual in its own region, by multiplying with proportion of habitat available
		fitness<-colSums(fitness*habitat.size[,as.integer(rownames(eco))])
  	##work out the survivors, up to population size of nindiv
		surv<-NULL
		for (i in (1:nreg)) {
			surv<-c(surv,sample(which(rownames(eco)==i),size=npop,prob=fitness[which(rownames(eco)==i)],replace=FALSE))}
	##WARNING - THE ORDER OF INDIVIDUALS AFTER HERE IS NOT RANDOM - HIGHER FITNESS MORE LIKELY WITH LOWER NUMBERS
	##RANDOM SAMPLING NEED TO USE RANDOM RATHER THAN PICKING FIRST 50, FOR EXAMPLE 
	##they will survive
		eco<-eco[surv,]
		repro<-repro[surv,]
	#	mtdna<-mtdna[surv,]

	##f) record stats
		##realized rate of dispersal
		eco.gen[gen,]<-rowSums(eco)
		repro.gen[gen,]<-rowSums(cbind(repro[,1:5],eco[,1:5],repro[,6:10],eco[,11:15]))
		fitness.gen[gen,]<-fitness[surv]
		effectsize.gen[gen,]<-colMeans(eco[,1:nloc.eco])
		effectsize.var.gen[gen,]<-apply(eco[,1:nloc.eco],2,var)
		
	}  ##end of gen loop
	  	
##save objects
file.name<-"ScenD"

save(eco,repro,file=paste(file.name,".rep",reps,".endmatrices",sep=""))
  
jpeg(file=paste(file.name,".rep",reps,".distribs.jpeg",sep=""),width=800,height=800,quality=100)
	  	
par(mfrow=c(2,2))
tmp<-sample(1:ncol(eco.gen),100)
matplot(eco.gen[,tmp],type="p",pch=20,col=rgb(0,0,0,0.01),cex=0.5,xlab="generation",ylab="eco trait")
matplot(repro.gen[,tmp],type="p",pch=20,col=rgb(0,0,0,0.01),cex=0.5,xlab="generation",ylab="repro trait")	  
matplot(fitness.gen[,tmp],type="p",pch=20,col=rgb(0,0,0,0.01),cex=0.5,xlab="generation",ylab="fitness",log="y")	  

dev.off()

jpeg(file=paste(file.name,".rep",reps,".sumstats.jpeg",sep=""),width=1200,height=900,quality=100)

par(mfrow=c(2,4))

plot(rowMeans(eco.gen),type="l",pch=20,col=1,cex=0.5,ylab="mean eco",xlab="generation")
plot(rowMeans(repro.gen),type="l",pch=20,col=1,cex=0.5,ylab="mean repro",xlab="generation")
plot(rowMeans(fitness.gen),type="l",pch=20,col=1,cex=0.5,ylab="mean fitness",xlab="generation")
matplot(effectsize.gen,type="l",lty=1,xlab="generation",ylab="mean effect per locus")
plot(apply(eco.gen,1,var),type="l",pch=20,col=1,cex=0.5,ylab="var eco",xlab="generation")
plot(apply(repro.gen,1,var),type="l",pch=20,col=1,cex=0.5,ylab="var repro",xlab="generation")
plot(apply(fitness.gen,1,var),type="l",pch=20,col=1,cex=0.5,ylab="var fitness",xlab="generation")
matplot(effectsize.var.gen,type="l",lty=1,xlab="generation",ylab="var effect per locus")
#plot((1:ngen)[!is.na(num.clust)],num.clust[!is.na(num.clust)],type="b",pch=20,col=1,cex=0.5,ylab="num clusters",xlab="generation")

dev.off()

save(eco.gen,repro.gen,file=paste(file.name,".rep",reps,".endmatrices.times",sep=""))

results.reps[reps,]<-c(rowMeans(eco.gen)[ngen],apply(eco.gen,1,var)[ngen],rowMeans(repro.gen)[ngen],apply(repro.gen,1,var)[ngen],rowMeans(fitness.gen)[ngen],apply(fitness.gen,1,var)[ngen],num.clust[gen])


}  #end of reps loop

Sys.time()-start.time


cols<-array(NA,100)
cols[pops[tmp]==1]<-rgb(0,0,0,0.01)
cols[pops[tmp]==2]<-rgb(0.5,0.5,0.5,0.01)
jpeg(file=paste(file.name,".rep",reps,".distribs.2tone.jpeg",sep=""),width=800,height=800,quality=100)
par(mfrow=c(2,2))
matplot(eco.gen[,tmp],type="p",pch=20,col=cols,cex=0.5,xlab="generation",ylab="eco trait")
matplot(repro.gen[,tmp],type="p",pch=20,col= cols,cex=0.5,xlab="generation",ylab="repro trait")	  
matplot(fitness.gen[,tmp],type="p",pch=20,col= cols,cex=0.5,xlab="generation",ylab="fitness",log="y")	  
dev.off()
