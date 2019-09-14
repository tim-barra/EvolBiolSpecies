###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##Figure 2.7, simulate a phylogeny, i.e. multiple speciation events

library(fields)
library(mclust)
library(splits)
library(lattice)

##1) SET GENERAL PARAMETERS

##NUMBER OF GENERATIONS TO RUN SIMULATION
	ngen<-2500
##NUMBER OF REPLICATES OF THE SIMULATION
	nreps<-1
##SEXUAL REPRODUCTION?
	sexual<-TRUE
##TOTAL NUMBER OF INDIVIDUALS
	nindiv<-1000
##NUMBER OF GEOGRAPHICAL REGIONS
	nreg<-1
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
##STRENGTH OF PURIFYING SELECTION
##IF BOTH VARIANCES = 0.1, MEAN S ~ 30%, IF ALPHASQ=0.01 AND OMEGASQ=0.1, MEAN S~5%
	omega.sq.eco<-0.003
##NUMBER OF HABITATS
	nopt<-1
##OPTIMUM PHENOTYPES IN EACH HABITAT - ROW, IN EACH REGION - COLUMNS
	opty<-0.2
	optimum.pheno<-matrix(opty,nrow=nopt,ncol=nreg) ##BEST TRAIT VALUE FOR EACH HABITAT
	colnames(optimum.pheno)<-paste("reg",1:nreg,sep=".")
	rownames(optimum.pheno)<-paste("opt",1:nopt,sep=".")
##RELATIVE AMOUNT OF EACH TYPE OF HABITAT IN EACH REGION
	habitat.size<-matrix(1,nrow=nopt,ncol=nreg)   
	colnames(habitat.size)<-paste("reg",1:nreg,sep=".")
	rownames(habitat.size)<-paste("opt",1:nopt,sep=".")
##PHENOTYPIC VARIANCE
	pheno.var.eco<-0.0015
	
	
##3) REPRODUCTIVE TRAIT PARAMETERS

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

genet.trait<-function(trait) {                     ##trait = eco or repro, matrix of allele values
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

#for (reps in (1:nreps)) {
	reps<-1

disp.gen<-array(0,ngen)
eco.gen<-matrix(NA,nrow=ngen,ncol=nindiv)
repro.gen<-matrix(NA,nrow=ngen,ncol=nindiv)
fitness.gen<-matrix(NA,nrow=ngen,ncol=nindiv)
effectsize.gen<-matrix(NA,nrow=ngen,ncol=nloc.eco)
effectsize.var.gen<-effectsize.gen
num.clust<-array(NA,ngen)

##FILL IN STARTING ARRAYS
load(paste("Scen0.rep",reps,".endmatrices",sep=""))
eco<-eco[1:nindiv,]
repro<-repro[1:nindiv,]
rownames(eco)<-pops
rownames(repro)<-pops

#mtdna<-matrix(sample(c("a","c","g","t"),lseq,replace=TRUE),nrow=npop*nreg,ncol=lseq,byrow=TRUE)
#rownames(mtdna)<-rep(1:nreg,each=npop)

##THIS MATRIX SPECIFIES THE INDIVIDUAL SLOTS THAT WILL BELONG TO EACH SPECIES AT THE END
##I.E. MODEL OF SUBDIVISION OF AN ANCESTRAL POPULATION, TOTAL NUMBER OF INDIVIDUALS
##REMAINS FIXED
split.pops<-matrix(1,nrow=5,ncol=nindiv)
split.pops[1,]<-1
split.pops[2:5,501:1000]<-2
split.pops[3:5,1:150]<-3
split.pops[4:5,151:300]<-4
split.pops[5,751:1000]<-5

##THESE ARE THE PHENOTYPIC OPTIMA FOR EACH FINAL SPECIES
opties<-c(0.2,0.5,0.1,0.35,0.45)


for (gen in (1:ngen)) {

		
##THIS SPECIFIES THE SPECIATION TIMES AND IMPLEMENTS THE NEXT SPLIT
	if (gen%in%c(500,1200,1500,2000)) {
	
	##NUMBER OF GEOGRAPHICAL REGIONS
	nreg<-nreg+1
	##NUMBER OF INDIVIDUALS PER SUB-POPULATION
	npop<-nindiv/nreg
	##INDEX OF WHICH POPULATION EACH INDIVIDUAL IS IN
	pops<-split.pops[nreg,]

	##DISPERSAL RATE PER INDIVIDUAL PER GENERATION
	disp.rate<-0.000
	##NUMBER OF HABITATS
	nopt<-1
	##OPTIMUM PHENOTYPES IN EACH HABITAT - ROW, IN EACH REGION - COLUMNS
	opty<-opties[1:nreg]
	optimum.pheno<-matrix(opty,nrow=nopt,ncol=nreg) ##BEST TRAIT VALUE FOR EACH HABITAT
	colnames(optimum.pheno)<-paste("reg",1:nreg,sep=".")
	rownames(optimum.pheno)<-paste("opt",1:nopt,sep=".")
	##RELATIVE AMOUNT OF EACH TYPE OF HABITAT IN EACH REGION
	habitat.size<-matrix(1,nrow=nopt,ncol=nreg)   
	colnames(habitat.size)<-paste("reg",1:nreg,sep=".")
	rownames(habitat.size)<-paste("opt",1:nopt,sep=".")
	
	rownames(eco)<-pops
	rownames(repro)<-pops
	
	}
	
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
		for (i in unique(rownames(eco))) {
			surv<-c(surv,sample(which(rownames(eco)==i),size=sum(rownames(eco)==i)/fecundity,prob=fitness[which(rownames(eco)==i)],replace=FALSE))}
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
file.name<-"Phylog"

save(eco,repro,file=paste(file.name,".rep",reps,".endmatrices",sep=""))
save(eco.gen,repro.gen,file=paste(file.name,".rep",reps,".endmatrices.times",sep=""))

#2D plot of eco trait over time
#matplot(eco.gen,type="p",pch=".",col=rgb(0,0,0,0.05),ylab="eco trait",xlab="Generations")	  	

##A - 3D plot of repro and eco variation over time

tmp<-NULL
for (i in (1:5)) {
tmp<-c(tmp,sample(which(rownames(eco)==i),size=50))}

jpeg(file=paste(file.name,".3Dplot.jpeg",sep=""),width=1200,height=600,quality=600)
cols<-eco.gen[,tmp]
	cols<-cols-min(cols)
	cols<-cols*1.75
	cloud(as.vector(repro.gen[,tmp])~rep(1:ngen,ncol(eco.gen))+as.vector(eco.gen[,tmp]),xlab="Generations",ylab="eco trait",zlab="repro trait",col=rgb(cols,cols,cols,0.1),pch=20,cex=0.5,screen=list(z=0),aspect=c(1,4),las=2)
dev.off()


##B  - network type representation 

jpeg(file=paste(file.name,".cobweb.jpeg",sep=""),width=1200,height=600,quality=600)


##calculate probability of reproductive isolation between all pairs of individuals, based on traits		
tmp2<-rdist(repro.gen[ngen,])	
p.repro<-exp(-((tmp2)^2)/(2*omega.sq.repro))
rownames(p.repro)<-rownames(repro)

##take a sample of 50 individuals of each species	
tmp<-NULL
samp<-50
for (i in (1:5)) {
tmp<-c(tmp,sample(which(rownames(eco)==i),size= samp))}

p.repro.tmp<-p.repro[tmp,tmp]

##plot individual location in phenotypic space
plot(eco.gen[ngen,tmp],repro.gen[ngen,tmp],pch=20,cex=0.4,ylab="repro trait",xlab="eco trait")

##between species - just add one line for each pair
eco.mean<-by(eco.gen[ngen,tmp],rownames(eco)[tmp],mean)
repro.mean<-by(repro.gen[ngen,tmp],rownames(repro)[tmp],mean)

	tmp2<-rdist(repro.mean)	
	##matrix of probabilities of interbreeding
	p.repro.between<-exp(-((tmp2)^2)/(2*omega.sq.repro))

	  x1<-matrix(eco.mean,nrow=length(tmp)/samp,,ncol=length(tmp)/samp,byrow=T)
	  x2<-matrix(eco.mean,nrow=length(tmp)/samp,ncol=length(tmp)/samp,byrow=F)
  	  y1<-matrix(repro.mean,nrow=length(tmp)/samp,ncol=length(tmp)/samp,byrow=T)
	  y2<-matrix( repro.mean,nrow=length(tmp)/samp,ncol=length(tmp)/samp,byrow=F)

cols<-p.repro.between
inc<-which(p.repro.between>0.01)

##plot links where p(interbreeding)>0.01
arrows(x1[inc],y1[inc],x2[inc],y2[inc],length=0,col=gray(1-cols[inc]),lwd=3)

whichy<-(cols<1)&(cols>0.05)

##annotate with probabilities
text(x=(x1[whichy]+x2[whichy])/2,y=(y1[whichy]+y2[whichy])/2,labels=round(cols[whichy],2))

##now add links for within-species connections
	 p.repro.tmp<-p.repro[tmp,tmp]
	 cols= p.repro.tmp[upper.tri(p.repro.tmp)]
	 x1<-p.repro.tmp
	 x2<-p.repro.tmp
	 y1<-p.repro.tmp
	 y2<-p.repro.tmp
	 x1[]<-eco.gen[ngen,tmp]
	 x2<-t(x1)	 
	 y1[]<-repro.gen[ngen,tmp]
	 y2<-t(y1)
	 x1<-x1[upper.tri(x1)]
	 	 x2<-x2[upper.tri(x2)]
	 y1<-y1[upper.tri(y1)]
	 y2<-y2[upper.tri(y2)]
##within species
withy<-as.matrix(dist(rownames(p.repro.tmp),upper=T,diag=T))
withy<-withy[upper.tri(withy)]

cols<-p.repro.tmp[upper.tri(p.repro.tmp)]

##add links, with transparency
arrows(x1[withy==0],y1[withy==0],x2[withy==0],y2[withy==0],length=0,col=gray(1-cols[withy==0],0.05))	
	  	
	  	
dev.off()
  	

