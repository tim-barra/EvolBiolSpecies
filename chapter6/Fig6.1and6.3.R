###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##Figures 6.1 and 6.3

##RUNS AN INDIVIDUAL-BASED SIMULATION OF ASEXUAL AND SEXUAL POPULATIONS
##FACED WITH ISOLATION/DIVERSIFYING SELECTION 

library(fields)

##1) SET GENERAL PARAMETERS

##36 seconds for 1000 generations, 200 individuals
##26 minutes for 10000 generations and nindiv 1000

##NUMBER OF GENERATIONS TO RUN SIMULATION
	ngen<-800
##NUMBER OF VERSIONS OF THE SIMULATION - 3 columns in fig 6.1, asex, sex easy, sex hard
	nreps<-3
##SEXUAL REPRODUCTION?
	sexual<-TRUE
##NUMBER OF GEOGRAPHICAL REGIONS
	nreg<-2
##TOTAL NUMBER OF INDIVIDUALS
	nindiv<-1000
##NUMBER OF INDIVIDUALS PER SUB-POPULATION
	npop<-nindiv/nreg
##INDEX OF WHICH POPULATION EACH INDIVIDUAL IS IN
	pops<-rep(1:nreg,each=npop)
##DISPERSAL RATE PER INDIVIDUAL PER GENERATION
	#disp.rate<-0.0000
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
	nopt<-2
##OPTIMUM PHENOTYPES IN EACH HABITAT - ROW, IN EACH REGION - COLUMNS
	optimum.pheno<-matrix(c(0.0,0.2,0.0,0.2),nrow=nopt,ncol=nreg) ##BEST TRAIT VALUE FOR EACH HABITAT
	colnames(optimum.pheno)<-paste("reg",1:nreg,sep=".")
	rownames(optimum.pheno)<-paste("opt",1:nopt,sep=".")
##RELATIVE AMOUNT OF EACH TYPE OF HABITAT IN EACH REGION
	habitat.size<-matrix(c(1.0,0.0,0.0,1.0),nrow=nopt,ncol=nreg)   
	colnames(habitat.size)<-paste("reg",1:nreg,sep=".")
	rownames(habitat.size)<-paste("opt",1:nopt,sep=".")
##PHENOTYPIC VARIANCE
	pheno.var.eco<-0.0015
	
	
##3) REPRODUCTIVE TRAIT PARAMETERS

##NUMBER OF LOCI CONTRIBUTING TO TRAIT DETERMINING FITNESS IN HABITATS
	nloc.repro<-10
##MUTATION RATE
	mu.repro<-0.0005
##VARIANCE OF EFFECT SIZES OF MUTATIONS 
	alpha.sq.repro<-0.0001
##STRENGTH OF ASSORTATIVE MATING
	omega.sq.repro<-0.003
##PHENOTYPIC VARIATION
	pheno.var.repro<-0.0015


##4) SEQUENCE DATA

##NUMBER OF SITES - ASSUME NEUTRAL, E.G. 3rd position sites = COULD HAVE SLOW AND FAST SITES TO BETTER RECOVER PHYLOGENY
	#lseq<-100
## SUBSTITUTION PARAMETERS FOR DNA
	#mu<-0.0002*4/3	 #scaling factor for simple mutation model


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

pick.males<-function(x) {sample(which(pops==pops[x]),fecundity,prob=p.repro[x,which(pops==pops[x])],replace=TRUE)}


##8) SIMULATION STARTS HERE

start.time<-Sys.time()

##range of dispersal rates to consider for fig 6.3
disp.rates<-rep(seq(0,1,0.1),each=5)
ndisp<-length(disp.rates)

##object to store results in across replicates
results.reps<-NULL

##loop through dispersal rates
for (j in (1:ndisp)) {

disp.rate<-disp.rates[j]

##loop through number of replicates for each dispersal rate
for (reps in (1:nreps)) {

##starting arrays for recording population metrics over time	
disp.gen<-array(0,ngen)
eco.gen<-matrix(NA,nrow=ngen,ncol=nindiv)
repro.gen<-matrix(NA,nrow=ngen,ncol=nindiv)
fst.gen<-array(NA,ngen)

##FILL IN STARTING ARRAYS FOR PHENOTYPES
eco<-matrix(0.0,ncol=nloc.eco*ncopy,nrow=nindiv) 
rownames(eco)<-rep(1:nreg,each=npop)               ##rows: all individuals in region 1 then region 2 etc.
colnames(eco)<-rep(1:ncopy,each=nloc.eco)          ##columns: allele 1 locus 1 to nloc then allele 2 locus 1 to nloc

repro<-matrix(0,ncol=nloc.repro*ncopy,nrow=nindiv) 
rownames(repro)<-rep(1:nreg,each=npop)               
colnames(repro)<-rep(1:ncopy,each=nloc.repro)             

#mitocondrial marker - not used here
#mtdna<-matrix(sample(c("a","c","g","t"),lseq,replace=TRUE),nrow=npop*nreg,ncol=lseq,byrow=TRUE)
#rownames(mtdna)<-rep(1:nreg,each=npop)

##loop through the number of generations
for (gen in (1:ngen)) {
	
	##a) random mating within each geographical region
	##Each individual lays fecundity eggs fathered by random other individual 
	females<-rep(1:nindiv,each=fecundity)
	if (reps==1) {   
		eco<-eco[females,]
		repro<-repro[females,]} else { 
	##pick random mate based on difference in the reproductive trait phenotype
		if (reps ==2) tmp<-rdist(rowSums(eco)+rnorm(nindiv,0,sqrt(pheno.var.repro)))		
		if (reps ==3) tmp<-rdist(rowSums(repro)+rnorm(nindiv,0,sqrt(pheno.var.repro)))	
	##or, can use difference across all loci - trait average same, but each locus different, will be incompatible still
	#tmp<-as.matrix(abs(dist(repro,upper=T,diag=T)))	
	##matrix of probabilities of interbreeding
		p.repro<-exp(-((tmp)^2)/(2*omega.sq.repro))
	##pick dads using pick.males function - replace=TRUE means polygyny, replace=FALSE = strict monogamy
		males<-as.vector(sapply(1:nindiv,pick.males))                    
		parents<-cbind(females,males)

	##b) Inheritance of traits and mtDNA
		eco<-genet.trait(eco)
		repro<-genet.trait(repro)}
	 	
	##c) Add mutation
	muts<-runif(nindiv*fecundity*nloc.eco*ncopy)<mu.eco
	eco[muts]<- rnorm(sum(muts),eco[muts],sqrt(alpha.sq.eco))
	muts<-runif(nindiv*fecundity*nloc.repro*ncopy)<mu.repro                   
	repro[muts]<-rnorm(sum(muts), repro[muts],sqrt(alpha.sq.repro))
	#muts<-runif(nindiv*fecundity*lseq)<mu                               
	#mtdna[muts]<-sample(c("a","c","g","t"),sum(muts),replace=TRUE)
	
	##e) Survival - selection happens at this step
	##calculate phenotype of ecological trait +/- phenotypic variation
		pheno<-rowSums(eco)+rnorm(nindiv,0,sqrt(pheno.var.eco)) 
	##fitness of each individual in each habitat based on exp(-((pheno-opt.pheno)^2)/(2*omega.sq.eco))
		fitness<-exp(-(sweep(optimum.pheno[,as.integer(rownames(eco))],2,pheno,"-")^2)/(2*omega.sq.eco))
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
		#mtdna<-mtdna[surv,]

	##d) Dispersal - disp.rate offspring enter disperser pool, then refil empty spaces left by them at random, i.e. some return to home region
	disps<-which(runif(nindiv)<disp.rate)
	##only bother if there's more than one disperser picked
	if (length(disps)>1) {
		new.place<-sample(disps)
		eco[disps,]<-eco[new.place,]
		repro[disps,]<-repro[new.place,]
		#mtdna[disps,]<-mtdna[new.place,]
		disp.gen[gen]<-sum(rownames(eco)[disps]!=rownames(eco)[new.place])}

	##f) record stats
		##realized rate of dispersal
		eco.gen[gen,]<-rowSums(eco)
		#repro.gen[gen,]<-rowSums(repro)
		#fitness.gen[gen,]<-fitness[surv]
		#effectsize.gen[gen,]<-colMeans(eco[,1:nloc.eco])
		#effectsize.var.gen[gen,]<-apply(eco[,1:nloc.eco],2,var)
		k<-kmeans(rowSums(eco),2)
		fst.gen[gen]<-k$betweenss/k$totss

	
	}  ##end of gen loop
	  	
##save objects
file.name<-c("Asexual","Sex.Easy","Sex.Hard")[reps]

##PLOT DISTRIBUTIONS OF ECOLOGICAL AND REPRODUCTIVE TRAIT - PANELS IN FIGURE 6.1
jpeg(file=paste(file.name,".disp.",disp.rates[j],".distribs.jpeg",sep=""),width=600,height=400,quality=100)
	  	
par(mfrow=c(1,1))
matplot(eco.gen,type="p",pch=20,col=rgb(0.0,0.0,0.0,0.01),cex=0.5,xlab="generation",ylab="eco trait",ylim=c(-0.1,0.3))
dev.off()

##work out some metrics at the end
k<-kmeans(rowSums(eco),2)
diff.mean<-abs(k$centers[1]-k$centers[2])
fst<-k$betweenss/k$totss
tmp<-summary(aov(rowSums(eco)~pops))
fst2<-tmp[[1]][1,2]/sum(tmp[[1]][1:2,2])

diff.mean2<-abs(diff(by(rowSums(eco),pops,mean)))
time.50<-which(abs(rowMeans(eco.gen[,(1+ nindiv/2):nindiv])-rowMeans(eco.gen[,1:(nindiv/2)]))>(diff.mean2/2))[1]
time.fst<-which(fst.gen>(fst))[1]

##store metrics
results.reps<-rbind(results.reps,c(mode = reps, disp.rate=disp.rate,diff.mean=diff.mean,fst=fst,diff.mean2=diff.mean2,time.50=time.50,time.fst =time.fst))


}  #end of reps loop

} ##end of j

Sys.time()-start.time

##save all results
write.table(results.reps,file="results.reps.txt",row.names=F)

##aggregate results for plotting fig. 6.3
agg.results<-aggregate(results.reps,list(results.reps[,1],results.reps[,2]),FUN=function(x) mean(x,na.rm=T))
se.results<-aggregate(results.reps,list(results.reps[,1],results.reps[,2]),FUN=function(x) sd(x)/sqrt(5))

##FIGURE 6.3

par(mfrow=c(1,2))

plot(agg.results[,6]~agg.results[,2],subset=agg.results[,1]==1,type="b",pch=16,ylab="Pst",xlab="Dispersal probability",ylim=range(agg.results[,6]))
lines(agg.results[,6]~agg.results[,2],subset=agg.results[,1]==2,type="b",pch=1)
lines(agg.results[,6]~agg.results[,2],subset=agg.results[,1]==3,type="b",pch=2)
#arrows(agg.results[agg.results[,1]==1,2],agg.results[agg.results[,1]==1,6]-se.results[agg.results[,1]==1,6],agg.results[agg.results[,1]==1,2],agg.results[agg.results[,1]==1,6]+se.results[agg.results[,1]==1,6],length=0)
#arrows(agg.results[agg.results[,1]==2,2],agg.results[agg.results[,1]==2,6]-se.results[agg.results[,1]==2,6],agg.results[agg.results[,1]==2,2],agg.results[agg.results[,1]==2,6]+se.results[agg.results[,1]==2,6],length=0)

plot(agg.results[,9]~agg.results[,2],subset=agg.results[,1]==1,type="b",pch=16,ylab="Time to 50% final divergence",xlab="Dispersal probability",ylim=range(agg.results[,9]))
lines(agg.results[,9]~agg.results[,2],subset=agg.results[,1]==2,type="b",pch=1)
lines(agg.results[,9]~agg.results[,2],subset=agg.results[,1]==3,type="b",pch=2)
