###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##FIGURE 3.8 Simulate metrics of isolation between two populations connected by dispersal

##NEED TO INSTALL THESE PACKAGES IF YOU DON'T HAVE THEM
library(ape)
library(fields)
library(phangorn)

##Trimmed down version of simulation model from chapter 2, only simulates mtDNA, i.e. single locus, no phenotypes

##1) SET GENERAL PARAMETERS

##NUMBER OF GENERATIONS TO RUN SIMULATION
	ngen<-1000
##NUMBER OF REPLICATES OF THE SIMULATION
	nreps<-100
##SEXUAL REPRODUCTION?
	sexual<-TRUE
##NUMBER OF GEOGRAPHICAL REGIONS
	nreg<-2
##TOTAL NUMBER OF INDIVIDUALS
	nindiv<-200
##NUMBER OF INDIVIDUALS PER SUB-POPULATION
	npop<-nindiv/nreg
##INDEX OF WHICH POPULATION EACH INDIVIDUAL IS IN
	pops<-rep(1:nreg,each=npop)

##DISPERSAL RATE PER INDIVIDUAL PER GENERATION
	#dispersal.rate<-0.0000
	for (dispersal.rate in c(0,0.002,0.01,0.02,0.1,0.2,1)) {

##4) SEQUENCE DATA

##NUMBER OF SITES - ASSUME NEUTRAL, E.G. 3rd position sites = COULD HAVE SLOW AND FAST SITES TO BETTER RECOVER PHYLOGENY
	lseq<-100
## SUBSTITUTION PARAMETERS FOR DNA
	mu<-0.0002*4/3	 #scaling factor for simple mutation model


##8) SIMULATION STARTS HERE

check.fixed<-function(x) {
	sum(seq.mat[1:10,x]%in%seq.mat[11:20,x])}

##sample every this number of generations
time.sample<-10

start.time<-Sys.time()

##array to store results
results.reps<-array(NA,dim=c(ngen/time.sample,12,nreps))

##length to run before isolation is imposed, i.e. burnin
pre.isolation<-nindiv*4


##loop through replicates of the simulation
for (reps in (1:nreps)) {

##FILL IN STARTING ARRAYS
mtdna<-matrix(sample(c("a","c","g","t"),lseq,replace=TRUE),nrow=npop*nreg,ncol=lseq,byrow=TRUE)
rownames(mtdna)<-rep(1:nreg,each=npop)
disp.gen<-array(0,ngen+pre.isolation)
#results<-matrix(NA,nrow=ngen,ncol=20)
results<-NULL

##loop through generations
for (gen in (1:(ngen+pre.isolation))) {
	
	##c) Add mutation
	muts<-runif(nindiv*lseq)<mu                               
	mtdna[muts]<-sample(c("a","c","g","t"),sum(muts),replace=TRUE)
		
	##e) Survival - selection happens at this step
	 	##work out the survivors, up to population size of nindiv
		surv<-NULL
		for (i in (1:nreg)) {
			surv<-c(surv,sample(which(rownames(mtdna)==i),size=npop,replace=TRUE))}
		##they will survive
			mtdna<-mtdna[surv,]

	##d) Dispersal - disp.rate offspring enter disperser pool, then refil empty spaces left by them at random, i.e. some return to home region
	if (gen<=pre.isolation) disp.rate<-1.0
	if (gen>pre.isolation) disp.rate<-dispersal.rate
	
	disps<-which(runif(nindiv)<disp.rate)
	##only bother if there's more than one disperser picked
	if (length(disps)>1) {
		new.place<-sample(disps)
		mtdna[disps,]<-mtdna[new.place,]
		disp.gen[gen]<-sum(rownames(mtdna)[disps]!=rownames(mtdna)[new.place])
		}
		

	##f) record stats
		##realized rate of dispersal
	
		if ((gen-pre.isolation)%in%(seq(1,ngen,time.sample))) {
	
		which.indiv<-NULL
		for (i in (1:nreg)) {
		which.indiv<-c(which.indiv,sample(which(rownames(mtdna)==i),size=10))}
		seq.mat<-as.DNAbin(mtdna[which.indiv,])
		rownames(seq.mat)<-paste(rownames(seq.mat),".",1:10,sep="")
		seq.dist<-dist.dna(seq.mat,model="raw",as.matrix=T)
	
		recon.tree<-upgma(seq.dist)
		seq.dist<-cophenetic(recon.tree)
		desc.nodes<-recon.tree$edge[recon.tree$edge[,1]==21,2]
		desc.a<-Descendants(recon.tree,desc.nodes[1])
		desc.b<-Descendants(recon.tree,desc.nodes[2])
		within.a<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)]]
		within.a<-within.a[upper.tri(within.a)]
		within.b<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)]]
		within.b<-within.b[upper.tri(within.b)]

		
		results<-rbind(results,c(gen-pre.isolation,disp.gen[gen],mean(within.a),mean(within.b),
						max(seq.dist),mean(c(within.a,within.b)),
						max(seq.dist)/(mean(c(within.a,within.b))),
						sum(lapply(1:ncol(seq.mat),check.fixed)==0),
						is.monophyletic(recon.tree,tips=paste("1",1:10,sep=".")),
						is.monophyletic(recon.tree,tips=paste("2",1:10,sep=".")),
						unlist(lapply(Descendants(recon.tree,recon.tree$edge[recon.tree$edge[,1]==21,2]),length))
						))
		}
		
	}  ##end of gen loop
	  	
	  	colnames(results)<-c("gen","disp.obs","mean.a","mean.b","max.tot","mean.ab","birky.index",
	  						"num.fixed","mono.a","mono.b","num.left","num.right")
	 
	results.reps[,,reps]<-results

	} 
	
		dimnames(results.reps)[[2]]<-c("gen","disp.obs","mean.a","mean.b","max.tot","mean.ab","birky.index",
	  						"num.fixed","mono.a","mono.b","num.left","num.right")
	 
	 	##save out array for each replicate
	 	save(results.reps,file=paste("Sim.mtdna.disp.rate",disp.rate,sep=""))

}  ##end of dispersal loop

Sys.time()-start.time



##RELOAD THE RESULTS FILES TO CONSTRUCT FIGURES - uses smoothing to iron out wrinkles,
##alternative would be to run more replicates above and take averages

file.name<-"Sim.mtdna.disp.rate"

ngen<-1000
params<-c(0,0.002,0.01,0.02,0.1,0.2,1)
time.sample<-10
	nreps<-100

all.results<-list()
mean.results<-list()
for (i in (1:7)) {
	load(paste(file.name,params[i],sep=""))
	new.results.reps<-array(NA,dim=c(ngen/time.sample,13,nreps))
	funky<-function(x) {
		return(apply(results.reps[,9:10,x],1,prod))}
	for (k in (1:100)) {
		new.results.reps[,,k]<-cbind(results.reps[,,k],prob.recip=funky(k))}
		dimnames(new.results.reps)[[2]]<-c("gen","disp.obs","mean.a","mean.b","max.tot","mean.ab","birky.index",
	  						"num.fixed","mono.a","mono.b","num.left","num.right","prob.recip")
	all.results[[i]]<-new.results.reps
	mean.results[[i]]<-apply(new.results.reps,c(1,2),mean)
	results.reps<-NULL
	new.results.reps<-NULL}
	
	
	j<-2
	plot(mean.results[[1]][,1],mean.results[[1]][,j],type="l")
for (i in (2:7)) {
	lines(mean.results[[i]][,1],mean.results[[i]][,j])}

	plot((mean.results[[1]][,16]/mean.results[[1]][,17]))
for (i in (2:5)) {
	lines((mean.results[[i]][,16]/mean.results[[i]][,17]))}

(mean.results[[1]][,16]/mean.results[[1]][,17])

##plot smoothed curves. N=200 in each pop (only females here), hence divide timescale by 200 to give in N generations


smoothy<-function(x,y) {
		tmp<-is.finite(y)
		x<-x[is.finite(y)]
		y<-smooth.spline(y[is.finite(y)])
		return(formula(y~x))}
par(mfrow=c(1,3))

	for (j in c(8,13,7)) {
	x<-mean.results[[1]][,1]
	y<-mean.results[[1]][,j]
	
	if (j==7) ylimm<-c( 3.5,11.5) else ylimm<-c(0,max(y))
	plot(x[is.finite(y)]/200,smooth.spline(y[is.finite(y)],df=6)$y,type="l",
		xlab="N generations", ylab=colnames(mean.results[[1]])[j],xlim=c(0,5.2))
text(max(x[is.finite(y)]/200),y[is.finite(y)][length(y[is.finite(y)])],params[1]/2)
for (i in (2:7)) {
	x<-mean.results[[i]][,1]
	y<-mean.results[[i]][,j]
	lines(x[is.finite(y)]/200,smooth.spline(y[is.finite(y)],df=6)$y)
	text(max(x[is.finite(y)]/200),smooth.spline(y[is.finite(y)],df=6)$y[length(y[is.finite(y)])],params[i]/2)
}}