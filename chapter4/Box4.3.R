###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##BOX 4.3 - running simulation code for table 4.1

##LOAD LIBRARIES
library(phybase)
library(phangorn)

##lists to store the results in
mt.results<-list()
nuc.results<-list()

##ratio of beta.between/beta.within as defined in box 4.3, page 74
##these are the four treatments used in table 4.1
bet.factors<-c(1,2,5,10)

##number of random replicates, set to 1000 to replicate table 4.1
num.reps<-10

##loop through the four treatments defined above 
for (treatments in (1:4)) {

##set up matrices to hold the results for each treatmen
mt.results[[treatments]]<-matrix(NA,nrow=num.reps,ncol=7)
nuc.results[[treatments]]<-matrix(NA,nrow=num.reps,ncol=7)

##now for each treatment loop through the number of replicates
for (rep in (1:num.reps)) {

##sample size per species
	num.indiv <- 10
##effective population size
	Ne <- 1
##divergence times of the two nodes on the species tree
	T<-c(15,20)
##species tree in text newick format
	sp.tree<-paste("((A:",T[1],",B:",T[1],"):",T[2]-T[1],",C:",T[2],");",sep="")
##species tree in phybase format for simulating gene trees
	species.tree<-read.tree.nodes(str=sp.tree)


##A) simulate mtDNA genealogy

ploidy <- 1 #1 for organelle, 4 for nuclear

##THIS SETS THE POPULATION SIZES FOR EACH BRANCH  
	species.tree$nodes[,5]<-Ne*ploidy
##THIS PART SIMULATES THE GENE TREE WITHIN THE SPECIES TREE
	mt.tree<-sim.coaltree.sp(rootnode=5,nodematrix=species.tree$nodes,nspecies=3,seq=rep(num.indiv,3),name=species.tree$names)
	mt.tree.phylo<-read.tree(text=mt.tree$gt)
	mt.tree.phylo$tip.label<-gsub("s","", mt.tree.phylo$tip.label)


##B) Nuclear locus

ploidy <- 4 #1 for organelle, 4 for nuclear

##THIS SETS THE POPULATION SIZES FOR EACH BRANCH
	species.tree$nodes[,5]<-Ne*ploidy
##THIS PART SIMULATES THE GENE TREE WITHIN THE SPECIES TREE
	nuc.tree<-sim.coaltree.sp(rootnode=5,nodematrix=species.tree$nodes,nspecies=3,seq=rep(num.indiv,3),name=species.tree$names)
	nuc.tree.phylo<-read.tree(text=nuc.tree$gt)
	nuc.tree.phylo$tip.label<-gsub("s","", nuc.tree.phylo$tip.label)


##C) Morphological trait

##trait value of ancestor, i.e. root node
	anc.x<-10
##mutational variation per Ne generations
	sigma2.x<-1
##Brownian motion - draw morphology of descendants from normal distribution with variance branch length X sigma2.x - Felsenstein 1985
	species.tree$nodes[nrow(species.tree$nodes),6]<-anc.x
	for (i in ((nrow(species.tree$nodes)-1):1)) {
		species.tree$nodes[i,6]<-rnorm(1,species.tree$nodes[species.tree$nodes[i,1],6],sd=sqrt(sigma2.x*bet.factors[treatments]*species.tree$nodes[i,4]))
	}
	anc.x<-c(anc.x,species.tree$nodes[4,6])

##vector of species means
	sp.x<-species.tree$nodes[species.tree$nodes[,2]==-9,6]
	
##now simulate neutral variation within each species - first assign species mean to each individual
	indiv.x<-rep(sp.x,each=num.indiv)
##then add random variate with variance 2*Ne*sigma2.x - Lynch and Hill 1985
	indiv.x<-rnorm(length(indiv.x),mean=indiv.x,sd=sqrt(2*Ne*sigma2.x))
	names(indiv.x)<-paste(rep(c("A","B","C"),each=num.indiv),rep(1:10,times=3),sep="")
	
##D) Convert genealogies to mt and nuc sequence trees - multiple by molecular mutation rate
		
##scale mt tree to be roughly magnitude assumed for mtDNA markers - theta = few %, i.e. Ne.u for mito = 1%
	mt.tree.phylo$edge.length<-mt.tree.phylo$edge.length/100
##scale nuc tree to be roughly magnitude assumed for nuclear markers - theta = few %, i.e. 4Ne.u for nuclear = 0.1%
	nuc.tree.phylo$edge.length<-nuc.tree.phylo$edge.length/4000

##and estimate mean theta and sd of mean for each species
	mt.dist<-cophenetic.phylo(mt.tree.phylo)
	mt.theta<-matrix(NA,nrow=4,ncol=2)
	colnames(mt.theta) <- c("mean","sd.mean")
	tmp.tally<-NULL
	for (i in (1:3)) {
	tmp<-mt.dist[grep(LETTERS[i],rownames(mt.dist)),grep(LETTERS[i],colnames(mt.dist))]
	tmp<-tmp[upper.tri(tmp)]
	tmp.tally<-c(tmp.tally,tmp)
	mt.theta[i,]<-c(mean(tmp),sqrt(var(tmp)/length(tmp)))}
	mt.theta[4,]<-c(mean(tmp.tally),sqrt(var(tmp.tally)/length(tmp.tally)))

##and for nuclear marker	
	nuc.dist<-cophenetic.phylo(nuc.tree.phylo)
	nuc.theta<-matrix(NA,nrow=4,ncol=2)
	colnames(nuc.theta) <- c("mean","sd.mean")
	tmp.tally<-NULL
	for (i in (1:3)) {
	tmp<-nuc.dist[grep(LETTERS[i],rownames(nuc.dist)),grep(LETTERS[i],colnames(nuc.dist))]
	tmp<-tmp[upper.tri(tmp)]
	tmp.tally<-c(tmp.tally,tmp)
	nuc.theta[i,]<-c(mean(tmp),sqrt(var(tmp)/length(tmp)))}
	nuc.theta[4,]<-c(mean(tmp.tally),sqrt(var(tmp.tally)/length(tmp.tally)))
	
##E) Calculate likelihoods - mtDNA

keep<-NULL
for (i in (1:3)) {
	##find mrca and descendents for each species
		mrcanc<-mrca.phylo(mt.tree.phylo ,grep(LETTERS[i],mt.tree.phylo$tip.label))
		all.desc<-Descendants(mt.tree.phylo,mrcanc,type="all")
		keep<-c(keep,min(all.desc))
		num.tips<-num.indiv*3}
	##mt species tree - remove theta/2 from terminal branches to account for ancestral polymorphism
		mt.sp.tree<-drop.tip(mt.tree.phylo,setdiff(1:num.tips,keep))
		mt.sp.tree$tip.label<-substring(mt.sp.tree$tip.label,1,1)
		terminals<-which(mt.sp.tree$edge[,2]<=3)
		mt.sp.tree$edge.length[terminals]<-mt.sp.tree$edge.length[terminals]-mt.theta[4,1]/2
	##pull out species means
		mt.sp.x<-array(NA,3)
		for (i in (1:3)) {
			mt.sp.x[i]<-mean(indiv.x[grep(mt.sp.tree$tip.label[i],names(indiv.x))])}
	##run ace
		#res<-ace(mt.sp.x,mt.sp.tree,type="continuous",method="ML")  ##ML
		res<-ace(mt.sp.x,mt.sp.tree,type="continuous",method="pic")  ##REML - fix ancestrals states as PIC values

	##make branch matrix
		mt.sp.matrix<-cbind(mt.sp.tree$edge,mt.sp.tree$edge.length,matrix(NA,nrow=4,ncol=2))
		colnames(mt.sp.matrix)<-c("anc","desc","length","x.anc","x.desc")
		mt.sp.matrix[terminals,5]<-mt.sp.x
		mt.sp.matrix[,4]<-res$ace[pmatch(mt.sp.matrix[,1],names(res$ace),duplicates.ok=T)]
		mt.sp.matrix[is.na(mt.sp.matrix[,5]),5]<-res$ace[pmatch(mt.sp.matrix[is.na(mt.sp.matrix[,5]),2],names(res$ace),duplicates.ok=T)]
		mt.sp.matrix<-data.frame(mt.sp.matrix)
		## likelihood function - based on Schluter et al. 1/(B^N) * exp(-Q/2B); that is N=number of nodes, but can substitute with N/2 if N is number of branches
		## https://www.statlect.com/fundamentals-of-statistics/normal-distribution-maximum-likelihood
 		## checked and consistent with ace in R: and correct with likelihood function for normal variate
		##function to calculate likelihood for within-species branches
		lik.between<-function(x) {
		return(sum(with(mt.sp.matrix,-0.5*log(x)-((x.anc-x.desc)^2)/(2*x*length))))}
		res.between<-optimise(lik.between,c(0.00001,100000),maximum=TRUE)
		log1<-res.between$objective
		sigma.1<-res.between$maximum

	##now make within-species matrix
		indiv.x<-indiv.x[pmatch(mt.tree.phylo$tip.label, names(indiv.x))]
		mt.within.matrix<-cbind(mt.tree.phylo$edge,matrix(NA,nrow=nrow(mt.tree.phylo$edge),ncol=3))
		colnames(mt.within.matrix)<-c("anc","desc","length","x.anc","x.desc")
		mt.within.matrix[mt.within.matrix[,2]<=num.tips,5]<-indiv.x
			for (i in (1:3)) {
			mrcanc<-mrca.phylo(mt.tree.phylo ,grep(LETTERS[i],mt.tree.phylo$tip.label))
			desc<-unlist(Descendants(mt.tree.phylo,mrcanc,type="all"))
			mt.within.matrix[mt.within.matrix[,2]%in%desc,4] <- mt.sp.x[pmatch(LETTERS[i],mt.sp.tree$tip.label)]
			mt.within.matrix[mt.within.matrix[,2]%in%desc,3]<-rnorm(sum(mt.within.matrix[,2]%in%desc),2*mt.theta[4,1],sd=2*mt.theta[4,2])
			}
		mt.within.matrix<-data.frame(na.omit(mt.within.matrix))

	##function to calculate likelihood for within-species branches
		lik.within<-function(x) {
		return(sum(with(mt.within.matrix,-0.5*log(x)-((x.anc-x.desc)^2)/(2*x*length))))}
	##optimise sigma2 for within species
		res.within<-optimise(lik.within,c(0.00001,100000),maximum=TRUE)
		log2<-res.within$objective
		sigma.2<-res.within$maximum
	##repeat to estimate single sigma2 for within and between
		mt.matrix<-rbind(mt.sp.matrix,mt.within.matrix)
		lik.all<-function(x) {
			return(sum(with(mt.matrix,-0.5*log(x)-((x.anc-x.desc)^2)/(2*x*length))))}
		res.all<-optimise(lik.all,c(0.00001,100000),maximum=TRUE)
		log.all<-res.all$objective
		sigma.all<-res.all$maximum
				
mt.results[[treatments]][rep,1:7]<-c(sigma.all,log.all,-2*log.all+2*1*(nrow(mt.matrix)/(nrow(mt.matrix)-1-1)),sigma.1,sigma.2,log1+log2,-2*(log1+log2)+2*2*(nrow(mt.matrix) / ( nrow(mt.matrix)-2-1)))

##nuclear likelihoods

keep<-NULL
for (i in (1:3)) {
	##find mrca and descendents for each species - made a fix here for when not monophyletic
		#mrcanc<-mrca.phylo(nuc.tree.phylo ,grep(LETTERS[i],nuc.tree.phylo$tip.label))
		#all.desc<-Descendants(nuc.tree.phylo,mrcanc,type="all")
		#keep<-c(keep,min(all.desc))
		num.tips<-num.indiv*3
		keep<-c(keep,grep(LETTERS[i],nuc.tree.phylo$tip.label)[1])
		}
	##nuc species tree - remove theta/2 from terminal branches to account for ancestral polymorphism
		nuc.sp.tree<-drop.tip(nuc.tree.phylo,setdiff(1:num.tips,keep))
		nuc.sp.tree$tip.label<-substring(nuc.sp.tree$tip.label,1,1)
		terminals<-which(nuc.sp.tree$edge[,2]<=3)
		nuc.sp.tree$edge.length[terminals]<-nuc.sp.tree$edge.length[terminals]-nuc.theta[4,1]/2
	##pull out species means
		nuc.sp.x<-array(NA,3)
		for (i in (1:3)) {
			nuc.sp.x[i]<-mean(indiv.x[grep(nuc.sp.tree$tip.label[i],names(indiv.x))])}
	##run ace
		#res<-ace(nuc.sp.x,nuc.sp.tree,type="continuous",method="ML")  ##ML
		res<-ace(nuc.sp.x,nuc.sp.tree,type="continuous",method="pic")  ##REML - fix ancestrals states as PIC values

	##make branch matrix
		nuc.sp.matrix<-cbind(nuc.sp.tree$edge,nuc.sp.tree$edge.length,matrix(NA,nrow=4,ncol=2))
		colnames(nuc.sp.matrix)<-c("anc","desc","length","x.anc","x.desc")
		nuc.sp.matrix[terminals,5]<-nuc.sp.x
		nuc.sp.matrix[,4]<-res$ace[pmatch(nuc.sp.matrix[,1],names(res$ace),duplicates.ok=T)]
		nuc.sp.matrix[is.na(nuc.sp.matrix[,5]),5]<-res$ace[pmatch(nuc.sp.matrix[is.na(nuc.sp.matrix[,5]),2],names(res$ace),duplicates.ok=T)]
		nuc.sp.matrix<-data.frame(nuc.sp.matrix)
		## likelihood function - based on Schluter et al. 1/(B^N) * exp(-Q/2B); that is N=number of nodes, but can substitute with N/2 if N is number of branches
		## https://www.statlect.com/fundamentals-of-statistics/normal-distribution-maximum-likelihood
 		## checked and consistent with ace in R: and correct with likelihood function for normal variate
		##function to calculate likelihood for within-species branches
		lik.between<-function(x) {
		return(sum(with(nuc.sp.matrix,-0.5*log(x)-((x.anc-x.desc)^2)/(2*x*length))))}
		res.between<-optimise(lik.between,c(0.00001,100000),maximum=TRUE)
		log1<-res.between$objective
		sigma.1<-res.between$maximum

	##now make within-species matrix
		indiv.x<-indiv.x[pmatch(nuc.tree.phylo$tip.label, names(indiv.x))]
		nuc.within.matrix<-cbind(nuc.tree.phylo$edge,matrix(NA,nrow=nrow(nuc.tree.phylo$edge),ncol=3))
		colnames(nuc.within.matrix)<-c("anc","desc","length","x.anc","x.desc")
		nuc.within.matrix[nuc.within.matrix[,2]<=num.tips,5]<-indiv.x
			for (i in (1:3)) {
			mrcanc<-mrca.phylo(nuc.tree.phylo ,grep(LETTERS[i],nuc.tree.phylo$tip.label))
			desc<-unlist(Descendants(nuc.tree.phylo,mrcanc,type="all"))
			nuc.within.matrix[nuc.within.matrix[,2]%in%desc,4] <- nuc.sp.x[pmatch(LETTERS[i],nuc.sp.tree$tip.label)]
			nuc.within.matrix[nuc.within.matrix[,2]%in%desc,3]<-rnorm(sum(nuc.within.matrix[,2]%in%desc),0.5*nuc.theta[4,1],sd=0.5*nuc.theta[4,2])
			}
		nuc.within.matrix<-data.frame(na.omit(nuc.within.matrix))

	##function to calculate likelihood for within-species branches
		lik.within<-function(x) {
		return(sum(with(nuc.within.matrix,-0.5*log(x)-((x.anc-x.desc)^2)/(2*x*length))))}
	##optimise sigma2 for within species
		res.within<-optimise(lik.within,c(0.00001,100000),maximum=TRUE)
		log2<-res.within$objective
		sigma.2<-res.within$maximum
	##repeat to estimate single sigma2 for within and between
		nuc.matrix<-rbind(nuc.sp.matrix,nuc.within.matrix)
		lik.all<-function(x) {
			return(sum(with(nuc.matrix,-0.5*log(x)-((x.anc-x.desc)^2)/(2*x*length))))}
		res.all<-optimise(lik.all,c(0.00001,100000),maximum=TRUE)
		log.all<-res.all$objective
		sigma.all<-res.all$maximum
				
	nuc.results[[treatments]][rep,1:7]<-c(sigma.all,log.all,-2*log.all+2*1*(nrow(nuc.matrix) / ( nrow(nuc.matrix) - 1 - 1)),sigma.1,sigma.2,log1+log2,-2*(log1+log2)+2*2*(nrow(nuc.matrix) / ( nrow(nuc.matrix) - 2 - 1)))


}  ##end rep loop

colnames(mt.results[[treatments]])<-c("sigma.joint","lik.joint","AICc.joint","sigma.between","sigma.within","lik.sep","AICc.sep")
colnames(nuc.results[[treatments]])<-c("sigma.joint","lik.joint","AICc.joint","sigma.between","sigma.within","lik.sep","AICc.sep")

}  ##end treatments loop


##Compile table 4.1, first mtDNA result
##mean ratio beta between/within,  CIs, mean deltaAIC, power 
mt.table<-NULL
	i<-1;null.rej<-quantile(mt.results[[i]][,3]-mt.results[[i]][,7],0.95)
print(null.rej)
for (i in (1:4)) {
	mt.table<-rbind(mt.table,c(beta.ratio=mean(mt.results[[i]][,4]/mt.results[[i]][,5]),
	quantile(mt.results[[i]][,4]/mt.results[[i]][,5],c(0.025,0.975)),
	mean.delta.AIC=mean(mt.results[[i]][,3]-mt.results[[i]][,7]),power=sum(mt.results[[i]][,3]-mt.results[[i]][,7]>null.rej)/num.reps
	))
}

##and nuclear result
nuc.table<-NULL
	i<-1;null.rej<-quantile(nuc.results[[i]][,3]-nuc.results[[i]][,7],0.95)
print(null.rej)
for (i in (1:4)) {
	nuc.table<-rbind(nuc.table,c(beta.ratio=mean(nuc.results[[i]][,4]/nuc.results[[i]][,5]),
	quantile(nuc.results[[i]][,4]/nuc.results[[i]][,5],c(0.025,0.975)),
	mean.delta.AIC=mean(nuc.results[[i]][,3]-nuc.results[[i]][,7]),power=sum(nuc.results[[i]][,3]-nuc.results[[i]][,7]>null.rej)/num.reps
	))
}

print(round(cbind(mt.table,nuc.table),2))


results3<-NULL
	i<-1;null.rej<-quantile(mt.results[[i]][,3]-mt.results[[i]][,7]+nuc.results[[i]][,3]-nuc.results[[i]][,7],0.95)
print(null.rej)
for (i in (1:4)) {
	results3<-rbind(results3,
	c(mean(mt.results[[i]][,3]-mt.results[[i]][,7]+nuc.results[[i]][,3]-nuc.results[[i]][,7]),
	sum(mt.results[[i]][,3]-mt.results[[i]][,7]+nuc.results[[i]][,3]-nuc.results[[i]][,7]>null.rej)/num.reps
	))
}

