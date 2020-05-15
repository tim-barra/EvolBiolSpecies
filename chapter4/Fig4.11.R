###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##Figure 4.11 


##LOAD LIBRARIES
library(phybase)
library(phangorn)


##matrix to store results across different divergence times
results.divs<-matrix(NA,nrow=40,ncol=16)

##loop to run simulations for each divergence time
for (divs in (1:40)) {

	##matrix to store metrics for properties of trees 
		alt.res<-matrix(NA,nrow=1000,ncol=9)
		colnames(alt.res)<-c("Birky.o","K.o","theta.o","Birky.n","K.n","theta.n","Birky.p","K.p","theta.p")
	
	##params for tree and trait simulation
		num.indiv <- 10
		Ne <- 1
		av.x<-10
		av.y<-15
		sigma2.x<-2
		sigma2.y<-3

	##divergence time
		div.time<-divs/2

##1000 replicates for these conditions
	for (i in (1:1000)) {

	##simulate species tree
	species.tree<-read.tree.nodes(str=paste("(A:",div.time,",B:",div.time,");",sep=""))

	##A) mtDNA

	ploidy <- 1 #1 for organelle, 4 for nuclear

	##THIS SETS THE POPULATION SIZES FOR EACH BRANCH - BY SETTING 4, THEN DIVERGENCE TIMES
	##ABOVE ARE IN UNITS OF N x number of generations (BECAUSE POPULATION SIZE IS GIVEN AS 4Ne)
	##set to 1 for haploid uniparental; 4 for nuclear; 
		species.tree$nodes[,5]<-Ne*ploidy
	##THIS PART SIMULATES THE GENE TREE WITHIN THE SPECIES TREE
	##10 INDIVIDUALS FOR EACH SPECIES
		gene.tree<-sim.coaltree.sp(rootnode=3,nodematrix=species.tree$nodes,nspecies=2,seq=rep(num.indiv,2),name=species.tree$names)
		gene.tree.phylo<-read.tree(text=gene.tree$gt)
		gene.tree.phylo$tip.label<-gsub("s","", gene.tree.phylo$tip.label)

	##extract metrics from the gene tree, variation within and between species
	seq.dist<-cophenetic(gene.tree.phylo)
	recon.tree<-gene.tree.phylo
	desc.nodes<-recon.tree$edge[recon.tree$edge[,1]==21,2]
	desc.a<-Descendants(recon.tree,desc.nodes[1])
	desc.b<-Descendants(recon.tree,desc.nodes[2])
	within.a<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)]]
	within.a<-within.a[upper.tri(within.a)]
	within.b<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)]]
	within.b<-within.b[upper.tri(within.b)]

	##store the relevant metrics
	alt.res[i,1]<-max(seq.dist)/(mean(c(within.a,within.b)))
	alt.res[i,2:3]<-c(max(seq.dist),mean(c(within.a,within.b)))

	##B) Nuclear locus

	ploidy <- 4 #1 for organelle, 4 for nuclear

	##THIS SETS THE POPULATION SIZES FOR EACH BRANCH - BY SETTING 4, THEN DIVERGENCE TIMES
	##ABOVE ARE IN UNITS OF N x number of generations (BECAUSE POPULATION SIZE IS GIVEN AS 4Ne)
	##set to 1 for haploid uniparental; 4 for nuclear; 
		species.tree$nodes[,5]<-Ne*ploidy
	##THIS PART SIMULATES THE GENE TREE WITHIN THE SPECIES TREE
	##10 INDIVIDUALS FOR EACH SPECIES
		gene.tree<-sim.coaltree.sp(rootnode=3,nodematrix=species.tree$nodes,nspecies=2,seq=rep(num.indiv,2),name=species.tree$names)
		gene.tree.phylo<-read.tree(text=gene.tree$gt)
		gene.tree.phylo$tip.label<-gsub("s","", gene.tree.phylo$tip.label)

	##extract metrics from the gene tree, variation within and between species
	seq.dist<-cophenetic(gene.tree.phylo)
	recon.tree<-gene.tree.phylo
	desc.nodes<-recon.tree$edge[recon.tree$edge[,1]==21,2]
	desc.a<-Descendants(recon.tree,desc.nodes[1])
	desc.b<-Descendants(recon.tree,desc.nodes[2])
	within.a<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)]]
	within.a<-within.a[upper.tri(within.a)]
	within.b<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)]]
	within.b<-within.b[upper.tri(within.b)]

	##store the relevant metrics
	alt.res[i,4]<-max(seq.dist)/(mean(c(within.a,within.b)))
	alt.res[i,5:6]<-c(max(seq.dist),mean(c(within.a,within.b)))


##C) morphological simulation
	mean.x<-rep(av.x,num.indiv*2)
	mean.y<-rep(av.y,num.indiv*2)
	##add neutral divergence between 2 species
	mean.x[1:num.indiv]<-rep(av.x,num.indiv)+rnorm(1,0,sqrt(sigma2.x * div.time * Ne * 2))
	mean.y[1:num.indiv]<-rep(av.y,num.indiv)+rnorm(1,0,sqrt(sigma2.y * div.time * Ne * 2))
	##add neutral variation within each species
	x.new<-rnorm(num.indiv*2,mean=mean.x,sd=sqrt(2*Ne*sigma2.x))
	y.new<-rnorm(num.indiv*2,mean=mean.y,sd=sqrt(2*Ne*sigma2.y))
		
##work out morphological Birky index
	morph.dist.x<-as.matrix(dist(x.new))
	within.a<-morph.dist.x[1:10,1:10]
	within.a<-within.a[upper.tri(within.a)]
	within.b<-morph.dist.x[11:20,11:20]
	within.b<-within.b[upper.tri(within.b)]
	between<-morph.dist.x[1:10,11:20]	

	##store metrics	
	alt.res[i,7]<-mean(between)/(mean(c(within.a,within.b)))
	alt.res[i,8:9]<-c(mean(between),mean(c(within.a,within.b)))
	
}  ##end of i loop

##record the median and 95% limits for this divergence time
results.divs[divs,]<-c(div.time,quantile(alt.res[,1],c(0.5,0.025,0.975)),quantile(alt.res[,4],c(0.5,0.025,0.975)),quantile(alt.res[,7],c(0.5,0.025,0.975)),quantile(alt.res[,7]/alt.res[,1],c(0.5,0.025,0.975)),quantile(alt.res[,7]/alt.res[,4],c(0.5,0.025,0.975)))

}  ##end of divs loop

##add column names - o = organelle, n = nuclear, p= phenotype
colnames(results.divs)<-c("div.time","Birky.o","Birky.o.lower","Birky.o.upper","Birky.n","Birky.n.lower","Birky.n.upper",
"Birky.p","Birky.p.lower","Birky.p.upper","Birky.p.o","Birky.p.o.lower","Birky.p.o.upper",
"Birky.p.n","Birky.p.n.lower","Birky.p.n.upper")

##PLOT FIGURE 4.11
par(mfrow=c(1,2))

plot(results.divs[,1],results.divs[,11],type="l",main="A)",ylim=c(0,1),xlab="Divergence time (Ne generations)",ylab="Birky index phenotype/Birky index organelle locus")
lines(results.divs[,1],results.divs[,12],lty=2)
lines(results.divs[,1],results.divs[,13],lty=2)

plot(results.divs[,1],results.divs[,14],type="l",main="B)",ylim=c(0,2),xlab="Divergence time (Ne generations)",ylab="Birky index phenotype/Birky index nuclear locus")
lines(results.divs[,1],results.divs[,15],lty=2)
lines(results.divs[,1],results.divs[,16],lty=2)