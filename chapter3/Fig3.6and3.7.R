###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##Box 3.1 Figures 3.6 and 3.7, the expected value of the Birky index (divergence/theta) 
##under a null model that the sample comes from a single population versus
##the alternative model that there are two isolated species

##LOAD LIBRARIES
library(phybase)
library(phangorn)

##First, work out expected stem lengths under null model 

##number of simulations
nsamp<-1000

##null distribution for stem lengths
null.res<-matrix(NA,nrow=nsamp,ncol=8)
for (i in (1:nsamp)) {
	tr<-sim.coaltree(20,2)
		tr<-read.tree(text=paste(tr,";",sep=""))
	seq.dist<-cophenetic(tr)
	recon.tree<-tr
	desc.nodes<-recon.tree$edge[recon.tree$edge[,1]==21,2]
	desc.a<-Descendants(recon.tree,desc.nodes[1])
	desc.b<-Descendants(recon.tree,desc.nodes[2])
	within.a<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)]]
	within.a<-within.a[upper.tri(within.a)]
	within.b<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)]]
	within.b<-within.b[upper.tri(within.b)]

null.res[i,1:4]<-c(max(seq.dist),mean(c(mean(within.a),mean(within.b))),mean(c(within.a,within.b)),
min(unlist(lapply(Descendants(recon.tree,recon.tree$edge[recon.tree$edge[,1]==21,2]),length))))
null.res[i,5:6]<-c(null.res[i,1]/null.res[i,2],null.res[i,1]/null.res[i,3])
null.res[i,7]<-sum(recon.tree$edge.length[recon.tree$edge[,1]==21])
null.res[i,8]<-null.res[i,7]/null.res[i,1]

				
}

##K = pairwise divergence across root; theta.1 = mean of thetas in each sub-clade; 
##theta.2 = mean of pairwise differences within subclade 1 and 2; number of individuals
##in smallest subclade; two measures of K/theta ratio using theta 1. and 2 in turn; length of stem branch to sub-clade, proportion as divergence across root.

colnames(null.res)<-c("K","theta.1","theta.2","min.samp","k.theta.1","k.theta.2","int.branch","int.branch.pro")


##NOW SIMULATE THE ALTERNATIVE MODEL FOR VARYING DIVERGENCE TIMES FORM 0 TO 12

##matrix to store summary stats for birky index at each divergent time
results.reps<-matrix(NA,nrow=13,ncol=7)

##number of simulations per divergence time
nsamp<-1000

##step through divergence times
for (rep in (0:12)) {

##alternative distribution for stem lengths
alt.res<-matrix(NA,nrow=nsamp,ncol=3)
colnames(alt.res)<-c("K.theta.2","K","theta.2")

div.time<-rep

for (i in (1:nsamp)) {
	species.tree<-read.tree.nodes(str=paste("(A:",div.time,",B:",div.time,");",sep=""))
##THIS SETS THE POPULATION SIZES FOR EACH BRANCH - BY SETTING 4, THEN DIVERGENCE TIMES
##ABOVE ARE IN UNITS OF N x number of generations (BECAUSE POPULATION SIZE IS GIVEN AS 4Ne)
##set to 1 for haploid uniparental; 4 for nuclear; 
species.tree$nodes[,5]<-1
##THIS PART SIMULATES THE GENE TREE WITHIN THE SPECIES TREE
##10 INDIVIDUALS FOR EACH SPECIES
gene.tree<-sim.coaltree.sp(rootnode=3,nodematrix=species.tree$nodes,nspecies=2,seq=rep(10,2),name=species.tree$names)
gene.tree.phylo<-read.tree(text=gene.tree$gt)
gene.tree.phylo$tip.label<-gsub("s","", gene.tree.phylo$tip.label)

	seq.dist<-cophenetic(gene.tree.phylo)
	recon.tree<-gene.tree.phylo
	desc.nodes<-recon.tree$edge[recon.tree$edge[,1]==21,2]
	desc.a<-Descendants(recon.tree,desc.nodes[1])
	desc.b<-Descendants(recon.tree,desc.nodes[2])
	within.a<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)]]
	within.a<-within.a[upper.tri(within.a)]
	within.b<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)]]
	within.b<-within.b[upper.tri(within.b)]

#alt.res[i,1:4]<-c(max(seq.dist),mean(c(mean(within.a),mean(within.b))),mean(c(within.a,within.b)),
#min(unlist(lapply(Descendants(recon.tree,recon.tree$edge[recon.tree$edge[,1]==21,2]),length))))
#alt.res[i,5:6]<-c(alt.res[i,1]/alt.res[i,2],alt.res[i,1]/alt.res[i,3])
#alt.res[i,7]<-sum(recon.tree$edge.length[recon.tree$edge[,1]==21])
#alt.res[i,8]<-alt.res[i,7]/alt.res[i,1]
alt.res[i,1]<-max(seq.dist)/(mean(c(within.a,within.b)))
alt.res[i,2:3]<-c(max(seq.dist),mean(c(within.a,within.b)))

				
}

results.reps[rep+1,]<-c(div.time,apply(alt.res[,2:3],2,median),quantile(alt.res[,1],c(0.5,0.025,0.975)),sum(alt.res[,1]>=quantile(null.res[,6],0.95,na.rm=T))/nsamp)

}

colnames(results.reps)<-c("div.time","median.K","median.theta","median.Ktheta","lower.Ktheta","higher.Ktheta","p.ktheta")

##PLOT FIGURE 3.6B Density of Birky index under null model
	plot(density(null.res[,6]),xlim=c(0,20),xlab="Birky.index",main="")
	crit<-quantile(null.res[,6],0.95)
	lines(c(crit,crit),c(0,0.5),lty=2)

##PLOT FIGURE 3.7 Birky index under alternative model and probability of exceeding 95% threshold of null
par(mfrow=c(1,2))
matplot(results.reps[,1],results.reps[,4:6],lty=c(1,3,3),col="black",type="l",ylab="Birky index",xlab="Time since isolation - N generations")
abline(quantile(null.res[,6],0.95),0,lty=2)
plot(results.reps[,1],results.reps[,7],ylab="P(Birky index > 10)", xlab="Time since isolation - N generations",type="l")


##results.reps -> estimate of K is actually K+theta, because of over-estimation of divergence time,
##while theta is underestimated because based on a sample. 
##plot(((0:12)+1)/(0.9),results.reps[,2])


