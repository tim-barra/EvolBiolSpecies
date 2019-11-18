###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##Figure 3.5

##LOAD LIBRARIES
library(phybase)
library(phangorn)


##FIGURE 3.5 top panel

div.times<-c(0.4,0.5,0.9,3)
plot.trees<-list()
par(mfrow=c(2,2))

for (i in (1:4)) {

div.time<-div.times[i]

##A) SIMULATE SINGLE GENE TREE FOR A 2-SPECIES DIVERGENCE
##THIS IS THE SPECIES TREE - DIVERGENCE TIME = 2, SPECIES A and SPECIES B
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

plot.trees[[i]]<-gene.tree.phylo
#plot(gene.tree.phylo)
}

##THIS NEEDS MANIPULATING A BIT TO BE IN THE SAME ORDER AS THE BOOK FIGURE 3.5
kronoviz(plot.trees,type="cladogram",horiz=F)



##FIGURE 3.5 middle and bottom panels

##Simulate probability of fixed difference, monophyly, reciprocal monophyly under
##the alternative model, of two species, each with 10 individuals sampled
##across a range of divergence times between the two species

results<-matrix(NA,nrow=4000,ncol=6)
for (i in (1:4000)) {

div.time<-i/1000

##A) SIMULATE SINGLE GENE TREE FOR A 2-SPECIES DIVERGENCE
##THIS IS THE SPECIES TREE - DIVERGENCE TIME = 2, SPECIES A and SPECIES B
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
#plot(gene.tree.phylo)

seqs<-simSeq(gene.tree.phylo,l=660,rate=0.05)
seq.mat<-as.character(seqs)[order(rownames(as.character(seqs))),]

check.fixed<-function(x) {
	sum(seq.mat[1:10,x]%in%seq.mat[11:20,x])}

##a) Know species groups
#seq.dist<-dist.dna(as.DNAbin(seq.mat),"raw", as.matrix=T)
#within.a<-seq.dist[1:10,1:10]
#within.a<-within.a[upper.tri(within.a)]
#within.b<-seq.dist[11:20,11:20]
#within.b<-within.b[upper.tri(within.b)]
#total<-seq.dist[upper.tri(seq.dist)]
#between<-seq.dist[1:10,11:20]

##b) Unknown species groups - look at basal sister clades
	seq.dist<-dist.dna(as.DNAbin(seq.mat),"raw", as.matrix=T)
	recon.tree<-upgma(seq.dist)
	seq.dist<-cophenetic(recon.tree)
	desc.nodes<-recon.tree$edge[recon.tree$edge[,1]==21,2]
	desc.a<-Descendants(recon.tree,desc.nodes[1])
	desc.b<-Descendants(recon.tree,desc.nodes[2])
	within.a<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.a)]]
	within.a<-within.a[upper.tri(within.a)]
	within.b<-seq.dist[rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)],rownames(seq.dist)%in%recon.tree$tip.label[unlist(desc.b)]]
	within.b<-within.b[upper.tri(within.b)]

#phi.stat<-(mean(total)-mean(c(mean(within.a),mean(within.b))))/mean(total)
#k.theta<-mean(between)/(mean(c(mean(within.a),mean(within.b))))

int.branches<-sum(recon.tree$edge.length[recon.tree$edge[,1]==21])
max.upgma<-branching.times(recon.tree)[1]*2
results[i,]<-c(div.time,sum(lapply(1:ncol(seq.mat),check.fixed)==0),
is.monophyletic(recon.tree,tips=paste("A",1:10,sep="")),
is.monophyletic(recon.tree,tips=paste("B",1:10,sep="")),max(seq.dist),mean(c(within.a,within.b)))

}

##quantities of interest
##any fixed differences
some.fixed<-results[,2]>0
##samples are reciprocally monophyletic between two species
recip.mono<-(rowSums(results[,3:4])==2)
##just one of the species samples is monophyletic, other paraphyletic
just.one.mono<-(rowSums(results[,3:4])>=1)
##p.ktheta<-(results[,5]/results[,6])>=quantile(null.res[,6],0.95)


##fit smoothed curves to the metrics to estimate how probs change with divergence time
m1<-glm(some.fixed~results[,1],family="binomial")
m2<-glm(just.one.mono~results[,1],family="binomial")
m3<-glm(recip.mono~results[,1],family="binomial")


##FIGURE 3.5 Plot middle and bottom panel

par(mfrow=c(2,1))

plot(predict(m1,type="response")~results[,1],type="l",col="black",ylab="Probability",xlab="Time since isolation - N generations")
lines(predict(m2,type="response")~results[,1],lty=2)
lines(predict(m3,type="response")~results[,1],lty=3)

legend(x=2.0,y=0.8,lty=c(1,2,3),cex=0.8,legend=c("At least 1 fixed difference","One species monophyletic","Reciprocal monophyly"))

plot(results[,5]~results[,1],pch=20,cex=0.5,col=gray(0.5,0.8),ylab="Average pairwise divergence",xlab="Time since isolation - N generations")
points(results[,6]~results[,1],pch=20,cex=0.5,col=gray(0,0.8))



