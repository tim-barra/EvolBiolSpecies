###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##Fig 4.1

##LOAD LIBRARIES
library(phybase)
library(phangorn)



##B) SIMULATE MULTIPLE LOCI FOR GIVEN SPECIES TREE

##OBJECT TO STORE THE TREES IN 
	tt<-NULL

##set up plotting window
	par(mfrow=c(2,4))

##The divergence times, T
	div.times<-c(0,0.5,5.5,6.7)

##loop through the four divergence times
for (plotty in (1:length(div.times))) {

##the species tree
	species.tree<-read.tree.nodes(str=paste("(A:",div.times[plotty],",B:",div.times[plotty],");",sep=""))
##set population sizes for each tip - expressed in 4XNe, so 4 means divergence times in
##units of Ne X number of generations
	species.tree$nodes[,5]<-4

##loop to simulate gene trees for 5 unlinked loci
	for (i in (1:5)) {

	##simulate gene tree
		gene.tree<-sim.coaltree.sp(rootnode=3,nodematrix=species.tree$nodes,nspecies=2,seq=rep(10,2),name=species.tree$names)
	##convert to phylo object
		gene.tree.phylo<-read.tree(text=gene.tree$gt)
	##trim names to look nicer
		gene.tree.phylo$tip.label<-gsub("s","", gene.tree.phylo$tip.label)
	##assign red and black as colours for each species
		tip.colors<-(substring(gene.tree.phylo$tip.label,1,1)=="A")+1
	##read tree
		tt[[i]]<-gene.tree.phylo
	} #end i loop

##make multi-tree object
	tt<-c(tt[[1]],tt[[2]],tt[[3]],tt[[4]],tt[[5]])
##calculate consensus
	sr<-consensus(tt)
	sr$tip.label<-sort(sr$tip.label)
##plot consensus of the 5 gene rees.
	par(mfg=c(1,plotty))
	plot(consensus(tt),direction="downwards")
##and densiTree
	par(mfg=c(2,plotty))
	densiTree(tt,consensus=consensus(tt),tip.color="red",col="black",width=1,direction = "downwards")

}

