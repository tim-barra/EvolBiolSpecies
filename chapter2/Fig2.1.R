###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##FIGURE 2.1 - this code uses simulations to make an illustrative figure, version in book included some manual steps

##load packages
library(phybase)

par(mfrow=c(2,1))

##A)

##make species tree
	tr<-read.tree.nodes(str="(((C:4,B:4):4,(E:5,D:5):3):2,A:10);")
##set population sizes for each tip -  branch lengths of species tree are in units of Ne
	tr$nodes[,5]<-1
##simulate gene tree within this species tree
	simmy<-sim.coaltree.sp(rootnode=9,nodematrix=tr$nodes,nspecies=5,seq=rep(5,5),name=tr$names)

##extract the tree and plot
tr.ape<-read.tree(text=simmy$gt)
plot(tr.ape,"cladogram",direction="downwards")

##to make the figure I saved 5 versions of this that had species A to E in the same order to overlay manually. Code in chapter 3 shows neater ways to do this.

##B)

##read in species tree in ape format instead of phybase
tr.ape<-read.tree(text="(((C:4,B:4):4,(E:5,D:5):3):2,A:10);")

##simulate two phenotypic traits onto the species tree assuming Brownian motion
x<-rTraitCont(tr.ape,sigma=0.05)
y<-rTraitCont(tr.ape,sigma=0.05)

##assume we have measured 1000 individuals of each species for phenotypic traits
x<-rep(x,each=1000)
y<-rep(y,each=1000)
##with normally distributed variation around the species mean
xx<-NULL
yy<-NULL
for (i in (1:5000)) {
xx<-c(xx,rnorm(n=1,mean=x[i],sd=0.02))
yy<-c(yy,rnorm(n=1,mean=y[i],sd=0.02))
}

##plot morphological traita
plot(xx,yy,xlab="trait1",ylab="trait2",pch=20,col=rgb(0,0,0,0.2))

