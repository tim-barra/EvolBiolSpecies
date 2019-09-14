###CODE FOR FIGURES IN 'THE EVOLUTIONARY BIOLOGY OF SPECIES' Oxford University Press, 2019.
##Author: Timothy G. Barraclough (t.barraclough@imperial.ac.uk)
##Citation: Timothy G. Barraclough 2019. The Evolutionary Biology of Species. Oxford Series in Ecology and Evolution. Oxford University Press, Oxford, UK.

##Figure 2.4, loss of diversity due to ecological drift

##population size
N<-1000
##deterministic 
curve(0.5*(1-sqrt(1-2*(2*0.1*0.9*exp(-x/N)))),from=0,to=100,ylab="Frequency of rarer species",xlab="Generations since contact",ylim=c(0,0.1))
curve(exp(-0.01*x)*0.1/(1-0.1+0.1*exp(-0.01*x)),from=0,to=100,add=T)
curve(exp(-0.1*x)*0.1/(1-0.1+0.1*exp(-0.1*x)),from=0,to=100,add=T)
curve(exp(-0.5*x)*0.1/(1-0.1+0.1*exp(-0.5*x)),from=0,to=100,add=T)

##This is the same as changes in frequency of two alleles in a population, hence used these expressions:
##deterministic change in p(allele) with selection, taken from http://www.genetics.org/content/186/1/295.full
##stochastic decline with drift, taken from http://www.tiem.utk.edu/~gavrila/583/C3.pdf