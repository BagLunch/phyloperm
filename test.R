#this will simulate traits along UCLN rates
require(NELSI)
require(geiger)
args <- commandArgs(trailingOnly = TRUE)
tree <- read.tree(toString(args[1]))
#mean.rate <- as.numeric(args[3])
#rate.sd <- as.numeric(args[4])

mol.tree = tree
mol.tree$edge.length = mol.tree$edge.length/18
ucln = simulate.uncor.lnorm(mol.tree,params = list(mean.log = 0.001,sd.log = 0.5))
rate = 0.15
q<-list(rbind(c(-rate, rate), c(rate, -rate)))
char = sim.char(ucln$phylogram,q,model="discrete",n=400)
char = data.frame(char)
outfl = paste(toString(args[1]),".tsv",sep="") 
write.table(char,sep = "\t", row.names=TRUE,col.names=FALSE,quote=FALSE,file = outfl)
