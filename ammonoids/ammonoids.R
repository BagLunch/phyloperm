#setwd("G:/My Drive/School/My papers/Arm number and biogeography/Phyloperm/ammonoids") #replace with your WD
#function for calculating phylo sig
moran.i <- function(dat,tr){W <- vcv(tr); diag(W) <- 0; moran.idx(dat,W)}

#reading in and exploring data
mcgowan <- read.csv("Mcgowan2004data.csv")
mcgowan$species <- as.factor(apply(cbind(as.character(mcgowan$Genus),as.character(mcgowan$species)),1,paste,collapse=" "))
mcgowan.cleaned <- mcgowan[!duplicated(mcgowan$species),] #lots of repeats for some reason
plot(mcgowan.cleaned$D,mcgowan.cleaned$W,pch=21,bg=ifelse(mcgowan.cleaned$FAMILY=="USSURITIDAE","black","white"))
lines(seq(0.001,0.6,0.001),1/seq(0.001,0.6,0.001))

#generating a taxonomy tree
frm <- ~SUPERFAMILY/FAMILY/species
at <- as.phylo(frm,mcgowan.cleaned,collapse=FALSE)
at$edge.length <- rep(1,length(at$edge.length))
write.tree(at,file="at.tre")
at <- read.tree(file="at.tre")

#uses signal-based permutations implemented in the phylo.permute function found in another folder
phylo.permute(at,mcgowan.cleaned$D,margin = 0.001)
phylo.permute(at,mcgowan.cleaned$W,margin = 0.001)

plot(phylo.permute(at,mcgowan.cleaned$D,margin = 0.001),phylo.permute(at,mcgowan.cleaned$W,margin = 0.001))
lines(seq(0.001,0.6,0.001),1/seq(0.001,0.6,0.001))

#ordinary permutations
ammo.op <- pbreplicate(1000,sum(sample(mcgowan.cleaned$W) > 1/mcgowan.cleaned$D))
#phylogenetic permutations
ammo.pp <- pbreplicate(1000,sum(mcgowan.cleaned$W > 1/phylo.permute(at,mcgowan.cleaned$D,margin = 0.001)))

#statistical significance
table(replicate(1000,moran.i(sample(mcgowan.cleaned$D),at)) > moran.i(mcgowan.cleaned$D,at))
table(replicate(1000,moran.i(sample(mcgowan.cleaned$W),at)) > moran.i(mcgowan.cleaned$W,at))

#reading in output from matlab functions dealing with minimum bounding triangles
ratios.op <- read.csv("ratios_operm.csv",header=FALSE)$V1
ratios.pp <- read.csv("ratios_phyloperm.csv",header=FALSE)$V1

#exploring morphospace occupation of individual families
cha<-function(x,y){
  chull(x,y)->i
  return(areapl(cbind(x[i],y[i])))}
sfs <- unique(mcgowan.cleaned$SUPERFAMILY)
supfam.chulls <- c();for(i in 1:78){
  supfam.chulls <- c(supfam.chulls,cha(mcgowan.cleaned$D[mcgowan.cleaned$FAMILY==unique(mcgowan.cleaned$FAMILY)[i]],
                                       mcgowan.cleaned$W[mcgowan.cleaned$FAMILY==unique(mcgowan.cleaned$FAMILY)[i]]))}
supfam.chulls <- supfam.chulls/cha(mcgowan.cleaned$D,mcgowan.cleaned$W)

#plotting
#A
#minimum bounding triangle
tri.x <- c(0,0,0.8115,0)
tri.y <- c(4.34,0.84,1.6414,4.34)
dev.off()
plot(mcgowan.cleaned$D,mcgowan.cleaned$W,pch=21,bg="darkorange",xlim=c(0,0.8115),ylim=c(0.84,4.34),xlab="D (I = 0.0717)",ylab="W (I = 0.0493)")
lines(seq(0.001,0.9,0.001),1/seq(0.001,0.9,0.001));
lines(mcgowan.cleaned$D[c(chull(mcgowan.cleaned$D,mcgowan.cleaned$W),316)],mcgowan.cleaned$W[c(chull(mcgowan.cleaned$D,mcgowan.cleaned$W),316)])
lines(tri.x,tri.y)

#B
par(mfrow=c(2,1))
hist(ammo.op,breaks=seq(6,50,2),border = FALSE,freq=FALSE,ylim=c(0,0.11),main="Genera over the W = 1/D line",xlab="")
hist(ammo.pp,breaks=seq(6,50,2),ylim=c(0,0.11),add=TRUE,col = alpha("black",0),freq=FALSE)
lines(c(6,6),c(0,0.1),col="darkorange2",lwd=4)

#C
hist(ratios.op,breaks=seq(0.5,0.9,0.02),border = FALSE,freq=FALSE,ylim=c(0,9),main="Triangularity",xlab="")
hist(ratios.pp,breaks=seq(0.5,0.9,0.02),ylim=c(0,9),add=TRUE,col = alpha("black",0),freq=FALSE)
lines(c(0.8535,0.8535),c(0,8),col="darkorange2",lwd=4)
