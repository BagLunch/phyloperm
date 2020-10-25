#setup ----
library(phytools);library(adephylo);library(combinat);library(pbapply);library(mvMORPH)
#setwd('G:/My Drive/School/My papers/Arm number and biogeography/Phyloperm') #set to your directory
#small trees of different sizes for playing with different methods (trees in directory methods/trees)
x <- read.tree("x.tre");w <- read.tree("w.tre");big <- read.tree("big.tre");fels <- read.tree("felsentree.tre")
y <- extract.clade(x,19);z <- extract.clade(x,18)

#functions ----
moran.i <- function(dat,tr){W <- vcv(tr); diag(W) <- 0; moran.idx(dat,W)}

#generates cyclic permutations of the data vec on a phylogeny tre
#rotates at nodes N times
cyclic.permutation <- function(tre,vec,N=sample((tre$Nnode*20):(tre$Nnode*40),1)){ #N is the number of permutations
  treperm <- tre
  for(i in 1:N){
    nod <- sample(unique(treperm$edge[,1]),1)
    Ndaughters <- sum(tre$edge[,1] == nod)
    treperm <- rotateNodes(treperm,nod,polytom=sort(sample(seq(Ndaughters),2)))}
  return(vec[match(treperm$tip.label,tre$tip.label)])}

#signal-based phylogenetic permutations, using moran's i
#returns rearrangements of vec whose phylogenetic signal on tre is within margin of the empirical signal
#first, permutes vec
#then, proposes swaps of randomly-selected pairs of observations and accepts if it brings phylogenetic signal closer to empirical signal
#if no swaps are accepted for stucklimit iterations, it restarts
phylo.permute <- function(tre,vec,margin=0.01,stucklimit=1000){
  physig.empirical <- moran.i(vec,tre)
  #print(paste("Empirical phylogenetic signal:",physig.empirical))
  vec.old <- sample(vec) #shuffled vector
  physig.old <- moran.i(vec.old,tre) #phylogenetic signal of shuffled vector
  #print(paste("Permuted phylogenetic signal:",physig.old))
  #initiate the main loop - break on condition that physig.old is sufficiently close to physig.empirical
  stuckcounter <- 0
  while(abs(physig.old - physig.empirical) > margin){
    if(stuckcounter>=stucklimit){ #if the algorithm got stuck
      #print("stuck!")
      stuckcounter <- 0
      vec.old <- sample(vec) #shuffled vector
      physig.old <- moran.i(vec.old,tre)}
    indices <- sample(seq(length(vec)),2)
    vec.new <- replace(vec.old,indices,vec.old[rev(indices)])
    physig.new <- moran.i(vec.new,tre)
    #if this configuration more closely approximates the physig of the empirical...
    if(abs(physig.new - physig.empirical) <= abs(physig.old - physig.empirical)){ 
      vec.old <- vec.new #replace the old value with the new one
      physig.old <- physig.new
      #print(paste("Permuted phylogenetic signal:",physig.new))
    }else{
      stuckcounter <- stuckcounter+1
    }}
  return(vec.old)}

#Lapointe-Garland permutations
#k is the "flattening" parameter from Lapointe & Garland 2001
LG.permute <- function(tre,vec,k){
  pmat <- (k-cophenetic.phylo(tre)/max(cophenetic.phylo(tre)))/sum(cophenetic.phylo(tre)[,1])
  out <- rep(NA,length(vec))
  for(i in seq(length(vec))){
    index <- sample(seq(length(vec)),1,prob=pmat[i,])
    out[index] <- vec[i]
    pmat[,index] <- 0}
  return(out)}

#false positives ----
x <- fastBM(big)
y <-fastBM(big)
replicate(cor.test(x,cyclic.permutation(big,y))$estimate)

ps.CP <- c();for(i in 1:1000){
  print(i)
  x <- fastBM(w)
  y <-fastBM(w)
  ps.CP <- c(ps.CP,sum(pbreplicate(500,abs(cor.test(x,cyclic.permutation(w,y))$estimate))>abs(cor.test(x,y)$estimate))/500)}
write.csv(ps.CP,file="CP_phyloperm_w_false_positives.csv")
ps.CP <- read.csv(file="CP_phyloperm_w_false_positives.csv",header = FALSE)$V1
ks.test(ps.CP,"punif",0,1)

ps.PP <- c();for(i in 1:1000){
  print(i)
  x <- fastBM(w)
  y <-fastBM(w)
  ps.PP <- c(ps.PP,sum(pbreplicate(500,abs(cor.test(x,phylo.permute(w,y,margin=0.001,stucklimit = 100))$estimate))>abs(cor.test(x,y)$estimate))/500)}
write.csv(ps.PP,file="PP_w_false_positives.csv")
ps.PP <- read.csv(file="PP_w_false_positives.csv",header = FALSE)$V1
ks.test(ps.PP,"punif",0,1)


ps.LGP <- c();for(i in 1:1000){
  print(i)
  x <- fastBM(w)
  y <-fastBM(w)
  ps.LGP <- c(ps.LGP,sum(abs(pbreplicate(500,abs(cor.test(x,LG.permute(w,y,1))$estimate)))>abs(cor.test(x,y)$estimate))/500)}
write.csv(ps.LGP,file="LG_phyloperm_w_false_positives.csv")
ps.LGP <- read.csv(file="LG_phyloperm_w_false_positives.csv",header=FALSE)$V1
ks.test(ps.LGP,"punif",0,1)

hist(ps.LGP,breaks=seq(0,1,0.05))
hist(ps.CP,breaks=seq(0,1,0.05))
hist(ps.PP,breaks=seq(0,1,0.05))


#false negatives ----
xy <- mvSIM(w,model="BM1",param = list(ntraits=2,sigma=matrix(c(1,0.75,0.75,1),nrow=2),theta=c(0,0)))
cor.test(xy[,1],xy[,2])$estimate
x1 <- pbreplicate(1000,cor.test(xy[,1],sample(xy[,2]))$estimate)
x2 <- pbreplicate(500,cor.test(xy[,1],cyclic.permutation(w,xy[,2]))$estimate)
x3 <- pbreplicate(1000,cor.test(xy[,1],LG.permute(w,xy[,2],1))$estimate)
table(x1 > cor.test(xy[,1],xy[,2])$estimate)

ps.FN.CP <- c();for(i in 1:1000){ #FN for false negative
  print(i)
  xy <- mvSIM(w,model="BM1",param = list(ntraits=2,sigma=matrix(c(1,0.75,0.75,1),nrow=2),theta=c(0,0)))
  ps.FN.CP <- c(ps.FN.CP,sum(pbreplicate(500,cor.test(xy[,1],cyclic.permutation(w,xy[,2]))$estimate)>cor.test(xy[,1],xy[,2])$estimate)/500)}
write.csv(ps.FN.CP,file="CP_phyloperm_w_false_negatives.csv")
ps.FN.CP <- read.csv(file="CP_phyloperm_w_false_negatives.csv",header = FALSE)$V1


ps.FN.PP <- c();for(i in 1:1000){
  print(i)
  xy <- mvSIM(w,model="BM1",param = list(ntraits=2,sigma=matrix(c(1,0.75,0.75,1),nrow=2),theta=c(0,0)))
  ps.FN.PP <- c(ps.FN.PP,sum(pbreplicate(500,cor.test(xy[,1],phylo.permute(w,xy[,2],margin=0.001,stucklimit = 100))$estimate)>cor.test(xy[,1],xy[,2])$estimate)/500)}
write.csv(ps.FN.PP,file="PP_w_false_negatives.csv")
ps.FN.PP <- read.csv(file="PP_w_false_negatives.csv",header = FALSE)$V1

ps.FN.LGP <- c();for(i in 1:1000){
  print(i)
  xy <- mvSIM(w,model="BM1",param = list(ntraits=2,sigma=matrix(c(1,0.75,0.75,1),nrow=2),theta=c(0,0)))
  ps.FN.LGP <- c(ps.FN.LGP,sum(pbreplicate(500,cor.test(xy[,1],LG.permute(w,xy[,2],1))$estimate)>cor.test(xy[,1],xy[,2])$estimate)/500)}
write.csv(ps.FN.LGP,file="LG_phyloperm_w_false_negatives.csv")
ps.FN.LGP <- read.csv(file="LG_phyloperm_w_false_negatives.csv",header = FALSE)$V1


#PICs
w.d <- multi2di(w,random=FALSE)
ps.PIC <- c();for(i in 1:1000){
  x <- fastBM(w.d)
  y <-fastBM(w.d)
  ps.PIC <- c(ps.PIC,cor.test(pic(x,w.d),pic(y,w.d))$p.value)}
ks.test(ps.PIC,"punif",0,1)

ps.FN.PIC <- c();for(i in 1:1000){
  print(i)
  xy <- mvSIM(w.d,model="BM1",param = list(ntraits=2,sigma=matrix(c(1,0.75,0.75,1),nrow=2),theta=c(0,0)))
  ps.FN.PIC <- c(ps.FN.PIC,cor.test(pic(xy[,1],w.d),pic(xy[,2],w.d))$p.value)}



#Uyeda's worst case----
up.pic <- c(); up.cp <- c(); up.sbp <- c(); for(i in 1:1000){
  print(i)
  x.fels <- fastBM(w.d); y.fels <- fastBM(w.d)
  x.fels <- x.fels + c(rep(0,4),rep(rnorm(1,sd=1000),4))
  y.fels <- y.fels + c(rep(0,4),rep(rnorm(1,sd=1000),4))
  up.pic <- c(up.pic,cor.test(pic(x.fels,w.d),pic(y.fels,w.d))$p.value)
  up.cp <- c(up.cp,sum(pbreplicate(500,cor.test(x.fels,cyclic.permutation(w.d,y.fels))$estimate)>
                         cor.test(x.fels,y.fels)$estimate)/500)
  up.sbp <- c(up.sbp,sum(pbreplicate(500,cor.test(x.fels,phylo.permute(w.d,y.fels))$estimate)>
                           cor.test(x.fels,y.fels)$estimate)/500)}
hist(up.pic,breaks=seq(0,1,0.05))
hist(up.cp,breaks=seq(0,1,0.05))
hist(up.sbp,breaks=seq(0,1,0.05))
write.csv(up.pic,file="Uyeda_pic.csv")
write.csv(up.cp,file="Uyeda_cp.csv")
write.csv(up.sbp,file="Uyeda_sbp.csv")



#Felsenstein's worst case----
felsentree <- read.tree('felsentree.tre')
x.fels <- fastBM(felsentree); y.fels <- fastBM(felsentree);plot(x.fels,y.fels)
#write.csv(cbind(x.fels,y.fels),"felsxy.csv")
xyfels <- read.csv("felsxy.csv")
orperms.fels <- pbreplicate(1000,cor.test(x.fels,sample(y.fels))$estimate)
pperms.mar2.fels <- pbreplicate(1000,cor.test(phylo.permute(tre=felsentree,vec=x.fels,margin=2),
                                              phylo.permute(tre=felsentree,vec=y.fels,margin=1))$estimate);beep()
pperms.mar05.fels <- pbreplicate(1000,cor.test(phylo.permute(tre=felsentree,vec=x.fels,margin=0.5),
                                               phylo.permute(tre=felsentree,vec=y.fels,margin=0.5))$estimate);beep()
pperms.mar03.fels <- pbreplicate(1000,cor.test(xyfels$x.fels,
                                               phylo.permute(tre=felsentree,vec=xyfels$y.fels,margin=0.3))$estimate);beep()
pperms.mar01.fels <- pbreplicate(1000,cor.test(xyfels$x.fels,
                                               phylo.permute(tre=felsentree,vec=xyfels$y.fels,margin=0.1))$estimate);beep()
pperms.mar001.fels <- pbreplicate(1000,cor.test(xyfels$x.fels,
                                                phylo.permute(tre=felsentree,vec=xyfels$y.fels,margin=0.01))$estimate);beep()
pperms.mar0005.fels <- pbreplicate(1000,cor.test(xyfels$x.fels,
                                                 phylo.permute(tre=felsentree,vec=y.fels,margin=0.005))$estimate);beep()
pperms.mar00001.fels <- pbreplicate(1000,cor.test(phylo.permute(tre=felsentree,vec=x.fels,margin=0.0001,stucklimit = 2000),
                                                  phylo.permute(tre=felsentree,vec=y.fels,margin=0.0001,stucklimit = 2000))$estimate);beep()
cperms.fels <- pbreplicate(1000,cor.test(xyfels$x.fels,
                                         cyclic.permutation(tr=felsentree,x=xyfels$y.fels,N=500))$estimate);beep()
lgperms.fels <- pbreplicate(1000,cor.test(xyfels$x.fels,LG.permute(tre=felsentree,vec=xyfels$y.fels,k=1))$estimate);beep()
hist(abs(orperms.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count")
hist(abs(pperms.mar2.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count")
hist(abs(pperms.mar05.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count")
hist(abs(pperms.mar03.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count")
hist(abs(pperms.mar01.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count")
hist(abs(pperms.mar001.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count")
hist(abs(pperms.mar00001.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count")
hist(abs(cperms.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count")
hist(abs(lgperms.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count")
# Figs --------------------------------------------------------------------
#Fig 3
#tree
plot(x.fels,y.fels,type="n",xaxt='n',yaxt='n',xlab="X (I = 0.419)",ylab="Y (I = 0.522)")
phylomorphospace(felsentree,cbind(x.fels,y.fels),node.size=c(0),label = 'off',add=TRUE)
points(x.fels,y.fels,pch=21,bg="darkorange2")
#hists
par(mfrow=c(6,1),mar=c(3,4.5,0,1))
hist(abs(orperms.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count",main="",xlab="")
hist(abs(pperms.mar05.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count",main="",xlab="")
hist(abs(pperms.mar03.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count",main="",xlab="")
hist(abs(pperms.mar01.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count",main="",xlab="")
hist(abs(pperms.mar001.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count",main="",xlab="")
hist(abs(cperms.fels),breaks=seq(0,1,0.02),ylim=c(0,110),ylab="Count",main="",xlab="")

#Fig 2
par(mfrow=c(3,2),mar=c(3,4.5,1,1))
hist(ps.LGP,breaks=seq(0,1,0.05),freq = FALSE,main="",xlab="");lines(c(0,1),c(1,1),col="red",lwd=2)
hist(ps.FN.LGP,breaks=seq(0,1,0.05),freq = FALSE,main="",xlab="",ylab="")
hist(ps.CP,breaks=seq(0,1,0.05),freq = FALSE,main="",xlab="");lines(c(0,1),c(1,1),col="red",lwd=2)
hist(ps.FN.CP,breaks=seq(0,1,0.05),freq = FALSE,main="",xlab="",ylab="")
hist(ps.PP,breaks=seq(0,1,0.05),freq = FALSE,main="",xlab="");lines(c(0,1),c(1,1),col="red",lwd=2)
hist(ps.FN.PP,breaks=seq(0,1,0.05),freq = FALSE,main="",xlab="",ylab="")
hist(ps.PIC,breaks=seq(0,1,0.05),freq = FALSE,main="",xlab="");lines(c(0,1),c(1,1),col="red",lwd=2)
hist(ps.FN.PIC,breaks=seq(0,1,0.05),freq = FALSE,main="",xlab="",ylab="",ylim=c(0,15))



