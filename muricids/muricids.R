#setup ----
require(phylolm)
require(pbapply)
require(phytools)
require(scales)
muricids <- read.csv('rangesmuri_phylogeneticAnalysis.csv')
muricids$latitude_lowest <- NA
for(i in 1:length(muricids$latitude_lowest)){
  muricids$latitude_lowest[i] <- ifelse(sign(muricids$south[i])!=sign(muricids$north[i]),0,
                                   ifelse(muricids$south[i]<0,abs(muricids$north[i]),muricids$south[i]))}
muricids$latitude_highest <- pmax(abs(muricids$south),abs(muricids$north))
s1 <- read.csv("Table_S1.csv")
s1$Species.name <- sapply(strsplit(s1$Species.name," "), "[[", 2)

#plotting/analysing functions and exercises ----
#functions
skyline <- function(xvals,yvals,color="black",switchXY=FALSE){ #xvals one longer than yvals
  for(i in 1:length(yvals)){   #horizontal lines
    if(switchXY){
      lines(c(yvals[i],yvals[i]),c(xvals[i],xvals[i+1]),col=color,lwd=2)
    }else{
      lines(c(xvals[i],xvals[i+1]),c(yvals[i],yvals[i]),col=color,lwd=2)}}
  for(i in 1:(length(yvals)-1)){ #vertical lines
    if(switchXY){
      lines(c(yvals[i],yvals[i+1]),c(xvals[i+1],xvals[i+1]),col=color,lwd=2)
    }else{
      lines(c(xvals[i+1],xvals[i+1]),c(yvals[i],yvals[i+1]),col=color,lwd=2)}}}

plotbylatmeans <- function(lowlats,highlats,trait,rangelimits=NA,xlims=c(-90,90),xlab="",ylab="",main="",switchXY=FALSE){
  if(is.na(rangelimits)){rangelimits <- sort(unique(c(lowlats,highlats)))}
  vals <- c()
  for(i in (1:(length(rangelimits)-1))){
    vals <- c(vals,mean(trait[highlats > rangelimits[i] & lowlats < rangelimits[i+1]],na.rm=TRUE))}
  #plot
  if(switchXY){
    plot(1,ylim=c(xlims[1],xlims[2]),xlim=c(min(vals,na.rm = TRUE),max(vals,na.rm = TRUE)),type='n',ylab=xlab,xlab=ylab,main=main)
  }else{
    plot(1,xlim=c(xlims[1],xlims[2]),ylim=c(min(vals,na.rm = TRUE),max(vals,na.rm = TRUE)),type='n',xlab=xlab,ylab=ylab,main=main)}
  skyline(rangelimits,vals,switchXY=switchXY)}

returnlatmeans <- function(lowlats,highlats,trait,rangelimits=NA){ #same thing as plotbylatmeans, but return instead of plot
  if(is.na(rangelimits)){rangelimits <- sort(unique(c(lowlats,highlats)))}
  vals <- c()
  for(i in (1:(length(rangelimits)-1))){
    vals <- c(vals,mean(trait[highlats > rangelimits[i] & lowlats < rangelimits[i+1]],na.rm=TRUE))}
  return(vals)}

#reading in phylogeny ----
mt <- read.tree("muricidae.dated.tre")
mt$tip.label <- sapply(strsplit(mt$tip.label,"_"), "[[", 2)
mt <- drop.tip(phy=mt,tip=mt$tip.label[!mt$tip.label %in% muricids$species])
muricids <- muricids[match(mt$tip.label,muricids$species),]

#logistic PGLS ----
mt.full <- read.tree("muricidae.dated.tre");mt.full$tip.label <- sapply(strsplit(mt.full$tip.label,"_"), "[[", 2)
s1 <- s1[match(mt.full$tip.label,s1$Species.name),]
temp <- s1$SST;names(temp) <- mt.full$tip.label
abslat <- abs(s1$Latitude...S.);names(abslat) <- mt.full$tip.label
pelagic <- as.numeric(s1$Pelagic=="yes");names(pelagic) <- mt.full$tip.label
feeding <- as.numeric(s1$Feeding=="yes");names(feeding) <- mt.full$tip.label
mtlmf <- phyloglm(feeding~abslat,phy = mt.full,method="logistic_IG10")
mtlmp <- phyloglm(pelagic~abslat,phy = mt.full,method="logistic_IG10")

#phylogenetic permutations----
#ordinary permutations
ops.1 <- pbreplicate(10000,cor.test(seq(61),returnlatmeans(muricids$latitude_lowest,muricids$latitude_highest,sample(muricids$MLD=="PV"),rangelimits=seq(0,61,1)),method="spearman")$estimate)
#phylogenetic permutations
pps.1 <- pbreplicate(1000,cor.test(seq(61),returnlatmeans(muricids$latitude_lowest,muricids$latitude_highest,phylo.permute(mt,as.numeric(muricids$MLD=="PV"),margin=0.01),rangelimits=seq(0,61,1)),method="spearman")$estimate)
#empirical correlation
cor.test(seq(61),returnlatmeans(muricids$latitude_lowest,muricids$latitude_highest,phylo.permute(mt,as.numeric(muricids$MLD=="PV"),margin=0.001),rangelimits=seq(0,61,1)),method="spearman")

#plotting ----
#A
dev.off();par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(c(1,length(muricids$species)),c(0,60),type="n",bty="n",xaxt="n",yaxt="n",ylab="Absolute latitude (\u00B0)",xlab="")
axis(side=2,at=seq(0,60,10),labels=seq(0,60,10))
wp <- wesanderson::wes_palette("Zissou1")
for(i in 1:length(muricids$species)){
  lines(c(i,i),c(muricids$latitude_lowest[i],muricids$latitude_highest[i]),col=ifelse(muricids$MLD[i]=="PV",wp[1],
                                                                                      ifelse(muricids$MLD[i]=="DD",wp[3],wp[5])),lwd=3)}

#B
plotbylatmeans(muricids$latitude_lowest,muricids$latitude_highest,muricids$MLD=="PV",rangelimits=seq(0,65,1),xlims = c(0,60),ylab="Proportion with planktonic larvae",switchXY = TRUE)

#C
plot.phylo(mt,direction="upwards")

#D
hist(abs(ops.1),breaks=seq(0,1,0.05),border = FALSE,freq=FALSE,ylim=c(0,8),
     main="Percent species with planktonic larvae against 1\u00B0 bins of absolute latitude",xlab="Absolute value of r")
hist(abs(pps.1),breaks=seq(0,1,0.05),col = alpha("black",0),freq=FALSE,ylim=c(0,8),add=TRUE)
lines(c(0.9270465,0.9270465),c(0,200),col="darkorange2",lwd=4)
