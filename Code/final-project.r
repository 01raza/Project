setwd('C:\\Users\\Raza\\Desktop\\Evolution\\Project\\Data') 
#dir()
data<- read.csv("brain.csv",stringsAsFactors=F)
data$brain.body.ratio<-(data$brain.mass/data$body.mass)
rownames(data) <- gsub(" ", "_", data[,2])
#head(data)


library(phytools)
Trees <- read.nexus("tree.nex")
Phy <- Trees[[1]]

Pruned <- data[Phy$tip.label,]
x <- Pruned$brain.body.ratio
names(x) <- rownames(Pruned)

pdf("Phyl_tree", height = 5, width = 10)
par(mfrow=c(1,2))
contMap(Phy, x, fsize = c(.1,1))
plotTree.wBars(Phy, x)
dev.off()

# This will take forever
phenogram(Phy, x, fsize=0.01)


### Analyze BAMM results
library(BAMMtools)
setwd('C:\\Users\\Raza\\Desktop\\bamm-2.5.0-Windows') 
Tree <- read.tree("phy.tre")
tree <- Tree[[1]]
eData <- getEventData(tree, "mar31_event_data.txt", burnin=0.5)
trait <- data[eData$tip.label, "brain.body.ratio"]
names(trait) <- eData$tip.label
plot(eData, spex="s") # set spex to e for extinction, netdiv for net diversification
plot(eData, spex="e")
plot(eData, spex="netdiv")

meanSpeciation <- eData$meanTipLambda
names(meanSpeciation) <- tree$tip.label

meanExtinction <- eData$meanTipMu
names(meanExtinction) <- tree$tip.label

meanDiversification<- meanSpeciation-meanExtinction
names(meanDiversification) <- tree$tip.label

meanSpeciation

plot(meanSpeciation, trait, xlab = "Mean Speciation", ylab = "Brain/Body Mass Ratio", col="orange",pch=1)
plot(meanExtinction, trait, xlab = "Mean Extinction", ylab = "Brain/Body Mass Ratio", col="pink",pch=1)
plot(meanDiversification, trait, xlab = "Mean Diversification", ylab = "Brain/Body Mass Ratio", main= "Mean diversification vs Trait with BAMM", col="brown",pch=1)
abline(lm(meanDiversification ~ trait))
## Maybe look into STRAPP
strappBBR <- traitDependentBAMM(eData, trait, reps=100)
strappBBR1<- traitDependentBAMM(eData, trait, reps=100,rate="net diversification")
strappBBR2<- traitDependentBAMM(eData, trait, reps=100,rate="extinction")
strappBBR
strappBBR1
strappBBR2



### DR Statistic
# DR metric / inverse equal splits
# WRITTEN BY PASCAL TITLE
DRstat <- function(tree) {
  
  spRate <- function(sp, tree) {
    #get branch lengths from root to tip
    edges <- vector()
    daughterNode <- match(sp, tree$tip.label)
    while (daughterNode != (length(tree$tip.label) + 1)) {
      parentNode <- tree$edge[which(tree$edge[,2] == daughterNode), 1]
      edges <- c(edges, tree$edge.length[which(tree$edge[,1] == parentNode & tree$edge[,2] == daughterNode)])
      daughterNode <- parentNode
    }
    
    res <- sum(sapply(1:length(edges), function(x) edges[x] * (1/(2 ^ (x-1)))))
    res <- res ^ (-1)
    
    return(res)
  }
  
  rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
  names(rates) <- tree$tip.label
  
  return(rates)
}

DRrates <- DRstat(tree)
plot(DRrates, trait,xlab = "DR rates", ylab = "Brain/Body Mass Ratio", main= "DR rates vs Trait with DR", col="10",pch=1)
abline(lm(DRrates ~ trait))
# this does NOT account for autoregression! (don't trust)
cor.test(DRrates, trait, method="spear")

