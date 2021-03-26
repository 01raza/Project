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

par(mfrow=c(1,2))
contMap(Phy, x, fsize = c(.1,1))
plotTree.wBars(Phy, x)

# This will take forever
phenogram(Phy, x, fsize=0.01)
