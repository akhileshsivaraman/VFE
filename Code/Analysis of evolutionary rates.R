##### Analysis of evolutionary rates #####
setwd("~/Documents/Imperial/UEBS/2.MS-Ecomorphology_lower_jaw/Data _ Analyses")

library(ape)
library(geiger)
library(phytools)

# load tree
tree <- read.nexus("Adam_ML_SharkTree.nex.txt")

# load pc data
data <- read.table("N=153 spp PC1 dietcode MCL averaged by spp NA diet pruned.txt", header=T,row.names=1)

# clean data
newtree <- treedata(tree, data,sort = T)
pc1 <- newtree$data[,2]
pc1 <- as.numeric(pc1)
names(pc1) <- rownames(newtree$data)


##### geiger rjmcmc #####
# bm
rjmcmc.bm(phy = newtree$phy, dat = pc1, type = "bm", ngen = 10e6, samp = 10000)
res <- load.rjmcmc("BM.result/")
plot(x = res, par = "shifts", legend = T, cex = 0.5, type = "fan") # colour error again
plot(x = res, par = "jumps", legend = T) # no jumps


# relaxed bm
rjmcmc.bm(phy = newtree$phy, dat = pc1, type = "rbm", ngen = 10000, samp = 10)
res1 <- load.rjmcmc("relaxedBM.result/")
plot(x = res1, par = "shifts", legend = T) # colour error
plot(x = res1, par = "jumps", legend = T) # no jumps

shifts <- res1$shifts
jumps <- res1$jumps
rates <- res1$rates
log <- res1$log

# the problem may lie in our data
# there's no shift to actually plot?

# jump-rbm
rjmcmc.bm(phy = newtree$phy, dat = pc1, type = "jump-rbm", ngen = 10000, samp = 10)
res2 <- load.rjmcmc("jump-relaxedBM.result/")
plot(x = res2, par = "shifts", legend = T) # colour error again
plot(x = res2, par = "jumps", legend = T) # some small jumps

# jump-bm
rjmcmc.bm(phy = newtree$phy, dat = pc1, type = "jump-bm", ngen = 10000, samp = 10)
res3 <- load.rjmcmc("jump-relaxedBM.result/")
plot(x = res3, par = "shifts", legend = T) # colour error again
plot(x = res3, par = "jumps", legend = T) # some small jumps

ratesres3 <- res3$rates

##### contmap #####
# estimated ancestral states
contMap(newtree$phy, pc1, fsize = 0.6, res = 100)


##### fitcontinuous #####
a <- fitContinuous(newtree$phy, pc1, model = "kappa")
b <- fitContinuous(newtree$phy, pc1, model = "lambda")
c <- fitContinuous(newtree$phy, pc1, model = "BM")
d <- fitContinuous(newtree$phy, pc1, model = "OU")
## cannot be plotted ffs
##### auteur - works for some reason (have a look at the data) #####
library(auteur)
rjmcmc.bm(newtree$phy, pc1, ngen = 10000, sample.freq = 100,
          fileBase = "auteur")
plot.new()
pdf("shark rates.pdf", 10, 20)
shifts.plot(newtree$phy, base.dir = "BM.auteur.parameters/", 
            burnin = 0.2, legend = F, edge.width = 2)
dev.off()
shifts.plot(newtree$phy, base.dir = "BM.auteur.parameters/", 
            burnin = 0.2, legend = T, edge.width = 2) # plots weirdly



##### using time-calibrated trees #####
setwd("~/Documents/Imperial/UEBS")
sharks <- read.nexus("shark.nex")
shark1 <- sharks[[1]]

# load pc data
newtree1 <- treedata(shark1, data)
pc1 <- newtree1$data[,2]
pc1 <- as.numeric(pc1)
names(pc1) <- rownames(newtree1$data)

geiger::rjmcmc.bm(phy = newtree1$phy, dat = pc1, 
                  ngen = 10000, samp = 100, type = "bm")
tcres <- load.rjmcmc("BM.result/")
plot(x = tcres, par = "shifts", legend = T) # colour error again

tclog <- tcres$log
