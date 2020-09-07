##### Applying shark MDS to phylogeny #####
setwd("~/Documents/Imperial/UEBS/2.MS-Ecomorphology_lower_jaw/Data _ Analyses")

library(ape)
library(geomorph)
library(dispRity)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(geiger)
library(phytools)
library(evobiR)
library(picante)

mds <- read.csv("Multidimensional scaling copy.csv", header = T)
rownames(mds) <- mds[,1] # this doesn't work because there are more than one lower jaw specimen for some of the species
# let's remove the duplicates
mds <- mds %>% distinct(specimen, .keep_all = T)
rownames(mds) <- mds[,1]
mds <- mds[,-1]

# read in trees
trees <- read.nexus("../../shark.nex")

# tip labels in the trees need to match the species present in the mds data
newtree1 <- treedata(trees[[1]], mds, sort = T)

##################
#### use apply ###
# or is there a way to take the average edge length
##################

#### Modelling evolution ####
# each PC axis is a continuous character for each species that we can map onto
# a phylogeny after estimating ancestral states
contMap(newtree1$phy, 
        newtree1$data[,1], # PC1
        fsize = c(0.6, 0.8), res = 100, type = "fan", legend = 250)
contMap(newtree1$phy, 
        newtree1$data[,2], # PC2
        fsize = c(0.6, 0.8), res = 100, type = "fan", legend = 250)
contMap(newtree1$phy, 
        newtree1$data[,3], # PC3
        fsize = c(0.6, 0.8), res = 100, type = "fan", legend = 250)

# more change along some lineages...
# rates of evolution
#### Analysis of evolutionary rates ####
pc1 <- newtree1$data[,1]
pc2 <- newtree1$data[,2]
pc3 <- newtree1$data[,3]

# set parameters
prop = 1.5        
iterations = 10e6 
sampling = 10000  
model = "rbm"

dhalfCauchy <- function(x, scale=1, log=F){
  if(any(x<0)) return(-Inf)
  density <- 2 /(pi * (1+(x/scale)^2))
  if(log==TRUE) density<-log(density/scale)
  return(density/scale)
}

# Define priors for the hyperparameter and measurment error (SE)
ratePrior <- function(x) dhalfCauchy(x, 25, log=TRUE)
sePrior <- function(x) dhalfCauchy(x, 25, log=TRUE)

# run the relaxed Brownian motion model
# tree1 and pc1
rjmcmc.bm(phy = newtree1$phy,
          dat = pc1,
          ngen = iterations,
          samp = sampling, 
          type = model,
          prop.width = prop,
          simple.start = T,
          dlnRATE = ratePrior,
          dlnSE = sePrior,
          filebase = "tree1-pc1-rjmcmc-result")
# load results
res <- load.rjmcmc("relaxedBM.tree1-pc1-rjmcmc-result/")
plot(res, par = "shifts", type = "fan", legend = T, label.offset = 4)

# tree1 and pc2
rjmcmc.bm(phy = newtree1$phy,
          dat = pc2,
          ngen = iterations,
          samp = sampling, 
          type = model,
          prop.width = prop,
          simple.start = T,
          dlnRATE = ratePrior,
          dlnSE = sePrior,
          filebase = "tree1-pc2-rjmcmc-result")
# load results
res <- load.rjmcmc("relaxedBM.tree1-pc2-rjmcmc-result/")
plot(res, par = "shifts", type = "fan", legend = T, label.offset = 4)



#### Simmap ####
data <- read.csv("Principal warps with diet.csv")
data <- data[,-1]
rownames(data) <- data[,1]
data <- data[,-1]

diettree <- treedata(trees[[1]], data, sort = T)
simmap <- make.simmap(diettree$phy, diettree$data[,43], nsim = 1)
# simmap can take a multiphylo object so try an apply function on treedata
plotSimmap(simmap, type = "fan", fsize = 0.5, offset = 1.5)
# add legend


#### Phylogenetic signal ####
sig <- phylosignal(x = newtree1$data[,1], phy = newtree1$phy, reps = 999)
