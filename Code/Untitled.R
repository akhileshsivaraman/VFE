setwd("~/Documents/Imperial/UEBS")

library(ape)
library(geiger)
library(dispRity)
library(vegan)
library(castor)
library(paleotree)

tree <- read.nexus("2.MS-Ecomorphology_lower_jaw/Data _ Analyses/Adam_ML_SharkTree.nex.txt")
# might need to time scale the tree as in Anderson et al (2013) - need to provide FADs for terminal taxa and a constraint for the age of the tree

# can try dating the tree based on evolutionary divergence
datedtree <- date_tree_red(tree, anchor_age = 1)

ages <- cbind(c("Hexanchus_griseus"), c(0))
datedtree1 <- setRootAge(tree, fixedAges = ages)
ages2 <- cbind(c("Glyphis_glyphis"), c(0))
datedtree2 <- setRootAge(tree, fixedAges = ages2) # same as dated tree 1

traits <- read.table("2.MS-Ecomorphology_lower_jaw/Data _ Analyses/N=153 spp PC1 dietcode MCL averaged by spp NA diet pruned.txt", header = T)
