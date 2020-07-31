##### Analysis of evolutionary rates #####
setwd("~/Documents/Imperial/UEBS/2.MS-Ecomorphology_lower_jaw/Data _ Analyses")

library(ape)
library(geomorph)

# load tree
tree <- read.nexus("Adam_ML_SharkTree.nex")

# load landmark data
landmarks <- readland.nts("Shark data with specimens.txt")

