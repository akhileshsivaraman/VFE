##### Part 1 Cookbook #####
# set working directory
# install and load packages
library(geomorph)
library(dispRity)
library(ggplot2)
library(ggrepel)
library(data.table)
library(ape)
library(tidyverse)
library(geiger)
library(phytools)
library(evobiR)
library(picante)

#### Loading landmark data ####
# nts or tps format
# description of the file format
landmarks <- readland.nts() # or readland.tps() which would be preferable for the addition step
dim(landmarks) # what are the dimensions of the file


#### Importing images and digistising landmarks ####
# digitising an image
digitize2d(# file containing the image,
           nlandmarks =, # number of landmarks in your scheme
           tpsfile = # file name to add your landmark data
           scale = 3,
           verbose = F)

# what to do when you have a lot of images you'd like to digitise in one go?
# firstly, make a list of the images
imagelist <- list.files(pattern = "IMG") # list.files() creates a list that includes all files in the folder with "IMG" in its name
digitize2d(imagelist,
           nlandmarks = , # number of landmarks in your scheme,
           tpsfile = , # file name to add your landmark data
           scale = 3,
           verbose = F)


#### Adding digitised landmark data to a file ####
# copy and paste the landmark data you generated into the file provided
# details depend on the file type (tbd)

# reload landmark data
landmarks <- readland.tps(file = ,
                          specID = "ID")

#### Procrustes transformation ####
# explain Procrustes transformation
proclandmarks <- gpagen(landmarks)
# plot all landmarks
plotAllSpecimens(proclandmarks$coords, # coordinates for each specimen are stored in the object coords
                 mean = F)

#### Thin plate splines ####
meanshape <- mshape(proclandmarks$coords) # find mean specimen shape
# we'll use the mean shape as a reference from which we determine the kinds of transformations involved in moving from the mean shape to a specimen's shape
# create list of specimens
plotRefToTarget(M1 = meanshape,
                M2 = specimenlist)

#### Principal warps and principal coordinate analysis ####
# plot specimens in a multivariate space based on the kinds of transformations involved in moving from the mean shape to a specimen's shape
pw <- geomorph.ordination(y)
pw <- as.data.frame(pw)
pw1 <- setDT(pw, keep.rownames = T) # move row names into a column
ggplot(pw1, aes(PC1, PC2)) + # plot specimens in a 2D space using the principal warps that account for most of the variation in shape
  geom_point(position = "jitter") + 
  geom_label_repel(aes(label = rn)) + # add labels to the plot
  theme_bw()




##### Part 2 Cookbook #####
### shark data ###
#### load landmark data ####
sharklndmks <- readland.nts("Shark data with specimens.txt")
# explore the data set
str(sharklndmks)
dim(sharklndmks)
dimnames(sharklndmks)
# what are the dimensions of the data set?
# how many specimens are there?
# how many landmarks are there per specimen?

#### Procrustes analysis ####
proclndmks <- gpagen(sharklndmks)

# plot all the specimens
plotAllSpecimens(proclndmks$coords, mean = T, label = T)

# visualise some of the specimens
dimnames(sharklndmks)[[3]][1]
Carcharhinus_dussumieri_Mandibula_L <- proclndmks$coords[,,1]
plot3d(x = Carcharhinus_dussumieri_Mandibula_L[,1],
       y = Carcharhinus_dussumieri_Mandibula_L[,2],
       z = Carcharhinus_dussumieri_Mandibula_L[,3],
       type = "p", size = 10)
dimnames(sharklndmks)[[3]][23]
Isurus_paucus_Mandibula_L <- proclndmks$coords[,,23]
plot3d(x = Isurus_paucus_Mandibula_L[,1],
       y = Isurus_paucus_Mandibula_L[,2],
       z = Isurus_paucus_Mandibula_L[,3],
       type = "p", size = 10)
# try plotting some other shark specimens in 3D

# thin plate splines
par(mar=c(0,0,0,0))
meanshape <- mshape(proclndmks$coords)
plotRefToTarget(M1 = meanshape,
                M2 = proclndmks$coords[,,1])
# try exploring thin plate splines of other specimens against the mean shape

#### Multidimensional scaling ####
# explain what's happening
mds <- geomorph.ordination(proclndmks)
# save mds to a csv file - we'll need this later
write.csv(mds, file = "Multidimensional scaling.csv")

# there are a lot of axes...
# how much variation does each axis account for?
# explain calculation
screedata <- apply(mds, 2, var) / sum(apply(mds, 2, var))*100
# scree plot
par(mar=c(5,5,2,2))
plot(screedata, type = "h",
     xlab = "PCA axis",
     ylab = "Percentage of variance accounted for")

# the first two axes of a PCA explain the most and second-most amount of variance
mds <- as.data.frame(mds)
mds <- setDT(mds, keep.rownames = T)
ggplot(mds, aes(PC1, PC2)) +
  geom_point(position = "jitter") +
  xlab(paste("PC1", "-", round(screedata[1], 2), "% variance")) +
  ylab(paste("PC2", "-", round(screedata[2],2), "% variance")) +
  theme_bw()

# there's one obvious outlier - let's see which specimen it is
ggplot(mds, aes(PC1, PC2)) +
  geom_text(aes(label = rn), size = 2) +
  theme_bw()


#### what does this mean for function? #####
# add diet data and use to colour points
mds <- read.csv("Multidimensional scaling copy.csv", header = T) # cleaned up the names of specimens, could ask students to do that?
mds <- mds %>% distinct(specimen, .keep_all = T)

diet <- read.table("N=153 PC1 scores MCL diet data averaged by spp NA diet pruned.txt", header = T)
colnames(diet)[1] <- "specimen"
diet <- diet[,-2] # if easier, have written to csv file so students don't need to clean their data but cleaning their data and understanding the structure of data could be a teaching point and is a transferrable skill

all <- merge(mds, diet) # have written to csv file too

ggplot(all, aes(PC1, PC2)) +
  geom_point(aes(colour = Prey), position = "jitter") +
  xlab(paste("PC1", "-", round(screedata[1], 2), "% variance")) +
  ylab(paste("PC2", "-", round(screedata[2], 2), "% variance")) +
  scale_colour_discrete(labels = c("Decapoda &
mollusca",
                                   "Fish &
cephalopoda",
                                   "Generalist",
                                   "Soft bodied
benthic invertebrates",
                                   "Zooplankton")) +
  theme(legend.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  theme_bw()


#### apply MDS to phylogeny/modelling evolution ####
rownames(mds) <- mds[,1]
mds <- mds[,-1]

# read in trees
trees <- read.nexus("../../shark.nex") # do we just give the students one of the trees?

# tip labels in the trees need to match the species present in the mds data
newtree1 <- treedata(trees[[1]], mds, sort = T)
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
# the amount of change that occurs across a lineage changes with each PC
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
# that's not particularly interesting...
# let's try the second PC axis
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

