##### Part 1 Cookbook #####
# set working directory
# install and load packages
library(geomorph)
library(dispRity)
library(ggplot2)
library(ggrepel)
library(data.table)

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