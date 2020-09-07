##### Generating PCA data from sharks #####
setwd("~/Documents/Imperial/UEBS/2.MS-Ecomorphology_lower_jaw/Data _ Analyses")

library(ape)
library(geomorph)
library(dispRity)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tidyverse)

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
mds <- read.csv("Multidimensional scaling copy.csv", header = T)
mds <- mds %>% distinct(specimen, .keep_all = T)

diet <- read.table("N=153 PC1 scores MCL diet data averaged by spp NA diet pruned.txt", header = T)
colnames(diet)[1] <- "specimen"
diet <- diet[,-2] # written to csv file for students

all <- merge(mds, diet) # written to csv

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



#### MANOVA ####
# does diet have an effect on shape
twoDproclndmks <- two.d.array(proclndmks$coords)
lm <- procD.lm(twoDproclndmks ~ dimnames(sharklndmks)[[3]], iter = 999)
lm$aov.table
# use diet (*diet or diet as the main variable)

# morpol.disparity

plot(lm,
     type = "regression",
     reg.type = "PredLine",
     predictor = #diet?
     )


