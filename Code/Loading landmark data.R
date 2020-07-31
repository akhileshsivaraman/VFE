##### Loading landmark data #####

library(geomorph)

setwd("~/Documents/Imperial/UEBS/geomorph")
a <- readland.nts("my.ply.nts.txt")
b <- as.data.frame.array(a) # not a 3D array


setwd("~/Documents/Imperial/UEBS/2.MS-Ecomorphology_lower_jaw/Data _ Analyses/")
c <- readland.nts("Shark data.txt") # removed specimen block
c[1,,] # first row of each triplet
c[,1,] # first six entries in the first column of a triplet
c[,,1] # all the data of the first triplet

d <- as.data.frame.array(c) # not a 3D array
d[1,,] # first row
d[,1,] # first column
d[,,1] # first two rows


e <- readland.nts("Shark data with specimens.txt") # L after total variables
e[1,,]
e[,,1]

f <- readland.nts("Shark data with specimens.txt") # L after total specimens
f[1,,]
f[,,1]

g <- readland.nts("Shark data with specimens.txt") # L after both
g[1,,]
g[,,1]

h <- readland.nts("Shark data with specimens.txt") # specimen block edited to remove spaces
h[1,,]
h[,,154] 
dim(h) # 3D array
h1 <- as.data.frame.array(h)
h2d <- two.d.array(h) # transform 3D array to 2D array

## adding data will require students to add the block of coordinates at the bottom of the matrix (maintain a gap)
## add the specimen name to the bottom of the specimen block and change the number of specimens
## will need to explain the structure of the data file then 