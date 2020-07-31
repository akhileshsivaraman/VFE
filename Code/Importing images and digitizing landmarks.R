##### Importing images and digitizing landmarks #####
setwd("~/Documents/Imperial/UEBS/cod")
library(geomorph)

# digitize one specimen
digitize2d("IMG_1016.JPG",
           nlandmarks = 14,
           tpsfile = "1016.tps",
           scale = 3,
           verbose = F)

# load saved tps file
a <- readland.tps("1016.tps.txt", specID = "ID")
dim(a) # 2D
a1 <- as.data.frame(a)



# workflow to digitize a bunch in one go
list <- list.files(pattern = "IMG")
digitize2d(list,
           nlandmarks = 14,
           tpsfile = "digitize.tps",
           scale = 3,
           verbose = F)

# load saved tps file
b <- readland.tps("digitize.tps.txt", specID = "ID")
dim(b) # 3D

