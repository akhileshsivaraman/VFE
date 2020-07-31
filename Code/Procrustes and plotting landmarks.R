##### Procrustes and plotting landmarks #####
library(geomorph)

### cod data
setwd("~/Documents/Imperial/UEBS/cod")
cod <- readland.tps("digitize.tps.txt", specID = "ID")
cod1 <- cod[,,1]

# procrustes transformation
codgpa <- gpagen(cod, Proj = T)
codgpa$coords
codgpa$Csize
par(mar = c(0,0,0,0))
plotAllSpecimens(codgpa$coords, mean = T, label = T)


# PCA
codpca <- gm.prcomp(codgpa$coords)


### shark data
setwd("~/Documents/Imperial/UEBS/2.MS-Ecomorphology_lower_jaw/Data _ Analyses/")
shark <- readland.nts("Shark data with specimens.txt")
shark1 <- shark[,,1]
plot3d(x = shark1[,1],
       y = shark1[,2],
       z = shark1[,3],
       size = 5)
lines3d(x = shark1[,1],
        y = shark1[,2],
        z = shark1[,3])
# doesn't look like anything

sharkgpa <- gm.prcomp(shark)
sharkgpa$A[,,1][,1]
plot3d(x = sharkgpa$A[,,1][,1],
       y = sharkgpa$A[,,1][,2],
       z = sharkgpa$A[,,1][,3])
lines3d(x = sharkgpa$A[,,1][,1],
       y = sharkgpa$A[,,1][,2],
       z = sharkgpa$A[,,1][,3])
# same as the above

# not useful to plot all as it's a mix of mandibles and meckel's
plotAllSpecimens(shark, mean = T, label = T)
plotAllSpecimens(sharkgpa$A, mean = T, label = T)



