##### Lab 2 in R #####
setwd("~/Documents/Imperial/UEBS/Lab 2 example")

#### Quantifying morphological diversity using discrete characters ####
library(ape)
a <- read.nexus.data("Mammal_data.nex")
b <- data.frame(t(sapply(a, c)))
c <- dist(b, method = "euclidean", diag = T, upper = T) # distance matrix
d <- hclust(c) # phenogram/hierarchical clustering | problem here
plot(d)


### NMDS
library(vegan)
library(ggplot2)
library(ggrepel)
e <- as.matrix(c) # generate matrix
f <- metaMDS(e, distance = "euclidean", trymax = 999) # NMDS
plot(f$points[,1],
     f$points[,2]) # simple plot
nmds <- as.data.frame(f$points)
nmds <- data.table::setDT(nmds, keep.rownames = T)
ggplot(nmds, aes(MDS1, MDS2)) + # plot with labels
  geom_point(position = "jitter") +
  geom_label_repel(aes(label = rn)) +
  theme_bw()


### character mapping
library(geiger)
library(phytools)
tree <- read.nexus("contree.tre")
newtree <- treedata(tree, b)
sim <- make.simmap(newtree$phy, newtree$data[,1], model = "ER", 1)
# doesn't work - because characters in the matrix are binary so we're not measuring variation of a phenotype


#### Examining shape variation using morphometry ####
library(geomorph)
z <- readland.tps("Cambrian_trilobites.TPS.txt", specID = "ID")
y <- gpagen(z) # Procrustes
plotAllSpecimens(y$coords, mean = F) # plot landmarks of all specimens
meanshape <- mshape(y$coords) # generate mean shape of all specimens
plotRefToTarget(M1 = meanshape, # thin plate splines
                M2 = y$coords[,,1]) # how to do this all in one go?
specimens <- y$coords 
# put each sheet into a list and apply plotRefToTarget to all elements in the list
# use a for loop?

par(mfrow=c(3,1))
par(mar=c(0,0,0,0))
plotRefToTarget(M1 = meanshape, 
                M2 = y$coords[,,1])
plotRefToTarget(M1 = meanshape, 
                M2 = y$coords[,,2])
plotRefToTarget(M1 = meanshape, 
                M2 = y$coords[,,3])


### dispRity - for ordination of procrustes points
library(dispRity)
pw <- geomorph.ordination(y) # principal warps
pw <- as.data.frame(pw)
pw <- data.table::setDT(pw, keep.rownames = T)
pw3 <- pw[,1:4] # take the number of axes that account for 95% of the variation?
ggplot(pw3, aes(PC1, PC2)) +
  geom_point(position = "jitter") +
  geom_label_repel(aes(label = rn)) +
  theme_bw()
# after plotting the students should look back at the specimens and hypothesise about how morphology changes along the axes
# evaluate differences in the morphospace with npMANOVA and/or ANOSIM


### PCA axes used with fitContinuous


# generate example tree via hierarchical clustering
pwnamed <- geomorph.ordination(y)
pcdist <- dist(pwnamed, method = "euclidean", diag = T, upper = T)
pcclustering <- hclust(pcdist)
plot(pcclustering)
pctree <- as.phylo(pcclustering)
par(mar = c(0,0,0,0))
plot(pctree)


# simulate evolutionary rates
traits <- c(1,2,3,1,2,3,2,3,4,1,1,
            2,3,2,2,1,4,1,2,1,4,2)
pwtraits <- as.data.frame(geomorph.ordination(y))
pwtraits <- cbind(pwtraits, traits)


charactertree <- treedata(pctree, pwtraits)
# plotting PCA to tree
simmap <- make.simmap(charactertree$phy, charactertree$data[,2], model = "ER", nsim = 1)
plotSimmap(simmap, fsize = 1) # PCA is continuous, simmap works with discrete characters

dis <- disparity(phy = charactertree$phy,
                 data = charactertree$data[,2]) # plot to a tree?

dis1 <- dispRity.through.time(data = pwnamed, 
                              tree = pctree,
                              time = 3)

# plotting characters to tree
simmap2 <- make.simmap(charactertree$phy, charactertree$data[,23], model = "ER", nsim = 1)
plotSimmap(simmap2)



##### compare.evol.rates with geomorph - separate by group (eg diet) #####
# measuring rates of trait evolutions as Borstein et al
traits1 <- factor(traits)
names(traits1) <- pctree$tip.label
g <- compare.evol.rates(A = specimens,
                        phy = pctree, 
                        gp = traits1, 
                        iter = 999)
summary(g) # evolutionary rate by group calculated
par(mar = c(5,4,2,2))
plot(g)


##### MCMC #####
# geiger rjmcmc
# continuous data - PC axes - stored in pwtraits
# only takes one trait at a time
pc1 <- pwtraits[,1]
names(pc1) <- pctree$tip.label
geiger::rjmcmc.bm(phy = pctree, dat = pc1, ngen = 10000, sample.freq = 10, type = "bm")
res <- load.rjmcmc("BM.result/")
plot(x = res, par = "shifts", burnin = 0.25, legend = T, show.tip = F, edge.width = 2) # plotting colour and shapes not working properly
plot(x = res, par = "jumps", burnin = 0.25, legend = T, show.tip = F, edge.width = 2) # no jumps to visualise here

plot(x = res, par = "shifts", burnin = 0.25, legend = T, show.tip = T)


## try with ggtree
library(ggplot2)
library(ggtree)
ggtree(res$phy, aes(colour = res$rates)) +
  scale_color_continuous(low = "blue", high = "red")
# error but on the right track



# evol.rate.mcmc in phytools
phyt <- evol.rate.mcmc(tree = pctree, x = pc1, ngen = 10000)
phytres <- summary(phyt)
plot(phytres, type = "min.split")
plot(phytres, type = "edge.prob") # nothing plotted



##### extra plots #####
## if we use something like diet to separate the species we can use ordiellipse and ordispider after NMDS
par(mar=c(2,2,2,2))
plot(x = pw$PC1,
     y = pw$PC2,
     col = c("limegreen", "firebrick", "steelblue", "gold")[pwtraits$traits],
     pch = 16)
ordiellipse(pwtraits,
            groups = traits,
            kind = "se", # not the correct kind in this case?
            col = c("limegreen", "firebrick", "steelblue", "gold"),
            draw = "polygon",
            border = 0)
ordispider(pwtraits,
           groups = traits,
           col = "grey") # could be messy if there's a lot of variation

## can use ordisurf with a continuous variable then apply contours?
plot(x = pw$PC1,
     y = pw$PC2,
     col = c("limegreen", "firebrick", "steelblue", "gold")[pwtraits$traits],
     pch = 16)

##### phylogenetic signal #####
library(picante)
sig <- phylosignal(pwtraits$traits, pctree, reps = 999) # tip labels in the same order
# K > 1 => covariance among species is stronger than expected under Browian motion evolution
# K > 1 => species tend to be independent with respect to their phylogenetic relationships

signal <- physignal(A = specimens, phy = pctree,iter = 999)
signal
# phylogenetic signal relative to what is expected under Brownian motion
# p-value given too - is it significantly different to the epected Brownian motion model


### from morphometrics R textbook
# is morphological disparity related to time divergence?
# carry out a mantel test
# cophenetic.phylo() on a tree - computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths
# compare to euclidean distance matrix
# make sure the order of the matrix is the same as the order of the tip labels - use match()
# to determine if there is a relationship between time of divergence and divergence in morphology we need to infer ancestral character states
# first need to determine variance-covariance among terminal taxa using vcv.phylo()
# vcv among terminal taxa corresponds to the shared history among taxa
