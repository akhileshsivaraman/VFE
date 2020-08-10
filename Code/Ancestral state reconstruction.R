library(ape)
library(phytools)
library(diversitree)

##### continuous characters #####
anoletree <- read.tree("http://www.phytools.org/eqg2015/data/anole.tre")
plotTree(anoletree)
svl <- read.csv("http://www.phytools.org/eqg2015/data/svl.csv",
                row.names=1) # data frame with a continuous vairable assigned to all species
svl <- as.matrix(svl)[,1] # convert column with continuous data to a matrix of numbers with names
fit <- fastAnc(anoletree, svl, vars = T, CI = T) # estimate ancestral states
obj <- contMap(anoletree, svl) # maps continuous character to a tree



##### discrete characters #####
data("anoletree")
x <- getStates(anoletree, "tips")
tree <- anoletree
rm(anoletree)

par(mar=c(0,0,0,0))
par(mfcol = c(1,1))
plot.phylo(tree, label.offset = 0.4, cex = 0.8)
cols <- setNames(palette()[1:length(unique(x))], sort(unique(x)))
tiplabels(pie = to.matrix(x, sort(unique(x))), piecol = cols, cex=0.3)
add.simmap.legend(colors = cols, prompt = F, 
                  x = 0.9*par()$usr[1],
                  y = -max(nodeHeights(tree)), fsize = 0.8)
add.simmap.legend(colors = cols)


# ER model (evolutionary rate?)
# marginal ancestral states
fitER <- ace(x, tree, model = "ER", type = "discrete")
fitER
round(fitER$lik.anc, 3)
plot.phylo(tree, label.offset = 0.4, cex = 0.8, type = "cladogram", use.edge.length = F)
nodelabels(node = 1:tree$Nnode+Ntip(tree),
           pie = fitER$lik.anc, piecol = cols, cex = 0.5)
tiplabels(pie = to.matrix(x, sort(unique(x))), piecol = cols, cex = 0.3)

# MCMC model
# posterior probabilites
mtree <- make.simmap(tree, x, model = "ER")
plot.phylo(mtree, edge.color = cols, type = "phylogram")

mtrees <- make.simmap(tree, x, model = "ER", nsim = 100)
par(mfrow = c(10,10))
null <- sapply(mtrees, plot, col = cols)

par(mfrow = c(1,1))
pd <- summary(mtrees, plot = F)
pd
plot(pd, fsize = 0.6, ftype = "i")

plot(mtrees[[1]], col = cols)
nodelabels(pie = pd$ace, piecol = cols, cex = 0.5)
