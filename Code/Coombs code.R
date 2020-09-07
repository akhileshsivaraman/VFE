##### using code from Ellen Coombs #####
# https://github.com/EllenJCoombs/Asymmetry-evolution-cetaceans/blob/master/Code%20for%20analyses/Clavel-models-shifts-jumps.R
library(ape)
library(geiger)
library(phytools)
library(evobiR)

setwd("~/Documents/Imperial/UEBS/2.MS-Ecomorphology_lower_jaw/Data _ Analyses")

tree <- read.nexus("Adam_ML_SharkTree.nex.txt")
data <- read.table("N=153 MORPHOJ PC1-3 scores MCL averaged by spp.txt", header = T, row.names = 1)
newtree <- treedata(tree, data)
newtree1 <- treedata(tree, data, sort = T)
tree <- newtree$phy
data <- as.data.frame(newtree$data)

reordered <- ReorderData(tree, data, taxa.names = "row names")
data <- reordered

# this creates a named number vector
subtree <- multi2di(tree, random = TRUE)
subtree_unroot <- unroot(subtree)
subtree_unroot1 <- multi2di(subtree_unroot)
data <- data[subtree_unroot1$tip.label, ] #this takes the label from the dataset
data <- as.vector(data$PC1) #pulls out a vector (not matrix) using the $ sign - just the radii here 
names(data) = subtree_unroot1$tip.label

#### rbm ####
#Set parameters 
prop = 1.5        # tuning parameter for the proposal window (it's just a starting value it will be improved in the model fit I think)
iterations = 10e6  # number of iterations for the chain (5 million)
sampling = 10000   # sampling frequency of the parameters (i.e. when you store the parameter along the chain)
model = "rbm"     # relaxed brownian motion (the model used in Eastman et al. 2011. More or less similar to BayesTrait)
filename = paste("testmcmc-rjmcmcREARRANGED.log", sep="", collapse="")

#Half-Cauchy distribution (see also Gelman 2006) - I used it for the prior density of the rate scalar instead of the default exponential distribution
dhalfCauchy <- function(x, scale=1, log=F){
  if(any(x<0)) return(-Inf)
  density <- 2 /(pi * (1+(x/scale)^2))
  if(log==TRUE) density<-log(density/scale)
  return(density/scale)
}

#Define priors for the hyperparameter and measurment error (SE)
ratePrior <- function(x) dhalfCauchy(x, 25, log=TRUE)
sePrior <- function(x) dhalfCauchy(x, 25, log=TRUE)



#run the rjmcmc
rjmcmc.bm(tree, data, prop.width=prop, ngen=iterations, samp=sampling, filebase=filename,
          simple.start=TRUE, type=model, dlnRATE=ratePrior, dlnSE=sePrior) # note: here I assume you're estimating a nuisance parameter as well

#retrieve the chain
chain_rearranged <- load.rjmcmc("relaxedBM.testmcmc-rjmcmcREARRANGED.log")

result=plot(chain_rearranged, par="shift", type= "fan", legend = T, cex = 0.6)
resulta=plot(chain_rearranged, par="jumps", type= "fan", legend = F, cex = 0.6)

# plotshift
library(fields)
plotShifts <-
  function(phylo, chain, burnin=1000, ...){
    #require(fields)
    args <- list(...)
    # options
    if(is.null(args[["fun"]])) args$fun <- mean
    if(is.null(args[["show.tip.label"]])) args$show.tip.label <- TRUE
    if(is.null(args[["horizontal"]])) args$horizontal <- TRUE
    if(is.null(args[["color"]])) args$color <- c("blue", "red")
    if(is.null(args[["scale"]])) args$scale <- FALSE
    if(is.null(args[["log"]])) args$log <- FALSE
    if(is.null(args[["palette"]])) args$palette <- FALSE
    if(is.null(args[["main"]])) args$main <- NULL
    if(is.null(args[["cex"]])) args$cex <- 0.8
    if(is.null(args[["width"]])) args$width <- 1
    if(inherits(chain, "mcmc")){
      tot <- nrow(chain)
      if(burnin>tot) stop("Error! the burnin value is higher than the chain length")
      chain <- chain[c(burnin:tot),-1]
      meanRate <- apply(chain, 2, args$fun)
      if(args$log==TRUE) meanRate <- log(meanRate)
    }else{
      meanRate <- chain
      if(args$log==TRUE) meanRate <- log(chain)
    }
    #check the order of the tree; prunning algorithm use "postorder"
    #if(attr(phylo,"order")!="postorder") phylo <- reorder.phylo(phylo, "postorder")
    #colors mapping
    if(any(args$palette==FALSE)){
      Colors = colorRampPalette(args$color)( 100 )
    }else{
      Colors = args$palette
    }
    #0 index induce error I scale it between 1 and 100
    linScale <- function(x, from, to) round( (x - min(x)) / max(x - min(x)) * (to - from) + from)
    col <- linScale(meanRate, from=1, to=100)
    if(args$scale==TRUE){
      phylo$edge.length <- phylo$edge.length*meanRate
    }
    plot(phylo, edge.color = Colors[col], show.tip.label = args$show.tip.label, main = args$main, cex = args$cex, edge.width = args$width)
    image.plot(z = as.matrix(meanRate),col = Colors,
               legend.only = T, horizontal = args$horizontal)
  }
plotShifts(reorder(subtree_unroot1,"cladewise") , chain=result$median.rates, color=c("blue","light blue", "green", "yellow", "orange", "dark orange", "red"))

#### jump-rbm ####
#Some parameters
prop = 1.5        # tuning parameter for the proposal window (it's just a starting value it will be improved in the model fit I think)
iterations = 10e6  # number of iterations for the chain
sampling = 10000 
model = "jump-rbm"     # relaxed brownian motion (the model used in Eastman et al. 2011. More or less similar to BayesTrait)
filename2 = paste("relaxedBM.testmcmc-jumprjmcmcREARRANGED.log", sep="", collapse="")

#run the rjmcmc
rjmcmc.bm(tree, data, prop.width=prop, ngen=iterations, samp=sampling, filebase=filename2,
          simple.start=TRUE, type=model) # note: here I assume you're estimating a nuisance parameter as well
chain_rearranged2 <- load.rjmcmc("relaxedBM.testmcmc-rjmcmcREARRANGED.log")

result2=plot(chain_rearranged2, par="jumps", type= "fan", legend = F, cex = 0.6)
result2a=plot(chain_rearranged2, par="shifts", type= "fan", legend = F, cex = 0.6)