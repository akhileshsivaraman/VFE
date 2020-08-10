##### rjmcmcm.bm example #####
# https://rdrr.io/cran/geiger/man/rjmcmc.bm.html
library(phytools)
library(geiger)
library(ape)
library(coda)
setwd("~/Documents/Imperial/UEBS/Code")

# create mock tree
phy <- ladderize(sim.bdtree(n=200), right = F)
r <- paste(sample(letters, 9, replace = T), collapse = "")
defpar <- par(no.readonly = T)

tmp <- ex.jumpsimulator(phy, jumps = 10)
dat <- tmp$dat
hist <- tmp$hist

ex.traitgram(phy, hist, alpha = 0) # plot history of trait change

## run analysis
rjmcmc.bm(phy, dat, ngen = 20000, samp = 500, filebase = r, type = "jump-bm")
outdir <- paste("jump-BM", r, sep = ".")
ps <- load.rjmcmc(outdir)
plot(x = ps, par = "jumps", burnin = 0.25)

## use coda to explore mcmc run
autocorr.plot(ps$log, ask = dev.interactive())
plot(ps$log, ask=dev.interactive())


# another example
scl <- ex.ratesimulator(phy, min=12, show.tip=FALSE)
dat <- rTraitCont(scl) # simulates continuous character evolution - check methods
rjmcmc.bm(phy, dat, ngen = 20000, samp = 500, type = "rbm")
ps <- load.rjmcmc("relaxedBM.result/")
plot(x = ps, par = "shifts", burnin = 0.25, legend = T, show.tip = T)


shifts <- ps$shifts
jumps <- ps$jumps
rates <- ps$rates
log <- ps$log

