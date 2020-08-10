##### auteur example #####
# install.packages("Imperial/UEBS/Papers/auteur/", repos = NULL, type = "source")
library(auteur)
n = 24
while (1) {
  phy=prunelastsplit(birthdeath.tree(b=1,d=0,taxa.stop=n+1))
  phy$tip.label=paste("sp",1:n,sep="")
  rphy=reorder(phy,"pruningwise")
  
  # find an internal edge
  anc=get.desc.of.node(Ntip(phy)+1,phy)
  branches=phy$edge[,2]
  branches=branches[branches>Ntip(phy) & branches!=anc]
  branch=branches[sample(1:length(branches),1)]
  desc=get.descendants.of.node(branch,phy)
  if(length(desc)>=4) break()
}

rphy = phy

rphy$edge.length[match(desc,phy$edge[,2])]=phy$edge.length[match(desc,phy$edge[,2])]*64

e=numeric(nrow(phy$edge))
e[match(c(branch,desc),phy$edge[,2])]=1
cols=c("red","gray")
plot(phy,edge.col=ifelse(e==1,cols[1],cols[2]), edge.width=2)
mtext("expected pattern of rates")

## simulate data on the 'rate-shifted' tree
dat=rTraitCont(phy=rphy, model="BM", sigma=sqrt(0.1)) # use to simulate PCA data

## run two short reversible-jump Markov chains
r=paste(sample(letters,9,replace=TRUE),collapse="")
lapply(1:2, function(x) rjmcmc.bm(phy=phy, dat=dat, ngen=10000, sample.freq=10, prob.mergesplit=0.1, simplestart=TRUE, prop.width=1, fileBase=paste(r,x,sep=".")))

# collect directories
dirs=dir("./",pattern=paste("BM",r,sep=".")) # rjmcmc results have been stored in documents (the wd in this case)
pool.rjmcmcsamples(base.dirs=dirs, lab=r)

# view rda file content
load(paste(paste(r,"combined.rjmcmc",sep="."),paste(r,"posteriorsamples.rda",sep="."),sep="/"))
print(head(posteriorsamples$rates))
print(head(posteriorsamples$rate.shifts))

## plot Markov sampled rates
plot.new()
shifts.plot(phy=phy, base.dir=paste(r,"combined.rjmcmc",sep="."), burnin=0.5, legend=TRUE, edge.width=2)

# clean-up: unlink those directories
unlink(dir(pattern=paste(r)),recursive=TRUE)



### phenogram
phenogram(phy = phy, model = "bm", plot = T)
# try this with shark data
shark <- read.nexus("../2.MS-Ecomorphology_lower_jaw/Data _ Analyses/Adam_ML_SharkTree.nex.txt")
phenogram(shark, model = "bm", plot = T) # doesn't work

##### use geiger #####
library(geiger)
geiger::rjmcmc.bm(phy = phy, dat = dat, ngen = 1000, samp = 10, type = "bm")
res <- load.rjmcmc(x = "BM.result/")

plot(x = res, par = "shifts", burnin = 0.25, legend = T, show.tip = F, edge.width = 2) # colour error
plot(x = res, par = "jumps", burnin = 0.25, legend = T, show.tip = F, edge.width = 2) # no error
