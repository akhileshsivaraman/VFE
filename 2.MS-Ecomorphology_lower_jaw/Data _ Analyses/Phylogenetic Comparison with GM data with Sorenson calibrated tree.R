##set working directory
setwd("/Users/akhileshsivaraman/Documents/Imperial/UEBS/2.MS-Ecomorphology_lower_jaw/Data _ Analyses")

##load libraries for analysis
library(ape)
library(geiger)
library(OUwie)
library(ouch)

###### alignment phylogenetic data with morphological data ######

##load mandible 3D data (MS-DOS).txt with PC scores of interest.
mc<-read.table("N=153 MORPHOJ PC1-3 scores MCL averaged by spp.txt",header=T,row.names=1)

##load body shape 2D data (MS-DOS).txt with PC scores of interest
#body<-read.table("filename",header=T,row.names=1)

##load newick/tre formatted tree: convert .newick file, -> .tre file
#tree<-read.tree(file = "Sorenson_tree.tre", text = NULL, tree.names = NULL, skip = 0, comment.char = "#", keep.multi = FALSE)
##load nexus formatted tree
tree<-read.nexus(file="Adam_ML_SharkTree.nex",tree.names = NULL, force.multi = FALSE)

##compare taxa in data and tree and drop tips from the tree that aren't included in dataset
newtree<-treedata(tree,mc)
newtree
plot(newtree$phy, cex=0.45)

#bodytree<-treedata(tree,body)
#bodytree
#plot(bodytree$phy)

###### fitcontinuous mandible 3D ######

##model fitting of modes of evolution

brownfit<-fitContinuous(newtree$phy,newtree$data[,1],model="BM")
brownfit
ebfit<-fitContinuous(newtree$phy,newtree$data[,1],model="EB")
ebfit
oufit<-fitContinuous(newtree$phy,newtree$data[,1],model="OU")
oufit
aic<-c(brownfit$opt$aic,ebfit$opt$aic,oufit$opt$aic)
aic

names(aic) <- c("BM","EB","OU")
aic
source("AKWeights.R")
TableAIC<-AkaikeWeights(aic)
write.table(TableAIC, "N=153 TableAICmodels_MCL.txt", quote = FALSE, sep = "\t")

###### fitcontinuous bodyshap 2D ######

#brownfit<-fitContinuous(bodytree$phy,bodytree$data[,1],model="BM")
#brownfit
#ebfit<-fitContinuous(bodytree$phy,bodytree$data[,1],model="EB")
#ebfit
#oufit<-fitContinuous(bodytree$phy,bodytree$data[,1],model="OU")
#oufit
#aic<-c(brownfit$opt$aic,ebfit$opt$aic,oufit$opt$aic)
#aic

#names(aic) <- c("BM","EB","OU")
#aic
#source("AKWeights.R")
#TableAIC<-AkaikeWeights(aic)
#write.table(TableAIC, "TableAICmodels_mand.txt", quote = FALSE, sep = "\t")

###### DTT (Disparity Through Time) ######

dispmc<-disparity(phy=NULL,newtree$data,index=c("avg.sq"))
dispmc

dttmand<-dtt(newtree$phy,newtree$data,index=c("avg.sq"),nsim=1000)
dttmand
