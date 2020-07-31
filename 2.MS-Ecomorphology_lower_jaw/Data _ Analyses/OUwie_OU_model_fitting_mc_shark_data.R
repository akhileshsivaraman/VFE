## script template provided by Bruno Frederich

###### Model fitting Analysis for phenotypic data grouped by DIET regime ######

## set working directory
setwd("~/Documents/Imperial/UEBS/2.MS-Ecomorphology_lower_jaw/Data _ Analyses")
library(geiger)
library(OUwie)
library(phytools)

## to load a newick formatted phylogeny file use read.tree("filename.newick")
## load a nexus formatted phylogeny file.
tree<-read.nexus(file="Adam_ML_SharkTree.nex",tree.names = NULL, force.multi = FALSE)

## diet data to run models on diet regimes: 1=fish_&_cephalopoda, 2=Soft-bodied_benthic_invertebrates, 3=decapoda_&_mollusca, 4=zooplankton, 5=generalist
diet<-read.table("N=153 Classifiers diet MC averaged by spp (sp, order, preycode) NA diets pruned.txt",header=T,row.names=1)

## prune tree and diet data set so the species congrue between both
newtree<-treedata(tree,diet)

## create a vector for the trait data
vec_trait<-newtree$data[,2]
names(vec_trait)<-rownames(newtree$data)

## map discrete character on your tree: model=ER when you test with equal rates of change of 1 state into another, model=ARD when you test your data with the assumptioon that the rate of change is different for every branch. Especially important if the rate of changing back to an ancestral state is less likely. 
simmap <- make.simmap(newtree$phy, vec_trait, model="ER", nsim=1)
plotSimmap(simmap, fsize=0.63)

## dataframe morphological and dietary data formatted in columns as: species;PC1; dietcode
data1<-read.table("N=153 spp PC1 dietcode MCL averaged by spp NA diet pruned.txt", header=T,row.names=1)

## single-rate Brownian motion
bm1<-OUwie(simmap,data1,model="BM1", simmap.tree=TRUE)
## Brownian motion with different rate parameters for each state on a tree
#bms<-OUwie(simmap,data=data,model="BMS", simmap.tree=T)
## Ornstein-Uhlenbeck model with a single optimum for all species
ou1<-OUwie(simmap,data=data1,model="OU1", simmap.tree=T)
## Ornstein-Uhlenbeck model with different state means and a single alpha and sigma^2 acting all selective regimes
oum<-OUwie(simmap,data=data1,model="OUM", simmap.tree=T)
## Ornstein-Uhlenbeck model that assumes different state means as well as multiple sigma^2 (=standard deviations)
#oumv<-OUwie(simmap,data=data1,model="OUMV", simmap.tree=T)
## Ornstein-Uhlenbeck model that assumes different state means as well as multiple alpha (=rates of diversification), this a constrained model, only useful if you can explain constrains morphologically 
ouma<-OUwie(simmap,data=data1,model="OUMA", simmap.tree=T)
## Ornstein-Uhlenbeck model that assumes different state means as well as multiple alpha and sigma^2 per selective regime
#oumva<-OUwie(simmap,data=data1,model="OUMVA", simmap.tree=T)
## fit continuous function to include Early Burst (EB) model
newtree2<-treedata(tree,data1)
vec_trait2<-as.numeric(newtree2$data[,"PC1"])
names(vec_trait2)<-rownames(newtree2$data)
ebfit<-fitContinuous(newtree2$phy,vec_trait2,model="EB")

plotSimmap(bm1$phy)
plotSimmap(ou1$phy)
plotSimmap(oum$phy)
plotSimmap(ouma$phy)
# all trees are the same

###### AIC comparison ######
aicScores<-c(bm1$AICc,ou1$AICc,oum$AICc,ouma$AICc,ebfit$opt$aicc)
names(aicScores)=c("BM1","OU1","OUM","OUMA","EB")
#min(aicScores)
#delta.aic <- aicScores - min(aicScores)
source("AKWeights.R")
TableAIC <- AkaikeWeights(aicScores)
write.table(TableAIC, "N=153 TableAICmodels OUwie MCL diet with simmap (BM1,OU1,OUM,OUMA,EB).txt", quote = FALSE, sep = "\t")

###### AIC comparison without complex OUMA model ######
aicScores2<-c(bm1$AICc,ou1$AICc,oum$AICc,ebfit$opt$aicc)
names(aicScores2)=c("BM1","OU1","OUM","EB")
#min(aicScores2)
#delta.aic <- aicScores2 - min(aicScores2)
source("AKWeights.R")
Table2AIC <- AkaikeWeights(aicScores2)
write.table(Table2AIC, "N=153 TableAICmodels OUwie MCL diet with simmap (BM1,OU1,OUM,EB).txt", quote = FALSE, sep = "\t")


