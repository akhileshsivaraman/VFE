##this takes a vector of AIC scores and spits out a table with delta AICs and Akaike weights, in descending order of fit.

AkaikeWeights<-function(AIC){
if(is.null(names(AIC)))
return("AIC scores have no associated names - you are a fucktard")

best<-min(AIC)
deltaAIC<-AIC-best
sumDeltaAIC<-sum(exp(-0.5 * deltaAIC))
Weights<-(exp(-0.5 * deltaAIC)/sumDeltaAIC) 

results<-data.frame("AIC"=AIC,"dAIC"=deltaAIC,"Weights"=Weights)

# original
#return(results[order(results$Weight,decreasing=T),])
#}
return(results[order(results$dAIC,decreasing=F),])
}
