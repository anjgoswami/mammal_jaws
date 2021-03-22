## Packages needed#
library(Morpho)
library(geomorph)
library(ape)
library(geiger)
library(phytools)
library(RRphylo)
library(dispRity)
library(devtools)
install_github("JClavel/mvMORPH", ref="devel_1.0.4")
library(mvMORPH)

#set the working directory#
setwd("C:/Users/Anne-Claire/Dropbox/all_extant_marsupial_crania/Mandible/Manuscript/Fabre_et_al_jaw_paper/Script/Data")

#Entire data set#
load(file="./means.Rdata")# procrustes superimposed coordinates for the entire data set#
load(file="./phy.Rdata")# slided coordinates#

#Extant species data set
load(file="./means_extant.Rdata")# procrustes superimposed coordinates for extant species data set#
load(file="./phy_extant.Rdata")# phylogeny pruned with extant species only#

#load the identification#
ident2<-read.csv(file="./Ident_tout_mean.csv", row.names = 1)


td<-treedata(phy,means) #compare taxa in data and tree#
phy=td$phy# pruned the tree by keeping taxa found in the data#

means=td$data # pruned the data by keeping taxa found in the phylogeny#

means<-means[match(phy$tip.label,rownames(means)),]#ordering data in the same order than phylogeny#
ident2<-ident2[match(phy$tip.label,rownames(ident2)),]#ordering identification in the same order than phylogeny#
 
 Infraclass<-ident2$Infraclass
 type<-ident2$type
 diet<-ident2$Diet3
 
 
 
 Y<-arrayspecs(means,114,3)#convert entire data set to an array#
 


########################################################################## 
#############Principal component analyses on entire data set#
##########################################################################
 PCA<-gm.prcomp(Y)
 
 PCscores<-PCA$x
 PCscores<- PCscores[match(phy$tip.label,rownames( PCscores)),]
 x <- PCscores[,1]
 y <- PCscores[,2]
 z <- PCscores[,3]
 d<-PCscores[,4]
 p13<-cbind(PCscores[,1],PCscores[,3])
 p14<-cbind(PCscores[,1],PCscores[,4])
 p23<-cbind(PCscores[,2],PCscores[,3])
 
 
 #Phylomorphospace

  color2<-c("red","gold" )
  
##############Axis 1,2 (you can explore the data depending on diet using the diet category rather than the infraclass###
 phylomorphospace(phy,PCscores[,1:2],a=NULL,label = "off",xlab="Principal component 1 (32.14%)",ylab="Principal component 2 (15.06%)", node.size=0)
 points(x,y,col=color2[Infraclass],pch=c(16, 17)[Infraclass],cex=2)
 
 phylomorphospace(phy,PCscores[,1:2],a=NULL,label = "horizontal",xlab="Principal component 1 (32.14%)",ylab="Principal component 2 (15.06%)", node.size=0)
 points(x,y,col=color2[Infraclass],pch=c(16, 17)[Infraclass],cex=2)
 
 

 #############################################################################################
 ###############################extracting diet and infraclass categories for extant data set only
 # Defining groups and taxa.
 gA <- ident2$type # contains data on fossil/extant
 taxaA <- rownames(ident2) # species names

 #group definition to select only extant taxa in identification and defining categories for extant data set
 testtaxa <- rownames(ident2[gA=="Fossil",])
 testtaxan <- row(ident2)[gA=="Fossil",1]
 trainingtaxa <- rownames(ident2[-testtaxan,]) # creating a dataframe that only contains taxa with known group affiliation.
 
 
 g <- gA[-testtaxan]#vector without fossil
 names(Infraclass)<- rownames(ident2)
 Infraclass_extant<-Infraclass[-testtaxan]
 Diet3<-ident2$Diet3
 names(Diet3)<- rownames(ident2)
 diet3_extant<-Diet3[-testtaxan]
 

 
 ##########################################################################################################################
 #######################Phylogenetic principal component analyses, preparing the data that will be input in Bayestrait#
 ########################################################################################################################
 ########For the entire data set
 ppca<-phyl.pca( phy, means, method="BM", mode="cov")
 Ppcscores <-ppca$S
 Ppcscores <- Ppcscores[match(phy_g$tip.label, rownames(Ppcscores)),]
 PC95<-Ppcscores[,1:32]#can be use as input in bayesian analyses
 ######For the extant data
 ppca<-phyl.pca( phy_extant, means_extant, method="BM", mode="cov")
 Ppcscores <-ppca$S
 Ppcscores <- Ppcscores[match(phy_extant$tip.label, rownames(Ppcscores)),]
 PC95<-Ppcscores[,1:32]]#can be use as input in bayesian analyses

 
 
 ###########################################################################
###############################DISPARITY 
##################################################################


####Disparity on entire data set#############################
#########################################################################################
############Geomorph disparity on entire data set

 gdf2<- geomorph.data.frame(coords=Y, Infraclass = Infraclass, diet=diet)#geomorph dataframe#

disparity_entire<-morphol.disparity(coords~Infraclass, groups = gdf2$Infraclass, data = gdf2,iter = 999)
 summary(disparity_entire)
 morphol.disparity(f1 = coords ~ Infraclass, groups = gdf2$Infraclass,iter = 999, data = gdf2) 
 
############dispRity disparity on entire data set

gdf2<- geomorph.data.frame(coords=means, Infraclass = Infraclass, diet=diet)#geomorph dataframe#
 
disparity<- dispRity.per.group(means,list( Eutheria=c(1:75,129:148), Metatheria=c(76:128,149:151)), metric=function(X) return(sum(X^2)/nrow(X)))
 
summary(disparity);
 
 
 
###Disparity extant taxa#############################
#########################################################################################
############Geomorph disparity on extant data set
 Y_extant<-arrayspecs(means_extant,114,3)#convert means_extant to an array#
 gdf2extant<- geomorph.data.frame(coords2= Y_extant, Infraclass_extant = Infraclass_extant)
 
 
 disparity_extant<-morphol.disparity(coords2~Infraclass_extant, groups =  gdf2extant$Infraclass_extant, data = gdf2extant,iter = 999)
 summary(disparity_extant)
 
###############################DispRity disparity  on extant data set#
 disparity<- dispRity.per.group(means_extant,list( Eutheria=c(1:75), Metatheria=c(76:128)), metric=function(X) return(sum(X^2)/nrow(X)))
 
 summary(disparity);
 
 
 ##########################################################################################
 #################################TEST OF CONVERGENCE (RRPHYLO)##############################
 ##########################################################################################
PCscores<- PCscores[,1:24]#95% of the overall variation#


Lingual<-ident2$Lingual
names(Lingual)<-rownames(ident2)

Carnivorous<-ident2$Carnivorous
names(Carnivorous)<-rownames(ident2)

Omnivorous<-ident2$Omnivorous
names(Omnivorous)<-rownames(ident2)

Folivorousbrowser<-ident2$Folivorousbrowser
names(Folivorousbrowser)<-rownames(ident2)

Mixed_feeder<-ident2$Mixed_feeder
names(Mixed_feeder)<-rownames(ident2)

Tuberivorous<-ident2$Tuberivorous
names(Tuberivorous)<-rownames(ident2)

Grazer<-ident2$Grazer
names(Grazer)<-rownames(ident2)

Frugivorous<-ident2$Frugivorous
names(Frugivorous)<-rownames(ident2)


Insectivorous<-ident2$Insectivorous
names(Insectivorous)<-rownames(ident2)



 RRphylo(phy,PCscores)->RRfel 
 
 
 ## Case 2. searching convergence within a single state
 search.conv(tree=phy, y=PCscores, state=Lingual,foldername = getwd())->SC.Lingual

 search.conv(tree=phy, y=PCscores, state=Carnivorous,foldername = getwd())->SC.Carnivorous

 search.conv(tree=phy, y=PCscores, state=Omnivorous,foldername = getwd())->SC.Omnivorous
 
 search.conv(tree=phy, y=PCscores, state=Mixed_feeder,declust=TRUE,foldername = getwd())->SC.Mixed_feeder

 search.conv(tree=phy, y=PCscores, state=Tuberivorous,declust=TRUE,foldername = getwd())->SC.Tuberivorous
 
 search.conv(tree=phy, y=PCscores, state=Folivorousbrowser,declust=TRUE, foldername = getwd())->SC.Folivorousbrowser

 search.conv(tree=phy, y=PCscores, state=Grazer,declust=TRUE, foldername = getwd())->SC.Grazer

 search.conv(tree=phy, y=PCscores, state= Frugivorous,declust=TRUE,foldername = getwd())->SC.Frugivorous

 search.conv(tree=phy, y=PCscores, state= Insectivorous,foldername = getwd())->SC.Insectivorous

resultppca<- rbind(SC.Lingual,SC.Carnivorous,SC.Omnivorous,SC.Folivorousbrowser,SC.Mixed_feeder,SC.Tuberivorous,SC.Grazer,SC.Frugivorous,SC.Insectivorous)
write.csv(resultppca,file="./resultppca.csv") 



##########################################################################################
##################################PHYLOGENETIC MANOVA#
##########################################################################################

Diet3<-ident2$Diet3
names(Diet3)<-rownames(ident2)

datas=list(shape=means,Diet3=Diet3,Infraclass=Infraclass)

###testing for jaw shape difference depending on diet on the entire dataset
##################################################################################
fit<-mvgls(shape~Diet3, data=datas,  phy, model="lambda", method=c("PL-LOOCV"))

#performing multivariate tests on generalized least squares linear model#
(multivariate_test <- manova.gls(fit, nperm=999, type="II", test="Pillai"))

#testing for jaw shape difference depending on infraclass (reproductive strategy) on the entire dataset
####################################################################################
fit1<-mvgls(shape~Infraclass, data=datas,  phy, model="lambda", method=c("PL-LOOCV"))

#performing multivariate tests on generalized least squares linear model#
(multivariate_test2 <- manova.gls(fit1, nperm=999, type="II", test="Pillai"))



###phylogenetic MANOVA on extant data set#
###################################################################################################
names(Infraclass_extant)<-rownames(means_extant)
names(diet3_extant)<-rownames(means_extant)

datae=list(shape=means_extant,Diet3=diet3_extant,Infraclass=Infraclass_extant)


####testing for jaw shape difference depending on infraclass (reproductive strategy) on the extant dataset
##################################################################################
fit2<-mvgls(shape~Infraclass, data=datae,  phy_extant, model="lambda", method=c("PL-LOOCV"))

#performing multivariate tests on generalized least squares linear model#
(multivariate_test3 <- manova.gls(fit2, nperm=999, type="II", test="Pillai"))

####testing for jaw shape difference depending on diet on the extant dataset
##################################################################################
fit3<-mvgls(shape~Diet3, data=datae,  phy_extant, model="lambda", method=c("PL-LOOCV"))

#performing multivariate tests on generalized least squares linear model#
(multivariate_test3 <- manova.gls(fit3, nperm=999, type="II", test="Pillai"))



########################################################################################################
###################################RATE OF EVOLUTION DEPENDING ON INFRACLASS (REPRODUCTIVE STRATEGIE)
########################################################################################################

###Rates depending on infraclass on the entire data set###########

# Making the simmap tree with mapped states
treec<-make.simmap(phy, Infraclass , model="ER", nsim=1)
col<-c("red","gold" ); names(col)<-c("Eutheria","Metatheria")

# Model fit
data <- list(Y=means)
#fit <- mvgls(Y~1, data=data, tree=tree, model="BMM", method = "H&L", error=TRUE) #


tree = phy
ancestral = treec
tips = Infraclass
paintAllTree <- function(tree, ancestral, tips){  
  
  if(inherits(ancestral, "describe.simmap")){
    names_grps <- colnames(ancestral$ace)
    statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
  }else{
    names_grps <- colnames(ancestral$lik.anc)
    statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
  }  
  
  combined = as.character(c(tips, statesNodes))
  treebis=tree
  for(i in sort(tree$edge[,2])){
    treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
  }  
  return(treebis)
}


# Estimation in ML and use these reconstructions to make the "simmap' trees ( as used in mvMORPH)#
ace_habitat <- ace(Infraclass, phy, type = "discrete", model="ARD")

# plot the reconstructions
plot(phy)
nodelabels(pie=ace_habitat$lik.anc)

# Convert ML reconstructions ML in SIMMAP trees
simmap_ace <- paintAllTree(tree, ace_habitat, as.character(Infraclass))
data <- list(Y=means)


## collection of SIMMAP trees
nsim = 100 # number of stochastic mapping
my_trees<-make.simmap(tree, Infraclass , model="ARD", nsim=nsim)
distrib <-summary(my_trees,plot=TRUE)


simulations <- sapply(1:nsim, function(x) {
  # model fit
  fit <- mvgls(Y~1, data=data, tree=my_trees[[x]], model="BMM", method = "H&L", error=TRUE)
  fit$param # estimate rates
})

boxplot(t(simulations))
rowMeans(simulations)

###Rates depending on infraclass on the extant data set###########
###################################################################

# Making the simmap tree with mapped states
treec<-make.simmap(phy_extant, Infraclass_extant , model="ER", nsim=1)
col<-c("red","gold" ); names(col)<-c("Eutheria","Metatheria")

# Model fit
data <- list(Y=means_extant)
#fit <- mvgls(Y~1, data=data, tree=tree, model="BMM", method = "H&L", error=TRUE) #


tree = phy_extant
ancestral = treec
tips = Infraclass_extant


# On estime en ML et on utilise les reconstructions pour faire un arbre "simmap", le format utilisé par mvMORPH...
ace_habitat <- ace(Infraclass_extant, phy_extant, type = "discrete", model="ARD")

# plot les reconstructions
plot(phy_extant)
nodelabels(pie=ace_habitat$lik.anc)


data <- list(Y=means_extant)


nsim = 100 # nombre de mapping stochastiques
my_trees<-make.simmap(phy_extant, Infraclass_extant , model="ARD", nsim=nsim)
distrib <-summary(my_trees,plot=TRUE)

simmap_summary <- paintAllTree(phy_extant, distrib, as.character(Infraclass_extant))


simulations <- sapply(1:nsim, function(x) {
  # model fit
  fit <- mvgls(Y~1, data=data, tree=my_trees[[x]], model="BMM", method = "H&L", error=TRUE)
  fit$param # Les taux estimés
})

boxplot(t(simulations))
rowMeans(simulations)





###############################################################################################
######################################Branch-specific rate reconstructions in BayesTraits 
###############################################################################################
##################################Checking trace#################################
#loading packages
library(BTRTools)
library(phytools)
library(coda)



tracePlots <- function(file, burnin=0, thinning=1, plot=TRUE, display=c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")){
  require(BTRTools)
  require(coda)
  
  rjout <- loadRJ(file, burnin = burnin, thinning = thinning)
  chain_out <- type.convert(rjout$rj_output)
  rownames(chain_out) = chain_out[,"It"]
  chain_out = chain_out[,-1]
  # Retrieve numerical
  index <- sapply(chain_out,function(x) is.numeric(x))
  chain <- mcmc(chain_out[,index])
  
  # plot the trace
  if(plot){
    plot(chain[,display])
  }
  
  # Just compute some statistics (autocorrelation...)
  cat("Effective sample size:","\n")
  print(effectiveSize(chain[,display]))
  
  # return results
  invisible(chain)
  
}

#(This will return the ESS (effective sample size) and trace plot)#
#If you need to compare several independent runs (e.g. using general mcmc diagnostic functions in "coda" package)#
###############################################
test = tracePlots(file=paste("./run 1/Pc95_tout_BM.txt.VarRates.txt"),burnin=4000)
test2 = tracePlots(file=paste("./run 2/Pc95_tout_BM.txt.VarRates.txt"),burnin=4000) #another dataset
#####other runs can be add in the same way ...#


my_list_of_chains = mcmc.list(list(test[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")], test2[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")]))



#Then simply use any functions from "coda" that works with objects of class "mcmc.list"
plot(my_list_of_chains)

#Gelman and Rubin's convergence diagnostic
gelman.diag(my_list_of_chains, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)
			
			
#Plotting rate onto tree#
source("./mytreerateplotter.R")#	available on Github here: https://github.com/anjgoswami/salamanders/blob/master/mytreerateplotter.R#



{
  ## load tree
  tree=phy
  
  ## EVOLUTIONARY MODELS ##
  # BAYES TRAITS- JAW ##
  bt.tree <- tree
  bt.tree <- reorder(bt.tree, "cladewise")
  bt.tree.node.no <- makeNodeLabel(bt.tree)
  rate.heatmap <- colorRampPalette(c("blue", "lightblue", "yellow", "orange", "red"), bias=1)
  plot(bt.tree, use.edge.length = F)
  nodelabels(text=bt.tree.node.no$node.label, frame="none", cex=0.5)
}



L_trad_tree_BTraits_tout<-BTRTools::rjpp(rjlog = paste("./Pc95_tout_BM.txt.VarRates.txt", sep=""),
                                         rjtrees = paste("./Pc95_tout_BM.txt.Output.trees", sep=""),
                                         tree = bt.tree, thinning = 1, burnin=4000)

cophylotrees <- cophylo(bt.tree,L_trad_tree_BTraits_tout$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
tree_trad_tout.a <- cophylotrees$trees[[1]]


#Categories to plot on tip of tree, not mandatory, here we will plot the infraclass#
Cat<-ident2
td <- treedata(tree_trad_tout.a, Cat);
Cat <- td$data;
Cat<- Cat[match(tree_trad_tout.a$tip.label, rownames(Cat)),]


Infraclass<-Cat[,1]

Infraclassbis<-as.factor(Infraclass)

cols<-setNames(c("red","gold"), levels(Infraclassbis))
col_tip = cols[Infraclassbis[tree_trad_tout.a$tip.label]]
zzz<-mytreebybranch(tree = tree_trad_tout.a, x=log(L_trad_tree_BTraits_tout$data$meanRate)[-1], 
               mode="edges",type="fan",
               palette = color3,cex=.3,
               edge.width=3,
               show.tip.label= TRUE, 
               tip.label.offset=4,
               tip.pch=20, 
               tip.cex=1, 
               tip.cols=col_tip,
               tip.offset=4)

##########################plot shift at the node

n<-BTRTools::plotShifts(L_trad_tree_BTraits_tout,type="fan",scalar="node", tips=TRUE, cex=0.3,shp = 24, colour = "grey")

####################function to obtain post probabilies
return_pprob <- function(PP, threshold=0.5, cl = "nOrgnNRate"){
  nodes <- PP$data$descNode[which((PP$data[ , cl] / PP$niter) >= threshold)]
  pprobs <- PP$data[which((PP$data[ , cl] / PP$niter) >= threshold) , cl] / PP$niter
  
  result = list(nodes=nodes, pprobs=pprobs)
  return(result)
}

#Plotting post probabilies on of the shift at the node
return_pprob(L_trad_tree_BTraits_tout, threshold = 0)

pp = return_pprob(L_trad_tree_BTraits_tout, threshold = 0)

fact=2
index_nodes = which(pp$nodes>Ntip(new.bt.tree))
pp$pprobs[index_nodes]
pp$nodes[pp$nodes>Ntip(new.bt.tree)]
plot(new.bt.tree, type="fan", show.tip.label = F)
nodelabels(cex=pp$pprobs[index_nodes]*fact, pch=24)

plot(new.bt.tree, type="fan", show.tip.label = F)
nodelabels(round(pp$pprobs[index_nodes], digits=3),cex=0.5)
 



