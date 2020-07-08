#############################################################################################
## Compare community modelling methods for bioregionalisation of simulated species dataset
## March 2018. N.Hill with input from S. Woolley
#############################################################################################

# Modelling process
# 1) RUN MODELS, DIAGNOSTICS, DETERMINE NUMBER OF GROUPS, OR SET GROUP=3, PREDICT GROUPS ACROSS SIMULATED REGION
# 2) Plot predicted distribution of groups across simulated region
# 3) Describe content of groups 
# 4) Describe environment of groups

# MODELS
# A) Cluster environment only
# 2 Stage Models- cluster then predict:
# B) cluster bioloigcal data, predict clusters with random forests
# 2 Stage Models- predict then cluster
# C) predict species using random forests, then cluster predictions
# D) predict species using Mistnet, then cluster predictions
# E) predict species using HMSC, then cluster predictions
# F) predict dissimilarities using GDM, cluster predicted dissimilarities
# G) predict biologically transformed environment, cluster predictions
# 1 Stage mdoel-based: cluster and predict
# H) Species Archetype Models
# I) Regions of Common Profile

#######################
## Set up---
#######################
# Get required libraries
library(vegan)          #Jaccard Index
library(cluster)        #Hierarchical clustering
library(fpc)            #Cluster statistics
library(extendedForest) # Random forests with conditional variable importance
#library(randomForest)  #Random Forests
#library(bbgdm)         #Bayesian Bootstrap Generalised Dissimilarity Models (and naive GDM)
library(gdm)            #Generalised Dissimilarity Models (GDM)
library(RCPmod)         #Regions of Common Profile Models v2.188
#devtools::install_github('skiptoniam/ecomix@dev')  #Species Archetype Models (SAMs)
library(ecomix) 
library(HMSC)           #Hierarchical Modelling of Species Communities
library(gradientForest) #Gradient Forests (an extension of Random Forests)
#devtools::install_github('davharris/mistnet')  #Mistnet stochastic, multispecies neural network
library(mistnet)
library(RColorBrewer)   #colour palettes
library(rasterVis)      #plotting rasters

setwd("C:\\Users\\hillna\\UTAS_work\\Projects\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\Simulation\\")
source("Simulation_Additional_Funcs.R")

#Load required files
#files in "simulate_communities" folder on dropbox

load("Sim_Setup/Many_covars_sim_fin.RData") 
#load("Sim_Setup/Many_covars_sim.RData") 
#sim_data= simulated species probabilities, occurrences, species' groups for entire region
#sp_200 = matrix of occurrences of 30 species at 200 sites to use as species dataset for analysis
#env_200= matrix of corresponding environmental conditions at 200 site to use as environmental dataset for analysis (scaled and centred)

load("Sim_Setup/sim_env_070518.RData") 
# env= raster brick of environmental data for entire region
# env_dat= raster brick of environmental data for entire region converted to matrix

species<-paste0("Sp", 1:30)
env_vars<-dimnames(env_dat[,3:10])[[2]]
dimnames(sp_200)[[2]]<-species

# Used scaled and centred environmental data for model-based analyses 
#i.e. sim_dat (entire grid) and env_200 or sample sites

#######################################################
## 1) Run models and any diagnostics and 
##  a) determine optimal number of groups
##  b) set number of groups =3
##  c) generate predictions of clusters for sim region
#######################################################

## -----------------------------------
## A) ENVIRONMENT ONLY CLUSTERING ----
## -----------------------------------

# Cluster euclidean distance matrix using Ward's SS
env_clust<-hclust(dist(sim_dat, method="euclidean"), method="ward.D2")
plot(env_clust) 

#determine number of groups
grp_stats(min_grp=2, max_grp=15, tree_clust= env_clust, dissim=dist(sim_dat, method="euclidean"), method="ward.D2")

#sil width, bootstrap suggests 11, CH suggests 4.

#cut at 11 
env11_clust<-cutree(env_clust, k=11)

#Set number of groups to 3
env3_clust<-cutree(env_clust, k=3)

rm(env_clust)

## ------------------------------------------------------------------------------------------------
## B) 2 STAGE: CLUSTER BIOLOGICAL DATA; PREDICT GROUPS USING RANDOM FORESTS & ENVIRONMENTAL DATA---
## ------------------------------------------------------------------------------------------------

## Cluster Jaccard distance matrix using Ward's SS
jac_dist<-vegdist(sp_200, method="jaccard",binary=TRUE)
dat_clust<-hclust(jac_dist, method="ward.D2")
plot(dat_clust,  labels=FALSE, xlab="", sub="")

#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= dat_clust, dissim=jac_dist, method="ward.D2")
#All suggests 2 groups

#set 2 groups
bio2_clust<-cutree(dat_clust,2)
#Set 3 groups
bio3_clust.2<-cutree(dat_clust,3)


## Use Random Forest to match environment and predict clusters across sim region
# 2 groups
bio2_rf<-randomForest(y= as.factor(bio2_clust), x=env_200[,2:9], importance=TRUE, corr.threshold = 0.65)
bio2_rf$confusion 

bio2_rf_pred<- predict(bio2_rf, sim_dat, type=c("prob"))

# 3 groups
bio3_rf<-randomForest(y= as.factor(bio3_clust), x=env_200[,2:9], importance=TRUE)
bio3_rf$confusion #more confusion amongst groups
bio3_rf_pred<- predict(bio3_rf, sim_dat, type=c("prob"))

rm(jac_dist, dat_clust)

## ----------------------------------------------------------------
## C) 2 STAGE: RANDOM FOREST PREDICT EACH SPECIES THEN CLUSTER ----
## ----------------------------------------------------------------

# Loop through running RF for each species (as classification problem)
rf_sp_mods<-list()
for(i in 1:dim(sp_200)[2]){
  #rf_sp_mods[[i]]<-randomForest(y=as.factor(sp_200[,i]), x=env_200[, env_vars], importance=TRUE)
  rf_sp_mods[[i]]<-randomForest(y=as.factor(sp_200[,i]), x=env_200[, env_vars], importance=TRUE, corr.threshold = 0.65)
}
names(rf_sp_mods)<-paste0("Sp",1:30)

# Loop through predictions for each species, generating probabilty of ocurrence output
rf_sp_pred<-data.frame(sim_dat[,1:2], 
                       matrix(nrow=dim(sim_dat)[1], ncol=dim(sp_200)[2], 
                              dimnames=list(NULL, paste0("Sp", 1:30))))

for(i in 1:dim(sp_200)[2]){                       
  rf_sp_pred[,(i+2)]<-predict(rf_sp_mods[[i]], sim_dat, type="prob")[,2]
}

# Cluster predictions

ssdm_clust<-hclust(dist(rf_sp_pred[,species], method="euclidean"), method="ward.D2")

#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= ssdm_clust, dissim=dist(rf_sp_pred[,species], method="euclidean"), method="ward.D2")
#suggests 2 groups

# set 2 groups
ssdm2_clust<-cutree(ssdm_clust, k=2)
# set 3 groups
ssdm3_clust<-cutree(ssdm_clust, k=3)

rm(ssdm_clust)


## ----------------------------------------------------
## D) 2 STAGE: MISTNET & CLUSTER SPECIES' PREDICTIONS ----
## ----------------------------------------------------
# Optimising mistnet parameters with 1 hidden layer using slightly modified code from 
# https://github.com/davharris/mistnet/extras/BBS-analysis/mistnet_cross-validation.R

set.seed(1)

# Choose hyperparams ------------------------------------------------------
## parameters to optimise:
# number of latent variables (2-5)
# number of neurons in 1 hidden layer (8-24, in intervals of 2)
# number of importance samples (n.importance.samples: use default)
# number of routes for gradient descent (n.minibatch: Use default)

hyperparams = data.frame(
  n.minibatch = 25, #default value
  sampler.size = rep(2:5,9),
  n.importance.samples = 25, #default value
  n.layer1 = rep(c(8,10,12,14,16,18,20,22,24),each=4),
  learning.rate = 0.1,
  fit.seconds = 90
)


# Fitting code ------------------------------------------------------------

fit = function(x, y, hyperparams, i){
  MNet_mod = mistnet(
    x = x,
    y = y,
    layer.definitions = list(
      defineLayer(
        nonlinearity = rectify.nonlinearity(),
        size = hyperparams$n.layer1[i],
        prior = gaussian.prior(mean = 0, sd = .1)
      ),
      
      defineLayer(
        nonlinearity = sigmoid.nonlinearity(),
        size = ncol(y),
        prior = gaussian.prior(mean = 0, sd = .1)
      )
    ),
    loss=bernoulliLoss(),
    updater = adagrad.updater(learning.rate = hyperparams$learning.rate[i]),
    sampler = gaussian.sampler(ncol = hyperparams$sampler.size[i], sd = 1),
    n.importance.samples = hyperparams$n.importance.samples[i],
    n.minibatch = hyperparams$n.minibatch[i],
    training.iterations = 0,
    initialize.biases = TRUE,
    initialize.weights = TRUE
  )
  MNet_mod$layers[[1]]$biases[] = 1 # First layer biases equal 1
  
  start.time = Sys.time()
  while(
    difftime(Sys.time(), start.time, units = "secs") < hyperparams$fit.seconds[i]
  ){
    MNet_mod$fit(100)
    cat(".")
    # Update prior variance
    for(layer in MNet_mod$layers){
      layer$prior$update(
        layer$weights, 
        update.mean = FALSE, 
        update.sd = TRUE,
        min.sd = .01
      )
    }
    # Update mean for final layer
      MNet_mod$layers[[2]]$prior$update(
      layer$weights, 
      update.mean = TRUE, 
      update.sd = FALSE,
      min.sd = .01
    )
  } # End while
  
  MNet_mod
}

# Cross-validation --------------------------------------------------------
fold.ids<-sample( x=1:5, size= nrow(sp_200),replace=TRUE)

out = list()

for(i in 1:dim(hyperparams)[[1]]){
  cat(paste0("Starting iteration ", i, "\n"))
  for(fold.id in 1:max(fold.ids)){
    cat(paste0(" Starting fold ", fold.id, "\n  "))
    in.val=fold.ids==fold.id   #added this line
    in.train = fold.ids != fold.id
    MNet_mod = fit(
      x=as.matrix(env_200[in.train,2:9]),
      y=sp_200[in.train,],
      hyperparams = hyperparams,
      i = i
    )
    
    cat("\n evaluating")
    
    #calculate log likelihood of predictions (replaced with my code)
    log.lik = sum(
      dbinom(
        sp_200[in.val,],
        size = 1, 
        prob = predict(MNet_mod, as.matrix(env_200[in.val,2:9]),n.importance.samples=dim(sp_200[in.val,])[[1]]),
        log = TRUE
      )
    )
    
    cat("\n")
    out[[length(out) + 1]] = c(
      iteration = i, 
      fold = fold.id, 
      seconds = hyperparams$fit.seconds[i],
      loglik = mean(log.lik)
    )
    
  } # End fold
} # End iteration


# Save CV results ---------------------------------------------------------

mistnet.results = merge(
  x = as.data.frame(do.call(rbind, out)),
  y = cbind(iteration = 1:nrow(hyperparams), hyperparams)
)
save(mistnet.results, file="Results/mistnet.results.Rdata")

# fit final model ---------------------------------------------------------

logliks<-aggregate(mistnet.results, by=list(mistnet.results$iteration), mean)
#indicates 2 latent variables and hidden layer with 8 neurons

MNet_mod = fit(
  x = as.matrix(env_200[,2:9]),
  y = sp_200, 
  hyperparams, 
  which.max(logliks[,5])
)

save(MNet_mod, file = "Results/mistnet.model.Rdata")

#generate predicitons for each species
MNet_pred<-predict(MNet_mod, as.matrix(sim_dat[,2:9]), n.importance.samples=25)
MNet_pred<-apply(MNet_pred,1:2, mean)

#cluster predictions
MNet_clust<-hclust(dist(MNet_pred, method="euclidean"), method="ward.D2")
plot(MNet_clust) 

#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= MNet_clust, dissim=dist(MNet_pred, method="euclidean"), method="ward.D2")

#sil suggests 5 grps, ch 8 groups, boot 2 groups!
MNet5_clust<-cutree(MNet_clust,5)

#choose 3 groups
MNet3_clust<-cutree(MNet_clust,3)

rm(MNet_clust)

## -------------------------------------------------------------------------------------------------
## E) 2 STAGE: HEIRARCHICAL MODELLING OF SPECIES COMMUNITIES (HMSC) & CLUSTER SPECIES' PREDICTIONS ----
## -------------------------------------------------------------------------------------------------

#Auto = spatial corordinates to calculate spatial autocorrelation structure
#get site row numbers (from entire grid) for 200 analysis sites

auto<-data.frame(sampling_unit=as.factor(sites), env_dat[sites,c("x","y")])

#define input data include sample level random effects that have a spatial structure
dat_hmsc<-as.HMSCdata(Y=sp_200, X=env_200[,2:9], 
                      Auto=auto, scaleX = FALSE)

#run model
set.seed(6)
mod_hmsc <- hmsc(dat_hmsc, family = "probit", niter = 10000, nburn = 1000, thin = 10)

#check mixing
par(mfrow=c(3,3), mar=c(2,2,2,2))
plot(as.mcmc(mod_hmsc, parameters = "paramX"))
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
#mixing looks OK

#check model fit for each species and overall
#plot fit against prevalence.
plot(apply(mod_hmsc$data$Y, 2, mean),  Rsquared(mod_hmsc, averageSp=FALSE), pch=19)
#looks Ok for most species

#R squared average for all species
Rsquared(mod_hmsc, averageSp=TRUE) #0.32


#generate point prediction for probaility of occurrence for each species
hmsc_pred_data<-as.HMSCdata(X= sim_dat[,2:9], 
                            Auto=data.frame(sampling_unit=as.factor(1:nrow(sim_dat)), env_dat[,c("x","y")]), scaleX=FALSE)

hmsc_pred<-predict(mod_hmsc, hmsc_pred_data)

#out of interest  plot predicted distributions for each species...
levelplot(rasterize(env_dat[,1:2], env,field= as.data.frame(hmsc_pred)), margin=FALSE)
#looks reasonable

#cluster predictions of species' probability of occurrence
hmsc_clust<-hclust(dist(hmsc_pred, method="euclidean"), method="ward.D2")
 
#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= hmsc_clust, dissim=dist(hmsc_pred, method="euclidean"), method="ward.D2")

#suggests 3 grps
hmsc3_clust<-cutree(hmsc_clust,3)

rm(auto, dat_hmsc, hmsc_clust)

## -----------------------------------------------------------
## F) 2 STAGE: GENERALISED DISSIMILARITY MODELS & CLUSTER ----
## ----------------------------------------------------------

### Set up data
## Using gdm package and bbgdm wrapper 

#format data as paired sitewise differences to input into gdm and bbgdm model
env_200a <- as.matrix(data.frame(site=1:200,env_dat[sites,1:2],env_200[,env_vars]))
gdmTab <- formatsitepair(bioData = data.frame(site=1:200,sp_200), bioFormat=1,
                         siteColumn = 'site', predData=env_200a, XColumn = 'x',YColumn = 'y')



#format prediction space environmental data as paired sites in order to generate predictions
#note- creating dummy bio data as gdm predict function requires formatsitepair input but ignores biodata
pred_envTab <- formatsitepair(bioData = data.frame(site=1:dim(sim_dat)[1],A=0, B=1), 
                              bioFormat=1, siteColumn = 'site', 
                              predData=as.matrix(data.frame(site=1:dim(sim_dat)[1],env_dat[,1:2],sim_dat[,env_vars])), 
                              XColumn = 'x',YColumn = 'y')


### naive gdm----
#model
non_bbgdm_mod <- gdm(gdmTab, geo = FALSE)

#predictions
non_bbgdm_pred <- predict_gdm(non_bbgdm_mod, pred_envTab, bbgdm=FALSE)
#format output vector as distance matrix
class(non_bbgdm_pred)='dist'
attr(non_bbgdm_pred,"Size")<-dim(env_dat)[1]
non_bbgdm_pred<- as.dist(non_bbgdm_pred)


# A) cluster dissimilarities directly
non_bbgdm_clust<-hclust(non_bbgdm_pred, method="ward.D2")
grp_stats(min_grp=2, max_grp=10, tree_clust= non_bbgdm_clust, dissim=non_bbgdm_pred, method="ward.D2")
# suggests 2 groups

#set 2 groups
non_bbgdm2_clust<-cutree(non_bbgdm_clust,2)
#set 3 groups
non_bbgdm3_clust<-cutree(non_bbgdm_clust,3)


# B) cluster tranformed environmental variables
non_bbgdm_env_trans <- gdm.transform(non_bbgdm_mod,sim_dat[,env_vars])
non_bbgdm_clust_env_trans<-hclust(dist(non_bbgdm_env_trans, method="euclidean"), method="ward.D2")
grp_stats(min_grp=2, max_grp=10, tree_clust= non_bbgdm_clust_env_trans, dissim=dist(non_bbgdm_env_trans, method="euclidean"), method="ward.D2")
#Suggests 2 clusters

non_bbgdm_clust2_env_trans<-cutree(non_bbgdm_clust_env_trans,2)
non_bbgdm_clust3_env_trans<-cutree(non_bbgdm_clust_env_trans,3)

### Bootstrapped gdm ----
#model
bbgdm_mod <- bbgdm(gdmTab, geo=FALSE, bootstraps=1000, ncores=1)

#predictions
bbgdm_pred <- predict_gdm(bbgdm_mod,pred_envTab, bbgdm=TRUE)
#format output vector as distance matrix
class(bbgdm_pred)='dist'
attr(bbgdm_pred,"Size")<-dim(sim_dat)[1]
bbgdm_pred<-as.dist(bbgdm_pred)


# A) cluster dissimilarities directly
bbgdm_clust<-hclust(bbgdm_pred, method="ward.D2")
grp_stats(min_grp=2, max_grp=10, tree_clust= bbgdm_clust, dissim=bbgdm_pred, method="ward.D2")
#Suggests 2 groups

#set 2 groups
bbgdm2_clust<-cutree(bbgdm_clust,2)
#set 3 groups
bbgdm3_clust<-cutree(bbgdm_clust,3)


# B) cluster tranformed environmental variables
bbgdm_env_trans <- bbgdm.transform (bbgdm_mod,sim_dat[,env_vars])
bbgdm_clust_env_trans<-hclust(dist(non_bbgdm_env_trans, method="euclidean"), method="ward.D2")
grp_stats(min_grp=2, max_grp=10, tree_clust= non_bbgdm_clust_env_trans, dissim=dist(non_bbgdm_env_trans_preds, method="euclidean"), method="ward.D2")
#Suggests 2 clusters

bbgdm2_clust_env_trans<-cutree(bbgdm_clust_env_trans,2)
bbgdm3_clust_env_trans<-cutree(bbgdm_clust_env_trans,3)


mods<-c( "GDM_Dissim2_HC","GDM_Dissim3_HC",
         "GDM_TransEnv2_HC","GDM_TransEnv3_HC",  
         "bbGDM_Dissim2_HC"  , "bbGDM_Dissim3_HC",
         "bbGDM_TransEnv2_HC", "bbGDM_TransEnv3_HC")

test<-as.data.frame(cbind(non_bbgdm2_clust,non_bbgdm3_clust,
  non_bbgdm_clust2_env_trans,non_bbgdm_clust3_env_trans,
  bbgdm2_clust, bbgdm3_clust,
  bbgdm2_clust_env_trans,  bbgdm3_clust_env_trans))
#names(test)<-mods

  
clust2<-stack()
# for (i in 1:length(mods)){#
  for (i in 1:ncol(test)){  
#  hc_rast<-rasterize(env_dat[,1:2], env, field=test[,mods[i]])
    hc_rast<-rasterize(env_dat[,1:2], env, field=test[,i])
    clust2<-stack(clust2, hc_rast)
}
#names(clust2)<-mods
names(clust2)<-names(test)



## ------------------------------------------------------------------
## G) 2 Stage: GRADIENT FOREST THEN CLUSTER TRANSFORMED ENVIRONMENT ----
## ------------------------------------------------------------------

#convert species'presence-absence to factor to treat as classification problem
sp_200_fac<-data.frame(lapply(as.data.frame(sp_200),factor))
names(sp_200_fac)<-species

#run gradient forest model
GF_mod<-gradientForest(cbind(env_200[,2:9], sp_200_fac), 
                       predictor.vars=env_vars, response.vars=species,
                       ntree=500, nbin=201)

# check performance for each species 
plot(GF_mod, plot.type = "P", show.names = T, horizontal = F,
      cex.axis = 1, cex.labels = 0.7, line = 2.5)
#looks OK for most species

#generate predictions (i..e biologically transformed environmental space)
GF_pred<-predict(GF_mod, sim_dat[,2:9])

#cluster biologically transformed environment (NOTE: GF package example clusters PCA loadings)
GF_clust<-hclust(dist(GF_pred, method="euclidean"), method="ward.D2")
plot(GF_clust) 

#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= GF_clust, dissim=dist(GF_pred, method="euclidean"), method="ward.D2")

#sil suggests 2 grps
GF2_clust<-cutree(GF_clust,2)

GF3_clust<-cutree(GF_clust,3)

rm(sp_200_fac, GF_clust)



## ------------------------------------------------
## H) 1 STAGE: SPECIES ARCHETYPE MODELS (SAM)----
## -------------------------------------------------

#format data
samdat <- make_mixture_data(sp_200, env_200)

# formula for archetypes GLMs
form<-as.formula(paste0(paste0('cbind(',paste(species,collapse = ", "),") ~  ", paste(env_vars, collapse= "+"))))

#run multiple groups and choose best using BIC and other diagnostics as per (ref)

test_mods <- list()
for(i in 1:6){
   test_mods[[i]] <- try(species_mix(archetype_formula = form,
                                     species_formula = ~1,
                                     data = samdat,
                                     n_mixtures = i+1,
                                     distribution = 'bernoulli',
                                     standardise = FALSE,
                                     control = species_mix.control(em_refit = 3, em_steps = 5,
                                                                   init_method = 'kmeans')))
 }

# look at change between consecutive archetypes
BICs <- sapply( 1:6, function(x) test_mods[[x]]$BIC)
plot(2:7, BICs, xlab="Archetypes", ylab="BIC") # indicates 3 groups

# look at change between consecutive archetypes
 plot(3:7, diff(BICs), xlab="Archetypes", ylab=expression(paste, Delta , "BIC"))
 abline(h = 0, col="red", lty=2)

# # check that min pi is >1/S (0.05)
test_pi<-sapply( test_mods, function(x) x$pi)
min_test_pi<-lapply( test_pi, function(x) min(x))
plot(2:7, unlist(min_test_pi), xlab="Archetypes", ylab="1/S")
abline(h = 0.05, col="red", lty=2)
#5 archetypes cutoff for 1/S

#choose 3 groups 
sam3_mod<-species_mix(archetype_formula = form,
                      species_formula = ~1,
                      data = samdat,
                      n_mixtures = 3,
                      distribution = 'bernoulli',
                      standardise = FALSE,
                      control = species_mix.control(em_refit = 3, em_steps = 5,
                                                    init_method = 'kmeans'))

#sam3_mod$vcov <- vcov(sam3_mod)
sam3_boot<-species_mix.bootstrap(sam3_mod,nboot=500, type="BayesBoot")

#generate predictions (gives predictions at both species and )
sam3_pred<-predict(sam3_mod, sam3_boot, newdata= sim_dat[,-1])

rm(form, BICs, test_pi, min_test_pi, test_tau)


## --------------------------------------------
## I) 1 STAGE: REGIONS OF COMMON PROFILE (RCP)----
## --------------------------------------------

# Conduct forward selection with linear terms only to select best model (scaled environmental covariates and number of RCPS)
# Forward selection takes some time and was run using the following code on a virtual machine with 10 cores. 
# Results of this process have been imported here


load("Results/RCP_fwd_sel_res.RData")
lin_dat<-as.data.frame(cbind(sp_200, env_200[,2:9]))

# NULL MODEL----
env_vars<-names(lin_dat)[31:38]


# NULL MODEL----
null_mod<-regimix(form.RCP=paste("cbind(",paste(species, collapse=", "),")~1"),
                  nRCP=1, data=lin_dat, dist="Bernoulli")
null_mod$BIC #[1] 7375.015

# STEP 1---
lin_step1<-fwd_step_linear(start_vars="", start_BIC=null_mod$BIC,
                           add_vars=env_vars,
                           species=species,  
                           data=lin_dat, nstarts=500, min.nRCP=2, max.nRCP =7, mc.cores=10)

matplot(t(lin_step1[,2:7]), type="l",lty=1, ylab="BIC", xlab="nRCP (-1)")
legend("topleft", legend=env_vars, lty=1, cex=0.7, col=1:7,bty="n")
rm(null_mod)
#2 RCPs- keep temp BIC=7172.7

# STEP 2---
lin_step2<-fwd_step_linear(start_vars="temp", start_BIC=7172.7,
                           add_vars=env_vars[!env_vars %in% "temp"],
                           species=species,  
                           data=lin_dat, nstarts=500, min.nRCP=2, max.nRCP =7, mc.cores=10)

matplot(t(lin_step2[,2:7]), type="l",lty=1, ylab="BIC", xlab="nRCP (-1)")
legend("topleft", legend=env_vars[!env_vars %in% "temp"], lty=1, cex=0.7, col=1:7,bty="n")
abline(h=7172, col="red", lty=2)
#3 RCPs, keep O2, BIC=7160.7



# STEP 3----
lin_step3<-fwd_step_linear(start_vars=c("temp", "O2"), start_BIC=7160.7,
                           add_vars=env_vars[!env_vars %in% c("temp", "O2")],
                           species=species,  
                           data=lin_dat, nstarts=500, min.nRCP=2, max.nRCP =7, mc.cores=10)

matplot(t(lin_step3[,2:7]), type="l",lty=1, ylab="BIC", xlab="nRCP (-1)")
legend("topleft", legend=env_vars[!env_vars %in% c("temp", "O2")], lty=1, cex=0.7, col=1:7,bty="n")
abline(h=7161, col="red", lty=2)
#sal 3 RCPs, BIC 7160 (keep going?)


# STEP 4----
lin_step4<-fwd_step_linear(start_vars=c("temp", "O2", "sal"), start_BIC=7159.5,
                           add_vars=env_vars[!env_vars %in% c("temp", "O2", "sal")],
                           species=species,  
                           data=lin_dat, nstarts=500, min.nRCP=2, max.nRCP =7, mc.cores=10)

matplot(t(lin_step4[,2:7]), type="l",lty=1, ylab="BIC", xlab="nRCP (-1)")
legend("topleft", legend=env_vars[!env_vars %in% c("temp", "O2", "sal")], lty=1, cex=0.7, col=1:7,bty="n")
abline(h=7160, col="red", lty=2)
# Definitely stop here

### Re-run multifit, find and get best model----
best_form<-as.formula(paste("cbind(",paste(species, collapse=", "),")~ 1+",paste(c("temp", "O2", "sal"), collapse="+")))

best_mod_multi<-RCPmod::regimix.multifit(form.RCP=best_form, data=lin_dat, nRCP=3, 
                                 inits="random2", nstart=500, dist="Bernoulli",mc.cores=1)

# Remove iterations where any RCP contains small number of sites (misfits)
BICs <- sapply( best_mod_multi, function(x) x$BIC)
minPosteriorSites <- sapply( best_mod_multi, function(x) min( colSums( x$postProbs)))
ObviouslyBad <- minPosteriorSites < 2
BICs[ObviouslyBad] <- NA

x<-which.min(BICs)
goodun <- best_mod_multi[[x]]

#get details (tidbits) associated with best model
rcp3_mod<-regimix(form.RCP=best_form, nRCP=3, 
                  data=lin_dat, dist="Bernoulli", inits = unlist(goodun$coef))

#generate bootstrap estimates of model parameters
rcp3_boot<-regiboot(rcp3_mod, nboot = 500,type = "BayesBoot", mc.cores = 10)


## Examine model diagnostics
plot.regimix(rcp3_mod, type="RQR") 
#looks fine

## generate predictions for sim region
rcp3_pred<-predict(rcp3_mod, rcp3_boot, newdata=sim_dat)


#save relevant outputs and tidy up
#fwd_sel_res=list(step1=lin_step1, step2= lin_step2,step3= lin_step3, step3=lin_step4)
#save(fwd_sel_res, rcp3_mod, rcp3_boot, file="RCP_fwd_sel_res.RData")
rm(lin_step1, lin_step2, lin_step3, lin_step_4, best_mod_multi, mod3_multi, BICs, minPosteriorSites, 
  ObviouslyBad, x,goodun, best_form)


##############################################
## COMPILE MODELS, CLUSTERS, PREDICTIONS ----
##############################################
save(rf_sp_mods, bio2_rf, bio3_rf, bbgdm_mod, nonbb_gdm_mod,
     GF_mod, MNet_mod, mod_hmsc, rcp3_mod, rcp3_boot,
     sam3_mod,    
     file="Results/models.RData")

save(rf_sp_pred, bbgdm_mod, nonbb_gdm_pred, GF_pred,
     MNet_pred, hmsc_pred, rcp3_pred,
     sam3_pred, 
     file="Results/preds.RData")

#compile hard clusters from 2 stage methods
hard_cluster2<-data.frame(env_dat[,c("x","y")], 
                          env=env11_clust, 
                          Sp_RF=ssdm2_clust, 
                          bbGDM=bbgdm2_clust, nonbbGDM=non_bbgdm2_clust,
                          GF=GF2_clust, MNet=MNet2_clust,MNet5_clust,
                          HMSC=hmsc3_clust)

hard_cluster3<-data.frame(env_dat[,c("x","y")], 
                          env=env3_clust,
                          Sp_RF=ssdm3_clust, 
                          bbGDM=bbgdm3_clust, nonbbGDM=non_bbgdm3_clust,
                          GF=GF3_clust, MNet=MNet3_clust,
                          HMSC=hmsc3_clust)
save(hard_cluster2,hard_cluster3, 
     bio2_clust, bio3_clust,
     bio2_rf_pred, bio3_rf_pred,
     rcp3_pred, 
     sam3_pred,
     file="Results/pred_clusters.RData")

#rm(list=ls(pattern="2_clust"))
#rm(list=ls(pattern="3_clust"))
