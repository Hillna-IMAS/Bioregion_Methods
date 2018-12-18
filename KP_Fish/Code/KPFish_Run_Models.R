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
# G) predict biologically transformed environment from GF, cluster predictions
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
library(bbgdm)          #Bayesian Bootstrap Generalised Dissimilarity Models (and naive GDM)
library(RCPmod)         #Regions of Common Profile Models v2.188
#devtools::install_github('skiptoniam/ecomix@dev')  #Species Archetype Models (SAMs)
library(ecomix) 
library(HMSC)           #Hierarchical Modelling of Species Communities
library(gradientForest) #Gradient Forests (an extension of Random Forests)
#devtools::install_github('davharris/mistnet')  #Mistnet stochastic, multispecies neural network
library(mistnet)
library(RColorBrewer)   #colour palettes
library(rasterVis)      #plotting rasters

setwd("C:\\Users\\hillna\\UTAS_work\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\KP_Fish\\")
source("Code/Additional_Funcs.R")

#Load required files
dat<- readRDS("Data/KP_Fish_Env.RDS")
# contains both presence-absence of 20 fish species at 524 sites as well as matching environmental variables 

load("Data/Prediction_space.Rda")
# env_raster is raster of all environmental variables for entire Kerguelen Plateau
# pred_sp is dataframe version of above

# Environmental variables to test
env_vars<-c("log_slope", "Av_depth" ,
            "sea.floor.temperature", 
            "log_current"   ,   "no3_bot_mean" ,            
            "T_mean"  , "chla_yearly_sd"  ,  "ssha_variability")

# Species to model
species=c("Antimora.rostrata" ,"Bathydraco.antarcticus", "Bathyraja.eatonii",                
          "Bathyraja.irrasa"  ,  "Bathyraja.murrayi",  "Champsocephalus.gunnari" , "Dissostichus.eleginoides" ,        
          "Etmopterus.viator"  ,  "Gobionotothen.acuta" , "Lepidonotothen.mizops"   ,  "Lepidonotothen.squamifrons",       
          "Lycodapus.antarcticus" ,"Macrourus.spp", "Mancopsetta.maculata"  , "Muraenolepis.spp",                
          "Notothenia.rossii"   , "Paradiplospinus.gracilis" ,  "Paraliparis.spp." ,"Zanclorhynchus.spinifer" ,         
          "Channichthys.rhinoceratus.velifer")

#isolate trawl x species matrix
fish_dat<-dat[,species]


#######################################################
## 1) Run models and any diagnostics and 
##  a) determine optimal number of groups
##  b) set number of groups to 4 or 5
##  c) generate predictions of clusters for Kerguelen Plateau region
#######################################################

## -----------------------------------
## A) ENVIRONMENT ONLY CLUSTERING ----
## -----------------------------------

# Cluster euclidean distance matrix using Ward's SS
env_clust<-hclust(dist(na.omit(pred_sp), method="euclidean"), method="ward.D2")


#determine number of groups
grp_stats(min_grp=2, max_grp=15, tree_clust= env_clust, dissim=dist(na.omit(pred_sp), method="euclidean"), method="ward.D2")

#sil width, bootstrap suggests 4, CH suggests 15.

#cut at 4 
env4_clust.2<-cutree(env_clust, k=4)

rm(env_clust)



## ------------------------------------------------------------------------------------------------
## B) 2 STAGE: CLUSTER BIOLOGICAL DATA; PREDICT GROUPS USING RANDOM FORESTS & ENVIRONMENTAL DATA----
## ------------------------------------------------------------------------------------------------

#create Jaccard distance matrix for presence-absence data 
# using sites that have more than one species. 5
jac_dist<-vegdist(fish_dat[rowSums(fish_dat)>1,], method="jaccard",binary=TRUE)

# Cluster distance matrix using Ward's distance
dat_clust<-hclust(jac_dist, method="ward.D2")

#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= dat_clust, dissim=jac_dist, method="ward.D2")
#All suggests 2 groups

#set 2 groups
bio2_clust<-cutree(dat_clust,2)
#Set 4 groups
bio4_clust<-cutree(dat_clust,4)


## Use Random Forest to match environment and predict clusters across sim region
# 2 groups
bio2_rf<-randomForest(y= as.factor(bio2_clust), x=dat[rowSums(dat[,species])>1,env_vars], importance=TRUE, corr.threshold = 0.65)
bio2_rf$confusion 

bio2_rf_pred<- predict(bio2_rf, na.omit(pred_sp), type=c("prob"))

# 4 groups
bio4_rf<-randomForest(y= as.factor(bio4_clust), x=dat[rowSums(dat[,species])>1,env_vars], importance=TRUE)
bio4_rf$confusion #more confusion amongst groups
bio4_rf_pred<- predict(bio4_rf, na.omit(pred_sp), type=c("prob"))

rm(jac_dist, dat_clust)


## ----------------------------------------------------------------
## C) 2 STAGE: RANDOM FOREST PREDICT EACH SPECIES THEN CLUSTER ----
## ----------------------------------------------------------------

# Loop through running RF for each species (as classification problem)
rf_sp_mods<-list()
for(i in 1:length(species)){
  rf_sp_mods[[i]]<-randomForest(y=as.factor(dat[,species[i]]), x=dat[, env_vars], importance=TRUE, corr.threshold = 0.65)
}
names(rf_sp_mods)<-species

# Loop through predictions for each species, generating probabilty of ocurrence output
rf_sp_pred<-data.frame(matrix(nrow=dim(na.omit(pred_sp))[1], ncol=length(species), 
                              dimnames=list(NULL, species)))

for(i in 1:length(species)){                       
  rf_sp_pred[,(i)]<-predict(rf_sp_mods[[i]], na.omit(pred_sp), type="prob")[,2]
}

# Cluster predictions

ssdm_clust<-hclust(dist(rf_sp_pred, method="euclidean"), method="ward.D2")

#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= ssdm_clust, dissim=dist(rf_sp_pred[,species], method="euclidean"), method="ward.D2")
#sil suggests 2 groups

# set 2 groups
ssdm2_clust<-cutree(ssdm_clust, k=2)
# set 4 groups
ssdm4_clust<-cutree(ssdm_clust, k=4)

rm(ssdm_clust)


## ----------------------------------------------------
## D) 2 STAGE: MISTNET & CLUSTER SPECIES' PREDICTIONS ----
## ----------------------------------------------------
# Optimising mistnet parameters with 1 hidden layer using slightly modified code from 
# https://github.com/davharris/mistnet/extras/BBS-analysis/mistnet_cross-validation.R

set.seed(6)

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
fold.ids<-sample( x=1:5, size= nrow(dat),replace=TRUE)

out = list()

for(i in 1:dim(hyperparams)[[1]]){
  cat(paste0("Starting iteration ", i, "\n"))
  for(fold.id in 1:max(fold.ids)){
    cat(paste0(" Starting fold ", fold.id, "\n  "))
    in.val=fold.ids==fold.id   #added this line
    in.train = fold.ids != fold.id
    MNet_mod = fit(
      x=as.matrix(dat[in.train,env_vars]),
      y=as.matrix(dat[in.train,species]),
      hyperparams = hyperparams,
      i = i
    )
    
    cat("\n evaluating")
    
    #calculate log likelihood of predictions (replaced with my code)
    log.lik = sum(
      dbinom(
        as.matrix(dat[in.val,species]),
        size = 1, 
        prob = predict(MNet_mod, as.matrix(dat[in.val,env_vars]),n.importance.samples=dim(dat[in.val,])[[1]]),
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
#indicates 4 latent variables and hidden layer with 12 neurons

MNet_mod = fit(
  x = as.matrix(dat[,env_vars]),
  y = as.matrix(dat[,species]), 
  hyperparams, 
  which.max(logliks[,5])
)

save(MNet_mod, file = "mistnet.model.RData")

#generate predictions for each species
MNet_pred<-predict(MNet_mod, as.matrix(na.omit(pred_sp[,env_vars])), n.importance.samples=25)
MNet_pred<-apply(MNet_pred,1:2, mean)

#cluster predictions
MNet_clust<-hclust(dist(MNet_pred, method="euclidean"), method="ward.D2")

#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= MNet_clust, dissim=dist(MNet_pred, method="euclidean"), method="ward.D2")

#sil suggests 4 grps, ch 10 groups, boot 4 groups!
MNet4_clust<-cutree(MNet_clust,4)

#choose 5 groups
MNet5_clust<-cutree(MNet_clust,5)

rm(MNet_clust)


## -------------------------------------------------------------------------------------------------
## E) 2 STAGE: HEIRARCHICAL MODELLING OF SPECIES COMMUNITIES (HMSC) & CLUSTER SPECIES' PREDICTIONS ----
## -------------------------------------------------------------------------------------------------

# Create orthogonal quadratics of environmental variables for HMSC, SAM and RCP
# Apply same transofrmations to prediction space
quad_dat<-poly_data(env_vars, degree=c(rep(2, length(env_vars))), 
                     id_vars = c("sample_ID", "Year", "Survey", "Lat", "Long"), species=species, data=dat)

env2_vars<-paste0(rep(env_vars, each=2), 1:2)
hmsc_dat<-quad_dat$rcp_data[,env2_vars]

trans_pred_space<-poly_pred_space(pred_space=na.omit(pred_sp), poly_output= quad_dat$poly_output, 
                                  vars=env_vars, sampling_levels=levels(quad_dat$rcp_data$Survey), sampling_factor="Survey")

# Fromat HMSC model inputs
#Y data needs to be 'numeric'
Y<-as.matrix(sapply(dat[,species], as.numeric))  
dimnames(Y)[[1]]<-rownames(dat)

#Auto = spatial coordinates to calculate spatial autocorrelation structure 
auto<-data.frame(sampling_unit=dat$sample_ID, dat[,c("Long","Lat")])
#one duplicated co-ordinate which causes problems for distance function. Remove row 86
auto<-auto[!duplicated(auto[,2:3]),]


#define input data include sample level random effects that have a spatial structure
dat_hmsc<-as.HMSCdata(Y=Y[rownames(auto),], X=hmsc_dat[rownames(auto),], 
                      Auto=droplevels(auto), scaleX = FALSE)

#run model
set.seed(6)
mod_hmsc <- hmsc(dat_hmsc, family = "probit", niter = 10000, nburn = 1000, thin = 10)

#check mixing
par(mfrow=c(3,3), mar=c(2,2,2,2))
plot(as.mcmc(mod_hmsc, parameters = "paramX"))
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
#mixing looks OK

#check model fit for each species and overall---
#plot fit against prevalence.
plot(apply(mod_hmsc$data$Y, 2, mean),  Rsquared(mod_hmsc, averageSp=FALSE), pch=19)
#looks Ok for most species

#R squared average for all species
Rsquared(mod_hmsc, averageSp=TRUE) #0.29

#look at AUC for each species
library(PresenceAbsence)
mod_fit<-predict(mod_hmsc)

sp_auc<-data.frame(Prev=round(apply(Y, 2, mean), 2),
                   AUC=NA, SD=NA)

for(i in 1:length(species)){
  temp<-auc(DATA=cbind(rownames(auto),Y[rownames(auto),species[i]], mod_fit[, species[i]]), na.rm=TRUE)
  sp_auc$AUC[i]<-round(temp$AUC,2)
  sp_auc$SD[i]<-round(temp$AUC.sd,2)
}
sp_auc[ order(sp_auc$AUC, decreasing=TRUE),]
#AUC looks reasonable for most species

#generate point prediction for probaility of occurrence for each species----
hmsc_pred_data<-as.HMSCdata(X= trans_pred_space[,env2_vars], 
                            Auto=data.frame(sampling_unit=as.factor(1:nrow(trans_pred_space)), na.omit(pred_sp)[,c("x","y")]), 
                            scaleX=FALSE)

hmsc_pred<-predict(mod_hmsc , hmsc_pred_data)

#cluster predictions of species' probability of occurrence----
hmsc_clust<-hclust(dist(hmsc_pred, method="euclidean"), method="ward.D2")

#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= hmsc_clust, dissim=dist(hmsc_pred, method="euclidean"), method="ward.D2")

#sil and boot suggests 2 grps
hmsc2_clust<-cutree(hmsc_clust,2)

#set 4 groups
hmsc4_clust<-cutree(hmsc_clust,4)

rm(auto, dat_hmsc, hmsc_clust)


## -----------------------------------------------------------
## F) 2 STAGE: GENERALISED DISSIMILARITY MODELS & CLUSTER (GDM_HC and bbGDM_HC) ----
## ----------------------------------------------------------

# formula for GDM and bbGDM
form<- paste("~1 +", paste(env_vars, collapse= "+"))


#Run bbGDM using sites with more than 1 species for calculating disimilarity matrix (same as hierarchical clustering on species' data)
gdm_mod<-bbgdm(form, dat[rowSums(dat[,species])>1, species], dat[rowSums(dat[,species])>1, env_vars], family="binomial", link="logit",
               dism_metric="number_non_shared", nboot= 100, geo=FALSE)  

rm(form)

save(gdm_mod, file=paste0(path,"Results/bbGDM/bbGM_mod.RData"))

#check diagnostics
resids<-diagnostics(gdm_mod)

pdf(file=paste0(path,"Results/bbGDM/bbGDM_diagnostics.pdf"), height=6, width=6)
par(mfrow=c(2,2))
plot(resids)
dev.off()
#model looks appropriate, but not much match between predicted and observed dissmilarities!

## predict turnover across region----
# Predicts pairwise dissimilarity between cells which are then clustered.
bbgdm_pred<-pred_gdm_dissim(gdm_mod, as.matrix(na.omit(pred_sp[,3:10])))

#cluster predictions
bbgdm_clust<-hclust(bbgdm_pred, method="ward.D2")
grp_stats(min_grp=2, max_grp=10, tree_clust= bbgdm_clust, dissim=dist(bbgdm_pred, method="euclidean"), method="ward.D2")


#sil width and boot suggests 2 grps
bbgdm2_clust<-cutree(bbgdm_clust,2)

#set at 4 groups
bbgdm4_clust<-cutree(bbgdm_clust,4)

##get naive GDM, predict and cluster----
nonbb_gdm_mod<-gdm_mod$starting_gdm
nonbb_gdm_pred<-pred_gdm_dissim(bbgdm_mod = gdm_mod, naive=TRUE, env_data = as.matrix(na.omit(pred_sp[,3:10])))
nonbb_gdm_clust<-hclust(nonbb_gdm_pred, method="ward.D2")

grp_stats(min_grp=2, max_grp=10, tree_clust= nonbb_gdm_clust, dissim=nonbb_gdm_pred, method="ward.D2")

#sil width and ch suggests 2 grps
nonbb_gdm2_clust<-cutree(nonbb_gdm_clust,2)

#set at 4 groups
nonbb_gdm4_clust<-cutree(nonbb_gdm_clust,4)



## ------------------------------------------------------------------
## G) 2 Stage: GRADIENT FOREST THEN CLUSTER TRANSFORMED ENVIRONMENT ----
## ------------------------------------------------------------------

#convert species'presence-absence to factor to treat as classification problem
sp_fac<-data.frame(lapply(as.data.frame(dat[,species]),factor))

#run gradient forest model
GF_mod<-gradientForest(cbind(dat[,env_vars], sp_fac), 
                       predictor.vars=env_vars, response.vars=species,
                       ntree=500, nbin=201)

# check performance for each species 
plot(GF_mod, plot.type = "P", show.names = T, horizontal = F,
     cex.axis = 1, cex.labels = 0.7, line = 2.5)
#full spectrum of model fits for each species

#generate predictions (i..e biologically transformed environmental space)
GF_pred<-predict(GF_mod, na.omit(pred_sp[,3:10]))

#cluster biologically transformed environment (NOTE: GF package example clusters PCA loadings)
GF_clust<-hclust(dist(GF_pred, method="euclidean"), method="ward.D2")

#determine number of groups
grp_stats(min_grp=2, max_grp=10, tree_clust= GF_clust, dissim=dist(GF_pred, method="euclidean"), method="ward.D2")

#sil and boot suggests 2 grps
GF2_clust<-cutree(GF_clust,2)

#set at 4 groups
GF4_clust<-cutree(GF_clust,4)

rm(sp_fac, GF_clust)

## ------------------------------------------------
## H) 1 STAGE: SPECIES ARCHETYPE MODELS (SAM)----
## -------------------------------------------------

# Format Model Inputs
sam_dat <- make_mixture_data(quad_dat$rcp_data[,species], quad_dat$rcp_data[, env2_vars])
form <- as.formula(paste0(paste0('cbind(',paste(species,collapse = ", "),") ~ 1 + ", paste(env2_vars, collapse= "+"))))

#run multiple groups and choose beset number of groups
test_mods <- list()
 for(i in 1:5){
   test_mods[[i]] <- try(species_mix(archetype_formula = form,
                                     species_formula = ~1,
                                     data = sam_dat,
                                     n_mixtures = i+1,
                                     distribution = 'bernoulli',
                                     standardise = FALSE,
                                     control = species_mix.control(em_prefit=TRUE, em_refit = 8, em_steps = 5,
                                                                   init_method = 'kmeans')))
 }


# look at change between consecutive archetypes
BICs <- sapply( 1:5, function(x) test_mods[[x]]$BIC)
plot(2:6, BICs, xlab="Archetypes", ylab="BIC") # indicates 4 groups

# look at change between consecutive archetypes
plot(3:6, diff(BICs), xlab="Archetypes", ylab=expression(paste, Delta , "BIC"))
abline(h = 0, col="red", lty=2)

# # check that min pi is >1/S (0.05)
test_pi<-sapply( test_mods, function(x) x$pi)
min_test_pi<-lapply( test_pi, function(x) min(x))
plot(2:6, unlist(min_test_pi), xlab="Archetypes", ylab="1/S")
abline(h = 0.05, col="red", lty=2)
#4 archetypes cutoff for 1/S

#Choose 4 archetypes,run final model
sam4_mod<-species_mix(archetype_formula = form,
                      species_formula = ~1,
                      data = sam_dat,
                      n_mixtures = 4,
                      distribution = 'bernoulli',
                      standardise = FALSE,
                      control = species_mix.control(em_prefit=TRUE, em_refit = 8, em_steps = 5,
                                                    init_method = 'kmeans'))

sam4_mod$vcov <- vcov(sam4_mod)

#generate predictions (gives predictions at both species and group level)
sam4_pred<-predict(sam4_mod, as.data.frame(trans_pred_space[env2_vars]))


## --------------------------------------------
## I) 1 STAGE: REGIONS OF COMMON PROFILE (RCP)----
## --------------------------------------------

# Conduct forward selection with orthogonal quadratic terms (both in or both out) to select best model
# Including survey (combination of year and season) as sampling factor
# Forward selection takes some time and was run using the following code on a virtual machine with 10 cores. 
# Results of this process have been imported here

rcp_dat<-quad_dat$rcp_data
rcp_vars<-env2_vars

##NULL MODEL---
null_mod<-regimix(form.RCP=as.formula(paste("cbind(",paste(species, collapse=", "),")~1")),form.spp= ~Survey,
                  nRCP=1, data=rcp_dat, dist="Bernoulli")
null_mod$BIC
#[1] 9251.175


#STEP 1----
ptm <- proc.time()
fwd_step1<-fwd_step(start_vars="", start_BIC=null_mod$BIC,
                    add_vars=rcp_vars,
                    species=species,  form.spp= ~Survey, 
                    data=rcp_dat, nstarts=500, min.nRCP=2, max.nRCP = 7, mc.cores=8)
proc.time() - ptm 

matplot(t(fwd_step1[,2:6]), type="l",lty=1, ylab="BIC", xlab="nRCP (-1)")
legend("topleft", legend=fwd_step1$Var, lty=1, cex=0.7, col=1:7,bty="n")
#keep depth, 4 RCPS, BIC=7012.4
rm(null_mod)

##STEP 2----
step2_vars<-c("Av_depth1","Av_depth2")

ptm <- proc.time()
fwd_step2<-fwd_step(start_vars=step2_vars, start_BIC=7012.4,
                    add_vars= rcp_vars[!rcp_vars %in% step2_vars],
                    species=species,  form.spp= ~Survey, 
                    data=rcp_dat, nstarts=500, min.nRCP=2, max.nRCP = 8, mc.cores=8)
proc.time() - ptm #7.2 hours

matplot(t(fwd_step2[,2:8]), type="l",lty=1, ylab="BIC", xlab="nRCP (-1)")
legend("topleft", legend=fwd_step2$Var, lty=1, cex=0.7, col=1:7,bty="n")
abline(h=7012.4, col="red", lty=2)
#keep T_mean, 5 RCPS, BIC= 6901.2

## STEP 3----
step3_vars<-c("Av_depth1","Av_depth2", "T_mean1", "T_mean2")

ptm <- proc.time()
fwd_step3<-fwd_step(start_vars=step3_vars, start_BIC=6901.2,
                    add_vars= rcp_vars[!rcp_vars %in% step3_vars],
                    species=species,  form.spp= ~Survey, 
                    data=rcp_dat, nstarts=500, min.nRCP=2, max.nRCP = 8, mc.cores=8)
proc.time() - ptm 

matplot(t(fwd_step3[,2:8]), type="l",lty=1, ylab="BIC", xlab="nRCP (-1)")
legend("topleft", legend=fwd_step3$Var, lty=1, cex=0.7, col=1:7,bty="n")
abline(h=6901.2, col="red", lty=2)

#STOP HERE. Best model includes depth and temp and 5 RCPs

## Rerun multifit, find and get best model----
best_form<-as.formula(paste("cbind(",paste(species, collapse=", "),")~",paste(step3_vars, collapse="+")))

best_mod_multi<-regimix.multifit(form.RCP=best_form, form.spp= ~ Survey, data=rcp_dat, nRCP=5, 
                                 inits="random2", nstart=500, dist="Bernoulli",mc.cores=8)

# Remove iterations where any RCP contains small number of sites (misfits)
BICs <- sapply( best_mod_multi, function(x) x$BIC)
minPosteriorSites <- sapply( best_mod_multi, function(x) min( colSums( x$postProbs)))
ObviouslyBad <- minPosteriorSites < 2
BICs[ObviouslyBad] <- NA

x<-which.min(BICs)
goodun <- best_mod_multi[[x]]

#get details (tidbits) associated with best model
rcp5_mod<-regimix(form.RCP=best_form,form.spp= ~Survey, nRCP=5, 
                  data=rcp_dat, dist="Bernoulli", inits = unlist(goodun$coef))

#generate bootstrap estimates of model parameters
rcp5_boot<- regiboot(rcp5_mod, nboot = 500,type = "BayesBoot", mc.cores = 8)

#save relevant outputs and tidy up
fwd_sel_res=list(step1=fwd_step1, step2= fwd_step2,step3= fwd_step3)
save(fwd_sel_res, rcp5_mod, rcp5_boot, file="Results/RCP/fwd_sel_res.Rda")


### DIAGNOSTICS- RESDIUAL PLOTS
plot.regimix(best_mod, type="RQR")
# Looks fine

# generate RCP predictions (mean and 95% CIs)
rcp5_spat_preds<-predict.regimix(rcp5_mod, rcp5_boot, newdata=trans_pred_space, mc.cores=4)

## Run RCP with four groups
best4_mod_multi<-regimix.multifit(form.RCP=best_form, form.spp= ~ Survey, data=rcp_dat, nRCP=4, 
                                 inits="random2", nstart=500, dist="Bernoulli",mc.cores=1)

# Remove iterations where any RCP contains small number of sites (misfits)
BICs <- sapply( best4_mod_multi, function(x) x$BIC)
minPosteriorSites <- sapply( best4_mod_multi, function(x) min( colSums( x$postProbs)))
ObviouslyBad <- minPosteriorSites < 2
BICs[ObviouslyBad] <- NA

x<-which.min(BICs)
goodun <- best4_mod_multi[[x]]

#get details (tidbits) associated with best model
rcp4_mod<-regimix(form.RCP=best_form,form.spp= ~Survey, nRCP=4, 
                  data=rcp_dat, dist="Bernoulli", inits = unlist(goodun$coef))

rcp4_boot<- regiboot(rcp4_mod, nboot = 500,type = "BayesBoot", mc.cores = 8)

rcp4_spat_preds<-predict.regimix(rcp4_mod, rcp4_boot, newdata=trans_pred_space, mc.cores=4)

rm(best_mod_multi, BICs, minPosteriorSites, 
   ObviouslyBad, x,goodun, step3_vars)

################################################################################
### Compile models, predictions, etc for further analysis and plotting----
################################################################################

save(rf_sp_mods, bio2_rf, bio4_rf, gdm_mod, nonbb_gdm_mod,
     GF_mod, MNet_mod, mod_hmsc, sam4_mod, rcp5_mod, rcp5_boot, 
     rcp4_mod, rcp4_boot,
     file="Results/models.RData")

save(rf_sp_pred,bio2_rf_pred, bio4_rf_pred, bbgdm_pred, nonbb_gdm_pred, GF_pred,
     MNet_pred, hmsc_pred, rcp4_spat_preds,rcp5_spat_preds,
     sam4_pred, 
     file="Results/preds.RData")

#compile hard clusters from 2 stage methods
hard_cluster_opt<-data.frame(na.omit(pred_sp)[,c("x","y")], 
                             Env_Only=env4_clust, 
                             SpRF_HC=ssdm2_clust, 
                             bbGDM_HC=bbgdm2_clust, GDM_HC=nonbb_gdm2_clust,
                             GF_HC=GF2_clust, MNet_HC=MNet4_clust,
                             HMSC_HC=hmsc2_clust)

hard_cluster4<-data.frame(na.omit(pred_sp)[,c("x","y")], 
                          Env_Only=env4_clust,
                          SpRF_HC=ssdm4_clust, 
                          bbGDM_HC=bbgdm4_clust, GDM_HC=nonbb_gdm4_clust,
                          GF_HC=GF4_clust, MNet_HC=MNet4_clust,
                          HMSC_HC=hmsc4_clust)

save(hard_cluster_opt,hard_cluster4, 
     bio2_clust, bio4_clust,
     bio2_rf_pred, bio4_rf_pred,
     rcp4_spat_preds,rcp5_spat_preds, 
     sam4_pred,
     file="Results/pred_clusters.RData")
