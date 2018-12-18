#############################################################################################
## Compare community modelling methods for bioregionalisation of simulated species dataset
## March 2018. N.Hill with input from S. Woolley
#############################################################################################

# Modelling process
# 1) Run models,diagnostics, determine number of groups or set number of groups to 3
# 2) Plot Predicted Distribution of groups across simulation region
# 3) Describe contents of groups
# 4) Describe environment of groups
# 5) PREDICTOR IMPORTANCE

#######################
## Set up---
#######################
# Get required libraries
#library(tidyr)         #data manipulation
library(ggplot2)       #plotting
#library(raster)        #spatial data
#library(RColorBrewer)  #colour palettes
#library(rasterVis)     #plotting rasters
library(RCPmod)         #Regions of Common Profile
library(bbgdm)          #Bayesian Bootstrap Generalised Dissimilarity Models (and naive GDM)
library(HMSC)           #Heirarchical Joint species Distribution modelling
library(mistnet)        #Multivariate neural networks
library(gradientForest) #Gradient Forests (an extension of Random Forests)
library(ecomix)         #Species Archetype Models (SAMs)
library(randomForest)   #Random Forests

setwd("C:\\Users\\hillna\\UTAS_work\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\Simulation\\")
source("Simulation_Additional_Funcs.R")

#Load required files
#files in "simulate_communities" folder on dropbox

load("Sim_setup/Many_covars_sim.RData") 
#sim_data= simulated species probabilities, occurrences, species' groups for entire region
#sp_200 = matrix of occurrences of 30 species at 200 sites to use as species dataset for analysis
#env_dat= matrix of corresponding environmental conditions at 200 site to use as environmental dataset for analysis


#load("sim_env.RData") 
# env= raster brick of environmental data for entire region
# env_dat= raster brick of environmental data for entire region converted to matrix

load("Results/models.RData")
#load("pred_clusters.RData")

species<-paste0("Sp", 1:30)

env_vars<-dimnames(sim_dat)[[2]][2:9]


############################################################################################################
## Take note of label switching from plot_pred_clusters
## Tabulate or get directly from models contents of groups- Limit to case where groups =3
## Generate dotchart of average and se of species' occurence in each group
## Generate partial plots or equivalent of the group response to the environment for GDM, GF, SAM & RCP
############################################################################################################


# A) Cluster environment only- can't do
# 2 Stage Models- cluster then predict:
# B) cluster biological data, predict clusters with random forests- from RF model
# 2 Stage Models- predict then cluster
# C) predict species using random forests, then cluster predictions- average of ind. species RFs
# D) predict species using Mistnet, then cluster predictions- calculated using network connection weightings and averaged over all species
# E) predict species using HMSC, then cluster predictions- based on variance decomposition
# F) predict dissimilarities using GDM, cluster predicted dissimilarities- Wald's test
# G) predict biologically transformed environment (GF)- from cumulaitve splits across all species
# 1 Stage model-based: cluster and predict
# H) Species Archetype Models- BIC reduction in model selection
# I) Regions of Common Profile- BIC reduction in model selection



## B) Cluster biology then RF predict----
tiff(file="Results/Plots/BioHC_RF_import.tiff",
     width=7, height=10, units="cm", compression = "lzw", res=1000)
varImpPlot(bio3_rf, type=1,main="", cex=0.9, cex.main=0.8)
dev.off()

## C) Sp_RF----

#calculate average importance over all species
rf_sp_varimp<-array(dim=c(length(env_vars),2,length(species)), dimnames=list(env_vars, c("MeanDecreaseAccuracy", "MeanDecreaseGini"), species))
for( i in 1:length(species)){
  rf_sp_varimp[, , i]<- rf_sp_mods[[i]]$importance[,3:4]
}

tiff(file="Results/Plots/spRF_Overall_imp.tiff", width=7, height=10, units="cm",
       compression = "lzw", res=1000)
imp<-apply(rf_sp_varimp, c(1,2), mean)
# modified from VarImpPlot function in RandomForest package
#op <- par(mfrow = c(1, 1), mar = c(4, 5, 4, 1), mgp = c(2,0.8, 0), 
#          oma = c(0, 0, 2, 0), no.readonly = TRUE)
#on.exit(par(op))

#for(i in 1:2){
for (i in 1:1) {
  ord <- rev(order(imp[, i], decreasing = TRUE)[1:length(env_vars)])
  #xmin <- if (colnames(imp)[i] %in% c("IncNodePurity", "MeanDecreaseGini")) 
  xmin <- if (colnames(imp)[i] %in% c("IncNodePurity")) 
    0
  dotchart(imp[ord, i], xlab = colnames(imp)[i], ylab = "", 
           main = "", cex=0.8)
}
mtext(outer = TRUE, side = 1, text = "Average Importance", cex = 0.5)
dev.off()

## D) HMSC----
# calculated for each environmental variable forspecies, then averaged across all species and plotted
hmsc_sp_var<-variPart(mod_hmsc, groupX=dimnames(mod_hmsc$data$X)[[2]])
hmsc_var<-data.frame(variable=c(env_vars, "latent"), mean=apply(hmsc_sp_var[,-1], 2, mean), se=apply(hmsc_sp_var[,-1], 2, function(x) sd(x)/sqrt(30)))
hmsc_var<-hmsc_var[order(hmsc_var$mean, decreasing=FALSE),]

tiff(file="Results/Plots/hmsc_imp.tiff", width=6, height=9, units="cm",
     compression = "lzw", res=1000)
cent<-barplot(hmsc_var$mean, beside=TRUE, space=0, names.arg=hmsc_var$variable,horiz=TRUE, las=1,xlim=c(0,0.25),
              cex.axis=0.7, cex.names=0.7)
arrows(y0=cent, x0=hmsc_var$mean + hmsc_var$se, x1=hmsc_var$mean - hmsc_var$se,  angle=90, code=3, length=0.05)
dev.off()

##E) Misnet ----
MNet_imp<-olden_mistnet(MNet_mod)
MNet_overall_imp<-MNet_imp$`Overall relative importance`[order(MNet_imp$`Overall relative importance`$mean, decreasing=FALSE),]/100
#MNet_overall_imp<-MNet_overall_imp[3:10,]


tiff(file="Results/Plots/MNet_imp.tiff", width=7, height=10, units="cm",
     compression = "lzw", res=1000)
cent<-barplot(MNet_overall_imp$mean, beside=TRUE, space=0, names.arg=row.names(MNet_overall_imp),horiz=TRUE, las=1,
              cex.names=0.7, cex.axis=0.7, xlim=c(0,0.35))
arrows(y0=cent, x0=MNet_overall_imp$mean + MNet_overall_imp$sd/sqrt(30), x1=MNet_overall_imp$mean - MNet_overall_imp$sd/sqrt(30),  angle=90, code=3, length=0.05)
dev.off()


tiff(file="Results/Plots/HMSC_MNet_imp.tiff", width=12, height=9, units="cm",
     compression = "lzw", res=1000)
par(mfrow=c(1,2))
cent<-barplot(hmsc_var$mean, beside=TRUE, space=0, names.arg=hmsc_var$variable,horiz=TRUE, las=1,xlim=c(0,0.25),
              cex.axis=0.7, cex.names=0.7)
arrows(y0=cent, x0=hmsc_var$mean + hmsc_var$se, x1=hmsc_var$mean - hmsc_var$se,  angle=90, code=3, length=0.05)

cent2<-barplot(MNet_overall_imp$mean, beside=TRUE, space=0, names.arg=row.names(MNet_overall_imp),horiz=TRUE, las=1,
              cex.names=0.7, cex.axis=0.7, xlim=c(0,0.35))
arrows(y0=cent2, x0=MNet_overall_imp$mean + MNet_overall_imp$sd/sqrt(30), x1=MNet_overall_imp$mean - MNet_overall_imp$sd/sqrt(30),  angle=90, code=3, length=0.05)
dev.off()

## F) bbGDM----
write.csv(round(bbgdm.wald.test(bbgdm_mod, gdm=TRUE),3), 
          file="Results/gdm_VarImp.csv")


## G) Gradient Forests----
tiff(file="Results/Plots/GF_HC_import.tiff",width=14, height=10, units="cm",
     compression = "lzw", res=1000)
plot(GF_mod, plot.type="Overall.Importance",plot.args=list(cex.axis=0.9, cex.names=0.9))
dev.off()

## H) SAM ----
## Run variable selection wiht 3 archetypes

#set up results
BICs<-as.data.frame(setNames(replicate((max.nRCP-min.nRCP+3),numeric(0), simplify = F), c("Var",paste0("RCP", rep(min.nRCP:max.nRCP)),"Start_BIC")))

#loop through adding variables  
sams_fwd_step<-function(start_vars, add_vars, n_mixtures, data){
  test_mods<-list()
  
  for(j in 1:length(add_vars)){
  add<-add_vars[j]
  temp_form<-as.formula(paste("cbind(",paste(species, collapse=", "),")~ 1+",
                              paste0(paste(start_vars, collapse="+"), "+", paste(add, collapse = "+"))))

  #run SAMs
  test_mods[[j]] <- try(species_mix(archetype_formula = temp_form,
                                    species_formula = ~1,
                                    data = data,
                                    n_mixtures = n_mixtures,
                                    distribution = 'bernoulli',
                                    standardise = FALSE,
                                    control = species_mix.control(em_refit = 3, em_steps = 5,
                                                                  init_method = 'kmeans')))
  }
  names(test_mods)<-add_vars
  return(test_mods)
  
}

samdat<-make_mixture_data(sp_200, env_200)
names(samdat)[1:30]<-species

#Null
sam_null<-species_mix(archetype_formula = as.formula(paste("cbind(",paste(species, collapse=", "),")~ 1")),
                      species_formula = ~1,
                      data = samdat,
                      n_mixtures = 1,
                      distribution = 'bernoulli',
                      standardise = FALSE,
                      control = species_mix.control(em_refit = 3, em_steps = 5,init_method = 'kmeans'))
#null model doesn't work:
#Error in stats::kmeans(beta, centers = G, nstart = 100) : 
#  more cluster centers than distinct data points.

#step1
sam_step1<-sams_fwd_step(start_vars="", add_vars = env_vars, n_mixtures=3, data= samdat)

BICs<-data.frame(var=env_vars, BIC=sapply( 1:8, function(x) sam_step1[[x]]$BIC))
#keep temp 7076

#step2
sam_step2<-sams_fwd_step(start_vars="temp", add_vars = env_vars[!env_vars %in% "temp"], n_mixtures=3, data= samdat)

BICs<-data.frame(var= env_vars[!env_vars %in% "temp"], BIC=sapply( 1:7, function(x) sam_step2[[x]]$BIC))
#keep O2 6957

#step3
sam_step3<-sams_fwd_step(start_vars=c("temp", "O2"), add_vars = env_vars[!env_vars %in% c("temp", "O2")], 
                         n_mixtures=3, data= samdat)

BICs<-data.frame(var= env_vars[!env_vars %in% c("temp", "O2")], 
                 BIC=sapply( 1:6, function(x) sam_step3[[x]]$BIC))
#NO3 6957. Same as step2- stop here.

