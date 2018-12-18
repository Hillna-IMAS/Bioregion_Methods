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

setwd("C:\\Users\\hillna\\UTAS_work\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\KP_Fish\\")
source("Code/Additional_Funcs.R")

#Load required files
load("Results/models.RData")
#load("pred_clusters.RData")
load("Data/Prediction_space.Rda")
dat<- readRDS("Data/KP_Fish_Env.RDS")


# Environmental variables to test
env_vars<-c("log_slope", "Av_depth" ,
            "sea.floor.temperature", 
            "log_current"   ,   "no3_bot_mean" ,            
            "T_mean"  , "chla_yearly_sd"  ,  "ssha_variability")

pretty_env<-c("Slope", "Depth" ,
            "Floor_temp", "Current" ,"NO3_mean" ,            
              "Surface_temp"  , "Chla_SD"  ,  "ssha_SD")

# Species to model
species=c("Antimora.rostrata" ,"Bathydraco.antarcticus", "Bathyraja.eatonii",                
          "Bathyraja.irrasa"  ,  "Bathyraja.murrayi",  "Champsocephalus.gunnari" , "Dissostichus.eleginoides" ,        
          "Etmopterus.viator"  ,  "Gobionotothen.acuta" , "Lepidonotothen.mizops"   ,  "Lepidonotothen.squamifrons",       
          "Lycodapus.antarcticus" ,"Macrourus.spp", "Mancopsetta.maculata"  , "Muraenolepis.spp",                
          "Notothenia.rossii"   , "Paradiplospinus.gracilis" ,  "Paraliparis.spp." ,"Zanclorhynchus.spinifer" ,         
          "Channichthys.rhinoceratus.velifer")

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
     width=7, height=8, units="cm", compression = "lzw", res=1000)

imp<-as.data.frame(importance(bio4_rf, type=1))
imp$var<-pretty_env
ord <- rev(order(imp$MeanDecreaseAccuracy, decreasing = TRUE))
par(mar=c(4,4,1,1))
dotchart(imp[ord,1], xlab = "", labels=imp$var[ord], ylab = "", 
         main = "", cex=0.8)

#mtext(outer = TRUE, side = 1, text = "Average Importance", cex = 0.5)
#varImpPlot(bio4_rf, type=1,main="", cex=0.9, cex.main=0.8)
dev.off()

## C) Sp_RF----

#calculate average importance over all species
rf_sp_varimp<-array(dim=c(length(env_vars),2,length(species)), dimnames=list(env_vars, c("MeanDecreaseAccuracy", "MeanDecreaseGini"), species))
for( i in 1:length(species)){
  rf_sp_varimp[, , i]<- rf_sp_mods[[i]]$importance[,3:4]
}

tiff(file="Results/Plots/spRF_Overall_imp.tiff", width=7, height=8, units="cm",
       compression = "lzw", res=1000)
imp<-as.data.frame(apply(rf_sp_varimp, c(1,2), mean))
imp$var<-pretty_env
# modified from VarImpPlot function in RandomForest package
#op <- par(mfrow = c(1, 1), mar = c(4, 5, 4, 1), mgp = c(2,0.8, 0), 
#          oma = c(0, 0, 2, 0), no.readonly = TRUE)
#on.exit(par(op))
par(mar=c(4,4,1,1))
#for(i in 1:2){
for (i in 1:1) {
  ord <- rev(order(imp[, i], decreasing = TRUE)[1:length(env_vars)])
  #xmin <- if (colnames(imp)[i] %in% c("IncNodePurity", "MeanDecreaseGini")) 
  xmin <- if (colnames(imp)[i] %in% c("IncNodePurity")) 
    0
  dotchart(imp[ord, i], xlab = "", labels=imp$var[ord], ylab = "", 
           main = "", cex=0.8)
}
#mtext(outer = TRUE, side = 1, text = "Average Importance", cex = 0.5)
dev.off()

## D) HMSC----
# calculated for each environmental variable forspecies, then averaged across all species and plotted
hmsc_sp_var<-variPart(mod_hmsc, groupX=rep(env_vars,each=2))
hmsc_var<-data.frame(variable=pretty_env, mean=apply(hmsc_sp_var[,-9], 2, mean), se=apply(hmsc_sp_var[,-9], 2, function(x) sd(x)/sqrt(30)))
hmsc_var<-hmsc_var[order(hmsc_var$mean, decreasing=FALSE),]

tiff(file="Results/Plots/hmsc_imp.tiff", width=7, height=10, units="cm",
     compression = "lzw", res=1000)
par(mar=c(4,6,1,1))
cent<-barplot(hmsc_var$mean, beside=TRUE, space=0, names.arg=hmsc_var$variable,horiz=TRUE, las=1,xlim=c(0,0.7),
              cex.axis=0.9, cex.names=0.9)
arrows(y0=cent, x0=hmsc_var$mean + hmsc_var$se, x1=hmsc_var$mean - hmsc_var$se,  angle=90, code=3, length=0.05)
dev.off()

##E) Misnet ----
MNet_imp<-olden_mistnet(MNet_mod)
MNet_overall_imp<-MNet_imp$`Overall relative importance`/100
MNet_overall_imp$var<-c(pretty_env, paste0("latent", 1:5))
MNet_overall_imp<-MNet_overall_imp[order(MNet_overall_imp$mean, decreasing=FALSE),]


tiff(file="Results/Plots/MNet_imp.tiff", width=7, height=10, units="cm",
     compression = "lzw", res=1000)
par(mar=c(4,6,1,1))
cent<-barplot(MNet_overall_imp$mean, beside=TRUE, space=0, names.arg=MNet_overall_imp$var,horiz=TRUE, las=1,
              cex.names=0.9, cex.axis=0.9, xlim=c(0,0.5))
arrows(y0=cent, x0=MNet_overall_imp$mean + MNet_overall_imp$sd/sqrt(30), x1=MNet_overall_imp$mean - MNet_overall_imp$sd/sqrt(30),  angle=90, code=3, length=0.05)
dev.off()


## F) bbGDM----
write.csv(round(bbgdm.wald.test(gdm_mod, gdm=TRUE),3), 
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
  
  for(j in 1:(length(add_vars)/2)){
    add<-add_vars[(j*2-1):(j*2)]
    temp_form<-as.formula(paste("cbind(",paste(species, collapse=", "),")~1 +",
                                paste0(paste(start_vars, collapse="+"), "+", paste(add, collapse = "+"))))
  #run SAMs
  test_mods[[j]] <- try(species_mix(archetype_formula = temp_form,
                                    species_formula = ~1,
                                    data = data,
                                    n_mixtures = n_mixtures,
                                    distribution = 'bernoulli',
                                    standardise = FALSE,
                                    control = species_mix.control(em_prefit=TRUE, em_refit = 5, em_steps = 8,
                                                                  init_method = 'kmeans')))
  }
  names(test_mods)<-add[1]
  return(test_mods)
  
}

quad_dat<-poly_data(env_vars, degree=c(rep(2, length(env_vars))), 
                    id_vars = c("sample_ID", "Year", "Survey", "Lat", "Long"), species=species, data=dat)

env2_vars<-paste0(rep(env_vars, each=2), 1:2)
sam_dat <- make_mixture_data(quad_dat$rcp_data[,species], quad_dat$rcp_data[, env2_vars])

#Null
sam_null<-species_mix(archetype_formula = as.formula(paste("cbind(",paste(species, collapse=", "),")~ 1")),
                      species_formula = ~1,
                      data = sam_dat,
                      n_mixtures = 1,
                      distribution = 'bernoulli',
                      standardise = FALSE,
                      control = species_mix.control(em_prefit=TRUE, em_refit = 5, 
                                                    em_steps = 8,init_method = 'kmeans'))
#null model doesn't work:
#Error in stats::kmeans(beta, centers = G, nstart = 100) : 
#  more cluster centers than distinct data points.

#step1
sam_step1<-sams_fwd_step(start_vars="", add_vars = env2_vars, n_mixtures=4, data= sam_dat)

BICs<-data.frame(var=env_vars, BIC=sapply( 1:8, function(x) sam_step1[[x]]$BIC))
#keep Av_depth 7080

#step2
sam_step2<-sams_fwd_step(start_vars=c("Av_depth1","Av_depth2") ,
                        add_vars = env2_vars[!env2_vars %in% c("Av_depth1","Av_depth2")],
                         n_mixtures=3, data= sam_dat)

BICs<-data.frame(var= env_vars[!env_vars %in% "Av_depth"], 
                 BIC=sapply( 1:7, function(x) sam_step2[[x]]$BIC))
#T_mean= 7121

#suggests Sam with just depth
sam_best<-species_mix(archetype_formula = as.formula(paste("cbind(",paste(species, collapse=", "),")~ 1 +Av_depth1 + Av_depth2 ")),
                      species_formula = ~1,
                      data = sam_dat,
                      n_mixtures = 4,
                      distribution = 'bernoulli',
                      standardise = FALSE,
                      control = species_mix.control(em_prefit=TRUE, em_refit = 5, 
                                                    em_steps = 8,init_method = 'kmeans'))
sam_best$vcov<-vcov(sam_best)

#generate prediction
trans_pred_space<-poly_pred_space(pred_space=na.omit(pred_sp), poly_output= quad_dat$poly_output, 
                                  vars=env_vars, sampling_levels=levels(quad_dat$rcp_data$Survey), sampling_factor="Survey")


pred_sam_best<-predict(sam_best, as.data.frame(trans_pred_space[,c("Av_depth1", "Av_depth2")]))

test<-rasterize(na.omit(pred_sp[,1:2]), env_raster, pred_sam_best$fit[,1])
