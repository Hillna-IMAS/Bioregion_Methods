#############################################################################################
## Compare community modelling methods for bioregionalisation of simulated species dataset
## March 2018. N.Hill with input from S. Woolley
#############################################################################################

# Modelling process
# 1) Run models,diagnostics, determine number of groups or set number of groups to 3
# 2) Plot Predicted Distribution of groups across simulation region
# 3) Describe contents of groups
# 4) DESCRIBE ENVIRONMENT OF GROUPS


#######################
## Set up---
#######################
# Get required libraries
library(tidyr)           #data manipulation
library(plyr)             #data manipulation
library(ggplot2)         #plotting
#library(raster)         #spatial data
#library(RColorBrewer)   #colour palettes
#library(rasterVis)      #plotting rasters
library(RCPmod)          #Regions of Common Profile
#library(bbgdm)          #Bayesian Bootstrap Generalised Dissimilarity Models (and naive GDM)
library(gdm)             # Generalised Dissimilarity Models
library(gradientForest)  #Gradient Forests (an extension of Random Forests)
library(ecomix)          #Species Archetype Models (SAMs)
library(SDMTools)        #weighted means and sd

setwd("C:\\Users\\hillna\\UTAS_work\\Projects\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\KP_Fish\\")
source("Code/Additional_Funcs.R")

#Load required files
load("Data/Prediction_space.Rda") 
load("Results/models.RData")
load("Results/pred_clusters.RData")
dat<- readRDS("Data/KP_Fish_Env.RDS")

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
## Tabulate or get directly from models contents of groups- Limit to case where groups =4
## Generate dotchart of average and se of species' occurence in each group
## Generate partial plots or equivalent of the group response to the environment for SAM & RCP
############################################################################################################


# A) Cluster environment only- directly from cluster results for entire region
# 2 Stage Models- cluster then predict:
# B) cluster bioloigcal data, predict clusters with random forests- tabulate from bio clusters
# 2 Stage Models- predict then cluster
# C) predict species using random forests, then cluster predictions- spatially match cluster prediction to sample & tabulate environment
# D) predict species using Mistnet, then cluster predictions- as above
# E) predict species using HMSC, then cluster predictions- as above
# F) predict dissimilarities using GDM, cluster predicted dissimilarities- as above, plus function plots
# G) predict biologically transformed environment (GF), cluster predictions- as above
# 1 Stage mdoel-based: cluster and predict
# H) Species Archetype Models- for ease of comparison, convert results to hard classes and use same process as for 2-stage methods, plus partial plots
# I) Regions of Common Profile- for ease of comparison, convert results to hard classes and use same process as for 2-stage methods, plus partial plots


## Account for label switching---
#convert probabalistic methods to give hard classes
hard_cluster4$BioHC_RF<-apply(bio4_rf_pred,1,which.max)
hard_cluster4$SAM<-apply(sam4_boot_pred$ptPreds,1,which.max)
hard_cluster4$RCP<-apply(rcp4_spat_preds[["ptPreds"]],1,which.max)

#fix label switching
hard_cluster4$SAM<-mapvalues(hard_cluster4$SAM ,from=c(1,2,3,4), to=c(2,4,1,3))
hard_cluster4$RCP<-mapvalues(hard_cluster4$RCP ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$SpRF_HC<-mapvalues(hard_cluster4$SpRF_HC ,from=c(1,2,3,4), to=c(2,4,1,3))
hard_cluster4$GDM_Dissim_HC<-mapvalues(hard_cluster4$GDM_Dissim_HC ,from=c(1,2,3,4), to=c(2,4,3,1))
hard_cluster4$GDM_TransEnv_HC<-mapvalues(hard_cluster4$GDM_TransEnv_HC ,from=c(1,2,3,4), to=c(2,3,1,4))
hard_cluster4$bbGDM_Dissim_HC<-mapvalues(hard_cluster4$bbGDM_Dissim_HC ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$bbGDM_TransEnv_HC<-mapvalues(hard_cluster4$bbGDM_TransEnv_HC, from=c(1,2,3,4), to=c(2,3,1,4))
hard_cluster4$GF_HC<-mapvalues(hard_cluster4$GF_HC ,from=c(1,2,3,4), to=c(3,2,4,1))
hard_cluster4$MNet_HC<-mapvalues(hard_cluster4$MNet_HC ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$HMSC_HC<-mapvalues(hard_cluster4$HMSC_HC ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$BioHC_RF<-mapvalues(hard_cluster4$BioHC_RF ,from=c(1,2,3,4), to=c(1,4,2,3))


#create raster stack of predictions and extract survey site values
mods2<-c("Env_Only", "BioHC_RF", "SpRF_HC","HMSC_HC" , "MNet_HC" , 
         "GDM_Dissim_HC","GDM_TransEnv_HC",  "bbGDM_Dissim_HC"  , "bbGDM_TransEnv_HC",
         "GF_HC", "SAM", "RCP")


##plot groups
rat2<-data.frame(ID=1:4, Group=paste0("Group", 1:4))

clust4<-stack()

for (i in 1: length(mods2)){
  
  hc_rast<-rasterize(na.omit(pred_sp)[,1:2], env_raster, field=hard_cluster4[,mods2[i]])
  hc_rast<-as.factor(hc_rast)
  levels(hc_rast) <- rat2
  
  clust4<-stack(clust4, hc_rast)
}

names(clust4)<-mods2

x11()
plot(clust4)


site_classes<-as.data.frame(raster::extract(clust4, dat[,5:6]))

### Extract environmental values at hard clasessed sites

#clusts<-names(site_classes)[1:10]
clusts<-names(site_classes)[1:12]
hclust_contents_SD<-list()
hclust_contents_SE<-list()

for( i in seq_along(clusts)){
  clust_vals<-env_match_vals( site_data= dat[,env_vars], 
                              pred_cluster_vals=site_classes[,clusts[i]] ,
                              site_index=1:nrow(dat))
  
  hclust_contents_SD[[i]]<-dotplot_env_tab(mean_df = clust_vals$mean,
                                       error_df = clust_vals$sd,
                                       nGrp=4, env=env_vars, method=clusts[i])
  
  hclust_contents_SE[[i]]<-dotplot_env_tab(mean_df = clust_vals$mean,
                                           error_df = clust_vals$se,
                                           nGrp=4, env=env_vars, method=clusts[i])
}
 
## Alternative method for RCP and SAMs that takes into account the probabalistic nature of bioregion predictions.
##RCP----
#using all posterior site probabilities
#rcp_calc_postprob<-list(as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=4)),
#                        as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=4)))

#for (i in 1:length(env_vars)){
#  rcp_calc_postprob[[1]][i,1]<-wt.mean(dat[,env_vars[i]],rcp4_mod$pis[,1])
#  rcp_calc_postprob[[1]][i,2]<-wt.mean(dat[,env_vars[i]],rcp4_mod$pis[,2])
#  rcp_calc_postprob[[1]][i,3]<-wt.mean(dat[,env_vars[i]],rcp4_mod$pis[,3])
#  rcp_calc_postprob[[1]][i,4]<-wt.mean(dat[,env_vars[i]],rcp4_mod$pis[,4])
  
#  rcp_calc_postprob[[2]][i,1]<-wt.sd(dat[,env_vars[i]],rcp4_mod$pis[,1])
#  rcp_calc_postprob[[2]][i,2]<-wt.sd(dat[,env_vars[i]],rcp4_mod$pis[,2])
#  rcp_calc_postprob[[2]][i,3]<-wt.sd(dat[,env_vars[i]],rcp4_mod$pis[,3])
#  rcp_calc_postprob[[2]][i,4]<-wt.sd(dat[,env_vars[i]],rcp4_mod$pis[,4])
#}


#rcp_calc_postprob[[3]]<-rcp_calc_postprob[[2]]/sqrt(nrow(dat))
#names(rcp_calc_postprob)<-c("mean", "sd", "se")

# rcp_post_contents_SD<-dotplot_env_tab(mean_df = rcp_calc_postprob$mean,
#                                     error_df = rcp_calc_postprob$sd,
#                                     nGrp=4, env=env_vars, method="RCP")

#rcp_post_contents_SE<-dotplot_env_tab(mean_df = rcp_calc_postprob$mean,
#                                     error_df = rcp_calc_postprob$se,
#                                     nGrp=4, env=env_vars, method="RCP")

#fix label switching issue
#rcp_post_contents_SD$Group<-mapvalues(rcp_post_contents_SD$Group ,from=c(1,2,3,4), to=c(2,3,4,1))
#rcp_post_contents_SE$Group<-mapvalues(rcp_post_contents_SE$Group ,from=c(1,2,3,4), to=c(2,3,4,1))


##SAMs----
# get spatial group predictions and match initial sites
#sam_rast<-rasterize(na.omit(pred_sp)[,1:2], env_raster, sam4_boot_pred$ptPreds)

#sam_probs<-raster::extract(sam_rast, dat[,c("Long", "Lat")])

#sam_calc_postprob<-list(as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=4)),
#                        as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=4)))

#for (i in 1:length(env_vars)){
#  sam_calc_postprob[[1]][i,1]<-wt.mean(dat[, env_vars[i]],sam_probs[,1])
#  sam_calc_postprob[[1]][i,2]<-wt.mean(dat[, env_vars[i]],sam_probs[,2])
#  sam_calc_postprob[[1]][i,3]<-wt.mean(dat[, env_vars[i]],sam_probs[,3])
#  sam_calc_postprob[[1]][i,4]<-wt.mean(dat[, env_vars[i]],sam_probs[,4])
  
#  sam_calc_postprob[[2]][i,1]<-wt.sd(dat[, env_vars[i]],sam_probs[,1])
#  sam_calc_postprob[[2]][i,2]<-wt.sd(dat[, env_vars[i]],sam_probs[,2])
#  sam_calc_postprob[[2]][i,3]<-wt.sd(dat[, env_vars[i]],sam_probs[,3])
#  sam_calc_postprob[[2]][i,4]<-wt.sd(dat[, env_vars[i]],sam_probs[,4])
#}

#sam_calc_postprob[[3]]<-sam_calc_postprob[[2]]/sqrt(nrow(dat))
#names(sam_calc_postprob)<-c("mean", "sd", "se")

#sam_post_contents_SD<-dotplot_env_tab(mean_df = sam_calc_postprob$mean,
#                                      error_df = sam_calc_postprob$sd,
#                                      nGrp=4, env=env_vars, method="SAM")

#sam_post_contents_SE<-dotplot_env_tab(mean_df = sam_calc_postprob$mean,
#                                      error_df = sam_calc_postprob$se,
#                                      nGrp=4, env=env_vars, method="SAM")

#fix label switching issue
#sam_post_contents_SD$Group<-mapvalues(sam_post_contents_SD$Group ,from=c(1,2,3,4), to=c(2,4,1,3))
#sam_post_contents_SE$Group<-mapvalues(sam_post_contents_SE$Group ,from=c(1,2,3,4), to=c(2,4,1,3))


#dotplot of contents ----
#create blank row/separator for each approach (clunky way of getting ggplot2 to have the legend I want!)
gtp<-dotplot_env_tab(mean_df = matrix(NA,nrow=8, ncol=4),
                    error_df = matrix(0,nrow=8, ncol=4),
                    nGrp=4, env=env_vars, method="Two-stage: Group then Predict")
ptg<-dotplot_env_tab(mean_df = matrix(NA,nrow=8, ncol=4),
                    error_df = matrix(0,nrow=8, ncol=4),
                    nGrp=4, env=env_vars, method="Two-stage: Predict then Group")
os<-dotplot_env_tab(mean_df = matrix(NA,nrow=8, ncol=4),
                   error_df = matrix(0,nrow=8, ncol=4),
                   nGrp=4, env=env_vars, method="One-stage")

#contents_SD<-rbind(do.call(rbind, hclust_contents_SD),
#                   rcp_post_contents_SD, sam_post_contents_SD, gtp, ptg, os)

contents_SD<-rbind(do.call(rbind, hclust_contents_SD),
                   gtp, ptg, os)

#contents_SD$Method<-mapvalues(contents_SD$Method ,from="BioHC_RF", to="BioHC_RF (Hard Class)")
contents_SD$Method<-mapvalues(contents_SD$Method ,from=c("BioHC_RF","SAM", "RCP"), to=c("BioHC_RF (Hard Class)", "SAM (Hard Class)", "RCP (Hard Class)"))


#contents_SE<-rbind(do.call(rbind, hclust_contents_SE),
#                   rcp_post_contents_SE, sam_post_contents_SE)

#set order of plotting
#contents_SD$Method<-contents_SD$Method<-factor(contents_SD$Method, 
#                        levels=c("RCP", "SAM",  "One-stage", 
#                              "GF_HC",  "bbGDM_Dissim_HC", "bbGDM_TransEnv_HC", "GDM_TransEnv_HC",  "GDM_Dissim_HC",
#                              "MNet_HC", "HMSC_HC", "SpRF_HC", "Two-stage: Predict then Group",
#                              "BioHC_RF (Hard Class)",  "Two-stage: Group then Predict", "Env_Only"))

contents_SD$Method<-contents_SD$Method<-factor(contents_SD$Method, 
                                               levels=c("RCP (Hard Class)", "SAM (Hard Class)",  "One-stage", 
                                                        "GF_HC",  "bbGDM_Dissim_HC", "bbGDM_TransEnv_HC", "GDM_TransEnv_HC",  "GDM_Dissim_HC",
                                                        "MNet_HC", "HMSC_HC", "SpRF_HC", "Two-stage: Predict then Group",
                                                        "BioHC_RF (Hard Class)",  "Two-stage: Group then Predict", "Env_Only"))

contents_SD$env<-mapvalues(contents_SD$env, from=env_vars, to=pretty_env)

contents_SD$env<-factor(contents_SD$env, levels=c("Depth", "Surface_temp",   "Floor_temp"  , "ssha_SD" ,  "Current", "NO3_mean", "Chla_SD" , "Slope" ))

names(contents_SD)[2]<-"Value"


p<-ggplot(data = contents_SD, 
          aes(x = Group, y = Value, ymin = lower, ymax = upper, colour = Method)) +
  geom_point(position = position_dodge(0.75), size=0.6) +
  geom_errorbar(position = position_dodge(0.75), width = 0.1) +
  coord_flip() +
  scale_colour_manual(name="Method", 
                      values=c( "darkblue","cyan",  "white", 
                                "darkslategray", "darkgreen", "darkseagreen2","chartreuse3", "chartreuse1",
                                "orange", "yellow"  , "pink", "white",
                                "purple" , "white", "grey"))+
  
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.y = element_text(face="italic"),
        axis.text.x = element_text(size=6),
        legend.key = element_blank()) +
  facet_wrap( ~env, ncol=4, scales="free") 

tiff(file="Results/Plots/Grp_env_SD.tiff", height=13, width=17, units="cm", res=1000)
p + 
scale_x_discrete(name="Bioregion", breaks=c(1,2,3,4), limits=c("1","2", "3", "4")) +
guides(color = guide_legend(reverse = TRUE))
dev.off()




#############################################################################
## Shape of environmental response of GROUPS -----
## A) bio_cluster then RF- RF partial plots 
## B) bbGDM- fitted function (but for 'turnover' not groups, so not presented)
## C) GF- distribution of split points (but for 'turnover' not for groups, so not presented)
## D) SAM- partial plot of species' group responses (prob group occurrence)
## E) RCP- partial plot of RCP group responses (prob group ocurrence)
############################################################################

## A) cluster bio then predict with RF----
#note label switching
#hard_cluster3$BioHC_RF<-mapvalues(hard_cluster3$BioHC_RF ,from=c(1,3), to=c(3,1))

tiff(file="Results/Plots/RF_partial_response.tiff", height=6, width =6, units="in", res=1000)
par(mfrow=c(3,3), mar=c(4,2,1,1), oma=c(3,3,1,1))

for(i in 1:length(env_vars)){
#for(i in 1:4){
partialPlot(bio4_rf, as.data.frame(dat[,env_vars]), env_vars[i], main="",
            xlab=pretty_env[i], ylab= "", col="darkolivegreen4", lty=3, lwd=1.5, ylim=c(-2, 2))

partialPlot(bio4_rf, as.data.frame(dat[,env_vars]), env_vars[i], add=TRUE,ylim=c(-2, 2),
            which.class="2", xlab=pretty_env[i], ylab= "", col="darkred", lty=1, lwd=1.5, main="")

partialPlot(bio4_rf, as.data.frame(dat[,env_vars]), env_vars[i], add=TRUE,ylim=c(-2, 2),
            which.class="3", xlab=pretty_env[i], ylab= "", col="grey", lty=2, lwd=1.5, main="")

partialPlot(bio4_rf, as.data.frame(dat[,env_vars]), env_vars[i], add=TRUE,ylim=c(-2, 2),
            which.class="4", xlab=pretty_env[i], ylab= "", col="orange1", lty=1, lwd=1.5, main="")

}
plot.new()

#manually match up group labels
#hard_cluster4$BioHC_RF<-mapvalues(hard_cluster4$BioHC_RF ,from=c(1,2,3,4), to=c(1,4,2,3))
legend("center", legend= c(paste0("Bioregion ", 1:4)), lty=c(3,2,1,1), xpd=NA, bty="n",
       col=c("darkolivegreen4", "grey", "orange1","darkred" ))
title(ylab="Response", outer=TRUE, line=1)
dev.off()



## B) bbGDM----
#bbgdm_response<-as.response(non_bbgdm_mod)

#tiff(file="Results/Plots/bbGDM_PPlots.tiff", height = 5, width=7, units="in", res=1000)
#par(mfrow=c(3,3), mar=c(4,4,2,2))

#  Splinesum <- Spline.05 <- Spline.95 <- NULL
#  Splinessum <- mapply(`%*%`, bbgdm_response$X, bbgdm_response$bspl)
#  Splines.05 <- mapply(`%*%`, bbgdm_response$X, bbgdm_response$bspl.05)
#  Splines.95 <- mapply(`%*%`, bbgdm_response$X, bbgdm_response$bspl.95)
#  for (i in 1:ncol(Splinessum)) {
#    plot(bbgdm_response$grid_real[, i], Splinessum[, i], type = "l", ylab = paste0("f(", 
#                                                                      pretty_env[i], ")"), xlab = pretty_env[i], 
#         ylim = range(c(Splinessum, Splines.05, Splines.95)))
#    polygon(c(bbgdm_response$grid_real[, i], rev(bbgdm_response$grid_real[, i])), c(Splines.05[, 
#                                                                     i], rev(Splines.95[, i])), col = "grey80", border = NA)
#    lines(bbgdm_response$grid_real[, i], Splinessum[, i], col = "black", 
#          type = "l", lwd = 2)
#    mtext(paste0("(", letters[i], ")"), adj = 0)
#  }
#plot(bbgdm_response)
#dev.off()


## C) SAMs----
#generate values for plotting (varying one environmental variable at a time, holding others at mean value)
part_vals<-partial_values(env_vars = env_vars, raw_data = dat[,env_vars])

#scale data using same scaling as used to build models
quad_dat<-poly_data(env_vars, degree=c(rep(2, length(env_vars))), 
                    id_vars = c("sample_ID", "Year", "Survey", "Lat", "Long"), species=species, data=dat)

part_sc<-lapply(part_vals, poly_pred_space, poly_output= quad_dat$poly_output, 
                vars=env_vars, sampling_levels=levels(quad_dat$rcp_data$Survey), sampling_factor="Survey")
env2_vars<-paste0(rep(env_vars, each=2), 1:2)

#predict to each set of variables and plot
tiff(file="Results/Plots/SAM_partial2.tiff", height=6, width=6, units="in", res=1000)

par(mfrow=c(3,3), mar=c(4,2,1,1), oma=c(4,5,1,1))
for(i in 1:length(part_vals)){
  temp<-as.data.frame(part_sc[[i]][,env2_vars])
  boots<- species_mix.bootstrap(sam4_mod, nboot=10, type="BayesBoot")
   pred<-predict(object=sam4_mod,boots, newdata=temp)
   #account for label switching
   dimnames(pred$ptPreds)[[2]]<- c(paste0("Bioregion", c(2,4,1,3)))
   pred$ptPreds<-as.data.frame(pred$ptPreds)
   pred$ptPreds<-pred$ptPreds[,c(paste0("Bioregion", c(1,2,3,4)))]
    matplot(x=as.vector(part_vals[[i]][,names(part_vals[[i]]) %in% env_vars[i]]), 
            y= pred$ptPreds, type='l', col=c("darkolivegreen4", "grey", "orange1", "darkred"), lty=c(3,2,1,1), xlab=pretty_env[i], ylab=NULL, lwd=1.5, ylim=c(0,1))

}
plot.new()
legend(x="center", legend=paste0("Bioregion ", 1:4), col=c("darkolivegreen4", "grey", "orange1", "darkred"),lty=c(3,2,1,1), 
       cex=1.2, bty="n")
mtext(text="Probability of Occurrence", side=2, line=1, outer=TRUE)
dev.off()

## E) RCPs----

#predict to each set of variables and plot
rcpmod_vars<-c("Av_depth", "T_mean")

tiff(file="Results/Plots/RCP_partial.tiff", height=2.5, width=6, units="in", res=1000)
par(mfrow=c(1,3), mar=c(4,2,1,1), oma=c(4,5,1,1))
for(i in 1:length(rcpmod_vars)){
  temp<-as.data.frame(part_sc[[rcpmod_vars[i]]])
  pred<-as.data.frame(predict(rcp4_mod, newdata=temp))
  #account for label switching
  names(pred)<- c(paste0("Bioregion", c(2,3,4,1)))
  pred<-pred[, c(paste0("Bioregion", c(1,2,3,4)))]
  pretty_lab<-c("Depth", "Surface_temp")
 # matplot(x=as.vector(part_vals[[i]][,names(part_vals[[i]]) %in% env_vars[i]]),
          matplot(x=as.vector(part_vals[[rcpmod_vars[i]]][,rcpmod_vars[i]]),
          y= pred, type='l', col=c( "darkolivegreen4", "grey", "orange1", "darkred"),lty=c(3,2,1,1), xlab=pretty_lab[i], ylab=NULL, lwd=1.5, ylim=c(0,1))
  
}
plot.new()
legend(x="center", legend=c(paste0("Bioregion ", c(1,2,3,4))), col=c("darkolivegreen4", "grey", "orange1", "darkred"), 
       lty=c(3,2,1,1),    cex=1.1, bty="n")
mtext(text="Probability of Occurrence", cex=0.7, side=2, line=1, outer=TRUE)
dev.off()


