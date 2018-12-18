#############################################################################################
## Compare community modelling methods for bioregionalisation of simulated species dataset
## March 2018. N.Hill with input from S. Woolley
#############################################################################################

# Modelling process
# 1) Run models,diagnostics, determine number of groups or set number of groups to 3
# 2) Plot Predicted Distribution of groups across simulation region
# 3) Describe contents of groups
# 4) DESCRIBE ENVIRONMENT OF GROUPS
# 5) Predictor importance

#######################
## Set up---
#######################
# Get required libraries
library(tidyr)           #data manipulation
library(ggplot2)        #plotting
#library(raster)         #spatial data
#library(RColorBrewer)   #colour palettes
#library(rasterVis)      #plotting rasters
library(RCPmod)         #Regions of Common Profile
library(bbgdm)          #Bayesian Bootstrap Generalised Dissimilarity Models (and naive GDM)
library(gradientForest) #Gradient Forests (an extension of Random Forests)
library(ecomix)     #Species Archetype Models (SAMs)
library(SDMTools)       #weighted means and sd

setwd("C:\\Users\\hillna\\UTAS_work\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\Simulation\\")
source("Simulation_Additional_Funcs.R")

#Load required files
#files in "simulate_communities" folder on dropbox

load("Sim_setup/Many_covars_sim.RData")
#sim_data= simulated species probabilities, occurrences, species' groups for entire region
#sp_200 = matrix of occurrences of 30 species at 200 sites to use as species dataset for analysis
#env_200= matrix of corresponding environmental conditions at 200 site to use as environmental dataset for analysis


load("Sim_setup/sim_env_070518.RData") 
# env= raster brick of environmental data for entire region
# env_dat= raster brick of environmental data for entire region converted to matrix

load("Results/models.RData")
load("Results/pred_clusters.RData")

env_vars<-dimnames(sim_dat)[[2]][2:9]

############################################################################################################
## Take note of label switching from plot_pred_clusters
## Tabulate or get directly from models contents of groups- Limit to case where groups =3
## Generate dotchart of average and se of species' occurence in each group
## Generate partial plots or equivalent of the group response to the environment for GDM, GF, SAM & RCP
############################################################################################################


# A) Cluster environment only- directly from cluster results for entire region
# 2 Stage Models- cluster then predict:
# B) cluster bioloigcal data, predict clusters with random forests- tabulate from bio clusters
# 2 Stage Models- predict then cluster
# C) predict species using random forests, then cluster predictions- spatially match cluster prediction to sample & tabulate environment
# D) predict species using Mistnet, then cluster predictions- as above
# E) predict species using HMSC, then cluster predictions- as above
# F) predict dissimilarities using GDM, cluster predicted dissimilarities- as above, plus function plots
# G) predict biologically transformed environment (GF), cluster predictions- as above, plus look at cumulative splits
# 1 Stage mdoel-based: cluster and predict
# H) Species Archetype Models- direct from model, plus partial plots
# I) Regions of Common Profile- direct from model, plus partial plots


## Account for label switching---
# correct previously identified label switching 
hard_cluster3$env<-mapvalues(hard_cluster3$env ,from=c(1,3), to=c(3,1))
#hard_cluster3$BioHC_RF<-mapvalues(hard_cluster3$BioHC_RF ,from=c(1,3), to=c(3,1))
hard_cluster3$Sp_RF<-mapvalues(hard_cluster3$Sp_RF ,from=c(1,2), to=c(2,1))
hard_cluster3$HMSC<-mapvalues(hard_cluster3$HMSC ,from=c(1,2,3), to=c(2,1,3))
hard_cluster3$MNet<-mapvalues(hard_cluster3$MNet ,from=c(1,3), to=c(3,1))
hard_cluster3$bbGDM<-mapvalues(hard_cluster3$bbGDM ,from=c(1,2), to=c(2,1))
hard_cluster3$GF<-mapvalues(hard_cluster3$GF ,from=c(1,2), to=c(2,1))




## Get rownumbers of sites used for building models to match environment at that site and predicted groups from heirarchical clustering
set.seed(66)       
site_ind<-sample(1:nrow(sim_dat), 200, replace=FALSE)

## 'True' Distribution
# from probability of group occurrence (using average prevalence of all species)
set.seed(42)
#simulation alphas 
betamean <- 0.3
betabeta <- 15
betaalpha <- betamean/(1-betamean) * betabeta
prevalences <- rbeta( 30, betaalpha, betabeta) 
alphas <- log( prevalences / ( 1-prevalences))  

# simulation betas
means<-as.matrix(data.frame(temp=c(0.75,0,-0.5), O2=c(0,-0.5,0), NO3=c(0,0,0), sal=c(0,0,0),
                            depth=c(0, 0,0), chla=c(0,0,0), ssh=c(0,0,0), curr=c(0,0,0)))

#calculate average probability of group occurrence for each cell across region
mean_alpha<-mean(alphas)
true_lps<-mean_alpha+ as.matrix(sim_dat[,2:9])%*% t(means)
true_grps<-exp(true_lps)/(1 +exp(true_lps))

# from probability of group occurrence extract probabilities at sampling sites
site_true_probs<-true_grps[site_ind,]
site_true_HC<-apply(site_true_probs, 1, which.max)

#true_prob<-list(as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=3)),
#                        as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=3)))

#for (i in 1:length(env_vars)){
#  true_prob[[1]][i,1]<-wt.mean(env_dat[site_ind,i+2],site_true_probs[,1])
#  true_prob[[1]][i,2]<-wt.mean(env_dat[site_ind,i+2],site_true_probs[,2])
#  true_prob[[1]][i,3]<-wt.mean(env_dat[site_ind,i+2],site_true_probs[,3])
  
#  true_prob[[2]][i,1]<-wt.sd(env_dat[site_ind,i+2],site_true_probs[,1])
#  true_prob[[2]][i,2]<-wt.sd(env_dat[site_ind,i+2],site_true_probs[,2])
#  true_prob[[2]][i,3]<-wt.sd(env_dat[site_ind,i+2],site_true_probs[,3])
#}


#true_prob[[3]]<-true_prob[[2]]/sqrt(nrow(sp_200))
#names(true_prob)<-c("mean", "sd", "se")

#trueprob_contents_SD<-dotplot_env_tab(mean_df = true_prob$mean,
#                                      error_df = true_prob$sd,
#                                      nGrp=3, env=env_vars, method="True")

#trueprob_contents_SE<-dotplot_env_tab(mean_df = true_prob$mean,
#                                      error_df = true_prob$se,
#                                      nGrp=3, env=env_vars, method="True")


## Hard classed version of "true"
true_HC_vals<-env_match_vals( site_data= env_dat[site_ind,3:10], 
                              pred_cluster_vals=site_true_HC ,
                              site_index=site_ind)
true_HC_contents_SD<-dotplot_env_tab(mean_df = true_HC_vals[[1]],
                                     error_df = true_HC_vals[[2]],
                                     nGrp=3, env=env_vars, method="True_Hard")

true_HC_contents_SE<-dotplot_env_tab(mean_df = true_HC_vals[[1]],
                                    error_df = true_HC_vals[[3]],
                                    nGrp=3, env=env_vars, method="True_Hard")


## Environment Only--
## Tabulate environment of grid cells belonging to cluster
env_clust_vals<-get_match_vals(site_data= env_dat[,3:10], 
                               pred_cluster_vals=hard_cluster3$env,
                               site_index = 1:nrow(env_dat))

env_clust_contents_SD<-dotplot_env_tab(mean_df = env_clust_vals$mean,
               error_df = env_clust_vals$sd,
               nGrp=3, env=env_vars, method="Env_Only")

env_clust_contents_SE<-dotplot_env_tab(mean_df = env_clust_vals$mean,
                                       error_df = env_clust_vals$se,
                                       nGrp=3, env=env_vars, method="Env_Only")

### 2 stage methods: cluster biology then predict.

#bio_clust_vals<-env_match_vals( site_data= env_dat[site_ind,3:10], 
#                                pred_cluster_vals=bio3_clust ,
#                                site_index=site_ind)

#bio_clust_contents_SD<-dotplot_env_tab(mean_df = bio_clust_vals$mean,
#                                   error_df = bio_clust_vals$sd,
#                                   nGrp=3, env=env_vars, method="BioHC_RF")

#bio_clust_contents_SE<-dotplot_env_tab(mean_df = bio_clust_vals$mean,
#                                       error_df = bio_clust_vals$se,
#                                       nGrp=3, env=env_vars, method="BioHC_RF")

## Hard Classed BioHC_RF
bio3_clust<-mapvalues(bio3_clust,from=c(1,3), to=c(3,1))

bio_clust_vals<-env_match_vals(site_data= env_dat[site_ind,3:10], 
                               pred_cluster_vals=bio3_clust ,
                               site_index=site_ind)

bio_clust_contents_SD<-dotplot_env_tab(mean_df = bio_clust_vals$mean,
                                   error_df = bio_clust_vals$sd,
                                   nGrp=3, env=env_vars, method="BioHC_RF")

bio_clust_contents_SE<-dotplot_env_tab(mean_df = bio_clust_vals$mean,
                                       error_df = bio_clust_vals$se,
                                       nGrp=3, env=env_vars, method="BioHC_RF")

### 2 stage methods: predict then heirarchical cluster

clusts<-names(hard_cluster3)[4:9]
hclust_contents_SD<-list()
hclust_contents_SE<-list()

for( i in seq_along(clusts)){
  clust_vals<-env_match_vals( site_data= env_dat[site_ind,3:10], 
                              pred_cluster_vals=hard_cluster3[site_ind,clusts[i]] ,
                              site_index=site_ind)
  
  hclust_contents_SD[[i]]<-dotplot_env_tab(mean_df = clust_vals$mean,
                                       error_df = clust_vals$sd,
                                       nGrp=3, env=env_vars, method=clusts[i])
  
  hclust_contents_SE[[i]]<-dotplot_env_tab(mean_df = clust_vals$mean,
                                           error_df = clust_vals$se,
                                           nGrp=3, env=env_vars, method=clusts[i])
}
 

##RCP----
#from posterior site probabilities
#rcp_calc_postprob<-list(as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=3)),
#                        as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=3)))

#for (i in 1:length(env_vars)){
#  rcp_calc_postprob[[1]][i,1]<-wt.mean(env_dat[site_ind,i+2],rcp3_mod$pis[,1])
#  rcp_calc_postprob[[1]][i,2]<-wt.mean(env_dat[site_ind,i+2],rcp3_mod$pis[,2])
#  rcp_calc_postprob[[1]][i,3]<-wt.mean(env_dat[site_ind,i+2],rcp3_mod$pis[,3])
  
#  rcp_calc_postprob[[2]][i,1]<-wt.sd(env_dat[site_ind,i+2],rcp3_mod$pis[,1])
#  rcp_calc_postprob[[2]][i,2]<-wt.sd(env_dat[site_ind,i+2],rcp3_mod$pis[,2])
#  rcp_calc_postprob[[2]][i,3]<-wt.sd(env_dat[site_ind,i+2],rcp3_mod$pis[,3])
#}


#rcp_calc_postprob[[3]]<-rcp_calc_postprob[[2]]/sqrt(nrow(sp_200))
#names(rcp_calc_postprob)<-c("mean", "sd", "se")

#rcp_post_contents_SD<-dotplot_env_tab(mean_df = rcp_calc_postprob$mean,
#                                     error_df = rcp_calc_postprob$sd,
#                                     nGrp=3, env=env_vars, method="RCP")

#rcp_post_contents_SE<-dotplot_env_tab(mean_df = rcp_calc_postprob$mean,
#                                     error_df = rcp_calc_postprob$se,
#                                     nGrp=3, env=env_vars, method="RCP")

#fix label switching issue
#rcp_post_contents_SD$Group<-mapvalues(rcp_post_contents_SD$Group ,from=c(1,2,3), to=c(3,1,2))
#rcp_post_contents_SE$Group<-mapvalues(rcp_post_contents_SE$Group ,from=c(1,2,3), to=c(3,1,2))

## Hard class RCPs
rcp_HC<-apply(rcp3_mod$postProbs, 1, which.max)
rcp_HC<-mapvalues(rcp_HC ,from=c(1,2,3), to=c(3,1,2))

rcp_HC_vals<-env_match_vals( site_data= env_dat[site_ind,3:10], 
                pred_cluster_vals=rcp_HC ,
                site_index=1:200)

rcp_HC_contents_SD<-dotplot_env_tab(mean_df = rcp_HC_vals$mean,
                                     error_df = rcp_HC_vals$sd,
                                     nGrp=3, env=env_vars, method="RCP_Hard")

rcp_HC_contents_SE<-dotplot_env_tab(mean_df = rcp_HC_vals$mean,
                                     error_df = rcp_HC_vals$se,
                                     nGrp=3, env=env_vars, method="RCP_Hard")

##SAMs----
# get spatial group predictions and match initial sites
sam_probs<-sam3_pred$fit[site_ind,]
sam_HC<-apply(sam_probs, 1, which.max)
sam_HC<-mapvalues(sam_HC ,from=c(1,2), to=c(2,1))

#sam_calc_postprob<-list(as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=3)),
#                        as.data.frame(matrix(data=NA, nrow=length(env_vars), ncol=3)))

#for (i in 1:length(env_vars)){
#  sam_calc_postprob[[1]][i,1]<-wt.mean(env_dat[site_ind,i+2],sam_probs[,1])
#  sam_calc_postprob[[1]][i,2]<-wt.mean(env_dat[site_ind,i+2],sam_probs[,2])
#  sam_calc_postprob[[1]][i,3]<-wt.mean(env_dat[site_ind,i+2],sam_probs[,3])
  
#  sam_calc_postprob[[2]][i,1]<-wt.sd(env_dat[site_ind,i+2],sam_probs[,1])
#  sam_calc_postprob[[2]][i,2]<-wt.sd(env_dat[site_ind,i+2],sam_probs[,2])
#  sam_calc_postprob[[2]][i,3]<-wt.sd(env_dat[site_ind,i+2],sam_probs[,3])
#}

#sam_calc_postprob[[3]]<-sam_calc_postprob[[2]]/sqrt(nrow(sp_200))
#names(sam_calc_postprob)<-c("mean", "sd", "se")

#sam_post_contents_SD<-dotplot_env_tab(mean_df = sam_calc_postprob$mean,
#                                      error_df = sam_calc_postprob$sd,
#                                      nGrp=3, env=env_vars, method="SAM")

#sam_post_contents_SE<-dotplot_env_tab(mean_df = sam_calc_postprob$mean,
#                                      error_df = sam_calc_postprob$se,
#                                      nGrp=3, env=env_vars, method="SAM")

#fix label switching issue
#sam_post_contents_SD$Group<-mapvalues(sam_post_contents_SD$Group ,from=c(1,2), to=c(2,1))
#sam_post_contents_SE$Group<-mapvalues(sam_post_contents_SE$Group ,from=c(1,2), to=c(2,1))

## Hard class SAM
sam_HC_vals<-env_match_vals( site_data= env_dat[site_ind,3:10], 
                pred_cluster_vals=sam_HC ,
                site_index=1:200)

sam_HC_contents_SD<-dotplot_env_tab(mean_df = sam_HC_vals$mean,
                                    error_df = sam_HC_vals$sd,
                                    nGrp=3, env=env_vars, method="SAM_Hard")

sam_HC_contents_SE<-dotplot_env_tab(mean_df = sam_HC_vals$mean,
                                      error_df = sam_HC_vals$se,
                                      nGrp=3, env=env_vars, method="SAM_Hard")

#dotplot of contents

contents_SD<-rbind(env_clust_contents_SD,
                   bio_clust_contents_SD, do.call(rbind, hclust_contents_SD),
                   rcp_HC_contents_SD, sam_HC_contents_SD, true_HC_contents_SD)

contents_SE<-rbind(env_clust_contents_SE,
                   bio_clust_contents_SE, do.call(rbind, hclust_contents_SE),
                   rcp_HC_contents_SE, sam_HC_contents_SE, true_HC_contents_SE)


#change to 'pretty' names
contents_SD$Method<-mapvalues(contents_SD$Method, from=c("Sp_RF", "nonbbGDM", "bbGDM","GF", "MNet", "HMSC"),
                              to=c("SpRF_HC", "GDM_HC","bbGDM_HC", "GF_HC", "MNet_HC", "HMSC_HC"))

contents_SD$Method<-contents_SD$Method<-factor(contents_SD$Method, 
                        levels=c("Env_Only", "BioHC_RF", "SpRF_HC", "HMSC_HC", "MNet_HC", "GDM_HC", 
                                 "bbGDM_HC", "GF_HC", "SAM_Hard", "RCP_Hard","True_Hard"))

contents_SD$env<-contents_SE$env<-factor(contents_SE$env, levels=c("temp", "O2",   "NO3"  , "sal" ,  "depth", "chla" , "ssh",   "curr" ))


contents_SE$Method<-mapvalues(contents_SE$Method, from=c("Sp_RF", "nonbbGDM", "bbGDM","GF", "MNet", "HMSC"),
                              to=c("SpRF_HC", "GDM_HC","bbGDM_HC", "GF_HC", "MNet_HC", "HMSC_HC"))


contents_SE$Method<-contents_SE$Method<-factor(contents_SE$Method, 
                                               levels=c("Env_Only", "BioHC_RF", "SpRF_HC", "HMSC_HC", "MNet_HC", "GDM_HC", 
                                                        "bbGDM_HC", "GF_HC", "SAM_Hard", "RCP_Hard", "True_Hard"))


p<-ggplot(data = contents_SD[contents_SD$env %in% c("temp", "O2", "NO3", "sal"),], 
#p<-ggplot(data = contents_SD,
          aes(x = Group, y = mean, ymin = lower, ymax = upper, colour = Method)) +
  geom_point(position = position_dodge(0.6), size=0.6) +
  geom_errorbar(position = position_dodge(0.6), width = 0.1) +
  coord_flip() +
  scale_colour_manual(name="Method", 
                      values = c("grey", "purple" ,"pink", "yellow"  , "orange", "green", "chartreuse4", "darkgreen",
                                 "cyan", "darkblue",  "red")) + #need to add enough colour for sampling levels here!
  #values=terrain.colors(14)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.y = element_text(face="italic"),
        legend.key = element_blank()) +
  facet_wrap( ~env, ncol=4, scales="free") 

tiff(file="Results/Plots/Grp_env_redSD.tiff", height=8, width=16, units="cm", res=1000)
p + scale_x_discrete(name="Group", breaks=c(1,2,3), limits=c("1","2", "3"))
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

tiff(file="Results/Plots/RF_partial_response.tiff", height=4, width =8, units="in", res=1000)
par(mfrow=c(2,4), mar=c(4,2,1,1), oma=c(4,5,1,1))

for(i in 1:length(env_vars)){
partialPlot(bio3_rf, as.data.frame(env_200[,2:9]), env_vars[i], main="",
            xlab=env_vars[i], ylab= "", col="orange1", lty=1, lwd=1.5, ylim=c(-0.9, 0.8))

partialPlot(bio3_rf, as.data.frame(env_200[,2:9]), env_vars[i], add=TRUE,ylim=c(-0.9, 0.8),
            which.class="2", xlab=env_vars[i], ylab= "", col="grey", lty=3, lwd=1.5, main="")

partialPlot(bio3_rf, as.data.frame(env_200[,2:9]), env_vars[i], add=TRUE,ylim=c(-0.9, 0.8),
            which.class="3", xlab=env_vars[i], ylab= "", col="darkolivegreen4", lty=2, lwd=1.5, main="")
}

title(ylab="Response", outer=TRUE, line=1)
dev.off()

## B) bbGDM----
bbgdm_response<-as.response(bbgdm_mod)

tiff(file="Results/Plots/bbGDM_PPlots.tiff", height = 4, width=4, units="in", res=1000)
par(mfrow=c(3,3), mar=c(4,4,2,2))
plot(bbgdm_response)
dev.off()

## C) Gradient Forest----
tiff(file="Results/Plots/GF_PPlots.tiff", height = 4, width=10, units="in", res=1000)
#par(oma=c(1,3,3,1))
plot(GF_mod, plot.type="Cumulative.Importance", plot.args=list(show.species=FALSE ))
#title(main="Gradient Forest", outer=TRUE)
dev.off()

## D) SAMs----
#generate values for plotting (varying one environmental variable at a time, holding others at mean value)
part_vals<-partial_values(env_vars = env_vars, raw_data = as.data.frame(env_dat[,3:10]))

#scale data usign same scaling as used to build models
scaling<-scale(env_dat[,3:10])

part_sc<-lapply(part_vals,scale, attr(scaling, "scaled:center"), attr(scaling, "scaled:scale"))

#predict to each set of variables and plot
tiff(file="Results/Plots/SAM_partial2.tiff", height=4, width=8, units="in", res=1000)

par(mfrow=c(2,4), mar=c(4,2,1,1), oma=c(4,5,1,1))
for(i in 1:length(part_vals)){
  temp<-as.data.frame(part_sc[[i]])
   pred<-predict(object=sam3_mod, temp)
   #account for label switching
   dimnames(pred$fit)[[2]]<- c(paste0("Group", c(2,1,3)))
    matplot(x=as.vector(part_vals[[i]][,names(part_vals[[i]]) %in% env_vars[i]]), 
            y= pred$fit, type='l', col=c("grey","darkolivegreen4", "orange1"), lty=c(3,2,1), xlab=env_vars[i], ylab=NULL, lwd=1.5, ylim=c(0,1))
    
    mtext(text="Probability of Occurrence", side=2, line=1, outer=TRUE)
    legend(x="topright", legend=dimnames(pred$fit)[[2]], col=c("grey","darkolivegreen4", "orange1"),lty=c(3,2,1), 
           cex=0.9, bty="n")
  }
dev.off()


dimnames(sam3_pred$fit)[[2]]<- c(paste0("Group", c(2,1,3)))

## E) RCPs----

#predict to each set of variables and plot
#rcpmod_vars<-c("temp", "O2", "sal")
tiff(file="Results/Plots/RCP_partial.tiff", height=4, width=8, units="in", res=1000)
par(mfrow=c(2,4), mar=c(4,2,1,1), oma=c(4,5,1,1))
for(i in 1:length(part_vals)){
  temp<-as.data.frame(part_sc[[i]])
  pred<-predict(rcp3_mod, newdata=temp)
  #account for label switching
  dimnames(pred)[[2]]<- c(paste0("Group", c(3,1,2)))
  matplot(x=as.vector(part_vals[[i]][,names(part_vals[[i]]) %in% env_vars[i]]), 
          y= pred, type='l', col=c( "orange1", "darkolivegreen4","grey"), lty=c(1,2,3), xlab=env_vars[i], ylab=NULL, lwd=1.5, ylim=c(0,1))
  
  #mtext(text="Probability of Occurrence", side=2, line=1, outer=TRUE)
  #legend(x="topright", legend=dimnames(pred)[[2]], col=c("orange1", "darkolivegreen4","grey"), lty=1, 
  #       cex=0.9, bty="n")
}
dev.off()


