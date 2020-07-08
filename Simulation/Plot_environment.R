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
library(tidyr)          #data manipulation
library(plyr)           #data manipulation
library(ggplot2)        #plotting
#library(raster)        #spatial data
#library(RColorBrewer)  #colour palettes
library(rasterVis)      #plotting rasters
library(RCPmod)         #Regions of Common Profile
#library(bbgdm)         #Bayesian Bootstrap Generalised Dissimilarity Models (and naive GDM)
library(gdm)            #Generalised Dissimilarity Modelling
library(gradientForest) #Gradient Forests (an extension of Random Forests)
library(ecomix)     #Species Archetype Models (SAMs)
library(SDMTools)       #weighted means and sd

setwd("C:\\Users\\hillna\\UTAS_work\\Projects\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\Simulation\\")
source("Simulation_Additional_Funcs.R")

#Load required files

load("Sim_setup/Many_covars_sim_fin.RData")
#load("Sim_setup/Many_covars_sim.RData")
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
## Generate partial plots or equivalent of the group response to the environment for RF, SAM & RCP
############################################################################################################


# A) Cluster environment only- directly from cluster results for entire region
# 2 Stage Models- cluster then predict:
# B) cluster bioloigcal data, predict clusters with random forests- tabulate from bio clusters
# 2 Stage Models- predict then cluster
# C) predict species using random forests, then cluster predictions- spatially match cluster prediction to sample & tabulate environment
# D) predict species using Mistnet, then cluster predictions- as above
# E) predict species using HMSC, then cluster predictions- as above
# F) i-predict dissimilarities using GDM, cluster predicted dissimilarities; ii- cluster biologically transformed environment
# G) predict biologically transformed environment (GF), cluster prediction
# 1 Stage mdoel-based: cluster and predict
# H) Species Archetype Models- for ease of comparison, convert results to hard classes and use same process as for 2-stage methods,  plus partial plots
# I) Regions of Common Profile-for ease of comparison, convert results to hard classes and use same process as for 2-stage methods, plus partial plots


## Account for label switching---
# correct previously identified label switching 
hard_cluster3$env<-mapvalues(hard_cluster3$env ,from=c(1,3), to=c(3,1))
hard_cluster3$Sp_RF<-mapvalues(hard_cluster3$Sp_RF ,from=c(1,2), to=c(2,1))
hard_cluster3$HMSC<-mapvalues(hard_cluster3$HMSC ,from=c(1,2,3), to=c(2,1,3))
hard_cluster3$MNet<-mapvalues(hard_cluster3$MNet ,from=c(1,3), to=c(3,1))
hard_cluster3$GDM_Dissim_HC<-mapvalues(hard_cluster3$GDM_Dissim_HC ,from=c(1,2), to=c(2,1))
hard_cluster3$GDM_TransEnv_HC<-mapvalues(hard_cluster3$GDM_TransEnv_HC ,from=c(1,2), to=c(2,1))
hard_cluster3$bbGDM_Dissim_HC<-mapvalues(hard_cluster3$bbGDM_Dissim_HC ,from=c(1,2), to=c(2,1))
hard_cluster3$bbGDM_TransEnv_HC<-mapvalues(hard_cluster3$bbGDM_TransEnv_HC ,from=c(1,2), to=c(2,1))
hard_cluster3$GF<-mapvalues(hard_cluster3$GF ,from=c(1,2), to=c(2,1))


# Quickly check labelling indeed corrected
class_pal<-c("darkolivegreen4","grey", "orange1")
clust3<-stack()
mods<-c( "Sp_RF","HMSC" , "MNet" , "GDM_Dissim_HC","GDM_TransEnv_HC",  
         "bbGDM_Dissim_HC"  , "bbGDM_TransEnv_HC",  "GF", "env")

#create factor attribute layer with as many levels as greatest number of clusters
rat<-data.frame(ID=1:11, Group=paste0("Bioregion", 1:11))
#ignore warning in following loop

for (i in 1:length(mods)){
  
  hc_rast<-rasterize(env_dat[,1:2], env, field=hard_cluster3[,mods[i]])
  hc_rast<-as.factor(hc_rast)
  levels(hc_rast) <- rat
  
  clust3<-stack(clust3, hc_rast)
}
names(clust3)<-mods
x11()
levelplot(clust3, col.regions=class_pal)

## 'True' Distribution ----
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
#sites= index of randomly selected sites for model building
site_true_probs<-true_grps[sites,]
site_true_HC<-apply(site_true_probs, 1, which.max)


## Hard classed version of "true"
true_HC_vals<-env_match_vals( site_data= env_dat[sites,3:10], 
                              pred_cluster_vals=site_true_HC ,
                              site_index=sites)
true_HC_contents_SD<-dotplot_env_tab(mean_df = true_HC_vals[[1]],
                                     error_df = true_HC_vals[[2]],
                                     nGrp=3, env=env_vars, method="Truth (Hard Class)")

true_HC_contents_SE<-dotplot_env_tab(mean_df = true_HC_vals[[1]],
                                    error_df = true_HC_vals[[3]],
                                    nGrp=3, env=env_vars, method="Truth(Hard Class)")


## Environment Only ----
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

### 2 stage methods: cluster biology then predict ----

## Hard Classed BioHC_RF
bio3_clust<-mapvalues(bio3_clust,from=c(1,3), to=c(3,1))

bio_clust_vals<-env_match_vals(site_data= env_dat[sites,3:10], 
                               pred_cluster_vals=bio3_clust ,
                               site_index=sites)

bio_clust_contents_SD<-dotplot_env_tab(mean_df = bio_clust_vals$mean,
                                   error_df = bio_clust_vals$sd,
                                   nGrp=3, env=env_vars, method="BioHC_RF (Hard Class)")

bio_clust_contents_SE<-dotplot_env_tab(mean_df = bio_clust_vals$mean,
                                       error_df = bio_clust_vals$se,
                                       nGrp=3, env=env_vars, method="BioHC_RF (Hard Class)")


### 2 stage methods: predict then heirarchical cluster----
clusts<-names(hard_cluster3)[4:11]
hclust_contents_SD<-list()
hclust_contents_SE<-list()

for( i in seq_along(clusts)){
  clust_vals<-env_match_vals( site_data= env_dat[sites,3:10], 
                              pred_cluster_vals=hard_cluster3[sites,clusts[i]] ,
                              site_index=sites)
  
  hclust_contents_SD[[i]]<-dotplot_env_tab(mean_df = clust_vals$mean,
                                       error_df = clust_vals$sd,
                                       nGrp=3, env=env_vars, method=clusts[i])
  
  hclust_contents_SE[[i]]<-dotplot_env_tab(mean_df = clust_vals$mean,
                                           error_df = clust_vals$se,
                                           nGrp=3, env=env_vars, method=clusts[i])
}
 

##RCP----

## Hard class RCPs
rcp_HC<-apply(rcp3_pred$ptPreds, 1, which.max)
rcp_HC<-mapvalues(rcp_HC ,from=c(1,2,3), to=c(3,1,2))


rcp_HC_vals<-env_match_vals( site_data= env_dat[sites,3:10], 
                pred_cluster_vals=rcp_HC[sites] ,
                site_index=1:200)

rcp_HC_contents_SD<-dotplot_env_tab(mean_df = rcp_HC_vals$mean,
                                     error_df = rcp_HC_vals$sd,
                                     nGrp=3, env=env_vars, method="RCP (Hard Class)")

rcp_HC_contents_SE<-dotplot_env_tab(mean_df = rcp_HC_vals$mean,
                                     error_df = rcp_HC_vals$se,
                                     nGrp=3, env=env_vars, method="RCP (Hard Class)")

##SAMs----
# get spatial group predictions and match initial sites
sam_probs<-sam3_pred$ptPreds[sites,]
sam_HC<-apply(sam_probs, 1, which.max)
sam_HC<-mapvalues(sam_HC,from=c(1,2,3), to=c(2,3,1))


## Hard class SAM
sam_HC_vals<-env_match_vals( site_data= env_dat[sites,3:10], 
                pred_cluster_vals=sam_HC ,
                site_index=1:200)

sam_HC_contents_SD<-dotplot_env_tab(mean_df = sam_HC_vals$mean,
                                    error_df = sam_HC_vals$sd,
                                    nGrp=3, env=env_vars, method="SAM (Hard Class)")

sam_HC_contents_SE<-dotplot_env_tab(mean_df = sam_HC_vals$mean,
                                      error_df = sam_HC_vals$se,
                                      nGrp=3, env=env_vars, method="SAM (Hard Class)")

## Dotplot of contents----
#create blank row/separator for each approach (clunky way of getting ggplot2 to have the legend I want!)
gtp<-dotplot_env_tab(mean_df = matrix(NA,nrow=8, ncol=3),
                    error_df = matrix(0,nrow=8, ncol=3),
                    nGrp=3, env=env_vars, method="Two-stage: Group then Predict")
ptg<-dotplot_env_tab(mean_df = matrix(NA,nrow=8, ncol=3),
                    error_df = matrix(0,nrow=8, ncol=3),
                    nGrp=3, env=env_vars, method="Two-stage: Predict then Group")
os<-dotplot_env_tab(mean_df = matrix(NA,nrow=8, ncol=3),
                   error_df = matrix(0,nrow=8, ncol=3),
                   nGrp=3, env=env_vars, method="One-stage")

contents_SD<-rbind(env_clust_contents_SD, gtp, ptg, os,
                   bio_clust_contents_SD, do.call(rbind, hclust_contents_SD),
                   rcp_HC_contents_SD, sam_HC_contents_SD, true_HC_contents_SD)

#contents_SE<-rbind(env_clust_contents_SE,
#                   bio_clust_contents_SE, do.call(rbind, hclust_contents_SE),
#                   rcp_HC_contents_SE, sam_HC_contents_SE, true_HC_contents_SE)


#change to 'pretty' names
contents_SD$Method<-mapvalues(contents_SD$Method, from=c("Sp_RF", "GF", "MNet", "HMSC" ),
                              to=c("SpRF_HC", "GF_HC", "MNet_HC", "HMSC_HC"))

contents_SD$Method<-contents_SD$Method<-factor(contents_SD$Method, 
                        levels=c("RCP (Hard Class)", "SAM (Hard Class)",  "One-stage", 
                                        "GF_HC",  "bbGDM_Dissim_HC", "bbGDM_TransEnv_HC", "GDM_TransEnv_HC",  "GDM_Dissim_HC",
                                        "MNet_HC", "HMSC_HC", "SpRF_HC", "Two-stage: Predict then Group",
                                        "BioHC_RF (Hard Class)",  "Two-stage: Group then Predict", "Env_Only",
                                        "Truth (Hard Class)"))

contents_SD$env<-factor(contents_SD$env, levels=c("temp", "O2",   "NO3"  , "sal" ,  "depth", "chla" , "ssh",   "curr" ))


#contents_SE$Method<-mapvalues(contents_SE$Method, from=c("Sp_RF", "GF", "MNet", "HMSC" ),
#                              to=c("SpRF_HC", "GF_HC", "MNet_HC", "HMSC_HC"))


#contents_SE$Method<-contents_SE$Method<-factor(contents_SE$Method, 
#                                               levels=c("RCP (Hard Class)", "SAM (Hard Class)",  "One-stage", 
#                                                        "GF_HC",  "bbGDM_Dissim_HC", "bbGDM_TransEnv_HC", "GDM_TransEnv_HC",  "GDM_Dissim_HC",
#                                                        "MNet_HC", "HMSC_HC", "SpRF_HC", "Two-stage: Predict then Group",
#                                                        "BioHC_RF (Hard Class)",  "Two-stage: Group then Predict", "Env_Only",
#                                                        "Truth (Hard Class)"))

names(contents_SD)[1:2]<-c("Bioregion","Value")


#p<-ggplot(data = contents_SD[contents_SD$env %in% c("temp", "O2", "NO3", "sal"),], 
p<-ggplot(data = contents_SD,
          aes(x = Bioregion, y = Value, ymin = lower, ymax = upper, colour = Method)) +
  geom_point(position = position_dodge(0.6), size=0.6) +
  geom_errorbar(position = position_dodge(0.6), width = 0.1) +
  coord_flip() +
  scale_colour_manual(name="Method", 
                          values=c( "darkblue","cyan",  "white", 
                         "darkslategray", "darkgreen", "darkseagreen2","chartreuse3", "chartreuse1",
                         "orange", "yellow"  , "pink", "white",
                         "purple" , "white", "grey", "red"))+
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.y = element_text(face="italic"),
        legend.key = element_blank()) +
  facet_wrap( ~env, ncol=4, scales="free") 

#tiff(file="Results/Plots/Grp_env_redSD.tiff", height=12, width=20, units="cm", res=1000)
tiff(file="Results/Plots/Grp_env_SD.tiff", height=16, width=20, units="cm", res=1000)
p + 
  scale_x_discrete(name="Bioregion", breaks=c(1,2,3), limits=c("1","2", "3")) +
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

tiff(file="Results/Plots/RF_partial_response.tiff", height=7, width =8, units="in", res=1000)
par(mfrow=c(3,3), mar=c(4,2,1,1), oma=c(4,5,1,1))

for(i in 1:length(env_vars)){
partialPlot(bio3_rf, as.data.frame(env_200[,2:9]), env_vars[i], main="",
            xlab=env_vars[i], ylab= "", col="orange1", lty=1, lwd=1.5, ylim=c(-0.9, 0.8))

partialPlot(bio3_rf, as.data.frame(env_200[,2:9]), env_vars[i], add=TRUE,ylim=c(-0.9, 0.8),
            which.class="2", xlab=env_vars[i], ylab= "", col="grey", lty=1, lwd=1.5, main="")

partialPlot(bio3_rf, as.data.frame(env_200[,2:9]), env_vars[i], add=TRUE,ylim=c(-0.9, 0.8),
            which.class="3", xlab=env_vars[i], ylab= "", col="darkolivegreen4", lty=1, lwd=1.5, main="")
}
plot.new()
legend(x="center", legend=paste0("Bioregion", 1:3), col=c("darkolivegreen4","grey", "orange1"), lty=1, 
       cex=1.2, bty="n")

title(ylab="Response", outer=TRUE, line=1)
dev.off()


## D) SAMs----
#generate values for plotting (varying one environmental variable at a time, holding others at mean value)
part_vals<-partial_values(env_vars = env_vars, raw_data = as.data.frame(env_dat[,3:10]))

#scale data usign same scaling as used to build models
scaling<-scale(env_dat[,3:10])

part_sc<-lapply(part_vals,scale, attr(scaling, "scaled:center"), attr(scaling, "scaled:scale"))

#predict to each set of variables and plot
#note label switching
#hard_cluster3$SAM<-mapvalues(hard_cluster3$SAM ,from=c(1,2,3), to=c(2,3,1))

tiff(file="Results/Plots/SAM_partial2.tiff", height=7, width=8, units="in", res=1000)

par(mfrow=c(3,3), mar=c(4,2,1,1), oma=c(4,5,1,1))
for(i in 1:length(part_vals)){
  temp<-as.data.frame(part_sc[[i]])
   pred<-predict(object=sam3_mod, newdata=temp)
   #account for label switching
   dimnames(pred)[[2]]<- c(paste0("Bioregion", c(2,3,1)))
   pred<-pred[,c(paste0("Bioregion", 1:3))]
   matplot(x=as.vector(part_vals[[i]][,names(part_vals[[i]]) %in% env_vars[i]]), 
           y= pred, type='l', col=c("darkolivegreen4", "grey", "orange1"), lty=c(1,1,1), xlab=env_vars[i], ylab=NULL, lwd=1.5, ylim=c(0,1))
}
plot.new()
mtext(text="Probability of Occurrence", side=2, line=1, outer=TRUE)
legend(x="center", legend=paste0("Bioregions", 1:3),  col=c("darkolivegreen4", "grey", "orange1"), lty=c(1,1,1), 
       cex=1.2, bty="n")
dev.off()

## E) RCPs----

#predict to each set of variables and plot
rcpmod_vars<-c("temp", "O2", "sal")
rcpmod_part_sc<-part_sc[names(part_sc)%in% rcpmod_vars]
rcpmod_part_vals<-part_vals[names(part_vals)%in% rcpmod_vars]
tiff(file="Results/Plots/RCP_partial.tiff", height=3, width=9, units="in", res=1000)
par(mfrow=c(1,4), mar=c(4,2,1,1), oma=c(4,5,1,1))
for(i in 1:length(rcpmod_vars)){
  temp<-as.data.frame(rcpmod_part_sc[[i]])
  pred<-predict(rcp3_mod, newdata=temp)
  #account for label switching
  dimnames(pred)[[2]]<- c(paste0("Bioregion", c(3,1,2)))
  xvals<-rcpmod_part_vals[[i]]
  matplot(x=as.vector(xvals[,names(xvals) %in% rcpmod_vars[i]]), 
         y= pred, type='l', col=c( "orange1", "darkolivegreen4","grey"), lty=c(1,1,1),
         xlab=rcpmod_vars[i], cex.lab=1.5, ylab=NULL, lwd=1.5, ylim=c(0,1))
}
plot.new()
legend(x="center", legend=dimnames(pred)[[2]], col=c("orange1", "darkolivegreen4","grey"), lty=1, 
       cex=1.2, bty="n")
mtext(text="Probability of Occurrence", side=2, line=0.9, outer=TRUE)
dev.off()

#plot RCP model pi's against environmental values.
plot_data<-as.data.frame(cbind(env_dat[sites, env_vars[1:4]],
                               rcp3_mod$pis))
plot_data<-as.data.frame(cbind(env_dat[sites, env_vars[1:4]],
                               rcp3_mod$postProbs))

names(plot_data)[5:7]<-paste0("RCP", c(3,1,2))
RCPS<-paste0("RCP", 1:3)
rcp_vars<-c( "temp", "O2", "sal")

colours=data.frame(color=c( "darkolivegreen4", "grey", "orange1"), poly=NA)
colours$poly<-adjustcolor(colours$color, alpha.f=0.35)

tiff(file="Results/Plots/RCP_pis.tiff", height=6, width=6, units="in", res=500)
par(mfcol=c(3,3), mar=c(1,1.5,2,2), oma=c(4,5,0,1))

for(i in 1:length(rcp_vars)){
  for (j in 1:length(RCPS)){
    plot(plot_data[,rcp_vars[i]], plot_data[,RCPS[j]], pch=21, ylim=c(0,1), cex=1.2, col=colours$colour[j], bg=colours$poly[j])
  }
}

#xlabels
mtext("Temp", at= 0.16, side=1, outer=TRUE, cex=0.9, line=2)
mtext("02", at= 0.5, side=1, outer=TRUE, cex=0.9, line=2)
mtext("sal", at= 0.82, side=1, outer=TRUE, cex=0.9, line=2)

#ylabels
mtext("Bioregion 1", at= 0.15, side=2, outer=TRUE, cex=0.9, line=2)
mtext("Bioregion 2", at= 0.48, side=2, outer=TRUE, cex=0.9, line=2)
mtext("Bioregion 3", at= 0.82, side=2, outer=TRUE, cex=0.9, line=2)
dev.off()


###plot environmental values and responses for temp, NO3, O2 and sal on one plot----
partial_preds <- list()
for (i in 1:4){
  temp<-as.data.frame(part_sc[[i]])
  part_preds_varible <- matrix(NA,nrow(temp),9)
  part_preds_varible[,1] <- partialPlot(bio3_rf, as.data.frame(env_200[,2:9]), env_vars[i], main="",which.class = '1',n.pt = 100,
                                        xlab=env_vars[i], ylab= "", col="orange1", lty=1, lwd=1.5, ylim=c(-0.9, 0.8))$y
  part_preds_varible[,2] <- partialPlot(bio3_rf, as.data.frame(env_200[,2:9]), env_vars[i], main="",which.class = '2',n.pt = 100,
                                        xlab=env_vars[i], ylab= "", col="orange1", lty=1, lwd=1.5, ylim=c(-0.9, 0.8))$y
  part_preds_varible[,3] <- partialPlot(bio3_rf, as.data.frame(env_200[,2:9]), env_vars[i], main="",which.class = '3',n.pt = 100,
                                        xlab=env_vars[i], ylab= "", col="orange1", lty=1, lwd=1.5, ylim=c(-0.9, 0.8))$y
  ## RF
  part_preds_varible[,4:6] <- predict(sam3_mod, newdata=temp) ## SAM
  part_preds_varible[,7:9] <- predict(rcp3_mod, newdata=temp)## RCP
  
  #fix label switchin here
  colnames(part_preds_varible) <- c(paste0("BioHC_RF_", c(3,2,1)),
                                    paste0("SAM_", c(2,3,1)),
                                    paste0("RCP_", c(3,1,2)))
  
  partial_preds[[i]] <- part_preds_varible 
}

library(reshape2)
part_preds_melted <- lapply(partial_preds,melt)
allpartials <- do.call('rbind',part_preds_melted)

allpartials$ID <- paste0(rep(c('BioHC_RF','SAM','RCP'),each=300),rep(env_vars[1:4], each = 900 ))
allpartials$ID <- factor(allpartials$ID, levels= c("BioHC_RFtemp","BioHC_RFO2", "BioHC_RFNO3",  "RFsal","SAMtemp","SAMO2","SAMNO3","SAMsal",
                                                   "RCPtemp","RCPO2","RCPNO3","RCPsal"))
#allpartials$grp_col <- rep(paste0('grp',rep(1:3,each=100)),12)
allpartials$grp_col <- sub(".*_", "", allpartials$Var2)
allpartials$env_id <- rep(env_vars[1:4], each = 900 )
allpartials$mod_id <- rep(rep(c('BioHC_RF','SAM','RCP'),each=300),4)
allpartials$env_id <- factor(allpartials$env_id,levels=env_vars[1:4])
allpartials$mod_id <- factor(allpartials$mod_id,levels=c('BioHC_RF','SAM','RCP'))


temp <- as.vector(part_vals[[1]][,names(part_vals[[1]]) %in% env_vars[1]])
O2 <- as.vector(part_vals[[2]][,names(part_vals[[2]]) %in% env_vars[2]])
NO3 <- as.vector(part_vals[[3]][,names(part_vals[[3]]) %in% env_vars[3]])
sal <- as.vector(part_vals[[4]][,names(part_vals[[4]]) %in% env_vars[4]])

allpartials$env <- c(rep(temp,9),rep(O2,9),rep(NO3,9),rep(sal,9))
 
head(allpartials)


library(ggplot2)
p2 <- ggplot(data=allpartials, aes(x=env, y = value, colour = grp_col)) + 
  geom_line() + 
  facet_grid(mod_id~env_id,scale='free',switch='x')+
  scale_color_manual(values=c("darkolivegreen4","grey","orange1"),
                     breaks=c("grp1", "grp2", "grp3"),
                     labels=c("Bioregion 1", "Bioregion 2", "Bioregion 3"))+
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(face = 'bold'), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = 'black',fill = NA))+
  labs(title = "B)", y = "Response", color = "")+
  expand_limits(y=c(0,1))


library(grid)
library(gridExtra)
tiff(file='Results/Plots/Fig7.tif',height = 11, width = 8,units= 'in', res = 100)
grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(8, 4)))
print(p,vp=vplayout(1:4,1:4))
print(p2,vp=vplayout(5:8,1:4))
dev.off()


