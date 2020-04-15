#############################################################################################
## Compare community modelling methods for bioregionalisation of simulated species dataset
## March 2018. N.Hill with input from S. Woolley
#############################################################################################

# Modelling process
# 1) Run models,diagnostics, determine number of groups or set number of groups to 3
# 2) Plot Predicted Distribution of groups across simulation region
# 3) DESCRIBE CONTENTS OF GROUPS
# 4) Describe environment of groups
# 5) Predictor importance

#######################
## Set up---
#######################
# Get required libraries
library(plyr)           #data manipulation
library(tidyr)          #data manipulation
library(ggplot2)        #plotting
library(SDMTools)       #weighted means and sd
library(rasterVis)      #plotting rasters

setwd("C:\\Users\\hillna\\UTAS_work\\Projects\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\Simulation\\")
source("Simulation_Additional_Funcs.R")

#Load required files

#load("Sim_Setup/Many_covars_sim.RData") 
load("Sim_Setup/Many_covars_sim_fin.RData") 
#sim_data= simulated species probabilities, occurrences, species' groups for entire region
#sp_200 = matrix of occurrences of 30 species at 200 sites to use as species dataset for analysis
#env_200= matrix of corresponding environmental conditions at 200 site to use as environmental dataset for analysis

load("Sim_setup/sim_env_070518.RData") 

#load("sim_env.RData") 
# env= raster brick of environmental data for entire region
# env_dat= raster brick of environmental data for entire region converted to matrix

load("Results/models.RData")
load("Results/pred_clusters.RData")
species=paste0("Sp", 1:30)
##########################################################################################
## Take note of label switching from plot_pred_clusters
## Tabulate or get directly from models contents of groups- Limit to case where groups =3
## generate dotchart of average and se of species' occurence in each group
##########################################################################################


# A) Cluster environment only- can't do
# 2 Stage Models- cluster then predict:
# B) cluster bioloigcal data, predict clusters with random forests- tabulate from bio clusters 
# 2 Stage Models- predict then cluster
# C) predict species using random forests, then cluster predictions- spatially match cluster prediction to sample & tabulate
# D) predict species using Mistnet, then cluster predictions- as above
# E) predict species using HMSC, then cluster predictions- as above
# F) predict dissimilarities using GDM, cluster predicted dissimilarities- as above
# G) predict biologically transformed environment with GF, cluster predictions- as above
# 1 Stage model-based: cluster and predict
# H) Species Archetype Models- direct from model
# I) Regions of Common Profile- direct from model


## Account for label switching---
# correct previously identified label switching in hard class methods

hard_cluster3$env<-mapvalues(hard_cluster3$env ,from=c(1,3), to=c(3,1))
hard_cluster3$Sp_RF<-mapvalues(hard_cluster3$Sp_RF ,from=c(1,2), to=c(2,1))
hard_cluster3$HMSC<-mapvalues(hard_cluster3$HMSC ,from=c(1,2,3), to=c(2,1,3))
hard_cluster3$MNet<-mapvalues(hard_cluster3$MNet ,from=c(1,3), to=c(3,1))
hard_cluster3$GDM_Dissim_HC<-mapvalues(hard_cluster3$GDM_Dissim_HC ,from=c(1,2), to=c(2,1))
hard_cluster3$GDM_TransEnv_HC<-mapvalues(hard_cluster3$GDM_TransEnv_HC ,from=c(1,2), to=c(2,1))
hard_cluster3$bbGDM_Dissim_HC<-mapvalues(hard_cluster3$bbGDM_Dissim_HC ,from=c(1,2), to=c(2,1))
hard_cluster3$bbGDM_TransEnv_HC<-mapvalues(hard_cluster3$bbGDM_TransEnv_HC ,from=c(1,2), to=c(2,1))
hard_cluster3$GF<-mapvalues(hard_cluster3$GF ,from=c(1,2), to=c(2,1))

#check labelling indeed fixed
class_pal<-c("darkolivegreen4","grey", "orange1")
clust3<-stack()
mods<-c( "Sp_RF","HMSC" , "MNet" , "GDM_Dissim_HC","GDM_TransEnv_HC",  
         "bbGDM_Dissim_HC"  , "bbGDM_TransEnv_HC",  "GF", "env")

#create factor attribute layer with as many levels as number of clusters
rat<-data.frame(ID=1:3, Group=paste0("Bioregion", 1:3))
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

## 'True' Distribution
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
#sites in simulation data= index of sites used in model building
site_true_probs<-true_grps[sites,]
site_true_HC<-apply(site_true_probs, 1, which.max)

## Hard class version of "true'
true_HC_vals<-get_match_vals( site_data= sp_200, 
                            pred_cluster_vals=site_true_HC ,
                            site_index=1:200)

true_HC_contents_SE<-dotplot_sp_tab(mean_df = true_HC_vals[[1]],
                                        error_df = true_HC_vals[[3]],
                                        nGrp=3, species=species, method="True_Hard")



### 2 stage methods: cluster biology then predict. Using RF probabalistic output

bio3_clust<-mapvalues(bio3_clust,from=c(1,3), to=c(3,1))

bio_clust_vals<-get_match_vals( site_data= sp_200, 
                          pred_cluster_vals=bio3_clust ,
                          site_index=1:200)


bio_clust_contents_SE<-dotplot_sp_tab(mean_df = bio_clust_vals[[1]],
                                      error_df = bio_clust_vals[[3]],
                                      nGrp=3, species=species, method="BioHC_RF")



### 2 stage methods: predict then heirarchical cluster
names(hard_cluster3)[4:7]<-c("SpRF_HC", "GF_HC", "MNet_HC", "HMSC_HC")
clusts<-names(hard_cluster3)[4:11]
hclust_contents_SE<-list()
 
for( i in seq_along(clusts)){
  clust_vals<-get_match_vals( site_data= sp_200, 
                              pred_cluster_vals=hard_cluster3[,clusts[i]] ,
                              site_index=sites)
  
  hclust_contents_SE[[i]]<-dotplot_sp_tab(mean_df = clust_vals[[1]],
                                          error_df = clust_vals[[3]],
                                          nGrp=3, species=species, method=clusts[i])
}


## Species Archetype Models----
# get spatial group predictions and match initial sites
sam_probs<-sam3_pred$ptPreds[sites,]
sam_HC <- apply(sam_probs, 1, which.max)


sam_HC_vals<-get_match_vals( site_data= sp_200, 
                            pred_cluster_vals=sam_HC ,
                            site_index=1:200)

sam_HC_contents_SE<-dotplot_sp_tab(mean_df = sam_HC_vals[[1]],
                                        error_df = sam_HC_vals[[3]],
                                        nGrp=3, species=species, method="SAM_Hard")
sam_HC_contents_SE$Group<-mapvalues(sam_HC_contents_SE$Group, from=c(1,2,3), to=c(2,3,1))

## Region of Common Profile Models----
#from model parameters
rcp_calc_contents<-calc_prev(boot_obj = rcp3_boot, mod_obj=rcp3_mod)

rcp_contents_SD<-dotplot_sp_tab (mean_df=t(rcp_calc_contents$mean),
                      error_df= t(rcp_calc_contents$sd),
                      nGrp=3, species= species, method="RCP_Coefs")


#take account of label switching
rcp_contents_SD$Group<-mapvalues(rcp_contents_SD$Group ,from=c(1,2,3), to=c(3,1,2))

rcp_contents_SE<-rcp_contents_SD

#hard classes
#rcp_HC<-apply(rcp3_pred$ptPreds[site_ind,], 1, which.max)
rcp_HC<-apply(rcp3_pred$ptPreds[sites,], 1, which.max)
rcp_HC<-mapvalues(rcp_HC ,from=c(1,2,3), to=c(3,1,2))

rcp_HC_vals<-get_match_vals( site_data= sp_200, 
                            pred_cluster_vals=rcp_HC ,
                            site_index=1:200)
rcp_HC_contents_SE<-dotplot_sp_tab(mean_df = rcp_HC_vals[[1]],
                                        error_df = rcp_HC_vals[[3]],
                                        nGrp=3, species=species, method="RCP_Hard")




#####################################################
## Generate dotplot to compare contents of groups----
#####################################################

#dotplot of contents
contents_SE<-rbind(true_HC_contents_SE, 
                   bio_clust_contents_SE, do.call(rbind, hclust_contents_SE),
                   rcp_contents_SE, rcp_HC_contents_SE, sam_HC_contents_SE)

#Tidy up names etc
contents_SE$Method<-factor(contents_SE$Method, 
                                               levels=c("BioHC_RF", "SpRF_HC", "HMSC_HC", "MNet_HC", 
                                                        "GDM_Dissim_HC", "GDM_TransEnv_HC", "bbGDM_Dissim_HC", "bbGDM_TransEnv_HC", 
                                                        "GF_HC", "SAM_Hard",  "RCP_Coefs","RCP_Hard", "True_Hard"))
contents_SE$species<-factor(contents_SE$species, 
                           levels=paste0("Sp", 1:30))

names(contents_SE)[2]<-"Prevalence"
names(contents_SE)[5]<-"Species"
contents_SE$Bioregion<-factor(contents_SE$Group, levels=c('1','2','3'),
                              labels=c("Bioregion_1", "Bioregion_2", "Bioregion_3"))


#truncate values to 0,1
contents_SE$lower<-ifelse(contents_SE$lower <0, 0,contents_SE$lower)
contents_SE$upper<-ifelse(contents_SE$upper >1, 1,contents_SE$upper)


p<-ggplot(data = contents_SE[contents_SE$Species %in% species[1:5],], 
#p<-ggplot(data = contents_SE,
          aes(x = Species, y = Prevalence, ymin = lower, ymax = upper, colour = Method)) +
  scale_y_continuous(limits = c(-0.1,1.1)) +
  geom_point(position = position_dodge(0.6), size=0.6) +
  geom_errorbar(position = position_dodge(0.6), width = 0.5) +
  coord_flip() +
  scale_colour_manual(name="Method", 
                      values = c("purple" ,"pink", "yellow"  , "orange", 
                                 "chartreuse1", "chartreuse3", "darkseagreen2", "darkgreen", "darkslategray",
                                 "cyan", "cornflowerblue", "darkblue", "red"))+ #need to add enough colour for sampling levels here!
  #values=terrain.colors(14)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.y = element_text(face="italic"),
        legend.key = element_blank()) +
  facet_wrap( ~Bioregion, ncol=4, scales="free") 

tiff(file="Results/Plots/Grp_sp_SE.tiff", height=10, width=18, units="cm", res=1000)
p 
dev.off()




