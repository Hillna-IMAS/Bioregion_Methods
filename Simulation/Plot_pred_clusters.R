#############################################################################################
## Compare community modelling methods for bioregionalisation of simulated species dataset
## March 2018. N.Hill with input from S. Woolley
#############################################################################################

# Modelling process
# 1) Run models,diagnostics, determine number of groups or set number of groups to 3
# 2) PLOT PREDICTED DISTRIBUTION OF GROUPS ACROSS SIMULATION REGION
# 3) Describe content of groups 
# 4) Describe environment of groups
# 5) Predictor importance

#######################
## Set up---
#######################
# Get required libraries
library(plyr)           #data manipulation
library(raster)         #spatial data
library(RColorBrewer)   #colour palettes
library(rasterVis)      #plotting rasters

#Load required files
# most files in "simulate_communities/Many_covars" folder on dropbox
setwd("C:\\Users\\hillna\\UTAS_work\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\Simulation\\")
source("Simulation_Additional_Funcs.R")

load("Sim_Setup/Many_covars_sim.RData")
load("Sim_setup/sim_env_070518.RData") 
# env= raster brick of environmental data for entire region
# env_dat= raster brick of environmental data for entire region converted to matrix

load("Results/pred_clusters.RData")


#####################################################################################
#1) plot predicted distribution of optimal number of groups determined by models----
#a) hard classes for all methods except SAMs
#b) probablistic for RF, RCP, SAMs
#####################################################################################

## plot "True' distribution of simulated groups (using mean prevalence from simulation parameters)
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

true_rast<-rasterize(env_dat[,1:2], env, true_grps)
names(true_rast)<-paste0("Group", 1:3)
save(true_rast, true_grps, file="Sim_setup/true_dist.RData")

# set probability colour pallette
prob_pal<-colorRampPalette(c("#FFFFCC", "#FFEDA0", "#FED967", "#FEB24C",  "#FD8D3C",
                             "#FC4E2A", "#E31A1C", "#BD0026", "#800026"), space = "rgb")
#plot
tiff(filename="Results/Plots/TrueProb.tiff", compression="lzw", 
     width=8, height=6, units="cm", res=1000)
levelplot(true_rast, col.regions=prob_pal(19),
          at=seq(0,1,length=18),
          colorkey=list(at=seq(0,1,length=18),col=prob_pal(19)),
          names.att=paste("Group", 1:3, sep=" "),
          scales=list(draw=FALSE), layout=c(3,1))
          #,colorkey=list(title="probabilty of occurrence", title.gpar=list(las=0)))
dev.off()


### Look for label switching of classes in hard class outputs
# (i.e. class1 in Sp_RF is class 2 in HMSC). Make all relative to Sp_RF

mods<-c("Sp_RF", "HMSC", "MNet5_clust", "bbGDM", "nonbbGDM", "GF", "env")

check_clust2<-list()
for (i in 1: length(mods)){
  check_clust2[[i]]<- table(hard_cluster2$Sp_RF, hard_cluster2[,mods[i]])
}
names(check_clust2)<-mods

# attempt to make labels match up to true classes
hard_cluster2$Sp_RF<-mapvalues(hard_cluster2$Sp_RF ,from=c(1,2), to=c(1,3))
hard_cluster2$bbGDM<-mapvalues(hard_cluster2$bbGDM ,from=c(1,2), to=c(2,1))
hard_cluster2$nonbbGDM<-mapvalues(hard_cluster2$nonbbGDM ,from=c(1,2), to=c(2,1))
hard_cluster2$GF<-mapvalues(hard_cluster2$GF, from=c(1,2), to=c(1,3)) 
hard_cluster2$HMSC<-mapvalues(hard_cluster2$HMSC, from=c(1,2,3), to=c(2,1,3)) 
hard_cluster2$MNet5_clust<-mapvalues(hard_cluster2$MNet5_clust ,from=c(1,3,5), to=c(5,1,3))
#hard_cluster2$env<-mapvalues(hard_cluster2$env ,from=c(1,2,3,4,5,7,9,10), to=c(9,3,2,7,10,4,1,5))

### Create raster stack with hard class groups for plotting
#set up
clust2<-stack()
pretty_names<-c("SpRF_HC", "HMSC_HC", "MNet_HC", "GDM_HC", "bbGDM_HC", "GF_HC", "Env_Only")

#create factor attribute layer with as many levels as greatest number of clusters
rat<-data.frame(ID=1:11, Group=paste0("Group", 1:11))
#ignore warning in following loop

for (i in 1:length(mods)){
  
  hc_rast<-rasterize(env_dat[,1:2], env, field=hard_cluster2[,mods[i]])
  hc_rast<-as.factor(hc_rast)
  levels(hc_rast) <- rat
  
  clust2<-stack(clust2, hc_rast)
}
names(clust2)<-pretty_names

##set up plot colours for hard classes
class_pal<-c("darkolivegreen4","grey", "orange1", "darkred", "cadetblue", "cornsilk", 
             "goldenrod", "darkolivegreen3", "darkorange2" ,"deeppink3","deepskyblue3")

#plot simulation groups

#plot hard cluster outputs
tiff(filename="Results/Plots/many_covar_HClusts.tiff", compression="lzw", 
     width=11, height=11, units="cm", res=1000)
levelplot(clust2, layout=c(3,3),scales=list(draw=FALSE), 
          col.regions=class_pal)
dev.off()

##plot probability of occurrence outputs
#BioHC_RF

dimnames(bio2_rf_pred)[[2]]<-c("Group3", "Group1")

tiff(filename="Results/Plots/BioHC_RF.tiff", compression="lzw", 
     width=6, height=6, units="cm", res=1000)
levelplot(rasterize(env_dat[,1:2], env,field= bio2_rf_pred[,c("Group1", "Group3")]), 
          par.settings=YlOrRdTheme, names.attr=c("", ""),
          scales=list(draw=FALSE))
dev.off()

#SAM
#account for label switching
dimnames(sam3_pred$fit)[[2]]<- c(paste0("Group", c(2,1,3)))
dimnames(sam3_pred$se.fit)[[2]]<- c(paste0("Group", c(2,1,3)))

tiff(filename="Results/Plots/SAM.tiff", compression="lzw", 
     width=8, height=8, units="cm", res=1000)
levelplot(rasterize(env_dat[,1:2], env, field=sam3_pred$fit[,paste0("Group", 1:3)]), 
          col.regions=prob_pal(19),
          at=seq(0,1,length=18),
          colorkey=list(at=seq(0,1,length=18),col=prob_pal(19)),
          names.attr=rep("",3),layout=c(3,1),
          scales=list(draw=FALSE))
dev.off()

#RCP
dimnames(rcp3_pred[["ptPreds"]])[[2]]<- c(paste0("Group", c(3,1,2)))
dimnames(rcp3_pred[["bootSEs"]])[[2]]<- c(paste0("Group", c(3,1,2)))

tiff(filename="Results/Plots/RCP.tiff", compression="lzw", 
     width=8, height=6, units="cm", res=1000)
levelplot(rasterize(env_dat[,1:2], env,field=rcp3_pred[["ptPreds"]][,paste0("Group", 1:3)]), 
          col.regions=prob_pal(19),
          at=seq(0,1,length=18),
          colorkey=list(at=seq(0,1,length=18),col=prob_pal(19)),
          names.attr=rep(" ",3),
          scales=list(draw=FALSE), layout=c(3,1))
dev.off()

# Final plot assembled in graphics software

#####################################################################################
#2) plot 3 groups as expected from simulation set  up ----
#a) hard classes for all methods
#####################################################################################

#convert probabalistic methods to give hard classes
hard_cluster3$True<-apply(true_grps,1,which.max)
hard_cluster3$BioHC_RF<-apply(bio3_rf_pred,1,which.max)
hard_cluster3$SAM<-apply(sam3_pred$fit,1,which.max)
hard_cluster3$RCP<-apply(rcp3_pred[["ptPreds"]],1,which.max)

#fix label switching
mods2<-c("True",  "BioHC_RF", "Sp_RF", "HMSC", "MNet", "nonbbGDM",
         "bbGDM", "GF", "SAM", "RCP", "env")
pretty_names2<-c("True", "BioHC_RF", "SpRF_HC", "HMSC_HC", 
                 "MNet_HC", "GDM_HC", "bbGDM_HC", "GF_HC","SAM", "RCP" ,"Env_Only")


check_clust3<-list()
for (i in 1: length(mods2)){
  check_clust3[[i]]<- table(hard_cluster3$True, hard_cluster3[,mods2[i]])
}
names(check_clust3)<-mods2

hard_cluster3$env<-mapvalues(hard_cluster3$env ,from=c(1,3), to=c(3,1))
hard_cluster3$BioHC_RF<-mapvalues(hard_cluster3$BioHC_RF ,from=c(1,3), to=c(3,1))
hard_cluster3$Sp_RF<-mapvalues(hard_cluster3$Sp_RF ,from=c(1,2), to=c(2,1))
hard_cluster3$HMSC<-mapvalues(hard_cluster3$HMSC ,from=c(1,2,3), to=c(2,1,3))
hard_cluster3$MNet<-mapvalues(hard_cluster3$MNet ,from=c(1,3), to=c(3,1))
hard_cluster3$bbGDM<-mapvalues(hard_cluster3$bbGDM ,from=c(1,2), to=c(2,1))
hard_cluster3$GF<-mapvalues(hard_cluster3$GF ,from=c(1,2), to=c(2,1))
hard_cluster3$SAM<-mapvalues(hard_cluster3$SAM ,from=c(1,2), to=c(2,1))
hard_cluster3$RCP<-mapvalues(hard_cluster3$RCP ,from=c(1,2,3), to=c(3,1,2))

##plot groups
rat2<-data.frame(ID=1:3, Group=paste0("Group", 1:3))

clust3<-stack()

for (i in 1: length(mods2)){
  
  hc_rast<-rasterize(env_dat[,1:2], env, field=hard_cluster3[,mods2[i]])
  hc_rast<-as.factor(hc_rast)
  levels(hc_rast) <- rat2
  
  clust3<-stack(clust3, hc_rast)
}

names(clust3)<-pretty_names2

#plot results

tiff(file="Results/Plots/Cluster3_HCpredictions.tiff",compression="lzw", 
     width=14, height=14, units="cm", res=1000)
levelplot(clust3,layout=c(4,3),scales=list(draw=FALSE), 
          col.regions=class_pal[1:3])
dev.off()

#####################################################
### Plot uncertainty for SAMs and RCPS

grey_pal<-rev(gray.colors(n=19, start=0.4, end=1))

tiff(file="Results/Plots/SAM_RCP_SEs.tiff",compression="lzw", 
     width=14, height=10, units="cm", res=1000)
plot(levelplot(rasterize(env_dat[,1:2], env,field=sam3_pred[["se.fit"]][,paste0("Group", 1:3)]), 
    layout=c(3,1),at=seq(0,0.4,length=18), col.regions=grey_pal,
    colorkey=list(at=seq(0,0.4,length=18),col=grey_pal),
    scales=list(draw=FALSE), main="SAM"), 
     split=c(1,1,1,2))

plot(levelplot(rasterize(env_dat[,1:2], env,field=rcp3_pred[["bootSEs"]][,paste0("Group", 1:3)]), 
     layout=c(3,1), at=seq(0,0.4,length=18),col.regions=grey_pal,
     colorkey=list(at=seq(0,0.4,length=18),col=grey_pal),
     scales=list(draw=FALSE), main="RCP"), 
     split=c(1,2,1,2),newpage=FALSE)
dev.off()
# uncertainty for SAMs looks suspiscously similar between groups- check with Skip


