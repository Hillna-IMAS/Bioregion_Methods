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

setwd("C:\\Users\\hillna\\UTAS_work\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\KP_Fish\\")
source("Code/Additional_Funcs.R")

#Load required files
load("Results/pred_clusters.RData")
load("Results/models.RData")
KP_Fish_Env<-readRDS("Data/KP_Fish_Env.RDS")
load("Data/Prediction_space.Rda")

species=c("Antimora.rostrata" ,"Bathydraco.antarcticus", "Bathyraja.eatonii",                
          "Bathyraja.irrasa"  ,  "Bathyraja.murrayi",  "Champsocephalus.gunnari" , "Dissostichus.eleginoides" ,        
          "Etmopterus.viator"  ,  "Gobionotothen.acuta" , "Lepidonotothen.mizops"   ,  "Lepidonotothen.squamifrons",       
          "Lycodapus.antarcticus" ,"Macrourus.spp", "Mancopsetta.maculata"  , "Muraenolepis.spp",                
          "Notothenia.rossii"   , "Paradiplospinus.gracilis" ,  "Paraliparis.spp." ,"Zanclorhynchus.spinifer" ,         
          "Channichthys.rhinoceratus.velifer")

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
#convert probabalistic methods to give hard classes
hard_cluster4$BioHC_RF<-apply(bio4_rf_pred,1,which.max)
hard_cluster4$RCP_Hard<-apply(rcp4_spat_preds[["ptPreds"]],1,which.max)
hard_cluster4$SAM_Hard<-apply(sam4_pred$fit,1,which.max)

#fix label switching
hard_cluster4$SAM_Hard<-mapvalues(hard_cluster4$SAM ,from=c(1,3), to=c(2,1))
hard_cluster4$RCP_Hard<-mapvalues(hard_cluster4$RCP ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$SpRF_HC<-mapvalues(hard_cluster4$SpRF_HC ,from=c(1,2,3,4), to=c(2,4,1,3))
hard_cluster4$GDM_HC<-mapvalues(hard_cluster4$GDM_HC ,from=c(1,2,3,4), to=c(1,3,2,4))
hard_cluster4$GF_HC<-mapvalues(hard_cluster4$GF_HC ,from=c(1,2,3,4), to=c(3,2,4,1))
hard_cluster4$MNet_HC<-mapvalues(hard_cluster4$MNet_HC ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$HMSC_HC<-mapvalues(hard_cluster4$HMSC_HC ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$BioHC_RF<-mapvalues(hard_cluster4$BioHC_RF ,from=c(1,2,3,4), to=c(1,4,2,3))

#create raster stack of predictions and extract survey site values
clusts<-names(hard_cluster4)[4:12]

##plot groups
rat2<-data.frame(ID=1:4, Group=paste0("Group", 1:4))

clust4<-stack()

for (i in 1: length(clusts)){
  
  hc_rast<-rasterize(na.omit(pred_sp)[,1:2], env_raster, field=hard_cluster4[,clusts[i]])
  hc_rast<-as.factor(hc_rast)
  levels(hc_rast) <- rat2
  
  clust4<-stack(clust4, hc_rast)
}

names(clust4)<-clusts

site_classes<-as.data.frame(raster::extract(clust4, KP_Fish_Env[,5:6]))



### 2 stage methods: predict then heirarchical cluster
hclust_contents_SE<-list()
 
for( i in 1:8){
  clust_vals<-get_match_vals( site_data= KP_Fish_Env[,species], 
                              pred_cluster_vals=site_classes[,clusts[i]] ,
                              site_index=1:nrow(KP_Fish_Env))
  
  hclust_contents_SE[[i]]<-dotplot_sp_tab(mean_df = clust_vals[[1]],
                                          error_df = clust_vals[[3]],
                                          nGrp=4, species=species, method=clusts[i])
}


## Region of Common Profile Models----
#from model parameters
rcp_calc_contents<-calc_prev(boot_obj = rcp4_boot, mod_obj = rcp4_mod, 
                             samp_fact = levels (as.factor(KP_Fish_Env$Survey)), 
                             calc_level = "overall")

rcp_contents_SD<-dotplot_sp_tab (mean_df=t(rcp_calc_contents$mean),
                      error_df= t(rcp_calc_contents$sd),
                      nGrp=4, species= species, method="RCP_coefs")

#take account of label switching
rcp_contents_SD$Group<-mapvalues(rcp_contents_SD$Group, from=c(1,2,3,4), to=c(2,3,4,1))

rcp_contents_SE<-rcp_contents_SD


#####################################################
## Generate dotplot to compare contents of groups----
#####################################################

#dotplot of contents
contents_SE<-rbind( do.call(rbind, hclust_contents_SE),
                   rcp_contents_SE)

contents_SE$Method<-factor(contents_SE$Method, 
                                               levels=c("BioHC_RF", "SpRF_HC", "HMSC_HC", "MNet_HC", "GDM_HC", "bbGDM_HC", "GF_HC", 
                                                        "RCP_Hard","RCP_coefs"))
#pretty species labels,
pretty_sp<-c("A.rostrata", "B.antarcticus"  , "B.eatonii" ,  "B.irrasa",                 
             "B.murrayi",   "C.gunnari", "D.eleginoides" , "E.viator"  ,"G.acuta" ,
             "L.mizops", "L.squamifrons",  "L.antarcticus", "Macrourus.spp" ,  "M.maculata" ,            
             "Muraenolepis.spp" , "N.rossii"  , "P.gracilis"  ,"Paraliparis.spp." ,                
             "Z.spinifer"  , "C.rhinoceratus")

contents_SE$species<-rep(pretty_sp, 36)

red_sp<-c( "B.antarcticus"  , "B.eatonii" ,                  
           "C.gunnari", "D.eleginoides" ,"G.acuta" ,
           "L.mizops", "L.squamifrons",   "Macrourus.spp" ,           
           "Muraenolepis.spp" , "Paraliparis.spp." ,                
           "C.rhinoceratus")


p<-ggplot(data = contents_SE[contents_SE$species %in% red_sp,], 
#p<-ggplot(data = contents_SE,
          aes(x = species, y = mean, ymin = lower, ymax = upper, colour = Method)) +
  scale_y_continuous(limits = c(-0.1,1.1)) +
  geom_point(position = position_dodge(0.6), size=0.6) +
  geom_errorbar(position = position_dodge(0.6), width = 0.5) +
  coord_flip() +
  scale_colour_manual(name="Method", 
                      values = c("purple" ,"pink", "yellow"  , "orange", "green", "chartreuse4", "darkgreen",
                                  "cornflowerblue", "darkblue")) + #need to add enough colour for sampling levels here!
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.y = element_text(face="italic"),
        legend.key = element_blank()) +
  facet_wrap( ~Group, ncol=4, scales="free_x") 

tiff(file="Results/Plots/Grp_sp_redSE.tiff", height=8, width=18, units="cm", res=1000)
p 
dev.off()



