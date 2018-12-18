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
setwd("C:\\Users\\hillna\\UTAS_work\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\KP_Fish\\")
source("Code/Additional_Funcs.R")

load("Data/Prediction_space.Rda") 
load("Results/pred_clusters.RData")

pred_sp<-na.omit(pred_sp)
#####################################################################################
#1) plot predicted distribution of optimal number of groups determined by models----
#a) hard classes for all methods except SAMs
#b) probablistic for RF, RCP, SAMs
#####################################################################################

### Look for label switching of classes in hard class outputs
# (i.e. class1 in Sp_RF is class 2 in HMSC). Make all relative to Sp_RF

mods<-c("Env_Only", "SpRF_HC","HMSC_HC" , "MNet_HC" , "GDM_HC", "bbGDM_HC"  ,  "GF_HC")
check_clust2<-list()
for (i in 1: length(mods)){
  check_clust2[[i]]<- table(hard_cluster2$Sp_RF, hard_cluster2[,mods[i]])
}
names(check_clust2)<-mods

# attempt to make labels match up to true classes
hard_cluster_opt$SpRF_HC<-mapvalues(hard_cluster_opt$SpRF_HC ,from=c(1,2), to=c(2,1))
hard_cluster_opt$GF_HC<-mapvalues(hard_cluster_opt$GF_HC, from=c(1,2), to=c(2,1)) 
hard_cluster_opt$HMSC_HC<-mapvalues(hard_cluster_opt$HMSC_HC, from=c(1,2), to=c(2,1)) 
hard_cluster_opt$MNet_HC<-mapvalues(hard_cluster_opt$MNet_HC ,from=c(1,2,3,4), to=c(2,3,4,1))


### Create raster stack with hard class groups for plotting
#set up
clust_opt<-stack()


#create factor attribute layer with as many levels as greatest number of clusters
rat<-data.frame(ID=1:5, Group=paste0("Group", 1:5))
#ignore warning in following loop

for (i in 1:length(mods)){
  
  hc_rast<-rasterize(pred_sp[,1:2], env_raster, field=hard_cluster_opt[,mods[i]])
  hc_rast<-as.factor(hc_rast)
  levels(hc_rast) <- rat
  
  clust_opt<-stack(clust_opt, hc_rast)
}
names(clust_opt)<-mods

##set up plot colours for hard classes
class_pal<-c("darkolivegreen4","grey", "orange1", "darkred", "cadetblue")

#plot simulation groups

#plot hard cluster outputs
tiff(filename="Results/Plots/Opt_HClusts.tiff", compression="lzw", 
     width=11, height=12, units="cm", res=1000)
levelplot(clust_opt, layout=c(2,4),scales=list(draw=FALSE), 
          col.regions=class_pal)
dev.off()

##plot probability of occurrence outputs
#BioHC_RF

tiff(filename="Results/Plots/BioHC_RF_prob.tiff", compression="lzw", 
     width=8, height=4, units="cm", res=1000)
levelplot(rasterize(pred_sp[,1:2], env_raster,field= bio2_rf_pred), 
         par.settings=YlOrRdTheme, names.attr=c("Group1", "Group2"),
          layout=c(2,1), scales=list(draw=FALSE))
dev.off()

#SAM
#account for label switching
dimnames(sam4_pred$fit)[[2]]<- c(paste0("Group", c(4,2,1,3)))
dimnames(sam4_pred$se.fit)[[2]]<- c(paste0("Group", c(4,2,1,3)))

prob_pal<-colorRampPalette(c("#FFFFCC", "#FFEDA0", "#FED967", "#FEB24C",  "#FD8D3C",
                             "#FC4E2A", "#E31A1C", "#BD0026", "#800026"), space = "rgb")

tiff(filename="Results/Plots/SAM_prob.tiff", compression="lzw", 
     width=8, height=8, units="cm", res=1000)
levelplot(rasterize(pred_sp[,1:2], env_raster, field=sam4_pred$fit[,paste0("Group", 1:4)]), 
          col.regions=prob_pal(19),
          at=seq(0,1,length=18),
          colorkey=list(at=seq(0,1,length=18),col=prob_pal(19)),
          names.attr=rep("",4),layout=c(2,2),
          scales=list(draw=FALSE))
dev.off()

#RCP
dimnames(rcp5_spat_preds[["ptPreds"]])[[2]]<- c(paste0("Group", c(1,2,4,5,3)))
dimnames(rcp5_spat_preds[["bootSEs"]])[[2]]<- c(paste0("Group", c(1,2,4,5,3)))

tiff(filename="Results/Plots/RCP_prob.tiff", compression="lzw", 
     width=10, height=8, units="cm", res=1000)
levelplot(rasterize(pred_sp[,1:2], env_raster,field=rcp5_spat_preds$ptPreds[,paste0("Group", 1:5)]), 
          col.regions=prob_pal(19),
          at=seq(0,1,length=18),
          colorkey=list(at=seq(0,1,length=18),col=prob_pal(19)),
          names.attr=rep(" ",5),
          scales=list(draw=FALSE), layout=c(3,2))
dev.off()

# Final plot assembled in graphics software

#####################################################################################
#2) plot 4 groups- hard classes to run further comparisons 
#####################################################################################

#convert probabalistic methods to give hard classes
hard_cluster4$BioHC_RF<-apply(bio4_rf_pred,1,which.max)
hard_cluster4$SAM<-apply(sam4_pred$fit,1,which.max)
hard_cluster4$RCP<-apply(rcp4_spat_preds[["ptPreds"]],1,which.max)

#fix label switching
hard_cluster4$SAM<-mapvalues(hard_cluster4$SAM ,from=c(1,3), to=c(2,1))
hard_cluster4$RCP<-mapvalues(hard_cluster4$RCP ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$SpRF_HC<-mapvalues(hard_cluster4$SpRF_HC ,from=c(1,2,3,4), to=c(2,4,1,3))
hard_cluster4$GDM_HC<-mapvalues(hard_cluster4$GDM_HC ,from=c(1,2,3,4), to=c(1,3,2,4))
hard_cluster4$GF_HC<-mapvalues(hard_cluster4$GF_HC ,from=c(1,2,3,4), to=c(3,2,4,1))
hard_cluster4$MNet_HC<-mapvalues(hard_cluster4$MNet_HC ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$HMSC_HC<-mapvalues(hard_cluster4$HMSC_HC ,from=c(1,2,3,4), to=c(2,3,4,1))
hard_cluster4$BioHC_RF<-mapvalues(hard_cluster4$BioHC_RF ,from=c(1,2,3,4), to=c(1,4,2,3))


mods2<-c("Env_Only", "BioHC_RF", "SpRF_HC","HMSC_HC" , "MNet_HC" , 
         "GDM_HC", "bbGDM_HC","GF_HC", "SAM", "RCP")

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

#plot results

tiff(file="Results/Plots/Cluster4_HCpredictions.tiff",compression="lzw", 
     width=14, height=16, units="cm", res=1000)
levelplot(clust4,layout=c(3,4),scales=list(draw=FALSE), 
          col.regions=class_pal[1:4])
dev.off()

#####################################################
### Plot uncertainty for SAMs and RCPS
hard_cluster4$RCP<-mapvalues(hard_cluster4$RCP ,from=c(1,2,3,4), to=c(2,3,4,1))
dimnames(rcp4_spat_preds$bootSEs)[[2]]<-paste0("Group", c(2,3,4,1))

grey_pal<-rev(gray.colors(n=19, start=0.4, end=1))

tiff(file="Results/Plots/SAM_RCP_SEs.tiff",compression="lzw", 
     width=14, height=10, units="cm", res=1000)
plot(levelplot(rasterize(pred_sp[,1:2], env_raster,field=sam4_pred[["se.fit"]][,paste0("Group", 1:4)]), 
    layout=c(4,1), col.regions=grey_pal,
    # at=seq(0,0.4,length=18),colorkey=list(at=seq(0,0.4,length=18),col=grey_pal),
    scales=list(draw=FALSE), main="SAM"), 
     split=c(1,1,1,2))

plot(levelplot(rasterize(pred_sp[,1:2], env_raster,field=rcp4_spat_preds[["bootSEs"]][,paste0("Group", 1:4)]), 
     layout=c(4,1), col.regions=grey_pal,
     #at=seq(0,0.4,length=18),colorkey=list(at=seq(0,0.4,length=18),col=grey_pal),
     scales=list(draw=FALSE), main="RCP"), 
     split=c(1,2,1,2),newpage=FALSE)
dev.off()
# uncertainty for SAMs has values of >1 for group 4

######################################################
### RCP vs SAM comparison ----
######################################################

### Plot probablistic output for RCP with 4 groups
tiff(file="Results/Plots/SAM_probs2.tiff",compression="lzw", 
     width=14, height=10, units="cm", res=1000)
levelplot(rasterize(pred_sp[,1:2], env_raster, field=sam4_pred$fit[,paste0("Group", 1:4)]), 
          col.regions=prob_pal(19),
          at=seq(0,1,length=18),
          colorkey=list(at=seq(0,1,length=18),col=prob_pal(19)),
          names.attr=paste("Archetype",1:4),layout=c(2,2),
          scales=list(draw=FALSE))
dev.off()

dimnames(rcp4_spat_preds[["ptPreds"]])[[2]]<- c(paste0("Group", c(2,3,4,1)))
dimnames(rcp4_spat_preds[["bootSEs"]])[[2]]<- c(paste0("Group", c(2,3,4,1)))

tiff(file="Results/Plots/RCP_prob2.tiff",compression="lzw", 
     width=14, height=10, units="cm", res=1000)
levelplot(rasterize(pred_sp[,1:2], env_raster,field=rcp4_spat_preds$ptPreds[,paste0("Group", 1:4)]), 
          col.regions=prob_pal(19),
          at=seq(0,1,length=18),
          colorkey=list(at=seq(0,1,length=18),col=prob_pal(19)),
          names.attr=paste0("RCP ", 1:4),
          scales=list(draw=FALSE), layout=c(2,2))
dev.off()

load("Results/models.RData")

## Tabulate SAM species composition
species=c("Antimora.rostrata" ,"Bathydraco.antarcticus", "Bathyraja.eatonii",                
          "Bathyraja.irrasa"  ,  "Bathyraja.murrayi",  "Champsocephalus.gunnari" , "Dissostichus.eleginoides" ,        
          "Etmopterus.viator"  ,  "Gobionotothen.acuta" , "Lepidonotothen.mizops"   ,  "Lepidonotothen.squamifrons",       
          "Lycodapus.antarcticus" ,"Macrourus.spp", "Mancopsetta.maculata"  , "Muraenolepis.spp",                
          "Notothenia.rossii"   , "Paradiplospinus.gracilis" ,  "Paraliparis.spp." ,"Zanclorhynchus.spinifer" ,         
          "Channichthys.rhinoceratus.velifer")

sam_sp_grps<-data.frame(species= sam4_mod$names$spp, group=apply(sam4_mod$taus, 1, which.max))
sam_sp_grps$group<-mapvalues(sam_sp_grps$group, from=c(1,2,3,4), to= c(4,2,1,3))
sam_sp_grps<-sam_sp_grps[order(sam_sp_grps$group),]
write.csv(sam_sp_grps, file="Results/sam_sp_grps.csv", row.names = FALSE)


## plot RCP species composition
KP_Fish_Env<-readRDS("Data/KP_Fish_Env.RDS")
rcp_contents<-calc_prev(boot_obj = rcp4_boot, mod_obj = rcp4_mod, 
          samp_fact = levels (as.factor(KP_Fish_Env$Survey)), 
          calc_level = "overall")

rcp_contents_SD<-dotplot_sp_tab (mean_df=t(rcp_contents$mean),
                                 error_df= t(rcp_contents$sd),
                                 nGrp=4, species= species, method="RCP")
#correct label switching
rcp_contents_SD$Group<-mapvalues(rcp_contents_SD$Group, from=c(1,2,3,4), to=c(2,3,4,1))

#pretty species labels,
pretty_sp<-c("A.rostrata", "B.antarcticus"  , "B.eatonii" ,  "B.irrasa",                 
"B.murrayi",   "C.gunnari", "D.eleginoides" , "E.viator"  ,"G.acuta" ,
"L.mizops", "L.squamifrons",  "L.antarcticus", "Macrourus.spp" ,  "M.maculata" ,            
"Muraenolepis.spp" , "N.rossii"  , "P.gracilis"  ,"Paraliparis.spp." ,                
"Z.spinifer"  , "C.rhinoceratus")

rcp_contents_SD$species<-rep(pretty_sp, 4)

red_sp<-c( "B.antarcticus"  , "B.eatonii" ,                  
              "C.gunnari", "D.eleginoides" ,"G.acuta" ,
             "L.mizops", "L.squamifrons",   "Macrourus.spp" ,           
             "Muraenolepis.spp" , "Paraliparis.spp." ,                
              "C.rhinoceratus")


#plot
#p<-ggplot(data = rcp_contents_SD,
p<-ggplot(data = rcp_contents_SD[rcp_contents_SD$species %in% red_sp,],
          aes(x = species, y = mean, ymin = lower, ymax = upper)) +
  scale_y_continuous(limits = c(-0.1,1.1)) +
  geom_point(position = position_dodge(0.6), size=0.8, color="blue") +
  geom_errorbar(position = position_dodge(0.6), width = 0.4, color="blue") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.y = element_text(face="italic"),
        legend.key = element_blank()) +
  facet_wrap( ~Group, ncol=4, scales="free_x") 

tiff(file="Results/Plots/RCP_contents_red.tiff", height=8, width=14, units="cm", res=1000)
p 
dev.off()

