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

######################################################
### RCP vs SAM comparison ----
######################################################

#SAM
#account for label switching
dimnames(sam4_pred$fit)[[2]]<- c(paste0("Group", c(4,2,1,3)))
dimnames(sam4_pred$se.fit)[[2]]<- c(paste0("Group", c(4,2,1,3)))

prob_pal<-colorRampPalette(c("#FFFFCC", "#FFEDA0", "#FED967", "#FEB24C",  "#FD8D3C",
                             "#FC4E2A", "#E31A1C", "#BD0026", "#800026"), space = "rgb")


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
