###############################################################################
### Generate Biological dataset to use for                                  ###
### Comparison of Community Modelling Methods for Bioregionalisation        ###
### Based on RSTS dataset used in Hill et al Diversity & Dist 2017          ###
###############################################################################


## Load in bio and env data used in D&D 2017, subset to RSTS data, remove commercially sensitive region
#DW: If you redact a box between longitude 74.47 and 74.8, and latitude -52.8 and -53.0 that should take care of them

path<-("C:\\Users\\hillna\\UTAS_work\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\")

bioenv<-readRDS(paste0(path, "Data\\comb_data.RDS"))
RSTS<-bioenv[bioenv$Dataset== "RSTS",]
RSTS<-subset(RSTS, Long <= 74.47 |  Long >= 74.8 |
                    Lat<= -53| Lat >= -52.8 )

dim(RSTS)
#[1] 524  61

table(RSTS$Survey)
#POKER_2006   POKER_2010   POKER_2013 SC Cruise 41 SC Cruise 57 SC Cruise 69 SD Cruise 59 
#    0            0            0          133          137          119          135 

## Re-examine dataset properties

species=c("Alepocephalus.antipodianus",  "Antimora.rostrata" ,"Bathydraco.antarcticus", "Bathyraja.eatonii",                
          "Bathyraja.irrasa"  ,  "Bathyraja.murrayi",  "Champsocephalus.gunnari" , "Dissostichus.eleginoides" ,        
          "Etmopterus.viator"  ,  "Gobionotothen.acuta" , "Lepidonotothen.mizops"   ,  "Lepidonotothen.squamifrons",       
          "Lycodapus.antarcticus" ,"Macrourus.spp", "Mancopsetta.maculata"  , "Muraenolepis.spp",                
          "Notothenia.rossii"   , "Paradiplospinus.gracilis" ,  "Paraliparis.spp." ,"Zanclorhynchus.spinifer" ,         
          "Channichthys.rhinoceratus.velifer")

env_vars<-c("log_slope", "Av_depth" ,
            "sea.floor.temperature", 
            "log_current"   ,   "no3_bot_mean" ,            
            "T_mean"  , "chla_yearly_mean"  ,  "ssha_variability")

# Species' prevalence
prev<-colSums(RSTS[,species]>0)

#min prevalence= 8 sites.

## Quickly checked species prevalence against original RSTS file before it was combined with POKER file
# one or 2 species with ~ 8 presences missed.
# cut new analysis file to species prevalent at > 10 sites (~2% prevalence).
# This now excludes Alepocephalus.antipodianus

RSTS<- subset(RSTS, select= -Alepocephalus.antipodianus)

# Correlation between environmental variables----

#vars used in DD paper
pairs(RSTS[,env_vars])
cor(RSTS[,env_vars])
#chl-a mean and temp mean now corelated at 0.8

# check all vars
pdf(file="all_env_corr.pdf", height=12, width=12)
pairs(RSTS[,38:57])
dev.off()

cor(RSTS[,38:57], use="pairwise.complete.obs")
# distance_max_ice_edge & chla_yearly_mean & T_var_ds & T_mean
# T_var & oxy_bot_SR  & oxy_bot_mean & salinity.at.the.sea.floor & bathymetry
# no3_bot_SR & oxy_bot_mean  & no3_bot_mean

temp_vars<-env_vars2<-c("bathymetry_slope", "caisom_floor_current_speed", "bathymetry" ,
                        "sea.floor.temperature", 
                        "no3_bot_mean" ,            
                        "T_mean"  , "chla_yearly_sd"  ,  "ssha_variability", "seaice_gt85")

# subset DD environmental space and look at maps of vars across region
library(raster)
env_rast<-brick(paste0(path, "Data/pred_masked"))
env_rast<-crop(env_rast, extent(70,78.3,-54.3, -50))

pdf(file="env_rast.pdf")
for (i in 4:length(temp_vars)){
  plot(env_rast, temp_vars[i])
}
dev.off()

#seaice looks weird. Remove.

env_vars2<-c("log_slope", "log_current", "Av_depth" ,
             "sea.floor.temperature", 
             "no3_bot_mean" ,            
             "T_mean"  , "chla_yearly_sd"  ,  "ssha_variability")


#plot GAMS for each species---
prelim_gams<-function(gam_var, species, data, family="gaussian",p_val=0.1, 
                      filename="GAM_Results", width=6, height=9, mfrow=c(4,3)){
  require(mgcv)
  
  #set up dataframe for GAM results
  res<-as.data.frame(matrix(NA,nrow=length(gam_var), ncol=length(species)),names=species,row.names=gam_var)
  names(res)<-species
  
  pdf(file=paste0(filename, ".pdf"), width=width, height=height)
  
  #for each species loop through GAM for each environmental variable, plot and record '1' in res if significant
  for(j in 1:length(species))
  {
    par(mfrow=mfrow,oma=c(0,0,1,0))
    
    for (i in 1:length(gam_var))
    { 
      print(paste(species[j], gam_var[i], sep="-"))
      vec<-paste0(species[j], "~s(",gam_var[i],")")
      
      test<-gam(as.formula(vec), family=family, data=data)
      res[i,j]<-ifelse(summary(test)$s.pv<p_val, 1,0)
      
      plot(test, main=gam_var[i], cex=0.8, xlab="", ylab="")
    }
    title(species[j],outer=TRUE)
  }
  dev.off()
  #calculate percentage of species for which each env variable is signficant
  var_imp<-rowSums(res)/length(species)*100
  
  return(var_imp)
}

gam_res<-prelim_gams(env_vars2, species[-1], data=RSTS, family="nb",p_val=0.1)

#  log_slope           log_current              Av_depth sea.floor.temperature 
#     45                    45                    90                    80 
#  no3_bot_mean            T_mean        chla_yearly_sd      ssha_variability 
#     55                    95                    45                    85 

#export data as RDS----
#Model_Building
#Do I need to rename surveys?
keep_dat<-c("Year","Survey", "sample_ID", "sample_date",  "Long", "Lat", "Av_depth", "Season",
                species[-1],
                env_vars2[-9])

saveRDS( droplevels(RSTS[, keep_dat]), "Data/KP_Fish_Env.RDS")

#Prediction environment
pred_sp<-as.data.frame(rasterToPoints(subset(env_rast, temp_vars)))
pred_sp$log_slope<- log(pred_sp$bathymetry_slope)
pred_sp$log_current<-log(pred_sp$caisom_floor_current_speed)
pred_sp$Av_depth<-pred_sp$bathymetry * -1
pred_sp<-pred_sp[,c("x", "y", env_vars2[-9])]

temp_rast<-subset(env_rast, temp_vars)
save(temp_rast, pred_sp, file="Prediction_space.Rda")

rm(keep_dat, temp_vars, env_vars, temp_rast)

