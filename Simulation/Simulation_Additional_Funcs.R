###################################################################################################################################
## Additional function for running, predicting, diagnosing and exploring community models and their outputs for bioregionalisation
## N.Hill March 2018 with input from S. Woolley, S. Foster
###################################################################################################################################

##################################
## Choosing Clusters ----
##################################

#Choose number of clusters based on average sihouette width (higher values better), Calinski and Harabasz index (CH; higher values better)
# and average stability of groups determined using bootstrap methods (higher values better)

grp_stats<- function(min_grp, max_grp,        # set values for number of groups to be tested
                    tree_clust,               # object output from hclust()
                    dissim,                   # either a dissimilarity/distance matrix of class 'dist' or call to dist() to produce one
                    method)                   # method used for hierarchical clustering (see hclust())
{
  require(fpc)
res<-data.frame(sil=rep(0,(max_grp-min_grp +1)),ch=rep(0,(max_grp-min_grp +1)),boot= rep (0,(max_grp-min_grp +1)), row.names = paste0("grp", min_grp:max_grp))
grps<-min_grp:max_grp

for (i in 1:length(grps)){
  test<-cutree(tree_clust,grps[i])
  test_stats<-cluster.stats(dissim, test)
  test_boot<-clusterboot(data=dissim, B=50, bootmethod="boot",clustermethod=disthclustCBI, method=method, k=grps[i])
  res[i,]<-round(c(test_stats$avg.silwidth, test_stats$ch, mean(test_boot$bootmean)),3)
}
return(res)
}



#################################
## bbGDM ----
#################################

#function to generate predicted pairwise dissimilarites for entire region based on environmental variables

pred_gdm_dissim<-function(bbgdm_mod,           # object output from bbgdm function
                          naive=FALSE,         # is object naive (nonbootstrap) gdm
                          env_data)            # marix or dataframe containing environmental variables to predict with (only include those used in bbgdm_mod)   
{
  # get model betas  
  if (naive== FALSE){
  betas<-bbgdm_mod$median.coefs.se
  }
 else{
   betas<-bbgdm_mod$starting_gdm$coef
 }
    
  
  #create table of prediction pairwise differences in environmental variables
  env_diffs <- bbgdm::diff_table_cpp(env_data)
  #remove intercept
  #env_diffs<-env_diffs[,-1]
  
  #transform using attributes of spline functions in explanaotry bbgdm model
  X <- lapply(1:ncol(env_diffs),function(x)bbgdm::spline_trans_for_pred(env_diffs[,x], attrib = bbgdm_mod$dissim_dat_params[[x]])) 
  Xall <- do.call(cbind,X)
  dim(Xall) # should be ndissimilarities * ncoefs from bbgdm model.
  
  # generate the linear predictor
  pred_dissim_linear_predictor <- cbind(1,Xall)%*%betas
  
  #convert to dissimilarities (assumes logit link was used in bbgdm_mod)
  predicted_dissimilarities <- plogis(pred_dissim_linear_predictor)
  
  #reformat into dissimilarity matrix
  dissim<-predicted_dissimilarities
  class(dissim)='dist'
  attr(dissim,"Size")<-nrow(env_data)
  dissim<-as.dist(dissim)
  
  return(dissim)
}

############################
## Mistnet ----
############################

#### Calculate variable importance using Olden method, adapted code from "olden" function in NeuralNetTools 
## Works for mistnet with one hidden layer
olden_mistnet<-function (mod_in){
  
  #get names of predictor variables and add latent factors
  x_names<-c(colnames(mod_in$x), paste0("latent", 1:mod_in$sampler$ncol))
  
  #get names of response variables
  y_names<-colnames(mod_in$y)
  
  out<-list()
  out[[1]]<-matrix(nrow=length(x_names), ncol=length(y_names))
  
  # get weights of input to hidden layer
  inp_hid<-mod_in$layers[[1]]$weights
  
  #get weights of hidden layer to output layer 
  hid_out<-mod_in$layers[[2]]$weights
  
  #Calculate sum(input-hidden * hidden-output) to get variable importance.
  #Loop through each output variable
  for (i in 1:length(y_names)){
    out[[1]][,i] <- inp_hid %*% matrix(hid_out[,i])
  }
  
  dimnames(out[[1]])<- list(x_names, y_names )
  
  out[[2]]<-apply(out[[1]], 2, function(x) abs(x)/sum(abs(x))*100 )
  
  out[[3]]<-data.frame (mean=apply(out[[2]], 1, mean), sd=apply(out[[2]], 1, sd))
  
  names(out)<-c("Variable influence", "Relative importance", "Overall relative importance")
  return(out)
}



############################
## RCP----
############################

#### Forward Selection-linear terms----
fwd_step_linear<-function(start_vars,         # names of variables to include in base model (usually from previous step)
                          start_BIC,          # value of BIC from best model in previous step
                          add_vars,           # names of variables to add (one at a time) in this step
                          species,            # names of species to model
                          dist= "Bernoulli",  # distribution for data model (see regimix())
                          nstarts=50,         # number of starts for regimix.multifit
                          form.spp=NULL,      # formula for species artifacts (see regimix())
                          data,               # dataframe that contains all the data to fit regimix model
                          min.nRCP=1,         # minimum number of RCPs to consider
                          max.nRCP,           # maximum number of RCPs to consider
                          mc.cores=1)         # number of cores if parallel processing
  {
  
  #set up results
  BICs<-as.data.frame(setNames(replicate((max.nRCP-min.nRCP+3),numeric(0), simplify = F), c("Var",paste0("RCP", rep(min.nRCP:max.nRCP)),"Start_BIC")))
  
  #loop through adding variables  
  for(j in 1:length(add_vars)){
    add<-add_vars[j]
    temp_form<-as.formula(paste("cbind(",paste(species, collapse=", "),")~ 1+",
                                paste0(paste(start_vars, collapse="+"), "+", paste(add, collapse = "+"))))
    #temp_form<-as.formula(paste("cbind(",paste(species, collapse=", "),")~",
    #                            paste0(paste(start_vars, collapse="+"), paste(add, collapse = "+"))))
    
    #run nRCPs  
    nRCPs_start<- list()
    for( ii in min.nRCP:max.nRCP)
      nRCPs_start[[ii-diff(c(1,min.nRCP))]] <- regimix.multifit(form.RCP=temp_form, form.spp= form.spp, data=data, nRCP=ii, 
                                                                inits="random2", nstart=nstarts, dist=dist, mc.cores=mc.cores)
    
    #get BICs
    RCP1_BICs <- sapply( nRCPs_start, function(x) sapply( x, function(y) y$BIC))
    #Are any RCPs consisting of a small number of sites?  (A posteriori) If so remove.
    RCP1_minPosteriorSites <- cbind( nrow(data), sapply( nRCPs_start[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))
    RCP1_ObviouslyBad <- RCP1_minPosteriorSites < 2
    RCP1_BICs[RCP1_ObviouslyBad] <- NA
    
    
    RCP1_minBICs <- apply( RCP1_BICs, 2, min, na.rm=TRUE)
    BICs[j,1:(ncol(BICs)-1)]<-c(paste0(add), round(RCP1_minBICs,1))
    
  }
  BICs$Start_BIC<-start_BIC
  return(BICs)
}


#### CALCULATE MEAN, SD AND CI OF EXPECTED PROBABILITY OF SPECIES' OCURRENCES IN EACH RCP----
# Mean, SDs and CIs calculate using bootstraps samples
# Option to calculate expected values for each botostrap or summaries for each sampling factor, or across all levels of the data
# currently only accomodates one sampling factor at a time



## Calculates average species prevalence and uncertainty estimates for each RCP using Bayesian bootstraps
## Accomodates single sampling factor
calc_prev<-function(boot_obj,                                              # Regiboot object 
                    mod_obj,                                               # Regimod  object
                    samp_fact=NULL,                                        # Vector of sampling factor names if present
                    calc_level=c("NULL", "bootstrap", "sample_fact", "overall"),   # Level at which to calculate expected prevalence
                    CI=c(0.025,0.975))                                     # Levels at which to calculate confidence intervals
{
  #require(tidyr)
  
  #set up coefficient extraction
  taus<-grepl("tau",dimnames(boot_obj)[[2]])
  alphas<-grepl("alpha",dimnames(boot_obj)[[2]])
  
  if (is.null(samp_fact)){
    
    res_all<-list()
    for(i in 1:dim(boot_obj)[1]){
      
      #extract and reformat coeficients 
      #alpha- OK as is
      temp_alphas<-boot_obj[i,alphas]
      
      #tau
      temp_tau <- boot_obj[i,taus]
      temp_tau <- matrix( temp_tau, nrow=length(mod_obj$names$RCPs)-1)
      tau_all <- rbind( temp_tau, -colSums( temp_tau))
      colnames( tau_all) <- mod_obj$names$spp 
      rownames( tau_all) <- mod_obj$names$RCPs
      
      #calculate values
      lps <- sweep( tau_all, 2, temp_alphas, "+") 
      res_all[[i]]<-as.matrix(round(exp( lps)/ (1+ exp(lps)),3))
    }
    
    overall_temp<-array(unlist(res_all), dim=c( length(mod_obj$names$RCPs),length(mod_obj$names$spp),nrow(boot_obj)))
    overall_res<-list( mean=round(apply(overall_temp, c(1,2), mean),3),
                       sd= round(apply(overall_temp, c(1,2), sd),3),
                       lower= round(apply(overall_temp, c(1,2), function(x) quantile(x, probs=CI[1])),3),
                       upper= round(apply(overall_temp, c(1,2), function(x) quantile(x, probs=CI[2])),3))
    
    dimnames(overall_res[[1]])<-dimnames(overall_res[[2]])<-dimnames(overall_res[[3]])<-dimnames(overall_res[[4]])<-list(mod_obj$names$RCPs, mod_obj$names$spp)
    return (overall_res)
  }
  
  if (! is.null(samp_fact)){
    #extract gammas and set up results frame 
    gammas<-grepl("gamma",dimnames(boot_obj)[[2]])
    res<-rep( list(list()), length(samp_fact)) 
    names(res)<-samp_fact
    
    for(i in 1:dim(boot_obj)[1]){
      print(i)
      
      #extract and reformat coeficients 
      #alpha- OK as is
      temp_alphas<-boot_obj[i,alphas]
      
      #tau
      temp_tau <- boot_obj[i,taus]
      temp_tau <- matrix( temp_tau, nrow=length(mod_obj$names$RCPs)-1)
      tau_all <- rbind( temp_tau, -colSums( temp_tau))
      colnames( tau_all) <- mod_obj$names$spp 
      rownames( tau_all) <- mod_obj$names$RCPs
      
      #gamma
      temp_gamma<-boot_obj[i, gammas]
      temp_gamma<-matrix(temp_gamma, nrow=length(mod_obj$names$spp))
      colnames(temp_gamma)<-mod_obj$names$Wvars
      rownames(temp_gamma)<-mod_obj$names$spp
      
      ## Level 1 of sampling factor (no gamma adjustment needed)
      lps <- sweep( tau_all, 2, temp_alphas, "+") 
      res[[1]][[i]]<-as.matrix(round(exp( lps)/ (1+ exp(lps)),3))
      
      ## other levels of sampling factor
      for(j in 1:length(mod_obj$names$Wvars)){
        lps_temp<-sweep( lps, 2, temp_gamma[,j], "+") 
        res[[j+1]][[i]]<- as.matrix(round(exp( lps_temp)/ (1+ exp(lps_temp)),3))
      }
    }
    
    if(calc_level=="bootstrap"){
      return(res)
    }
    
    
    #Compile list of summaries at the sampling factor level
    if(calc_level=="samp_fact"){
      samp_res<-rep( list(list()), length(samp_fact))
      names(samp_res)<-samp_fact
      
      for(k in 1: length(samp_fact)){
        samp_res[[k]]<-list(mean=round(apply(simplify2array(res[[k]]), c(1,2), mean),3),
                            sd=round(apply(simplify2array(res[[k]]), c(1,2), sd),3),
                            lower=round(apply(simplify2array(res[[k]]), c(1,2), function(x) quantile(x, probs=CI[1])),3),
                            upper=round(apply(simplify2array(res[[k]]), c(1,2), function(x) quantile(x, probs=CI[2])),3))
      }
      return(samp_res)
    }
    
    #compile summaries across each bootstrap for all sampling factors
    if(calc_level=="overall"){
      #extract ith bootstrap values for each sampling factor
      #average species values across sampling factor for each bootstrap
      #perform calculations across averaged bootstrap values
      overall_temp<-list()
      for(i in 1:dim(boot_obj)[1]){
        get_vals<-lapply(res, function(x) x[[i]])
        overall_temp[[i]]<-apply(simplify2array(get_vals), c(1,2), mean)
      }
      overall=list(mean=round(apply(simplify2array(overall_temp), c(1,2), mean),3),
                   sd= round(apply(simplify2array(overall_temp), c(1,2), sd),3),
                   lower= round(apply(simplify2array(overall_temp), c(1,2), function(x) quantile(x, probs=CI[1])),3),
                   upper= round(apply(simplify2array(overall_temp), c(1,2), function(x) quantile(x, probs=CI[2])),3))
      return(overall)
    }
  }
}






############################################################################################
## Tabulate species average and SD of occurrence or predicted probabilities 
## OR average and SD of environmental variables for each group ---
############################################################################################


#match species occurrences (or environment values) with predicted clusters from heirarchical clustering
get_match_vals<-function(site_data,          #dataframe or matrix containing data of species or environment (rows) at sites (cols)
                         pred_cluster_vals,  # vector of values of predicted cluster
                         site_index)         # index values to extract from pred_cluster_vals
{
  pred_match<-data.frame(site_data, clust=as.factor(pred_cluster_vals[site_index]))
  mean_df<-as.data.frame(t(round(aggregate(.~clust, data=pred_match, mean)[,-1],3)))
  sd_df<-as.data.frame(t(round(aggregate(.~clust, data=pred_match, sd)[,-1],3)))
  se_df<-sd_df/sqrt(length(site_index))
  res<-list(mean_df, sd_df, se_df)
  names(res)<-c("mean", "sd", "se")
  return(res)
}

env_match_vals<-function(site_data,          #dataframe or matrix containing data of species or environment (rows) at sites (cols)
                         pred_cluster_vals,  # vector of values of predicted cluster
                         site_index)         # index values to extract from pred_cluster_vals
{
  pred_match<-data.frame(site_data, clust=as.factor(pred_cluster_vals))
  mean_df<-as.data.frame(t(round(aggregate(.~clust, data=pred_match, mean)[,-1],3)))
  sd_df<-as.data.frame(t(round(aggregate(.~clust, data=pred_match, sd)[,-1],3)))
  se_df<-sd_df/sqrt(length(site_index))
  res<-list(mean_df, sd_df, se_df)
  names(res)<-c("mean", "sd", "se")
  return(res)
}


# format table of species means and SDs or SEs ready for compiling a dotplot.
dotplot_sp_tab<-function(mean_df,                     #data frame containing mean values. row=species, cols= groups
                         error_df,                    #data frame containing SD or SE values. row=species, cols= groups
                         nGrp,                        #number of groups
                         species,                     #vector of species names
                         method)                      #name of analysis method
{
  require(tidyr)
  
  temp<-gather(as.data.frame(mean_df), key= Group, value=mean)
  temp$Group<-rep(1:nGrp, each=length(species))
  
  temp_error<-gather(as.data.frame(error_df), key= Group, value=sd)[,2]
  temp$lower<-temp$mean- temp_error
  temp$upper<-temp$mean+ temp_error
  temp$species<-rep(species,nGrp)
  temp$Method<-as.character(method)
  
  return(temp)
}

# format table of environment means and sds ready for compiling a dotplot.
dotplot_env_tab<-function(mean_df,    #data frame containing mean values. row=species, cols= groups
                          error_df,   #data frame containing SD or SE values. row=species, cols= groups
                          nGrp,       #number of groups
                          env,        #vector of environmental variable names
                          method)     #name of analysis method
{
  require(tidyr)
  temp_error<-gather(as.data.frame(error_df), key= Group, value=error)[,2]
  
  temp<-gather(as.data.frame(mean_df), key= Group, value=mean)
  temp$Group<-rep(1:nGrp, each=length(env))
  temp$lower<-temp$mean- temp_error
  temp$upper<-temp$mean+ temp_error
  temp$env<-rep(env,nGrp)
  temp$Method<-as.character(method)
  
  return(temp)
}


###################################################################
## Partial response plots for SAMs and RCPs----
###################################################################

#### partial plots----
#1) generate sequence data to use for plotting (use function: partial_values)
#2) transform plotting data to same scale as used to build models (lapply on partial_values using scale function)
#3) run predictions and plot results (use function: plot_partial)

partial_values<-function(env_vars,       # env_vars= environmental variable names used to build models
                         raw_data,       # raw_data= original raw data dataframe (i.e. unscaled, untransformed)
                         rep_out=100)    # number of values along sequence for plotting
{
  
  #list to contain results
  dat_frames<-list()
  
  #for each variable in env_vars generate sequence from min to max 
  #and paste mean of all other variables
  for(i in 1:length(env_vars)){
    temp=as.data.frame(matrix(nrow=rep_out, ncol = length(env_vars)))
    dat_frames[[i]]<-as.data.frame(cbind(part= seq(from=min(raw_data[,names(raw_data) %in% env_vars[i]]), 
                                                   to=max(raw_data[,names(raw_data) %in% env_vars[i]]),
                                                   length=rep_out),
                                         do.call(rbind,
                                                 replicate(rep_out,colMeans(raw_data[,names(raw_data) %in% env_vars[-i]],na.rm=TRUE)
                                                           , simplify=FALSE))))
    #name variables
    names(dat_frames[[i]])[1]<-env_vars[i]
    dat_frames[[i]]<-dat_frames[[i]][,env_vars]
  }
  #name list element  -added
  names(dat_frames)<-env_vars
  
  dat_frames
}

### predict onto partial,scaled dataframe and plot partial plots. Point predictions only atm!
plot_partial<-function(env_vars,             # names of environmental variables used in model
                       partial_space_out,    # scaled values for prediction (i.e. output of step 2)
                       partial_vals_out,     # unscaled values for plotting (i.e. output of partial_values)
                       model=c("sam", "rcp"),# predict using sam or rcp
                       object,               # SAM or RCP model to use for predictions
                       mfrow=c(2,2))         # layout for plotting
                       #file_name=NULL,       # filename for saving pdf of plot
                       #width=8, height=8)    # size of pdf 
{
  #pdf(file= paste0(file_name, "PPlots.pdf"), height=height, width=width)
  #palette(cols)
  par(mfrow=mfrow, mar=c(5,4,2,2))
  for(i in 1:length(partial_vals_out)){
    temp<-as.data.frame(partial_space_out[[i]])
    
    if(model=="sam"){
    pred<-predict(object=object, new.obs=temp)
    matplot(x=as.vector(partial_vals_out[[i]][,names(partial_vals_out[[i]]) %in% env_vars[i]]), 
            y= pred$fit, type='l', 
            xlab=env_vars[i], ylab="Probability of Occurrence", ylim=c(0,1), lwd=3)
    
    legend(x="topright", legend=dimnames(pred$fit)[[2]], col=1:length(dimnames(pred$fit)[[2]]), lty=1, cex=0.9, bty="n")
    }
  
   if(model=="rcp"){
     pred<-predict(object=object, newdata=temp)
     matplot(x=as.vector(partial_vals_out[[i]][,names(partial_vals_out[[i]]) %in% env_vars[i]]), 
             y= pred, type='l', 
             xlab=env_vars[i], ylab="Probability of Occurrence", ylim=c(0,1), lwd=3)
     legend(x="topright", legend=dimnames(pred)[[2]], col=1:length(dimnames(pred)[[2]]), lty=1, cex=0.9, bty="n")
   }
}
  

}

### PLOT HARD CLASS VERSION OF PREDICTIONS----

hclass_maps<-function(   SpPreds,               # Output of predict.regimix 
                         coast =NULL,           # Object of class spatialPolygons representing coastline
                         pred_space,            # Untransformed pred space (or dataframe where columns 1:2 are coordinates for plotting)
                         rast,                  # Raster of same size and attributes as prediction space (used in call to 'rasterize')
                         cols,                  # Colours to use for plottigng classes
                         title_text=NULL)       # Title text
{
  hc_rast<-rasterize(pred_space[,1:2], rast, field=apply(SpPreds[[2]], 1, which.max))
  names(hc_rast)<-"RCP"
  
  #create table of factors for raster plotting
  hc_rast<-ratify(hc_rast)
  rat <- levels(hc_rast)[[1]]
  rat$RCP <- paste0( "RCP ", 1:hc_rast@data@max)
  rat$col <- cols
  levels(hc_rast) <- rat
  
  #plot results
  plot(hc_rast, colNA="white",col=rat$col, legend=FALSE)
  if(!is.null(coast)){
    plot(coast,add=TRUE, col= "grey")
  }
  legend("topright", legend=rat$RCP, fill=cols)
  title(title_text)
}


