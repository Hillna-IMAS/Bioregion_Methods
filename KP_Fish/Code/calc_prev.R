

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
