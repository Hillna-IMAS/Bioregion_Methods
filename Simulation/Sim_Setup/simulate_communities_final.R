
## Simulate communities for comparison of methods

#Using multivariate normal mixture model code adapted from Skipton Woolley (ref)
# Generated 3 groups of responses to hypothetical environmental layers for a region.
# Only 2 covariates affect the distribution of species- the remaining 6 have little or no effect.
# Species have differing prevalence (specific intercepts instead of group means).


#Load required libraries
library(mvtnorm)
library(MASS)
library(scatterplot3d)
library(raster)
library(rasterVis)

##Simulation function----
species_theta_generation_mixed_proportions <- function(alphas, means,variances,covariances,
                                                       nSp,dat,mix.prop,dist='negbin',
                                                       phi=NULL,plot=TRUE){
  nSigma <- dim(means)[1]
  dimSigma <- dim(means)[2]
  sigmas <- array(0L,c(dimSigma,dimSigma,nSigma))
  species_parameters <- array(0L,c(nSp,dimSigma +1,nSigma))
  for(i in 1:nSigma){
    tmp <-sigmas[,,i]
    tmp <- diag(variances[,i])
    tmp[row(tmp)!=col(tmp)] <- covariances[i]
    sigmas[,,i]<-tmp
    species_parameters[,,i] <- c(alphas, mvtnorm::rmvnorm(nSp, means[i,], sigmas[,,i]))
  }
  
  
  grp <- rmultinom(nSp, 1, mix.prop)
  row.names(grp)<-1:nrow(grp)
  grp <- apply(ifelse(grp==1,as.numeric(row.names(grp)),0),2,sum)-1
  # df <- data.frame(grp=grp,do.call(cbind,plyr::alply(species_parameters,3)))
  theta <- matrix(0L,nSp,dimSigma+1)
  for(i in 1:nSigma){
    row_tmp <- which(grp == c(i-1))
    theta[row_tmp,] <- as.matrix(species_parameters[row_tmp,,i])
  }
  out <- matrix(0, dim(dat)[1], nSp)
  probs<-matrix(0, dim(dat)[1], nSp)# NH added
  X <- as.matrix(dat)
  for (s in 1:nSp){
    if (dist == "bernoulli") {
      lgtp <- X %*% theta[s, ]
      p <- exp(lgtp)/(1 + exp(lgtp))
      probs[,s]<-p
      out[, s] <- rbinom(dim(X)[1], 1, p)
    }
    if (dist == "negbin") {
      if(is.null(phi))phi<-rep(1,nSp)
      tmp <- rep(1e+05, dim(X)[1])
      while (max(tmp, na.rm = T) > 50000 | sum(tmp) < 2) {
        lpd_sp <- X %*% theta[s,]
        p <- exp(lpd_sp)
        tmp<-rnbinom(dim(X)[1],mu=p,size=1/phi[s])
      }
      probs[,s] <-p
      out[,s] <- tmp
    }
    if (dist == "poisson") {
      tmp <- rep(1e+05, dim(X)[1])
      while (max(tmp, na.rm = T) > 50000 | sum(tmp) < 2) {
        lpd_sp <- X %*% theta[s,]
        p <- exp(lpd_sp)
        tmp<-rpois(dim(X)[1],lambda =p)
      }
      probs[,s]<-p
      out[,s] <- tmp
    }
  }
  if(plot)for(i in 1:dim(theta)[2]-1)filled.contour(MASS::kde2d(theta[,1], theta[,i+1]),main=paste0("covar",i))
  return(list(sp_probs= probs, sp_data=out, thetas=theta, mu=means,sigma=sigmas,group=grp))
}


### Simulation Settings ----

## Load up the simulated environmental data    
setwd("C:\\Users\\hillna\\UTAS_work\\Projects\\Antarctic_BioModelling\\Analysis\\Community_modelling\\Comm_Analysis_Methods\\Simulation\\Sim_setup\\")
load("./sim_env_070518.RData")

# generate realisation of data for entire survey region
sim_dat<-data.frame(intercept=1,env_dat[,3:ncol(env_dat)])
sim_dat[,-1] <-scale(sim_dat[,-1])

### Generate thetas
#Generate "species" with responses to covariates and mix it up (skew the distributions)
#Create intercept and covariate values.

# intercept values (overall prevelance of species), based on values observed for species in Kerguelen PLateau RCP models
set.seed(42)
#parameters for the beta distribution for the SAM alpha parameter
betamean <- 0.3
betabeta <- 15
betaalpha <- betamean/(1-betamean) * betabeta
prevalences <- rbeta( 30, betaalpha, betabeta) #prevalences with mean of betaalpha/(betaalpha+betabeta)
curve( dbeta( x, betaalpha, betabeta), from=0.001, to=0.999) #the distribution that it is drawn from
alphas <- log( prevalences / ( 1-prevalences))  #put them on the right scale-- but note that this is conditional on all covars being zero

#alphas <- rnorm(30, mean=-1.5, sd=1) 
plot(density(alphas), main="Alphas")
prev<-exp(alphas)/(1+ exp(alphas))
hist(prev, main="Prevalence")

# setup the betas - how species group should on average respond (means) to each covariate.
#columns= covariate, rows= groups
means<-as.matrix(data.frame(temp=c(0.75,0,-0.5), O2=c(0,-0.5,0), NO3=c(0,0,0), sal=c(0,0,0),
                            depth=c(0, 0,0), chla=c(0,0,0), ssh=c(0,0,0), curr=c(0,0,0)))
#variances - variation around mean response to covariate
variances <- matrix(c(rep(c(0.05,0.05),each=3), rep( rep( 0.01, 3), each=dim(means)[2]-2)),dim(means)[2],dim(means)[1],byrow=T)

#Zero covariance
covariances <- rep(0,dim(variances)[1])

#mixing proportion, the proportion species to each group.
mix_prop <- c(0.3, 0.4,0.3)
nSp <- 30  #lots of species (to see pattern in distribution)


### Generate simulated data----
set.seed(100)
sim_data <- species_theta_generation_mixed_proportions(alphas,means,variances,covariances,nSp,sim_dat,mix_prop,dist='bernoulli',plot=FALSE)

set.seed(66)
sites<-sample(1:nrow(sim_dat), 200, replace=FALSE)
sp_200<-sim_data$sp_data[sites,]
env_200<-sim_dat[sites,]
save(sim_data,sites, sp_200,env_200, sim_dat , file="Many_covars_sim.RData")


##Visualise environment, groups and distribution of species----

#plot environment layers
tiff(filename="Plots/Simulated_Truth/Sim_env.tiff", compression="lzw", 
     width=14, height=12, units="cm", res=1000)
plot(env, asp=1)
#title(main="Simulated Environment")
dev.off()

#plot betas
tiff(filename="Plots/Simulated_Truth/group_betas.tiff", compression="lzw", 
     width=10, height=10, units="cm", res=1000)
plot(x=sim_data$thetas[,2], y= sim_data$thetas[,3], col=sim_data$group+1,
     xlab="temp", ylab="O2",pch=20, main="Species' betas")
legend("bottomright", legend=paste0("Group ", 1:3),  col=1:length( unique( sim_data$group)), pch=20, cex=0.8)
dev.off()

#plot True distribution of groups
mean_alpha<-mean(alphas)
true_lps<-mean_alpha+ as.matrix(sim_dat[,2:9])%*% t(means)
true_grps<-exp(true_lps)/(1 +exp(true_lps))
grp_raster<-rasterize(env_dat[,1:2], env, field= true_grps)
names(grp_raster)<-paste0("Group", 1:3)

tiff(filename="Plots/Simulated_Truth/Prob_groups.tiff", compression="lzw", 
     width=10, height=5, units="cm", res=1000)
levelplot(grp_raster, par.settings=YlOrRdTheme, main="True Group Distribution")
dev.off()


#plot without intercept
true_lps2<-as.matrix(sim_dat[,2:9])%*% t(means)
true_grps2<-exp(true_lps2)/(1 +exp(true_lps2))
grp_raster2<-rasterize(env_dat[,1:2], env, field= true_grps2)
names(grp_raster2)<-paste0("Group", 1:3)
levelplot(grp_raster2, par.settings=YlOrRdTheme, main="True Group Distribution")

# hard class groups
#class_pal<-c("darkolivegreen4", "orange1","grey")
hc<-apply(true_grps,1, which.max)
hc_raster<-rasterize(env_dat[,1:2], env, field= hc)
hc_raster<-as.factor(hc_raster)
levels(hc_raster)<-data.frame(ID=1:3, Group=paste0("Group", 1:3))

tiff(filename="Plots/Simulated_Truth/HC_groups.tiff", compression="lzw", 
     width=10, height=10, units="cm", res=1000)
levelplot(hc_raster, par.settings=YlOrRdTheme, main="True Groups: Hard Class", margin=FALSE)
dev.off()

#plot species' occurrence probabilities across region
sp_prob_brick<- rasterize(x= SpatialPointsDataFrame(coords=env_dat[,1:2], data=as.data.frame(sim_data$sp_probs)), 
                          y=raster(ncols=51,nrows=51, extent(144.95, 150.05, -40.05, -34.95)))
sp_prob_brick<-dropLayer(sp_prob_brick,1)
names(sp_prob_brick)<-paste0("Sp", 1:30)

# set probability colour pallette
prob_pal<-colorRampPalette(c("#FFFFCC", "#FFEDA0", "#FED967", "#FEB24C",  "#FD8D3C",
                             "#FC4E2A", "#E31A1C", "#BD0026", "#800026"), space = "rgb")

tiff(filename=paste0("Plots/Simulated_Truth/Grp_1species.tiff"), compression="lzw", 
     width=10, height=12, units="cm", res=1000)
levelplot(sp_prob_brick, layers=which(sim_data$group==0), layout=c(4,4),  main= "Group 1",
          col.regions=prob_pal(19),
          at=seq(0,1,length=18),
          colorkey=list(at=seq(0,1,length=18),col=prob_pal(19)))
dev.off()

tiff(filename=paste0("Plots/Simulated_Truth/Grp_2species.tiff"), compression="lzw", 
     width=10, height=12, units="cm", res=1000)
levelplot(sp_prob_brick, layers=which(sim_data$group==1), layout=c(4,3),main= "Group 2",
          col.regions=prob_pal(19),
          at=seq(0,1,length=18), colorkey=FALSE)
dev.off()


tiff(filename=paste0("Plots/Simulated_Truth/Grp_3species.tiff"), compression="lzw", 
     width=10, height=12, units="cm", res=1000)
levelplot(sp_prob_brick, layers=which(sim_data$group==2), layout=c(4,3), main= "Group 3",
          col.regions=prob_pal(19),
          at=seq(0,1,length=18), colorkey=FALSE)
dev.off()


