### Create environmental data for simulation 
#code modified from example in https://github.com/skiptoniam/bbgdm


library(raster)
xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
d <- as.matrix(dist(xy))
w <- exp(-1/nrow(xy) * d)
w2<-exp(-1/nrow(xy) * d^2)
ww <- chol(w)


# Gradient variables
set.seed(6)
xy$temp <- t(w) %*% rnorm(nrow(xy),0, 0.1)
xy$temp <- scales::rescale(xy$temp,to=range(16:20))

set.seed(789)
#sim temp
xy$O2 <- t(w) %*% rnorm(nrow(xy),0, 0.5)
xy$O2 <- scales::rescale(xy$O2,to=range(0.5:5))

set.seed(555)
#sim NO3
xy$NO3 <- t(w) %*% rnorm(nrow(xy),0, 0.3)
xy$NO3 <- scales::rescale(xy$NO3,to=range(0.5:15))

set.seed(498)
#sim temp
xy$sal <- t(w) %*% rnorm(nrow(xy),0, 0.1)
xy$sal <- scales::rescale(xy$sal,to=range(32:36))


# Patchier variables
set.seed(66)
xy$depth <- t(ww) %*% rnorm(nrow(xy),0, 0.1)
xy$depth <- scales::rescale(xy$depth,to=range(10:100))

set.seed(123)
xy$chla <- t(ww) %*% rnorm(nrow(xy),0, 0.1)
xy$chla <- scales::rescale(xy$chla,to=range(0:3))

set.seed(3)
xy$ssh <- t(ww) %*% rnorm(nrow(xy),0, 0.2)
xy$ssh <- scales::rescale(xy$ssh,to=range(0:3))

set.seed(498)
xy$curr <- t(ww) %*% rnorm(nrow(xy),0, 1)
xy$curr <- scales::rescale(xy$curr,to=range(0.1:20))


# Make raster stack
coordinates(xy) <- ~x+y
env <- rasterize(xy, raster(points2grid(xy)), fields=c("temp", "O2","NO3", "sal", "depth","chla", "ssh", "curr"))
env<-dropLayer(env,1)
plot(env)

# matrix form
env_dat<-rasterToPoints(env)

#save outputs
save(env, env_dat, file="C:/Users/hillna/Dropbox/simulate_communities/Many_covars/sim_env_070518.RData")




