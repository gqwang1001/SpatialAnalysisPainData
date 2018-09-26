library(R.matlab)
library(rstan)

# read data
dat = readMat("./data/extdata.mat")
data = dat$mean.ext
data1 = data.frame(sx = rep(1:25, 16), sy = rep(1:16,each=25), y= as.vector(data[,,1]), subj = rep(1,25*16))
for(k in 2:20){
  data1 = rbind(data1,
                data.frame(sx = rep(1:25, 16), sy = rep(1:16,each=25), 
                           y= as.vector(data[,,k]), subj = rep(k,25*16)))
}

# MCMC implement using RStan
Nx = 25; Ny = 16; Nsubj = 20;
##### Prepare data
dist = array(dim = c(Nx*Ny,Nx*Ny, Nsubj))
coord = array(dim = c(Nx*Ny,2,Nsubj))
y = matrix(nrow = Nx*Ny, ncol = Nsubj)
for (i in 1:20){
  dist[,,i] = as.matrix(dist(data1[data1$subj==i,1:2]))
  y[,i] = data1[data1$subj==i,3]
  coord[,,i] = as.matrix(data1[data1$subj==i,1:2])
  }


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data.input = list(Nx = Nx, Ny = Ny, Nsubj = Nsubj, N = Nx*Ny,
                  y = y, coord = coord, dist = dist)


fit = stan(file = "spatial.stan", data = data.input, chains = 4, iter = 1e3)
saveRDS(fit, file = "spatialfit.rds")

fit = readRDS("spatialfit.rds")

rstan::traceplot(fit,pars = c("theta_x[2]","phi","psi","sigma[1]"))
rstan::traceplot(fit,pars = c("phi"))
rstan::traceplot(fit,pars = c("phi"))
