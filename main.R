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
dist = array(dim = c(Nsubj,Nx*Ny,Nx*Ny))
X = matrix(ncol = Nx*Ny, nrow = Nsubj)
coordx = array(dim = c(Nsubj,Nx*Ny))
coordy = array(dim = c(Nsubj,Nx*Ny))
for (i in 1:20){
  dist[i,,] = as.matrix(dist(data1[data1$subj==i,1:2]))
  X[i,] = t(data1[data1$subj==i,3])
  coordx[i,] = t(as.matrix(data1[data1$subj==i,1]))
  coordy[i,] = t(as.matrix(data1[data1$subj==i,2]))
  
  }
y = X[1, ]
dist.y = dist[1,,]

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data.input = list(Nx = Nx, Ny = Ny, Nsubj = Nsubj, N = Nx*Ny,
                  X = X, y = y, dist_y = dist.y, coordx = coordx,coordy = coordy, dist = dist, 
                  s_tx = 10, s_ty = 10, s_s = 1, apsi = 1, bpsi = 1, aphi = 1, bphi = 1)


fit = stan(file = "spatial1.stan", data = data.input, chains = 4, iter = 1e3, seed = 1,
           pars = c("theta_x","theta_y","sigma","phi","psi"))
           
saveRDS(fit, file = "spatialfit.rds")

# fit = readRDS("spatialfit.rds")
# 
# rstan::traceplot(fit,pars = c("theta_x[1]","phi","psi","sigma[1]"))
# rstan::traceplot(fit,pars = c("phi"))
# rstan::traceplot(fit,pars = c("phi"))
# rstan::traceplot(fit,pars = c("theta_x[2]"))
# rstan::traceplot(fit,pars = c("sigma[2]"))
