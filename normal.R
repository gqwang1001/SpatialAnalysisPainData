library(rstan)

model_file = 'normal.stan'
iterations = 500
N = 1000
mu = 100
sigma = 10
y = rnorm(N, mu,sigma) # simulate data

stan_data = list(N=N, y=y) # data passed to stan 
    # set up the model
stan_model = stan(model_file, data = stan_data, chains = 0)
stanfit = stan(fit = stan_model, data = stan_data, iter=iterations)
print(stanfit,digits=2)
