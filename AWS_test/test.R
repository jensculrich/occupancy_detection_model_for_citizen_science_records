# test model for AWS

N = 200
mu = 0
sigma = 1

set.seed(13)
y = rnorm(n=N, mean=mu, sd=sigma)

stan_data <- c("N", 
               "y") 

# Parameters monitored
params <- c(
  "mu",
  "sigma"
)

# MCMC settings
n_iterations <- 800
n_thin <- 1
n_burnin <- 400
n_chains <- 3
n_cores <- 4

stan_model <- "./AWS_test/test.stan"

## Call Stan from R
stan_out_sim <- stan(stan_model,
                     data = stan_data, 
                     #init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     open_progress = FALSE,
                     #save_warmup = FALSE,
                     cores = n_cores)

print(stan_out_sim, digits = 3)
