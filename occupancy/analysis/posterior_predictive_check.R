## --------------------------------------------------
### Posterior Predictive Check
# Tukey Freeman (species-binned) discrepancy 

# Read in a model file
stan_out <- readRDS("./occupancy/model_outputs/bombus_30km_300minpop3minpersp3_5.rds")

# Read in a model file
#stan_out <- readRDS("./occupancy/simulation/stan_sim_out.rds")


## --------------------------------------------------
### PPC

# as data frame
fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

list_of_draws <- as.data.frame(stan_out)

# Bayesian P-values (P_average should be > 0.10 otherwise the model fit would be highlu questionable)
P_average_citsci = vector(length = n_species)

# Note the column number is just the column number of the first species P in the df
for(i in 1:n_species){
  P_average_citsci[i] = fit_summary$summary[196+i,1]
}

# Bayesian P-values (P_average should be > 0.10 otherwise the model fit would be highlu questionable)
P_average_museum = vector(length = n_species)

# Note the column number is just the column number of the first species P in the df
for(i in 1:n_species){
  P_average_museum[i] = fit_summary$summary[262+i,1]
}

# Visual evaluation of fit
# Could plot by species or all species on same plot but coloured differently
plot(fit_summary$summary[262,1], fit_summary$summary[262,1], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 150),
     xlim = c(0, 150))
abline(0, 1, lwd = 2, col = "black")

# I binned the discrepancies by species, so the column numbers are just the first species
# Should be close to 1. 
# If the mean of actual data is greater (value > 1)
# then the model underpredicts the real variation in counts.
# If the mean of actual data is less (value < 1)
# then the model overpredicts the variation in counts.
mean(list_of_draws[,319]) / mean(list_of_draws[,269])

# Should be close to 50% or 0.5
# similarly the actual data should be further away from  
# the expected value of the count about half of the time,
# (versus a count generated using the abundance rate and detection rate)
mean(list_of_draws[,44] > list_of_draws[,19])

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out_sim, pars = c(
  "mu_eta_0",
  "eta_site_area",
  "mu_p_citsci_0",
  "mu_p_museum_0",

))

# pairs plot
pairs(stan_out_sim, pars = c(
  "mu_eta_0",
  "eta_site_area",
  "mu_p_citsci_0",
  "mu_p_museum_0",
  "phi",
  "gamma_0",
  "gamma_1"
))

# Posterior predictive check
list_of_draws <- as.data.frame(stan_out_sim)