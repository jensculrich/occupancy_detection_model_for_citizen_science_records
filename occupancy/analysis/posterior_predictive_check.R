## --------------------------------------------------
### Posterior Predictive Check
# Chi-squared discrepancy 

# Read in a model file
stan_out_sim <- readRDS("./occupancy/model_outputs/_.rds")

# as data frame
list_of_draws <- as.data.frame(stan_out_sim)

## --------------------------------------------------
# Abundance

# Evaluation of fit
plot(list_of_draws$fit_citsci, list_of_draws$fit_new_citsci, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(25000, 90000),
     xlim = c(25000, 90000))
abline(0, 1, lwd = 2, col = "black")

# Should be close to 1. 
# If the mean of actual data is greater (value > 1)
# then the model underpredicts the real variation in counts.
# If the mean of actual data is less (value < 1)
# then the model overpredicts the variation in counts.
mean(list_of_draws$fit_citsci) / mean(list_of_draws$fit_new_citsci)

# Should be close to 50% or 0.5
# similarly the actual data should be further away from  
# the expected value of the count about half of the time,
# (versus a count generated using the abundance rate and detection rate)
mean(list_of_draws$fit_new_citsci > list_of_draws$fit_citsci)

## --------------------------------------------------
# Occupancy

# Evaluation of fit
plot(list_of_draws$fit_museum, list_of_draws$fit_occupancy_new_museum, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 500),
     xlim = c(0, 500))
abline(0, 1, lwd = 2, col = "black")

# Should be close to 1. 
# If the mean of actual data is greater (value > 1)
# then the model underpredicts the real variation in counts.
# If the mean of actual data is less (value < 1)
# then the model overpredicts the variation in counts.
mean(list_of_draws$fit_museum) / mean(list_of_draws$fit_occupancy_new_museum)

# Should be close to 50% or 0.5
# similarly the actual data should be further away from  
# the expected value of the count about half of the time,
# (versus a count generated using the abundance rate and detection rate)
mean(list_of_draws$fit_occupancy_new_museum > list_of_draws$fit_museum)

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out_sim, pars = c(
  "mu_eta_0",
  "eta_site_area",
  "mu_p_citsci_0",
  "mu_p_museum_0",
  "phi",
  "gamma_0",
  "gamma_1"
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