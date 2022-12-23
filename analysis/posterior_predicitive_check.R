## --------------------------------------------------
### Posterior Predictive Check
# Chi-squared discrepancy 

# Read in a model file
stan_out_sim <- readRDS("./simulation/stan_out_sim_abundance.rds")

# should now also write a posterior predictive check into the model
list_of_draws <- as.data.frame(stan_out_sim)

# Evaluation of fit
plot(list_of_draws$fit, list_of_draws$fit_new, main = "", xlab =
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
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)

# Should be close to 50% or 0.5
# similarly the actual data should be further away from  
# the expected value of the count about half of the time,
# (versus a count generated using the abundance rate and detection rate)
mean(list_of_draws$fit_new > list_of_draws$fit)