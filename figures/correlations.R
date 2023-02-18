# File to make results figures
# jcu, started dec 5, 2022.

library(tidyverse)

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/bombus_25km_750minpop1minpersp4_3_sq.RDS")
species_names <- readRDS("./figures/bombus_names_20min.RDS")

list_of_draws <- as.data.frame(stan_out)

n_species <- length(species_names)

n_draws <- 2400
posterior_length <- nrow(list_of_draws)

corr <- vector(length = n_draws)
samples <- matrix(nrow = n_species, ncol = 2)

for(draw in 1:n_draws){

  for(species in 1:n_species){
    
    # for each species, draw a random sample from the posterior for:
    # psi0 urban areas (starts at column 21)
    samples[species,1] <- list_of_draws[sample.int(posterior_length, 1),
                                        21+species] 
    # and psi species natural habitat (starts at column 103)
    samples[species,2] <- list_of_draws[sample.int(posterior_length, 1),
                                        103+species]
    
    # calculate the correlation between the two terms across all species
    # and add to a vector of correlations, then
    # rerun the 
    corr[draw] <- cor(samples[,1], samples[,2])
    
  }

}

# mean correlations across n_species
mean <- mean(corr)
# sd of correlations across n_species
sd <- sd(corr)

# 95% confidence interval
upper95 <- mean + 1.95*sd
lower95 <- mean - 1.95*sd

# 90%
upper90 <- mean + qnorm(0.95)*sd
lower90 <- mean - qnorm(0.95)*sd

# 50%
upper50 <- mean + qnorm(0.75)*sd
lower50 <- mean - qnorm(0.75)*sd

hist(corr, breaks = 20, 
     main="Correlation between urban occupancy rate and \n 
     species-specific effect of natural habitat on urban occupancy")
abline(v = mean, col="black", lwd=3, lty=1)
abline(v = cbind(lower95, upper95), col="blue", lwd=3, lty=2)
abline(v = cbind(lower90, upper90), col="red", lwd=3, lty=2)
#abline(v = cbind(lower50, upper50), col="red", lwd=3, lty=2)
