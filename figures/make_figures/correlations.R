# File to make results figures
# jcu, started dec 5, 2022.

library(tidyverse)

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/bombus_15km_1200minpop10minpersp4_3_bbna.RDS")
stan_out2 <- readRDS("./occupancy/model_outputs/bombus_30km_1000minpop10minpersp4_3nonurban_bbna.RDS")

species_names <- readRDS("./figures/bombus_names_15km.RDS")
species_names2 <- readRDS("./figures/bombus_names_30km_nonurban.RDS")

list_of_draws <- as.data.frame(stan_out)
list_of_draws2 <- as.data.frame(stan_out2)

n_species <- length(species_names)
n_species2 <- length(species_names)

# remove species that didn't occur in the urban data set
# these occur in species_names2 but not in species_names
# 17 is the number of columns before psi species starts
check <- which(is.na(charmatch(species_names2, species_names)))
check_plus <- which(is.na(charmatch(species_names2, species_names))) + 17

species_names2 <- species_names2[-(c(check))]
# which(is.na(charmatch(species_names2, species_names)))
list_of_draws2 <- list_of_draws2[,-(c(check_plus))] 

n_draws <- 1000
posterior_length <- nrow(list_of_draws)

corr <- vector(length = n_draws)
samples <- matrix(nrow = n_species, ncol = 2)

for(draw in 1:n_draws){

  for(species in 1:n_species){
    
    # for each species, draw a random sample from the posterior for:
    # psi0 non-urban habitat (starts at column 18) from list_of_draws2
    samples[species,1] <- list_of_draws2[sample.int(posterior_length, 1),
                                        17+species] 
    # and psi species natural habitat (starts at column 104) 
    # from list_of_draws
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
     main="Correlation between range-wide occupancy rate and \n 
     species-specific effect of natural habitat on urban occupancy")
abline(v = mean, col="black", lwd=3, lty=1)
abline(v = cbind(lower95, upper95), col="blue", lwd=3, lty=2)
abline(v = cbind(lower90, upper90), col="red", lwd=3, lty=2)
abline(v = cbind(lower50, upper50), col="red", lwd=3, lty=2)

# there is a weak negative association between need for natural habitat and range-wide occupancy
# species that are at lower occupancy within their range are marginally more likely to 
# need lots of natural habitat to persist in an urban landscape
# Species that can do well in urban landscapes without large patches of natural habitat
# are more likely to be at high occupancy across their range.