# File to make results figures
# jcu, started dec 5, 2022.

library(tidyverse)

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/large_files/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits.rds")
stan_out2 <- readRDS("./occupancy/model_outputs/syrphidae/non_urban/syrphidae_40km_1000minpop_2minUniqueDetections_3ints_3visits_.rds")

species_names <- readRDS("./figures/species_names/syrphidae_names_10km_urban.RDS")
species_names2 <- readRDS("./figures/species_names/syrphidae_names_40km_nonurban.RDS")

list_of_draws <- as.data.frame(stan_out)
list_of_draws2 <- as.data.frame(stan_out2)

n_species <- length(species_names)
n_species2 <- length(species_names2)

# remove species that didn't occur in the nonurban data set
# these occur in species_names but not in species_names2
# 157 is the number of columns before psi_herb_shrub_forest species starts
check <- which(is.na(charmatch(species_names, species_names2)))
check_plus <- which(is.na(charmatch(species_names, species_names2))) + 157

species_names <- species_names[-(c(check))]
#which(is.na(charmatch(species_names2, species_names)))
list_of_draws <- list_of_draws[,-(c(check_plus))] 

n_draws <- 1000
posterior_length <- nrow(list_of_draws)
posterior_length2 <- nrow(list_of_draws2)

correlation <- vector(length = n_draws)
samples <- matrix(nrow = n_species, ncol = 2)

for(draw in 1:n_draws){

  for(species in 1:n_species){
    
    # for each species, draw a random sample from the posterior for:
    # psi0 non-urban habitat (starts at column 14) from list_of_draws2
    samples[species,1] <- list_of_draws2[sample.int(posterior_length2, 1),
                                        13+species] 
    # and psi species natural habitat (starts at column 158) 
    # from list_of_draws
    samples[species,2] <- list_of_draws[sample.int(posterior_length, 1),
                                        157+species]
    
    # calculate the correlation between the two terms across all species
    # and add to a vector of correlations, then
    # rerun the 
    correlation[draw] <- cor(samples[,1], samples[,2])
    
  }

}

# mean correlations across n_species
mean <- mean(correlation)
# sd of correlations across n_species
sd <- sd(correlation)

# 95% confidence interval
upper95 <- mean + 1.95*sd
lower95 <- mean - 1.95*sd

# 90%
upper90 <- mean + qnorm(0.95)*sd
lower90 <- mean - qnorm(0.95)*sd

# 50%
upper50 <- mean + qnorm(0.75)*sd
lower50 <- mean - qnorm(0.75)*sd

hist(correlation, breaks = 20, 
     main="Correlation between species-specific range-wide occupancy rate and \n 
     species-specific effect of natural habitat on urban occupancy (hoverflies)",
     xlim = c(-1, 1),
     xlab = "Pearson's correlation coefficient",
     ylab = paste0("Frequency in ", n_draws, " draws from posterior distributions"))
abline(v = mean, col="black", lwd=3, lty=1)
abline(v = cbind(lower95, upper95), col="blue", lwd=3, lty=2)
#abline(v = cbind(lower90, upper90), col="green", lwd=3, lty=2)
abline(v = cbind(lower50, upper50), col="red", lwd=3, lty=2)


