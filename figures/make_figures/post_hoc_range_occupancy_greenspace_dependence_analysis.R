# Post hoc analysis to test the association between rangewide occurrence rate
# and species-specific effect of natural greenspace area on urban occupancy.

# This file can be used to generate figure S19

# jcu, started dec 5, 2022.

library(rstan)
library(tidyverse)
library(rstanarm)

################################################################################
## bumble bees

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/large_files/bombus_10km_1200minpop_1minUniqueDetections_4ints_3visits_.rds")
stan_out2 <- readRDS("./occupancy/model_outputs/bombus/non_urban/bombus_40km_1000minpop_1minUniqueDetections_4ints_3visits_.RDS")

species_names <- readRDS("./figures/species_names/bombus_names_10km_urban.RDS")
species_names2 <- readRDS("./figures/species_names/bombus_names_40km_nonurban.RDS")

fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
fit_summary2 <- rstan::summary(stan_out2)
View(cbind(1:nrow(fit_summary2$summary), fit_summary2$summary)) # View to see which row corresponds to the parameter of interest


list_of_draws <- as.data.frame(stan_out)
list_of_draws2 <- as.data.frame(stan_out2)

n_species <- length(species_names)
n_species2 <- length(species_names)

# remove species that didn't occur in the urban data set
# these occur in species_names2 but not in species_names

check <- which(is.na(charmatch(species_names2, species_names))) 
# looks like it's species 18 (insularis) occurs in non-urban sample but not in urban
# we will need to remove it for purposes of this post hoc analysis
# 84 is the number of columns before psi species starts
check_plus <- which(is.na(charmatch(species_names2, species_names))) + 84

species_names2 <- species_names2[-(c(check))]
which(is.na(charmatch(species_names2, species_names))) # should now be 0
list_of_draws2 <- list_of_draws2[,-(c(check_plus))] 

n_draws <- 100
n_posterior_samples <- 500

posterior_length <- nrow(list_of_draws)
posterior_length2 <- nrow(list_of_draws2)

correlation <- vector(length = n_draws)
samples <- matrix(nrow = n_species, ncol = 2)

post_posterior <- vector(length = n_draws * n_posterior_samples)

for(draw in 1:n_draws){
  
  for(species in 1:n_species){
      
    # for each species, draw a random sample from the posterior for:
    # psi0 non-urban habitat (starts at column 85) from list_of_draws2
    samples[species,1] <- list_of_draws2[draw,84+species] 
    # and psi species natural habitat (starts at column 152) 
    # from list_of_draws
    samples[species,2] <- list_of_draws[draw,
                                        151+species]
  }
  
  # calculate the correlation between the two terms across all species
  # and add to a vector of correlations, then
  # rerun the 
  samples <- as.data.frame(samples)
  
  # regression on rangewide occupancy as a predictor of species-specific effect of 
  # natural greenspace area. More negative slope means that species with large
  # range occupancy are not very dependent on natural greenspace, whereas rarer species
  # tend to rely on natural greenspace more substantially.
  stan_fit_draw <- rstanarm::stan_glm(V2~V1, data=samples, family = "gaussian")
  list_out <- as.data.frame(stan_fit_draw)[,2]
  list_out <- sample(list_out, n_posterior_samples)
  
  post_posterior[(1+((draw-1)*n_posterior_samples)):(draw*n_posterior_samples)] <- list_out
   
}

mean(post_posterior)
#returns the 22 and 77th percentiles of the input values
quantiles <- quantile(post_posterior, probs = c(0.025,0.25, 0.5, 0.75, 0.975))

hist(post_posterior, breaks = 20, 
     main = "",
     xlim = c(-1, 1),
     xlab = "estimated association between rangewide occupancy and \nspecies-specific effect of natural greenspace area",
     ylab = paste0("frequency in ", n_draws*n_posterior_samples, " draws from posterior distributions"))
abline(v = quantiles[3], col="black", lwd=3, lty=1)
abline(v = cbind(quantiles[1], quantiles[5]), col="blue", lwd=3, lty=2)
#abline(v = cbind(lower90, upper90), col="green", lwd=3, lty=2)
abline(v = cbind(quantiles[2], quantiles[4]), col="red", lwd=3, lty=2)


################################################################################
## repeat for hoverflies

## --------------------------------------------------
## Read in model run results


stan_out <- readRDS("./occupancy/model_outputs/large_files/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits_.rds")
stan_out2 <- readRDS("./occupancy/model_outputs/syrphidae/non_urban/syrphidae_40km_1000minpop_2minUniqueDetections_3ints_3visits_.rds")

species_names <- readRDS("./figures/species_names/syrphidae_names_10km_urban.RDS")
species_names2 <- readRDS("./figures/species_names/syrphidae_names_40km_nonurban.RDS")

fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
fit_summary2 <- rstan::summary(stan_out2)
View(cbind(1:nrow(fit_summary2$summary), fit_summary2$summary)) # View to see which row corresponds to the parameter of interest

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

n_draws <- 100
n_posterior_samples <- 500

posterior_length <- nrow(list_of_draws)
posterior_length2 <- nrow(list_of_draws2)

correlation <- vector(length = n_draws)
samples <- matrix(nrow = n_species, ncol = 2)

post_posterior <- vector(length = n_draws * n_posterior_samples)

for(draw in 1:n_draws){
  
  for(species in 1:n_species){
    
    # for each species, draw a random sample from the posterior for:
    # psi0 non-urban habitat (starts at column 85) from list_of_draws2
    samples[species,1] <- list_of_draws2[draw,13+species] 
    # and psi species natural habitat (starts at column 152) 
    # from list_of_draws
    samples[species,2] <- list_of_draws[draw,
                                        157+species]
  }
  
  # calculate the correlation between the two terms across all species
  # and add to a vector of correlations, then
  # rerun the 
  samples <- as.data.frame(samples)
  
  # regression on rangewide occupancy as a predictor of species-specific effect of 
  # natural greenspace area. More negative slope means that species with large
  # range occupancy are not very dependent on natural greenspace, whereas rarer species
  # tend to rely on natural greenspace more substantially.
  stan_fit_draw <- rstanarm::stan_glm(V2~V1, data=samples, family = "gaussian")
  list_out <- as.data.frame(stan_fit_draw)[,2]
  list_out <- sample(list_out, n_posterior_samples)
  
  post_posterior[(1+((draw-1)*n_posterior_samples)):(draw*n_posterior_samples)] <- list_out
  
}

mean(post_posterior)
#returns the 22 and 77th percentiles of the input values
quantiles <- quantile(post_posterior, probs = c(0.025,0.25, 0.5, 0.75, 0.975))

hist(post_posterior, breaks = 20, 
     main = "",
     xlim = c(-1, 1),
     xlab = "estimated association between rangewide occupancy and \nspecies-specific effect of natural greenspace area",
     ylab = paste0("frequency in ", n_draws*n_posterior_samples, " draws from posterior distributions"))
abline(v = quantiles[3], col="black", lwd=3, lty=1)
abline(v = cbind(quantiles[1], quantiles[5]), col="blue", lwd=3, lty=2)
#abline(v = cbind(lower90, upper90), col="green", lwd=3, lty=2)
abline(v = cbind(quantiles[2], quantiles[4]), col="red", lwd=3, lty=2)