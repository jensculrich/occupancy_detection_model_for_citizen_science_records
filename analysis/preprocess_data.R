### Prepeare data to feed to model
# jcu; started oct 27, 2022

# will need to assign occupancy intervals and visit numbers
# for occupancy (as opposed to abundance), will need to 
# filter down to one unique occurrence per species*site*interval*year

## The data we will need to prepare to feed to the model:
# V <- V # an array of detection data
# n_species # number of species
# n_sites # number of sites
# n_intervals # number of occupancy intervals 
# n_visits # number of visits in each interval

library(tidyverse)

# spatially explicit occurrence data
df <- read.csv("./data/data_urban_occurrences.csv")

## --------------------------------------------------
# assign occupancy intervals

era_start <- 2000 # We previously filtered records to those that occur after this year
era_end <- 2022 # latest year of data
n_intervals <- 7

total_years <- era_end - era_start
remainder <- total_years %% 7
years_per_interval <- (total_years - remainder) / n_intervals

test <- df %>%
  # assign year as - year after era_start
  mutate(occ_year = (year - 2000)) %>%
  # remove data from years that are in the remainder
  # occupancy intervals have to be equal in length for the model to process
  # e.g. we can't have intervals of 3 years and then one interval with only 1 year
  filter(!occ_year %in% (remainder:1))
  # now assign the years into 1:n_intervals


