### Prepeare data to feed to model
# jcu; started oct 27, 2022

# will need to assign occupancy intervals and visit numbers
# for occupancy (as opposed to abundance), will need to 
# filter down to one unique occurrence per species*site*interval*year

library(tidyverse)

# spatially explicit occurrence data
df <- read.csv("./data/data_urban_occurrences.csv")

## --------------------------------------------------
# 

