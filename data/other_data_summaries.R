library(tidyverse)

## SYRPHIDAE
# how many different species per average museum sample (if a museum sample occurs)?
# how many different genera per average museum sample (if a museum sample occurs)?

# read either the syrphidae data 
df <- read.csv(paste0("./data/occurrence_data/", "syrphidae", "_data_all.csv"))

## --------------------------------------------------
# Prep the data

# make the df into a spatial file
df$decimalLongitude <- na_if(df$decimalLongitude, '')
df$decimalLatitude <- na_if(df$decimalLatitude, '')

df_museum_species <- df %>%
  
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude)) %>%
  
  # remove records (if any) missing species level identification
  filter(species != "") %>%
  
  # filter to citizen science data only
  filter(basisOfRecord == "research_collection") %>% 
  
  # determine whether a community sampling event occurred
  # using collector name might be overly conservative because for example
  # the data includes recordedBy == J. Fulmer *and* recordedBy J. W. Fulmer
  # instead grouping by collections housed in the same institution from the same year
  # within a site
  #group_by(institutionCode, year, grid_id) %>%
  #mutate(n_species_sampled = n_distinct(species)) %>%
  # filter(n_species_sampled >= min_species_for_community_sampling_event) %>%
  
  # one unique row per site*species*occ_interval*visit combination
  group_by(year, institutionCode) %>% 
  mutate(n_species = n_distinct(species)) %>%
  slice(1) %>%
  ungroup()

mean(as.numeric(df_museum_species$n_species))


df_museum_genera <- df %>%
  
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude)) %>%
  
  # remove records (if any) missing species level identification
  filter(species != "") %>%
  
  # filter to citizen science data only
  filter(basisOfRecord == "research_collection") %>% 
  
  # determine whether a community sampling event occurred
  # using collector name might be overly conservative because for example
  # the data includes recordedBy == J. Fulmer *and* recordedBy J. W. Fulmer
  # instead grouping by collections housed in the same institution from the same year
  # within a site
  #group_by(institutionCode, year, grid_id) %>%
  #mutate(n_species_sampled = n_distinct(species)) %>%
  # filter(n_species_sampled >= min_species_for_community_sampling_event) %>%
  
  # one unique row per site*species*occ_interval*visit combination
  group_by(year, institutionCode) %>% 
  mutate(n_genera = n_distinct(genus)) %>%
  slice(1) %>%
  ungroup()

mean(as.numeric(df_museum_genera$n_genera))

(ratio = 
    mean(as.numeric(df_museum_genera$n_genera)) / 
    mean(as.numeric(df_museum_species$n_species)))

## --------------------------------------------------
# How much uncertainty in coordinate estimates?

## --------------------------------------------------
# How much uncertainty in coordinate estimates given that they reported 
# and after we've filtered out > 10km

## SYRPHIDS
df <- read.csv(paste0("./data/occurrence_data/", "syrphidae", "_data_all.csv"))

temp1 <- df %>%
  filter(!is.na(coordinateUncertaintyInMeters))

nrow(temp1) / nrow(df)

temp2 <- df %>%
  filter(coordinateUncertaintyInMeters < 10000,
         institutionCode == "iNaturalist")

(mean(temp2$coordinateUncertaintyInMeters))
(sd(temp2$coordinateUncertaintyInMeters))
(median(temp2$coordinateUncertaintyInMeters))

temp3 <- df %>%
  filter(coordinateUncertaintyInMeters < 10000,
         institutionCode != "iNaturalist")

(median(temp3$coordinateUncertaintyInMeters))

rm(df, temp1, temp2, temp3)
gc()

## BUMBLE BEES
df <- read.csv(paste0("./data/occurrence_data/", "bbna_private/", "bbna_trimmed.csv"))


## --------------------------------------------------
# How many contributors???

## --------------------------------------------------

## SYRPHIDS
df <- read.csv(paste0("./data/occurrence_data/", "syrphidae", "_data_all.csv"))

df <- df %>%
  filter(institutionCode == "iNaturalist",
         year > 2000) 

# how many identification contributors for community science data
temp <- df %>%
  group_by(identifiedBy) %>%
  add_tally()

max(temp$n)
max(temp$n) / nrow(df)

temp <- df %>%
  group_by(identifiedBy) %>%
  add_tally() %>%
  slice(1) %>%
  ungroup()

temp2 <- temp %>%
  filter(n > 99)
