### Data Collection for pollinator occurrence data
# jcu; started oct 27, 2022

# Collect data of pollinator occurrence from GBIF
# data should be limited to a taxonomic group - Syrphidae
# to a time period - 2001-2022
# and to a geographic area - state of California or beyond

# sections inculde:
# 1: GBIF DOWNLOAD
# 2: POST-PROCESSING
# 3: EXPLORATORY MAPPING

library(tidyverse)
library(rgbif)
# library(plotly) # for plotting occurrence change through time
# library(leaflet) # for mapping GBIF data (exploratory map viewing)

## --------------------------------------------------
# 1) GBIF DOWNLOAD
## --------------------------------------------------


## --------------------------------------------------
# Enter GBIF credentials

# Enter your GBIF user info here before proceeding

# I will delete this before pushing this file to protect my account privacy
# must be rewritten every time running this script

user="jensj27" 
pwd="" 
email="@"

## --------------------------------------------------
# Syrphidae

# Syrphidae West
## --------------------------------------------------
# Enter download filters
taxonKey <- 6920 # Family - Syrphidae
# Administrative areas 
stateProvince <- c("USA.5_1", "USA.38_1", "USA.48_1", "USA.48_1", "USA.13_1", "USA.29_1",
                   "USA.3_1", "USA.32_1", "USA.45_1", "USA.6_1", "USA.44_1", "USA.51_1",
                   "USA.37_1", "USA.17_1", "USA.42_1", "USA.35_1", "USA.24_1", "USA.28_1",
                   "USA.19_1", "USA.26_1", "USA.4_1", "USA.16_1", "USA.50_1", "USA.23_1") 
basisOfRecord <- c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION') # museum specimens and human obs (i.e. inaturalist)
hasCoordinate <- TRUE # excludes records without coordinates
years <- seq(2000, 2022, 1)

## --------------------------------------------------
# Download data from GBIF 
# use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
down_code = occ_download(
  pred("taxonKey", taxonKey),
  pred_in("basisOfRecord", basisOfRecord),
  pred_in("gadm", stateProvince), 
  pred("hasCoordinate", hasCoordinate),
  pred_in("year", years),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email)

# now get the records from your GBIF downloads page
setwd("./data/occurrence_data/")
occurrences <- occ_download_get(down_code[1], overwrite = TRUE)
occurrences_df <- occ_download_import(occurrences)

# https://www.gbif.org/occurrence/download/0266686-220831081235567
# write.csv(occurrences_df, "syrphidae_data_west.csv")

# Syrphidae East
## --------------------------------------------------
# Enter download filters
taxonKey <- 6920 # Family - Syrphidae
# Administrative areas
stateProvince <- c( "USA.20_1", "USA.30_1", "USA.46_1", "USA.22_1", "USA.7_1", "USA.40_1",
                    "USA.33_1", "USA.39_1", "USA.36_1", "USA.31_1", "USA.8_1", "USA.21_1",
                    "USA.9_1", "USA.49_1", "USA.47_1", "USA.34_1", "USA.41_1", "USA.11_1",
                    "USA.10_1", "USA.1_1", "USA.25_1", "USA.43_1", "USA.18_1", "USA.15_1",
                    "USA.14_1") 
basisOfRecord <- c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION') # museum specimens and human obs (i.e. inaturalist)
hasCoordinate <- TRUE # excludes records without coordinates
years <- seq(2000, 2022, 1)

## --------------------------------------------------
# Download data from GBIF 
# use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
down_code = occ_download(
  pred("taxonKey", taxonKey),
  pred_in("basisOfRecord", basisOfRecord),
  pred_in("gadm", stateProvince), 
  pred("hasCoordinate", hasCoordinate),
  pred_in("year", years),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email)

# now get the records from your GBIF downloads page
setwd("./data/occurrence_data/")
occurrences <- occ_download_get(down_code[1], overwrite = TRUE)
occurrences_df <- occ_download_import(occurrences)

# https://www.gbif.org/occurrence/download/0266689-220831081235567

# write.csv(occurrences_df, "syrphidae_data_east.csv")

## --------------------------------------------------
# Combine data downloads

occurrences_df_all <- rbind(read.csv("./data/occurrence_data/syrphidae_data_east.csv"),
                             read.csv("./data/occurrence_data/syrphidae_data_west.csv"))

df_updated <- occurrences_df_all %>%
  # replace all Eumerus with Eumerus sp.
  mutate(species = ifelse(genus == "Eumerus", "Eumerus sp.", species)) %>%
  # replace all Chrysogaster with Chrysogaster sp.
  mutate(species = ifelse(genus == "Chrysogaster", "Chrysogaster sp.", species)) %>%
  # replace Eoseristalis (genus name) with Eristalis (genus name)
  mutate(species = gsub("Eoseristalis", "Eristalis", species)) %>%
  
  
  # update names to match BBNA dataset (community_science or research_collection)
  mutate(basisOfRecord = gsub("HUMAN_OBSERVATION", "community_science", basisOfRecord)) %>%
  mutate(basisOfRecord = gsub("PRESERVED_SPECIMEN", "research_collection", basisOfRecord)) %>%
  
  rename("id" = "X")

# write.csv(df_updated, "./data/occurrence_data/syrphidae_data_all.csv")


## --------------------------------------------------
# Need to add a unique id column

df <- read.csv("./data/occurrence_data/syrphidae_data_all.csv")

df <- df %>% mutate(id = row_number())

# write.csv(df, "./data/occurrence_data/syrphidae_data_all.csv")

df <- read.csv("./data/occurrence_data/bombus_data_all.csv")

df <- df %>% mutate(id = row_number())

# write.csv(df, "./data/occurrence_data/bombus_data_all.csv")


## --------------------------------------------------
# 3) DATA EXTRACTION FROM BBNA DATA - BOMBUS
## --------------------------------------------------

occurrences_df <- read.csv("./data/occurrence_data/bbna_private/all.shareable.bbna.12.07.2022.csv")

min_year <- 2000 # what year to cut data for the plot
min_number_records <- 0 # threshold records for plot

## --------------------------------------------------
# Let's see what issues have cropped up in the dataset

# How many records are missing species-level Identification?
species_names <- occurrences_df %>% select(species) %>% filter(species == "") %>% nrow()

# what institutions
institution_names <- occurrences_df %>% 
  dplyr::select(Institution) %>% 
  group_by(Institution) %>% 
  slice(1) %>%
  ungroup()

citsci_institutions <- c("BeeSpotter online citizen science portal: www.beespotter.org",
                         "bumblebeewatch.org", "Direct observation", 
                         "Flickr", "iNaturalist", "photo only", "photos", "photos/ observation only",
                         "www.beespotter.org", 
                         "The Xerces Society for Invertebrate Conservation; Bumble Bee Watch")

# and the research data only 
temp <- occurrences_df %>%
  filter(country == "USA")

# what data.sources
source_names <- temp %>% 
  dplyr::select(data.source) %>% 
  group_by(data.source) %>% 
  slice(1) %>%
  ungroup()

# www.bumblebeewatch.org 05-24-2017
# Bumble Bee Watch 10-26-2022
# BeeSpotter 10-03-2018
# iNaturalist 12-01-2022
# Flickr
# Xerces Society--dates various
# Xerces Society 11/14/2013
# Xerces Society 2012 Citizen Science
# Xerces Society B. franklini 02/26/2012

citsci_data_sources <- c("www.bumblebeewatch.org 05-24-2017", "Bumble Bee Watch 10-26-2022",
                         "BeeSpotter 10-03-2018", "iNaturalist 12-01-2022", "Flickr", 
                         "Xerces Society--dates various", "Xerces Society 11/14/2013",
                         "Xerces Society 2012 Citizen Science", "Xerces Society B. franklini 02/26/2012")


# and the citsci data only 
citsci_data <- temp %>%
  filter(data.source %in% citsci_data_sources) %>%
  mutate(basisOfRecord = "community_science")
# and the research data only 
collections_data <- temp %>%
  filter(!data.source %in% citsci_data_sources)  %>%
  mutate(basisOfRecord = "research_collection")

temp <- rbind(citsci_data, collections_data)

rm(citsci_data, collections_data, occurrences_df)
gc()

temp <- temp %>%
  filter(year >= min_year)
gc()

temp <- temp %>%
  mutate(species = str_replace(species, "vancouverensis", "bifarius"))


# trimmed to unly USA records, from 2000 or later, 
# and renamed any vancouverensis as bifarius
# later also manually trimmed any occurrences from Alaska
write.csv(temp, "./data/occurrence_data/bbna_private/bbna_trimmed.csv")
