### Data Collection for pollinator occurrence data
# jcu; started oct 27, 2022

# Collect data of pollinator occurrence from GBIF
# data should be limited to a taxonomic group - Syrphidae
# to a time period - 2001-2022
# and to a geographic area - state of California

# sections inculde:
# 1: GBIF DOWNLOAD
# 2: POST-PROCESSING
# 3: EXPLORATORY MAPPING

library(tidyverse)
library(rgbif)
library(plotly) # for plotting occurrence change through time
library(leaflet) # for mapping GBIF data (exploratory map viewing)

## --------------------------------------------------
# 1) GBIF DOWNLOAD
## --------------------------------------------------

## --------------------------------------------------
# Enter GBIF credentials

# Enter your GBIF user info here before proceeding

# I will delete this before pushing this file to protect my account privacy
# must be rewritten every time running this script

user="jensj27" 
pwd="Ceratina_1802" 
email="jensj27@gmail.com"

## --------------------------------------------------
# Enter download filters
taxonKey <- 6920 # Family - Syrphidae
# taxonKey <- 1540915 # test one species - Allograpta obliqua
stateProvince <- "USA.5_1" # Administrative area of California
basisOfRecord <- c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION') # museum specimens and human obs (i.e. inaturalist)
hasCoordinate <- TRUE # excludes records without coordinates
# hasGeospatialIssue <- FALSE  

## --------------------------------------------------
# Download data from GBIF 
# use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
down_code = occ_download(
  pred("taxonKey", taxonKey),
  pred_in("basisOfRecord", basisOfRecord),
  pred("gadm", stateProvince), 
  pred("hasCoordinate", hasCoordinate),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email)

# now get the records from your GBIF downloads page
setwd("./data")
occurrences <- occ_download_get(down_code[1], overwrite = TRUE)
occurrences_df <- occ_download_import(occurrences)

## --------------------------------------------------
# 2) POST-PROCESSING 
## --------------------------------------------------

# filter options
era_start <- 2000 # filter to records AFTER this year
# min_records_per_species <- 5 # filter to species with at least this many records

## --------------------------------------------------
# filter out specimens that weren't ID'd to species level

# remove all occurrences without any species designation assigned
occurrences_df <- occurrences_df %>%
  filter(species != "")

# are any records missing lat/long?
is.na(regular_contemporary_occurrences_df$decimalLatitude) %>% head(100)
# no! no occurrences missing lat/long 
# (we included "hasCoordinate" == TRUE in our query)
  
## How many species do we have records from in the download? 
## How many records per species?
species_counts <- occurrences_df %>% 
  count(species, sort = TRUE)

n_species <- nrow(species_counts)
head(species_counts)

## --------------------------------------------------
# Filter occurrences to a time era/s of interest

## Let's filter the occurrences to a time range that we are interested in 
# here we will use data from 2001 to present (22 years of data)
contemporary_occurrences_df <- occurrences_df %>%
  filter(year > era_start)

# how many occurrences did we exclude?
nrow(occurrences_df) - nrow(contemporary_occurrences_df)
# ~5,000 occurrences excluded
# what proportion of the total occurrences are from the contemporary era?
1 - ((nrow(occurrences_df) - nrow(contemporary_occurrences_df)) / nrow(occurrences_df))
# ~85% of data is from after the year 2000

## How many species do we have records from in the download IN THE CONTEMPORARY ERA? 
## How many records per species IN THE CONTEMPORARY ERA?
contemporary_species_counts <- contemporary_occurrences_df %>% 
  count(species, sort = TRUE)

n_species_contemporary <- nrow(contemporary_species_counts)
head(contemporary_species_counts)

# Let's say we want to exclude taxa with fewer than some arbitrary number of small 
# observations, min_records_per_species, which may be e.g. only transiently occurring in the region
# or not enough observations to indicate that this is a solid species concept
# we'll get a df of species that occur slightly more regularly
regular_contemporary_occurrences_df <- contemporary_occurrences_df %>%
  # add a count for each species
  group_by(species) %>%
  add_tally() %>%
  ungroup() # %>%
  # then filter out species who's count is less than min_records_per_species
  # we won't filter now, rather we'll filter after we fit the data to our geographic sites,
  # which likely will not encompass the entire geographic area of the downloads (California)
  # filter(n > min_records_per_species) 

regular_contemporary_occurrences_df %>% 
  count(species, sort = TRUE) %>% 
  drop_na(species) %>% 
  filter(n > 100) %>% 
  ggplot(aes(x = reorder(species, n), y = n, fill = species)) + 
  geom_bar(stat = "identity", how.legend = FALSE) + 
  labs(x = "Species with > 100 observations since 2000", y = "Number of Occurrence Records (observations)") + 
  coord_flip()

# where is the occurrence information coming from?
regular_contemporary_occurrences_df %>% count(institutionCode, sort = TRUE)
# most occurrences from iNaturalist
# the Natural History Museum of Los Angeles County (LACM) is the only other big contributor to these data
# Essig Museum of Entomology (EMEC) at UC Berkley is the next biggest contributor but my a wide difference

# how are the number of records changing through time
out <- regular_contemporary_occurrences_df %>% 
  count(year) %>% 
  ggplot(aes(x = year, y = n, group = 1)) + 
  xlim(2000, 2022) + 
  geom_line() +
  labs(y = "Total Number of Occurrence Records")
plotly::ggplotly(out)

regular_contemporary_occurrences_df %>% select(issue) %>% head()

# write.csv(occurrences_df, "data_unfiltered.csv")
# write.csv(regular_contemporary_occurrences_df, "data_filtered.csv")

## --------------------------------------------------
# 3) EXPLORATORY MAPPING
## --------------------------------------------------

# setwd("..")
df <- read.csv("./data/data_filtered.csv")

# Rename Latitude and Longitude
df <- dplyr::rename(df, lat = decimalLatitude, 
                                  long = decimalLongitude)

# Let's just map one 'typical' species for now
row_per_species <- df %>%
  group_by(species) %>%
  filter(row_number()==1)

median_n <- median(row_per_species$n)

df_median_species <- df %>%
  filter(n == median_n)

p <- leaflet::leaflet(df_median_species) %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(~long, ~lat, 
             popup = paste("Species: ", df_median_species$species, "<br>",
                           "Year: ", df_median_species$year))
p 

df_500 <- df %>%
  sample_n(500)

q <- leaflet::leaflet(df_500) %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(~long, ~lat, 
             popup = paste("Species: ", df_500$species, "<br>",
                           "Year: ", df_500$year, "<br>",
                           "Obs. per species: ", df_500$n))
q 
