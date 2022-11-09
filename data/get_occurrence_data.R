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

# write.csv(occurrences_df, "data_unfiltered.csv")

## --------------------------------------------------
# 2) POST-PROCESSING 
## --------------------------------------------------

occurrences_df <- read.csv("./data/data_unfiltered.csv")

## --------------------------------------------------
# Let's see what issues have cropped up in the dataset

# How many records are missing species-level Identification?
occurrences_df %>% select(species) %>% filter(species == "") %>% head()

# are there any records with listed issues
(occurrences_df %>% select(issue) %>% filter(issue != ""))
# many with 'taxon match fuzzy', 'rounded coordinate', 'recorded date unlikely'
# will eventually need to decide what to do with these issue points

# are any records missing lat/long?
is.na(occurrences_df$decimalLatitude) %>% head(100)
# no! no occurrences missing lat/long 
# (we included "hasCoordinate" == TRUE in our query)

## --------------------------------------------------
# where is the occurrence information coming from?
occurrences_df %>% count(institutionCode, sort = TRUE)
# most occurrences from iNaturalist
# the Natural History Museum of Los Angeles County (LACM) is the only other big contributor to these data
# Essig Museum of Entomology (EMEC) at UC Berkley is the next biggest contributor but my a wide difference

# the column basisOfRecord splits the records into either 
# 'HUMAN_OBSERVATION' - iNaturalist data
# or 'PRESERVED_SPECIMEN' - a museum record.
occurrences_df %>% count(basisOfRecord, sort = TRUE)
# these two different sources of data may have very different observation processes
# and therefore in the data analysis it would be wise to split the way these
# different types of data contribute to the likelihood.

## --------------------------------------------------
# Let's look at some basic summaries of the dataset

## How many species do we have records from in the download? 
## How many records per species?
species_counts <- occurrences_df %>% 
  count(species, sort = TRUE)

n_species <- nrow(species_counts)
head(species_counts)

# Let's add these total count of records for each species
occurrences_df <- occurrences_df %>%
  # add a count for each species
  group_by(species) %>%
  add_tally() %>%
  ungroup() %>%
  rename("total_count" = "n")

## --------------------------------------------------
# Filter occurrences to a time era/s of interest
# Plot number of occurrence records since some given time
min_year <- 2000 # what year to cut data for the plot
min_number_records <- 100 # threshold records for plot

contemporary_occurrences_df <- occurrences_df %>%
  filter(year > min_year)

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

occurrences_df %>% 
  # filter to after min year and count records since year
  filter(.$year > min_year) %>%
  count(species, sort = TRUE) %>% 
  drop_na(species) %>% 
  filter(n > min_number_records) %>% 
  
  # draw plot
  ggplot(aes(x = reorder(species, n), y = n, fill = species)) + 
  geom_bar(stat = "identity", show.legend = FALSE) + 
  labs(x = "Species with > 100 records since 2000", 
       y = "Number of occurrence records") + 
  coord_flip()

# how are the number of records changing through time
out <- occurrences_df %>% 
  
  # count for iNat vs Museum
  group_by(basisOfRecord) %>% 
  count(year) 
  
  # make a chrono plot
  ggplot(aes(x = year, y = n, col = as.factor(basisOfRecord))) + 
  xlim(1960, 2022) + # choose some years for the axes
  geom_line() +
  labs(y = "Total Number of Occurrence Records")  +
  scale_colour_manual(name = "Basis of Records", 
                        labels = c("iNat Observations", "Preserved Specimens"),
                        values=c("red","blue"))

plotly::ggplotly(out)

# Chronological record counts split by basis of record
ggplot(out, aes(x = year, y = n, col = as.factor(basisOfRecord))) + 
  xlim(1960, 2022) + # choose some years for the axes
  geom_line() +
  labs(y = "Number of Records")  +
  scale_colour_manual(name = "Basis of Records", 
                      labels = c("iNat Observations", "Preserved Specimens"),
                      values=c("red","blue")) +
  theme_bw() +
  xlab("Year") +
  ylab("Number of Records") +
  theme(legend.position = c(0.20, 0.8),
        legend.title = element_text(colour="black", size=14, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)
        )

## --------------------------------------------------
# 3) EXPLORATORY SPATIAL MAPPING
## --------------------------------------------------

# use data frame created above that includes only record after the year min_year
df <- contemporary_occurrences_df

# Rename Latitude and Longitude
df <- dplyr::rename(df, lat = decimalLatitude, 
                                  long = decimalLongitude) %>%
  # add a count of records per species
  group_by(species) %>%
  add_tally()

# Let's just map one 'typical' species for now and then one widely abundant species
row_per_species <- df %>%
  group_by(species) %>%
  filter(row_number()==1)

median_n <- median(row_per_species$n)
max_n <- max(row_per_species$n)

df_median_species <- df %>%
  # there are two species with the median value
  filter(species == "Dasysyrphus creper")

df_max_species <- df %>%
  # there are two species with the median value
  filter(n == max_n)

# plot records for a median species
p <- leaflet::leaflet(df_median_species) %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(~df_median_species$long, ~df_median_species$lat, 
             popup = paste("Species: ", df_median_species$species, "<br>",
                           "Year: ", df_median_species$year))
p 

# plot records for the max species
q <- leaflet::leaflet(df_max_species) %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(~df_max_species$long, ~df_max_species$lat, 
             popup = paste("Species: ", df_max_species$species, "<br>",
                           "Year: ", df_max_species$year))
q

# plot some random sampled reccords from all species
df_500 <- df %>%
  sample_n(500, replace = TRUE)

r <- leaflet::leaflet(df_500) %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(~long, ~lat, 
             popup = paste("Species: ", df_500$species, "<br>",
                           "Year: ", df_500$year, "<br>",
                           "Obs. per species: ", df_500$n))
r 
