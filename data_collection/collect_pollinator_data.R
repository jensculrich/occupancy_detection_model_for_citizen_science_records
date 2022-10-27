### Data Collection for pollinator occurrence data
# jcu; started oct 10, 2022

# Collect data of pollinator occurrence from GBIF
# data should be limited to a taxonomic group - Syrphidae
# to a time period - 2001-2022
# and to a geographic area - state of California

library(tidyverse)
library(rgbif)

## --------------------------------------------------
# Enter GBIF credentials

# Enter your GBIF user info here before proceeding
# user=user 
# pwd=pwd 
# email=email
# I will delete this before pushing this file to protect my account privacy
# must be rewritten every time running this script
user="jensj27" 
pwd="" 
email="jensj27@gmail.com"

## --------------------------------------------------
# Enter download filters
taxonKey <- 6920 # Family - Syrphidae
stateProvince <- "USA.5_1" # Administrative area of California
basisOfRecord <- c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION') # museum specimens and human obs (i.e. inaturalist)
hasCoordinate <- TRUE # excludes records without coordinates
# hasGeospatialIssue <- FALSE  

## --------------------------------------------------
# Download data from GBIF 
down_code = occ_download(
  pred_in("taxonKey", taxonKey),
  pred_in("basisOfRecord", basisOfRecord),
  pred("stateProvince", stateProvince), 
  pred("hasCoordinate", hasCoordinate),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email)


