---
title: "README.md"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Using NHC records to test whether urban land cover variables impact city-wide pollinator abundance and occupancy rate.  

The four high level folders include: the abundance-occupancy model, the more simple occupancy model, the data, and a methods outline.

### ./abundance-occupancy/
The model presented is an abundance-occupancy model that utilizes the link between abundance and occupancy to more precisely estimate the relative abundance of multiple species in pollinator taxonomic groups (separately for Syrphidae (hover fly family) and Bombus (bumble bee genus)) at large, city-wide scale. The relative abundances of species are compared among urban areas to test whether hypothesized urban habitat variables significantly effect pollinator populations at the city-wide scale. The abundance-occupancy models are held in the folder './abundance-occupancy/'. For more details on running the abundance-occupancy folder see below...

### ./occupancy/
Prior to constructing the abundance-occupancy model, I constructed a more straightforward integrated multi-species occupancy model. The occupancy model tests whether the hypothesized urban drivers impact presence or absence of species at the city-wide scale. Although more simple to run and more direct in its interpretation, the occupancy model provides a more coarse understanding of how these urban drivers might be affecting pollinator populations, as it may take many years or even decades for negative effects to cause species-level extinction at broad, city-wide spatial scales considered here. The occupancy models are held in the folder './occupancy/'. For more details on running the abundance-occupancy folder see below...

### ./data/
Includes the occurrence and spatial data used in the analysis. As of now, this file is ignored in git commits, because the spatial data files are too big to transfer back and forth between a local computer and github. Eventually I will prep all of the data into a zip file that could be included as a supplement to the manuscript that could be downloaded and then unzipped here so that others can run the analysis from remote computers.

#### get_occurrence_data.R
Prepare a request for NHC data from GBIF. The download will be stored in the data folder. To limit the requests from GBIF, this call was only used once, and is not included in the model building/running process. The file also enables visualization and summary of the downloaded data. This second half of the file is not called in the model building/running process and is purely for visualization.   

### ./occurrence_data/

#### ./syrphidae/

#### ./bombus/

### ./spatial_data/
The spatially explicit environmental data

#### ./impervious_surface/
#### 300m_aggregate_crop_nlcd_2016_impervious.tif
2016 urban imp surface cover from https://www.mrlc.gov/data/nlcd-2016-percent-developed-imperviousness-conus.
The data was originally at a 30m x 30m grain. I aggregated the dataset up to averages at 300m x 300m to speed downstream processing times.

#### ./population_density/
#### gpw_v4_population_density_rev11_2015_30_sec.tif
2015 population density from
MRLC data https://www.mrlc.gov/data
30 arc second grain, extent is across all of world

#### ./time_integrated_ndvi
#### TIN2014_wUSAeM250m_v2.tif
2014 NDVI metric from 
https://earthexplorer.usgs.gov/ 
> Vegetation Monitoring > Phenology > eMODIS Phenology > 250m res - 2014 - Western NA V2


### ./methods/
Includes a detailed word document ("methods.doc") outlining my methods for data collection and analysis.



# Model running details

## ./abundance-occupancy/
To fit the model to simulated data call the model using rstan::stan() from ./simulation/simulate_abundance.R.
To fit the model to real data call the custom 'prep_data()' function in ./analysis/run_analysis.R after defining study dimensions (temporal grain, spatial grain, minimum threshold for considering a site to be 'urban', minimum records per species, etc.), and then call the model from within this file using rstan::stan().
All other files are either called in a chain from the 'prep_data()' function, are auxiliary, or are used for post-processing or figure production.

### ./simulation/
#### simulate_abundance.R
Start with simulating a dataset that includes hypothetical abundance counts from a citizen science NHC dataset and hypothetical presence/absence records from the museum records NHC dataset generated from a known hypothetical underlying ecological state and set of observation processes. Fitting the model to the simulated data should accurately estimate the parameter choices used to generate the data.

#### _.RDS
A saved simulation output that can be assessed at a later time. Note, the parameter estimates will depend on the 'targets' specified in the specific run of ./simulation/simulate_abundance.R

### ./models/
#### model_abundance.stan
This file contains the model. Call rstan::stan() from ./simulation/simulate_abundance.R or from ./analysis/run_analysis.R to compile and fit the model. rstan::stan() will return a stanfit object.

### ./model_outputs/
#### _.RDS
A saved model fit output that can be assessed at a later time, or used for figures etc.

### ./data_prep/
#### explore_spatial_data.R
A data exploration file. Use to view plots of the spatial land cover data and the spatial NHC records overlayed on the landscape. This file is not called in the model building/running process and is purely for visualization.


#### get_spatial_data.R
Contains a function to prepare the spatial data based on the study design dimensions chosen in run_analysis.R. Function will also then match NHC records to spatial sites. 

#### get_species_ranges.R
Contains a function to prepare an indicator array that tell will tell stan which sites are within range for each species.

### ./analysis/
#### prep_data.R
Contains a function to prepare the data given the study design dimension choices. This function will call lower level functions in get_spatial_data.R and get_species_ranges.R to obtain the spatial data requested.

#### run_analysis.R
First choose study design dimensions. Then call 'prep_data()' from source ./analysis/prep_data.R based on these choices. The returned data can then be bundled into a form that can be read by stan. Call rstan::stan() with the abundance model in the directory and pointing to the bundled data to fit the model.

#### posterior_predictive_check.R
Conduct a posterior predictive check (chi-squared discrepancy test) on the stanfit object.


## ./occupancy/
Folder is laid out identically to the ./abundance-occupancy/ folder. However, the simulation, real data preparation, and model are all built to hand presence/absence NHC records as opposed to abundance counts.


## ./data/
Note, download the compressed folder and unzip it here at this directory location to access and manipulate the data using the data prep workflow held in ./abundance-occupancy/or ./occupancy/

### ./occurrence_data/
NHC occurrence data obtained from [GBIF](https://www.gbif.org/)

### ./syrphidae/

#### data_unfiltered_CA.csv
Syrphid NHC data from California.

### ./bombus/

### ./spatial_data/
### ./population_density/
#### gpw_v4_population_density_rev11_2015_30_sec.tif
Human population density. Dataset obtained from the [NASA Socioeconomic Data and Applications Center](https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download) was collected in 2015. Data are at a resolution of 30 seconds (about 500m by 500m), and are stored in the WGS84 latitude/longitude (EPSG:4326) geographic coordinate system. 

### ./impervious_surface/
#### nlcd_2016_impervious_aggregated_300m.tif
Impervious surface across the study extent. File was pre-cropped to match the extent and then aggregated to average value per 300m. The aggregation is to speed up the site variable extraction (and save energy and computing resources)  when the spatial and temporal grains are adjusted across a loop for the sensitivity analysis. The original dataset obtained from the [the USGS National Landcover Database](https://www.mrlc.gov/data/nlcd-2016-percent-developed-imperviousness-conus) has an extent across western United States. The Dataset was collected in 2016 and represents the percent of developed surface over every 30-meter pixel in the extent, and are stored in the Albers Conical Equal Area (EPSG:9822) geographic coordinate system.

### ./plant_cover/
#### .tif
_ plant cover across the study extent. File was pre-cropped to match the extent and then aggregated to average value per 300m. The aggregation is to speed up the site variable extraction (and save energy and computing resources)  when the spatial and temporal grains are adjusted across a loop for the sensitivity analysis. The original dataset obtained from the [the USGS National Landcover Database]() has an extent across western United States. The Dataset was collected in ... and represents the percent of developed surface over every 30-meter pixel in the extent, and are stored in the _ geographic coordinate system.

### ./california_urban_areas
Currently not using this data, it was originally used to map specific city names to sites.
#### _.shp
Administrative area place names



