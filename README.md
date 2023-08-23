---
title: "README"
output:
  pdf_document: default
  html_document: default
---

# Use NHC records to test whether urban landscape habitat quantity and quality predict urban pollinator occurrence.

(under construction)

To run an analysis, open ./occupancy/analysis/run_model.R.
you will be prompted to enter details that will specify the model run, including:
taxonomic group, spatial grain and temporal divisions of occupancy intervals.

This repository holds three public subdirectories: ./data/, ./occupancy/, and ./figures. 

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 20px;}
</style>
<div class = "blue">
## ./occupancy/
</div>
<br>
Contains simulation files (./simulation/), data prep(./data_prep/), models (./models/), model implementation (./analysis/) and model outputs (./model_outputs/).

The models used for the main analyses (held in the models subdirectory) are labelled "model_syrphidae.stan" and "model_bombus.stan". I also used a simplified model to estimate range wide occurrence rate ("model_taxon_simple.stan") to conduct a post hoc correlation analysis with species-specific effects of range wide occurrence (a proxy for habitat generality) and species-specific association with natural habitat area in urban landscapes. The folder also holds a reparameterized random effects model (non-centered random effects; "model_taxon_reparameterized_rand_effects.stan"). I explored whether this parameterization resulted in more stable scale estimates for site-specific random effects but since it did not, I do not use this reparametrization for the final analyses.

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 20px;}
</style>
<div class = "blue">
#### How to run a model:
</div>
<br>
Run a model (./analysis/run_model.R) by specifying the data level constraints (taxonomic group, spatial grain and temporal divisions of occupancy intervals) and by tweaking the HMC settings if desired.

After specifiying the data level constraints at the top of the file, you will be offered to prepare data for an analysis (prep_data() function). This preparation ~10 minutes to run, so if it's already been done (and the data has been saved), you can go ahead and skip down to load previously saved data and then enter the HMC setting before running the model. Prepared data .rds files are held in the ./analysis/prepped_data/ folder

There are a few quick model diagnostic check functions listed after the model run call so that you can conveniently check the model run properties in this same file.  

#### Data prep (./data_prep/)
The function prep_data() (called by run_model.R) is held in prep_data.R. After being called, this function will then communicate with get_spatial_data.R to define sites and gather covariate data before attributing detections to sites and inferring a sampling process. The function will also communicate with get_species_ranges.R to determine which sites are in range and which are out of range (should be treated as NA's) based on each species's distribution.

#### Model outputs (./model_outputs/)
These are the saved HMC runs.

#### Simulation (./simulation/)
Simulate a community of pollinators inferring the same ecological process and observation process. Also provides options to "break" some of the assumptions of our model and observe the consequences on the parameters. For example, our model ignores the possibility of research collections community surveys recovering zero species. We can see how consequential this is for the results by simulating some "hidden" community surveys that are overlooked. 

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 20px;}
</style>
<div class = "blue">
## ./data/
</div>
<br>
Includes both spatial data and occurrence records data (need to go through and see what is private v public)

### Spatial Data (./spatial_data/): 

#### Ecoregion 1 (./na_cec_eco_l1/)
Contains shapefile for ecoregion 1 in the north america (broadest spatial clustering unit for analysis).
add source

#### Ecoregion 3 (./NA_CEC_Eco_Level3/)
Contains shapefile for ecoregion 3 in the north america (intermediate spatial clustering unit for analysis).
add source

#### Metropolitan areas (./tl_2019_us_cbsa/)
Contains shapefile for metroplitan areas in the United States.
add source

#### Population Density (./population_density/)
Contains a raster with population density at _ resolution from _ year.
add source

#### Land cover (./land_cover/)
Contains a raster with land cover data (NLCD data) at 30m resolution from 2016. 
add source

### Occurrence Data (./occurrence_data/): 
Contains occurrence data for bumble bees (BBNA (private folder due to data rights and size)) and for hoverflies (from GBIF (private folder due to size))

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 20px;}
</style>
<div class = "blue">
## ./figures/
</div>
<br>
make figures for the manuscript and/or for presentations
