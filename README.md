---
title: "README.md"
output:
  html_document: default
---

# Use NHC records to test whether urban landscape habitat quantity and quality predict urban pollinator occurrence.

(under construction)

To run an analysis, open ./occupancy/analysis/run_model.R.
you will be prompted to enter details that will specify the model run, including:
taxonomic group, spatial grain and temporal divisions of occupancy intervals.

This repository holds three public subdirectories: ./data/, ./occupancy/, and ./figures. 

### ./occupancy/
Contains simulation files (./simulation/), data prep(./data_prep/), models (./models/), model implementation (./analysis/) and model outputs (./model_outputs/).

Run a model by specifying the data level constraints (taxonomic group, spatial grain and temporal divisions of occupancy intervals) and by tweaking the MCMC settings if desired.

The models used for the main analysis (held in the models subdirectory) are labelled "model_syrphidae.stan" and "model_bombus.stan". We also used a simplified model to estimate range wide occurrence rate to conduct a post hoc correlation analysis with species-specific effects of range wide occurrence (a proxy for habitat generality) and species-specific association with natural habitat area in urban landscapes. The folder also holds a reparameterized random effects model (non-centered random effects). We explored whether this parameterization resulted in more stable scale estimates for site-specific random effects but do not use this for our final analysis.




