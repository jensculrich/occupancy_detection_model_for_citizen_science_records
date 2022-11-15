---
title: "README.md"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## ./simulation/
simulate_model0.R - simulate some spatially-temporally heterogenous GBIF data for a taxonomic group of species.
The data is 'collected' over several occupancy intervals, with each site visited multiple times (repeat surveys) during the interval.
Occupancy intervals represent a set of years, and visits = the sum of all observations of the community per year.
Here, species are simulated to potentially occur at every site in the system - i.e., assuming no range restrictions;
and all species are potentially sampled at each site visit - i.e., assuming complete community sampling at each site visit.


## ./models/
model0.stan - model used to analyze the simulated data from simulate_simple.R. The uncertainty in the parameter estimates 
from the model output should encompass the point values used to simulate the simple data. We should also use this model to test how model
decisions, including priors, influence the models inference.
