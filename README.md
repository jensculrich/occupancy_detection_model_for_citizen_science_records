---
title: "README.md"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### Using NHC records to test whether impervious surface and plant cover in urban areas impact pollinator abundance and occupancy.  

## ./abundance-occupancy/
The model presented is an abundance-occupancy model that utilizes the link between abundance and occupancy to more precisely estimate the relative abundance of multiple species in pollinator taxonomic groups (separately for Syrphidae (hover fly family) and Bombus (bumble bee genus)). The relative abundances of species are compared among cities to test whether hypothesized urban habitat variables significantly effect pollinator populations at the city-wide scale. The abundance-occupancy models are held in the folder './abundance-occupancy/'. For more details on the abundance-occupancy folder see below...

## ./occupancy/
Prior to constructing the abundance-occupancy model, I constructed a more straightforward integrated multi-species occupancy model. The occupancy model tests whether the hypothesized urban drivers impact presence or absence of species at the city-wide scale. Although more simple to run and more direct in its interpretation, the occupancy model provides a more coarse understanding of how these urban drivers might be affecting pollinator populations, as it may take many years or even decades for negative effects to cause species-level extinction at broad, city-wide spatial scales considered here. The occupancy models are held in the folder './occupancy/'. For more details on the abundance-occupancy folder see below...

## ./methods/
Includes a detailed word document ("methods.doc") outlining my methods for data collection and analysis.


## ./data/
Includes the occurrence and spatial data used in the analysis. As of now, this file is ignored in git commits, because the spatial data files are too big to transfer back and forth between a local computer and github. Eventually I will prep all of the data into a zip file that could be included as a supplement to the manuscript that could be downloaded and then unzipped here so that others can run the analysis from remote computers.

## ./abundance-occupancy/
### ./simulation/
#### simulate_abundance.R
Simulate a dataset that includes abundance counts from a citizen science NHC dataset and presence/absence records from the museum records NHC dataset