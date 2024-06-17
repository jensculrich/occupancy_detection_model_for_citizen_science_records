Metadata for the data used in "Ulrich et al. (2023) Urban landscapes with more natural habitat support higher pollinator diversity. (manuscript submitted for consideration)".

#-------------------------------------------------------------------------------

syrphidae_data_all.csv - detections of hoverflies (syrphidae).

This dataset comprises 145,572 specimen records and spans 2000-2022. These records have been sourced from GBIF (https://www.gbif.org/).

GBIF.org. GBIF Occurrence Download. https://doi.org/10.15468/dl.nga26z (2023).
GBIF.org. GBIF Occurrence Download. https://doi.org/10.15468/dl.n5cmwv (2023).

The data is a csv file with 51 variables.
Please see the GBIF API for variable definitions: (https://www.gbif.org/developer/occurrence#parameters)

We made some changes to the data:
(1) "basisOfRecord"" for all iNaturalist data labelled as "community_science"
(2) "basisOfRecord"" for all records from digitized musuem collections data  labelled as "research_collection"
(3) Made the following taxonomic changes:
  - replace all Eumerus with Eumerus sp.
  - replace all Chrysogaster with Chrysogaster sp.
  - replace Eoseristalis (genus name) with Eristalis (genus name)
(4) add a unique row number - column 1, "id" (numeric)

#-------------------------------------------------------------------------------

syrphidae_nativity.csv - detections of hoverflies (syrphidae).

The data is a csv file with 3 variables: 

"species" species name (for row matching) (character)

"nativity" non-native to north america (== 0), all other species are considered native (numeric)
Nativity data obtained from https://bugguide.net/node/view/248891

"notes" notes if any (character)