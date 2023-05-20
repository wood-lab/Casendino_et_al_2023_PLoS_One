# Casendino_et_al_2023_PLoS_One

Data and code associated with Casendino et al. (2023) 
"Two decades of change in sea star abundance at a subtidal site in Puget Sound, Washington"

Data are available in this repository, and are also permanently archived in Dryad:

Casendino, Helen et al. (2023), Two decades of change in sea star abundance at a subtidal site in Puget Sound, Washington, Dryad, Dataset, https://doi.org/10.5061/dryad.cz8w9gj7q


## Files included

1. invertebrates_dataset.csv
2. invertebrates_dataset_1997.csv 
3. trawling_dataset.csv
4. trawling_dataset_1999.csv
5. WADeptEcology_CTD.csv
6. data_processing_function.R
7. modeling_analysis_and_figures.R
8. changepoint_analysis_and_figures.R

## Description of the data

invertebrates_dataset.csv - documents invertebrates caught in sampling bottom trawls in Puget Sound, WA, from 1999 to 2019. 
Includes trawl information such as year, date, time, and depth. 
Categorizes invertebrates by informal group (e.g., crab, sea star, or shrimp), common name (e.g., vermillion star, sand star, or spiny red star), genus species (e.g., Pycnopodia helianthoides) and the number caught during a particular trawl. 

invertebrates_dataset_1997.csv - documents invertebrates caught in sampling bottom trawls in Puget Sound, WA, in 1997. 
The row labelled “TRAWL” includes labels to indicate trawl time, where numbers 1, 2, 3, 4, and 5 correspond to afternoon, evening, night, early morning, and morning trawls. 
Depth is also listed. Sea stars collected during trawls are listed under the “echinoderms” heading by scientific name, and corresponding columns show how many individuals were caught during particular trawls.   

trawling_dataset.csv - includes trawling information from 2000 to 2019. 
The column labelled “Station Name” includes basic information such as the type of trawl employed and depth sampled, and “Date” includes the day, month and year of each trawl. 
The exact time and depth of the start and end of each trawl tow is also listed, as well as the Mean Lower Low Water depth. 
Coordinates of the start and end of trawl tows are listed in degrees decimal minutes form, as well as the distance, speed and direction of each trawl. 

trawling_dataset_1999.csv - reports trawling data from 1989 to 2000, including the same information presented in Dataset_trawlmaster.csv. 
From 1995 onwards, spatial coordinates are listed in degrees decimal minutes form.

WADeptEcology_CTD.csv – reports environmental data (including temperature data) collected by the Washington State Department of Ecology from 1999 to 2017 at a site in Puget Sound ~5.5 miles southeast of our trawling sites. 

data_processing_function.R - standardizes trawling data, environmental data, and sea star catch data into a single data frame to be called as a function in later analysis. Also used to calculate average sea star catch by species. 

modeling_analysis_and_figures.R - tests generalized linear mixed models of step changes and gradual trends in high-susceptibility and moderate-susceptibility sea star catch with different model specifications to identify which model is best supported for each susceptibility group based on AICc values. The best fit models are also used to assess the effect of temperature on high-susceptibility and moderate-susceptibility sea star catch over the sampling period. Code for the plots of model fit (Figs 3 and 4 in the manuscript) is included. Code for the plot of species abundance overtime (Fig 2) is also included. 

changepoint_analysis_and_figures.R - uses the “cpm” package in RStudio to look for step changes in high-susceptibility and moderate-susceptibility sea star catch at depth across the sampling period. Code for the plots of changepoints within susceptibility categories (Figs 5 and 6 in the manuscript) is included.
