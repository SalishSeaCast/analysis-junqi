# Weekly Meeting on Sep. 03, 2026

## Heat Flux in Fjords

In CaSR results, the ocean is losing more heat into the atmosphere compared to HRDPS result in fjords, causing the SST difference.

Wind Stress flux looks good.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/00Mar2023_flux/00Mar2023_SST.ipynb

## Nitrate and Diatom in CaSR and HRDPS Driven Results

No significant difference in surface and top 10 m nitrate or diatoms caused by SST difference. The main difference happens in the open ocean.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/00Mar2023_flux/00Mar2023_Biol.ipynb

## Weights Files

CaSR is using grid points involving land area in the fjords. 

Maybe use a remote point in the open ocean to balance it?

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_forcing_comparison/Analysis_weights/Fjords.ipynb


## To Do

Take a look at the grid points not used and see if it's more correlated with the observation data. 

Take a look at June, September and December. Change the `time namelist` and `yaml` file. `time 000` is labeled in `restart` file (plus 1 when using) and `time end` can be calculated by 2160 * Number of days. There are 2160 time steps in a day, 40 s step size. 
