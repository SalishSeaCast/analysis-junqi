# Weekly Meeting on Sep. 17, 2026

## HRDPS 10 (Subsampled) Results

1. Same SST problem in the fjords, but relatively better, not as serious as original CaSR outputs.

2. A significant surface velocity error (an error larger than CaSR) happens on Dec.31 in mid SoG, but the currents are similar between CaSR and HRDPS 10 outputs.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/00MMM2023/HRDPS_10km_2023_Results.ipynb

## HRDPS 1km Blow-up, again

Somehow the HRDPS 1km output has an extremely low SST in Mar, 2023. No explanations yet. Flux looks weird but most forcing variables look quite all right. Might need another look.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/00Mar2023_flux/00Mar2023_HRDPS1.ipynb 

## Corrected CaSR Results

Looks good. Still slightly different in the fjords but much better than uncorrected CaSR forcing field.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_weights/weights_CaSR_temp.ipynb

## Manman's Email?



## To Do

Take another look at HRDPS 1 Weights file and the forcing field through it.

SSS, halocline depth and strength (how deep and steep), RMSE velocity averaged over the surface (or in boxes), near surface nitrate, diatoms and flagellates. 

Long time scale. 1d looks good so far and we could have `1h` outputs anytime we like. 

2 ways of Temperature correlation: 1, find a sea point , or 2, height correction. Use height to correct the air temperature and see the bias compared to the observational data.

