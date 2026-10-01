# Weekly Meeting on Oct.01, 2026

## Correction (in the morning, nothing new)

Link to the paper: https://journals.ametsoc.org/view/journals/phoc/37/4/jpo2984.1.xml

1km HRDPS: https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/00Mar2023_flux/00Mar2023_HRDPS1.ipynb

Correction Compare: https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_forcing_vs_observation/Temp_Correction_vs_Model.ipynb

## HRDPS 1km Weights Files

Figured out the problem. 1km HRDPS covers the Salish Sea, but not all the bathy area.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_weights/HRDPS_1km_Weights/weights_HRDPS_1km_what's_wrong.ipynb

Regenerated.

Also used HRDPS 2.5 to interpolate the missing dates.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Data_Conversion/Data_HRDPS_1km_conversion/interpolate_missing_days.ipynb

## Physical and Biological Compare

Link: https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/00MMM2023/All_compare.ipynb

Probably not so reliable. The processing has been quite painful.

## To Do

### Wind correction

first compare the different correction plans. Maybe the correction at Pam Rocks will make wind worse (even weaker). 

Take a look at winds in the different models, to see if they are really so different.

### All compare debug

Correct the hailcline calculation, maybe get rid of the bottom. 

Why HRDPS 2.5 has low salinity at specific regions?

### CaSR Biological Response

Figure out why CaSR has low Nitrate and high diatom/flagellates. (They are not dino).

### Rerun the simulations 

Run everything with 6.5K corrections for all models and find out the difference. See what the next largest problem is.

