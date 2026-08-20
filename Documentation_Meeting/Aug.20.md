# Weekly Meeting on August 20, 2026

## SST

The `CaSR` results show lower SST in the fjords compared to default results, and the difference increases as time goes.

The profiles at the most extreme grid point indicates that it was probably not caused by mixing.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/Test_results_comparison/00mar23_CaSRvsHRDPS.ipynb

What about precipitation or mixing?

## Salinity Profile

The `CaSR` salinity somehow seemed to increase at the most extreme point? (forgot to chang the title)

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/Test_results_comparison/00mar23_CaSRvsHRDPS.ipynb

## Wind Stress Flux

The `CaSR` wind stress flux is lower than `HRDPS` results at fjords.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/Analysis_flux/01mar23_flux.ipynb

## Heat Flux (without solar radiation)

The `CaSR` downward heat flux is lower than that of `HRDPS`.

Maybe the SST decrease was caused by surface cooling? But the air temperature for `CaSR` is generally higher than `HRDPS`.

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_forcing_comparison/Analysis_spring/Spring_variables_inspector.ipynb

## Still working on HRDPS 1km

The conversion was incorrect. There were accumulative variables, which I thought was hourly.

## Output Variables

No `ptrc` outputs. I will change the configuration and output frequencies and run them again. (Go to biol for nitrate, not ptrc)


## To Do

Pull everything on Mar 2023 together and figure out the story. 



