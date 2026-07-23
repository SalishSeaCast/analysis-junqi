# Weekly Meeting on July 23

## Wind-driven Coastal Upwelling

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_atmo_events/Wind_event_26jul.ipynb

Strong wind from Northwest was observed in Jun and early July, but no evidence indicates that the wind in 2026 is significantly stronger or more persistent than that in other years.

The low-oxygen was probably caused by wind-driven coastal upwelling, but wind alone cannot explain even lower oxygen compared to other years.


## Comparison between Outputs form CaSR and HRDPS

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_model_comparison/Analysis_results_comparison/Results_CaSRvsHRDPS.ipynb

RMSE increases as time goes by, from 0.022 to 0.126. The error of low current velocities increased first and error of all velocities  increased later. 

It seems that there is an systemic error that the surface current velocity driven by CaSR is higher than that driven by HRDPS. 

Spatially, the current velocity increased in CaSR results in fjords or inlets while that decreased in the open ocean at 12:30. However, CaSR surface velocity is higher in all regions than HRDPS at 23:30.

## More Weighted Files

Still waiting for them?

## To Do

1. Rotate the grids from SalishSeaCast resylts. The u and v are along the grids instead of east-north. Also take a look at the staggered-grids.

2. Take a look at the temperature and salinity.

3. Run CaSR simulate for a month to see if the dfference of surface current was caused by pressure offsets.

4. Take a look at biological response in the 1 month simulation.
