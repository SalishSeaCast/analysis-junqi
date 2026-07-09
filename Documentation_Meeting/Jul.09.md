# Weekly Meeting on July 9, 2026

## Weighted Files

The corrected example of `CaSR` forcing filed has been saved to `/ocean/jqiu/CaSR/CaSR_y2023m03d01.nc`. Might take a while to generate the weighted file.

## Variables Impacting Spring Blooms

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_model_comparison/Analysis_spring/Spring_variables_inspector.ipynb

According to `Hindcast of the timing of the spring phytoplankton bloom in the Strait of Georgia, 1968–2010`, atmospheric variables including solar radiation, cloud, wind (through mixing) dominate the timing of spring blooms while temperature affects the frequency of spring blooms.

https://www.sciencedirect.com/science/article/pii/S0079661113000773

`CaSR` and `HRDPS` actually have some differences in those variables. May affect the timing of spring blooms.


## What Simulations to Do?

Besides `HRDPS (2.5km)` and `CaSR`, we also have `HRDPS (1km)`. We could (but will not) recalibrate (or tune) the model values if necessary.

It's easy to swap the variables just by not changing the namlist. 

## To Do

Get a Coarsened example `HRDPS` file (by just picking points) and 1 km, to generate weighted files.
