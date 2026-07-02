# Weekly Meeting on July 2, 2026

## SalishSeaCast on `nibi`

Example run was a success.

In https://ubc-moad-docs.readthedocs.io/en/latest/alliance-computing.html, `def-allen` doesn't have the files we need. I changed it to `export PROJECT=$HOME/projects/rrg-allen`.

## `CaSR` Forcing Fields

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Data_CaSR_Conversion/CaSR2opr.ipynb

The values match well. `x`, `y` have different formats in 2 forcing fields and time counters are defined differently.


## Monthly Solar Radiation

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_model_comparison/Radiation%20Comparison.ipynb

I believe that the monthly `CaSR` radiation at `Sandheads` is systematically higher that that of `HRDPS`. The daily values can be either higher or lower.

Different from the expectation that `CaSR` is less extreme.

## To Do

Set `CaSR` tim counter unlimited and seconds since 1970. "unlimited" variable time, to tell NEMO that which dimension is time. 

The format ifference between float 32/64 should be fine.

Weight file. So that to match the grids of the forcing field and the model. 

Try the correct `CaSR` at the same date with the example `SalishSeaCast` yaml first, to avoid further changing.

