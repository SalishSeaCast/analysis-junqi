# Weekly Meeting on August 8

## Surface velocity blow-up

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/Test_results_comparison/mar23_CaSRvs_HRDPS.ipynb


### What went wrong?

In the `mar23` run, the `SST` decreased drastically (probably since `time step = 18` in the North) and to an extreme extent.

`Surface Velocity` deviated from standard (`HRDPS`) results at approximately the same time starting from the north. The velocity exploded at `time step = 32`, when the highest surface velocity over 8 m/s, in the North, with the whole flow field totally different from the standard results.

### What is not the problem?

I believe that all the units are correct. See the format in `mar23_CaSRvs_HRDPS.ipynb`. 

The wind speed and solar radiation have some errors but not significant enough to cause all the troubles.

### CaSR wind only + all other HRDPS forcings

https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_results_comparison/Test_results_comparison/01mar23_CaSR_wind_only.ipynb

Looks good. No significant error of `SST` or `surface velocity` within the 24 hours. 


## To Do


Zoom in the fjords to see the forcing fields.

Swap more variables and run influx analysis.

