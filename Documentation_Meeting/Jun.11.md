# Weekly Meeting on April 11

## Solar radiation (downward short-wave radiation)

See https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_model_comparison/Radiation%20Comparison.ipynb

The radiation of `CaSR` is lightly lower than `HRDPS`. 

## Wind speed and wind stress curl

See https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Analysis_model_comparison/Wind%20Comparison.ipynb

The wind speeds between `HRDPS` and 2 derived low-resolution models match well, and there are systemic errors between `HRDPS` and `CaSR`.

However, the wind stress curl tells a different story. `HRDPS` and `CaSR` have similar wind stress curls (of the same trend), but w derived models exhibit systemic errors. Subsampling and meaning are introducing errors through the spatial non-linearity of wind speed.

## Variables of `HRDPS` and `CaSR`

See https://github.com/SalishSeaCast/analysis-junqi/blob/main/Analysis_Atmospheric_Forcing/Data_CaSR_Conversion/CaSR2opr.ipynb


1. `PercentCloud` is missing in `CaSR`.
2. How was `x` and `y` defined?