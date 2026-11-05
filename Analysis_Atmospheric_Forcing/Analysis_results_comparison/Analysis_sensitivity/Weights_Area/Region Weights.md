# Regional Volume Weights — README

## Description

This folder contains pre-calculated volume weights for seven regions in the Salish Sea, used to calculate regional volume-weighted mean values from model outputs.

Each region has a separate NetCDF file (`*_weights.nc`).

## File Structure

Each file contains:

- **Variable:** `weights`
- **Dimensions:** `(deptht: 40, y: 898, x: 398)`
- **Units:** dimensionless
- **Attribute:** `total_volume_m3` — total seawater volume of the region (m³)

Weights are calculated as:

\[
w(z,y,x) = \frac{e3t(z)\times area(y,x)\times tmask(z,y,x)}{V_{region}}
\]

where \(V_{region}\) is the sum of all valid grid-cell volumes within the region.

- Inside the region (ocean cells): normalized volume weights
- Outside the region and on land: `NaN`
- Sum of all valid weights in each file: **1.0**

Volumes are based on the reference vertical grid thickness (`e3t_1d`).

## Usage

The weights have the same spatial dimensions as the model output, so no regional slicing is required.

```python
import xarray as xr

# Load weights
weights = xr.open_dataset(
    "central_sog_weights.nc"
)["weights"]

# Model variable dimensions:
# (time_counter, deptht, y, x)

# Regional volume-weighted mean
regional_mean = (variable * weights).sum(
    dim=("deptht", "y", "x")
)
```

The result has dimensions `(time_counter,)`, representing the regional volume-weighted mean at each time step.

**Note:** This direct calculation assumes that the model variable has valid values at all ocean cells with non-NaN weights. If additional missing values exist, renormalize the weights over valid cells:

```python
regional_mean = (
    (variable * weights).sum(dim=("deptht", "y", "x"))
    /
    weights.where(variable.notnull()).sum(
        dim=("deptht", "y", "x")
    )
)
```

Ensure that the model output and weights use the same grid and vertical depth levels.