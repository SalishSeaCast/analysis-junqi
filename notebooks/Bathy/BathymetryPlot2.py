# This is a program that plots bathymetry data using a given dataset.
# The codes are adapted from the following website: https://github.com/SalishSeaCast/tools/blob/main/analysis_tools/Exploring%20netCDF%20Datasets%20Using%20xarray.ipynb
# This version is designed to test a high resolution map from Cartopy.
# Author: Junqi Qiu
# Date of creation: 2025-09-24

#-----------------------------------------------------
# Read the bathemetry data
import numpy as np
import xarray as xr

ds = xr.open_dataset('https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSnBathymetryV21-08')

# Alternatively, we can use a local dataset, which doesn't work for now. 
#lds = xr.open_dataset('../../NEMO-forcing/grid/bathy_meter_SalishSea2.nc')

# print(ds)

# Coordinates:
#  * gridY       (gridY) int32 4kB 0 1 2 3 4 5 6 ... 891 892 893 894 895 896 897
#  * gridX       (gridX) int16 796B 0 1 2 3 4 5 6 ... 391 392 393 394 395 396 397

# Data variables:
#    bathymetry  (gridY, gridX) float64 3MB ...
#    latitude    (gridY, gridX) float64 3MB ...
#    longitude   (gridY, gridX) float64 3MB ...

# print(ds.bathymetry.units, ds.bathymetry.long_name)
# metres Sea Floor Depth

#-----------------------------------------------------
# Plot the bathymetry data

import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from salishsea_tools import viz_tools

# Add coordinate and plot according to the coordinates
# lons = np.linspace(-127, -121, 400)
# lats = np.linspace(46, 52, 900)
lats = ds.variables['latitude']
lons = ds.variables['longitude']
bathy = ds.variables['bathymetry'][:]

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# Plot the bathymetry data using pcolormesh
# We must specify the transform, which is the coordinate system of the data
cmap = plt.get_cmap('winter_r')
cmap.set_bad('burlywood')
mesh = ax.pcolormesh(lons, lats, bathy, cmap=cmap, transform=ccrs.PlateCarree())


# Add map features for better context

ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.LAND, facecolor='burlywood')
ax.add_feature(cfeature.OCEAN, facecolor='skyblue')
ax.add_feature(cfeature.BORDERS, linestyle=':')

# Set the extent of the map to focus on your 6-degree area
# The order is [lon_min, lon_max, lat_min, lat_max]
ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()], crs=ccrs.PlateCarree())

# Add gridlines and labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False # Hide labels at the top
gl.ylabels_right = False # Hide labels on the right

# Add color bar
cbar = fig.colorbar(mesh, ax=ax, shrink=0.7)
cbar.set_label('Water Depth [meters]')

# Add axis labels
ax.set_xlabel('Longitude'.format(longitude=lons))
ax.set_ylabel('Latitude'.format(latitude=lats))
cbar.set_label('Water Depth [meters]')

# Save the plot to have a look, as X11 doesn't work
# plt.savefig('Bathy_test_2.png')
print("Plot saved to Bathy_test_2.png")
