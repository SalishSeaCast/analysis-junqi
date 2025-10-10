# This is a program that plots bathymetry data using a given dataset.
# The codes are adapted from the following website: https://github.com/SalishSeaCast/tools/blob/main/analysis_tools/Exploring%20netCDF%20Datasets%20Using%20xarray.ipynb
# Author: Junqi Qiu
# Date of creation: 2025-09-23

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

from salishsea_tools import viz_tools

bathy = ds.variables['bathymetry'][:]

# Create figures and axises
# fig, ax = plt.subplots(1, 1, figsize=(10, 8))

# Color mesh plots
# ax.pcolormesh(bathy)


# Adjust the size of the picture
# fig, ax = plt.subplots(1, 1, figsize=(10, 8))
# ax.set_aspect(5/4.4)
# ax.pcolormesh(bathy)

# Add a color bar
# fig, ax = plt.subplots(1, 1, figsize=(10, 8))
# viz_tools.set_aspect(ax)
# mesh = ax.pcolormesh(bathy)
# Change the style of the colormap which suits the bathymetry better
# mesh = ax.pcolormesh(bathy, cmap='winter_r')
# fig.colorbar(mesh)

# Add coordinate and plot according to the coordinates
lats = ds.variables['latitude']
lons = ds.variables['longitude']

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
viz_tools.set_aspect(ax, coords='map', lats=lats)
cmap = plt.get_cmap('winter_r')
cmap.set_bad('burlywood')
mesh = ax.pcolormesh(lons[:], lats[:], bathy[:], cmap=cmap)
cbar = fig.colorbar(mesh)
plt.axis((-123.5, -122.5, 48, 49))

# Add axis labels
ax.set_xlabel('Longitude'.format(longitude=lons))
ax.set_ylabel('Latitude'.format(latitude=lats))
cbar.set_label('Water Depth [meters]')

# Save the plot to have a look, as X11 doesn't work
plt.savefig('Bathy_test.png')
print("Plot saved to Bathy_test.png")