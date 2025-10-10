# This is a python file
# This program aims to generate velocity fields
# Following the example from https://github.com/SalishSeaCast/tools/blob/main/analysis_tools/Plotting%20Velocity%20Fields%20on%20Horizontal%20Planes.ipynb
# Written by Junqi Qiu (2025-10-01)

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from IPython.display import display, Math, Latex

from salishsea_tools import (
    nc_tools,
    viz_tools,
)

# Please cd to the dictionary "/home/jqiu/Programing/VelocityFields" before running the code

# Load the data
u_vel = nc.Dataset('SalishSea_1h_20131231_20131231_grid_U.nc')
v_vel = nc.Dataset('SalishSea_1h_20131231_20131231_grid_V.nc')
w_vel = nc.Dataset('SalishSea_1h_20131231_20131231_grid_W.nc')

# Check the variables and dimensions of the nc files
# nc_tools.show_dimensions(u_vel)
# nc_tools.show_variables(u_vel)

# 3 velocities are called 'vozocrtx', 'vomecrty', 'depthu' (I don't know why they are named like this)

# Extract volocities from the nc files
ugrid = u_vel.variables['vozocrtx']
vgrid = v_vel.variables['vomecrty']
zlevels = u_vel.variables['depthu']
timesteps = u_vel.variables['time_counter']

# Extract the grid information
t = nc.Dataset('SalishSea_1h_20131231_20131231_grid_T.nc')
lats = t.variables['nav_lat']
lons = t.variables['nav_lon']

t, zlevel = 4, 0
y_slice = np.arange(250, 370)
x_slice = np.arange(200, 320)

ugrid_tzyx = np.ma.masked_values(ugrid[t, zlevel, y_slice, x_slice], 0)
vgrid_tzyx = np.ma.masked_values(vgrid[t, zlevel, y_slice, x_slice], 0)
u_tzyx, v_tzyx = viz_tools.unstagger(ugrid_tzyx, vgrid_tzyx)
speeds = np.sqrt(np.square(u_tzyx) + np.square(v_tzyx))
max_speed = viz_tools.calc_abs_max(speeds)

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
viz_tools.set_aspect(ax)

ax.streamplot(
    x_slice[1:], y_slice[1:], u_tzyx, v_tzyx,
    linewidth=7*speeds/max_speed,
)
viz_tools.plot_land_mask(ax, 'bathymetry_201702.nc', xslice=x_slice, yslice=y_slice)

ax.set_xlim(x_slice[0], x_slice[-1])
ax.set_ylim(y_slice[0], y_slice[-1])
ax.grid()

ax.set_xlabel('x Index')
ax.set_ylabel('y Index')
ax.set_title(u'Date = 2013-12-31, depth \u2248 {d:.2f}{z.units}'.format(d=zlevels[zlevel], z=zlevels))

# Save the plot to have a look
plt.savefig('VelocityField2.png')
print("Plot saved to VelocityField2.png")