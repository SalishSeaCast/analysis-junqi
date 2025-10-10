# This is a Python program that analyzes wind data by generating stream plots
# (I hope it works)
# Author: Junqi Qiu
# Date: 2025-10-09

# ---------------------------------------
# The data comes from /results/forcing/atmospheric/continental2.5/nemo_forcing/fcst/
# This dataset is a forecast of boundary conditions for NEMO model
# ---------------------------------------

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from IPython.display import display, Math, Latex
# We will need to transfer t(seconds since 1970-01-01 00:00:00) to a real date)
import datetime

# Try the old way to plot the map
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from salishsea_tools import (
    nc_tools,
    viz_tools,
)

indir='/home/jqiu/Programing/Wind/'
indir_results='/home/jqiu/Programing/Wind/Wind_Speed_Plots/'
fnames=['hrdps_y2025m10d01.nc','hrdps_y2025m10d02.nc','hrdps_y2025m10d03.nc','hrdps_y2025m10d04.nc',
'hrdps_y2025m10d05.nc', 'hrdps_y2025m10d06.nc','hrdps_y2025m10d07.nc', 'hrdps_y2025m10d08.nc',
'hrdps_y2025m10d09.nc','hrdps_y2025m10d10.nc']

# i=1
# ds=nc.Dataset(indir+fnames[i])
# print(ds)

# ---------------------------------------
# Data Infirmation:
#  dimensions(sizes): time_counter(24), y(230), x(190)
#  variables(dimensions): float64 time_counter(time_counter), int64 y(y), int64 x(x),
#       float64 nav_lon(y, x), float64 nav_lat(y, x), float32 LHTFL_surface(time_counter, y, x), 
#       float32 PRATE_surface(time_counter, y, x), float32 RH_2maboveground(time_counter, y, x), 
#       float32 atmpres(time_counter, y, x), float32 precip(time_counter, y, x), 
#       float32 qair(time_counter, y, x), float32 solar(time_counter, y, x), 
#       float32 tair(time_counter, y, x), float32 therm_rad(time_counter, y, x), 
#       float64 u_wind(time_counter, y, x), float64 v_wind(time_counter, y, x)

#for i in num(fnames):
#    ds=nc.Dataset(indir+fnames[i])

# --------------------------------------- 
# See the information of the variables
# print("\n--- u_wind Variable Information ---")
# print(u_wind_dataarray)

# ---------------------------------------
# u_wind Variable Information:
# long_name: U-Component of Wind at 10m
#    units: m s-1
#    standard_name: eastward_wind
#    ioos_category: atmospheric
#    coordinates: nav_lat nav_lon
# unlimited dimensions: time_counter
# current shape = (24, 230, 190)
# ---------------------------------------

# Plot the wind field at the first time step,

#t_index=0
# print(t[t_index])

# Create empty lists to store the time series, mean wind speeds, and max wind speeds
datetime_seris=[]
mean_speeds=[]
max_speeds=[]


for day_index in range(len(fnames)):
    ds=nc.Dataset(indir+fnames[day_index])

    # Access variables u_wind and v_wind
    # 此时 ds['u_wind'] 是一个 xarray.DataArray 对象
    t=ds.variables['time_counter']
    u_wind_dataarray = ds.variables['u_wind']
    v_wind_dataarray = ds.variables['v_wind']

    if day_index==9:
        t_index_max=13
    else:
        t_index_max=24

    for t_index in range(t_index_max):
        
        u_wind=u_wind_dataarray[t_index,:,:]
        v_wind=v_wind_dataarray[t_index,:,:]

        # Calculate the real time
        real_time=datetime.datetime(1970,1,1)+datetime.timedelta(seconds=int(t[t_index]))-datetime.timedelta(hours=8)

        # Calculate the wind speeds
        speeds = np.sqrt(np.square(u_wind) + np.square(v_wind))[:]
        max_speed = viz_tools.calc_abs_max(speeds)
        mean_speed = np.nanmean(speeds)

        # Append the results to the lists
        datetime_seris.append(real_time)
        mean_speeds.append(mean_speed)
        max_speeds.append(max_speed)
        
# Plot the time series of mean wind speeds and max wind speeds

fig=plt.figure(figsize=(12, 6))
plt.plot(datetime_seris, mean_speeds, label='Mean Wind Speed', color='blue', marker='o')
plt.plot(datetime_seris, max_speeds, label='Max Wind Speed', color='red', marker='x')
plt.xlabel('Date and Time')
plt.ylabel('Wind Speed (m/s)')
plt.title('Time Series of Mean and Max Wind Speeds (Oct 1 - Oct 10, 2025)')
plt.xticks(rotation=45)
plt.legend()
plt.grid()
plt.tight_layout()

# Save the plot to have a look
plt.savefig(indir_results+'WindSpeed.png')
print("Plot saved to WindSpeed.png")
plt.close('all')
