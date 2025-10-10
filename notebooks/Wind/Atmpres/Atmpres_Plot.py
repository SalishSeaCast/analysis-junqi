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
indir_results='/home/jqiu/Programing/Wind/Atmpres/'
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
# Let's begin with the wind data on Oct. 9 in Vancouver, meaning that it is Oct. 10 in UTC
# i=9 # Remenber Python index starts from 0

# Oct. 8 in Vancouver
i=8 # Remenber Python index starts from 0
ds=nc.Dataset(indir+fnames[i])

# Access variables u_wind and v_wind
# 此时 ds['u_wind'] 是一个 xarray.DataArray 对象
t=ds.variables['time_counter']
atmpres_dataarray=ds.variables['atmpres']
u_wind_dataarray = ds.variables['u_wind']
v_wind_dataarray = ds.variables['v_wind']
lats=ds.variables['nav_lat']
lons=ds.variables['nav_lon']

lats1d=lats[:]
lons1d=lons[:]
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

for t_index in range(24):
    real_time=datetime.datetime(1970,1,1)+datetime.timedelta(seconds=int(t[t_index]))-datetime.timedelta(hours=8)
# print(real_time)

    u_wind=u_wind_dataarray[t_index,:,:]
    v_wind=v_wind_dataarray[t_index,:,:]

    # ---------------------------------------
    # Start Plotting
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    ax.set_extent([lons1d.min(), lons1d.max(), lats1d.min(), lats1d.max()], crs=ccrs.PlateCarree())
    min_pressure = 100000  # 设定气压范围
    max_pressure = 103000

    # Plot Atmospheric Ptressure first as the background
    atmpres=atmpres_dataarray[t_index,:,:]
    atmpres_cmap = plt.get_cmap('viridis')
    atmpres_cmap.set_bad('lightgray')
    
    

    # 1. 填充背景色 (pcolormesh)
    atmpres_mesh = ax.pcolormesh(
        lons1d, lats1d, atmpres[:], 
        cmap=atmpres_cmap, 
        alpha=0.6,
        vmin=min_pressure,  # 确保颜色映射范围固定
        vmax=max_pressure
    )

    cbar = fig.colorbar(atmpres_mesh, ax=ax, orientation='vertical', pad=0.02, aspect=16, shrink=0.8)
    cbar.set_label(atmpres_dataarray.units)
    
    # 2. 叠加等值线 (contour) - 推荐
    atmpres_contour = ax.contour(
        lons1d, lats1d, atmpres[:], 
        levels=np.arange(960, 1040, 4), # 等值线间隔通常取 4hPa 或 5hPa
        colors='k', 
        linewidths=0.7
    )
    ax.clabel(atmpres_contour, fmt='%d', fontsize=8, colors='k') 


    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    # ax.add_feature(cfeature.LAND, color='lightgray')

    # Add grids and labels
    gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False # 隐藏顶部标签
    gl.right_labels = False # 隐藏右侧标签

    speeds = np.sqrt(np.square(u_wind) + np.square(v_wind))[:]
    max_speed = viz_tools.calc_abs_max(speeds)

    # 目标投影是 PlateCarree
    target_proj = ax.projection # 获取当前ax使用的投影

    # 创建源坐标系
    source_crs = ccrs.PlateCarree()

    # 将经纬度(lons, lats)和零高度 转换到目标投影坐标
    # 注：必须先展平lons和lats，进行转换后再重新塑形
    transformed_points = target_proj.transform_points(source_crs, lons1d, lats1d)

    # 3. 提取转换后的 X 和 Y 坐标（transformed_points[:,:,0] 是 X，transformed_points[:,:,1] 是 Y）
    X_proj = transformed_points[..., 0]
    Y_proj = transformed_points[..., 1]

    # Try streamplot this time
    ax.streamplot(
        X_proj, Y_proj, u_wind, v_wind,
        linewidth=7*speeds/max_speed,
    )

    ax.set_title('Date: '+ str(real_time) +' Wind Field')

    # Save the plot to have a look
    plt.savefig(indir_results+str(real_time)+'AtmPresField.png')
    print("Plot saved to AtmPresField.png")
    plt.close('all')
