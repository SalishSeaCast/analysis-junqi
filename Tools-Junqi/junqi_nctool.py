
import xarray as xr
import pandas as pd
from IPython.display import display, Markdown

def nc_disp(file_path: str, detailed: bool = False) -> None:
    """
    Display the contents of a NetCDF file natively in Jupyter Notebook.
    
    Args:
        file_path (str): Path to the NetCDF file.
        detailed (bool): If True, displays long names, units, and global attributes. 
                         Defaults to False.
    """
    try:
        with xr.open_dataset(file_path) as ds:
            # 渲染主标题
            display(Markdown(f"### NetCDF Summary"))
            display(Markdown(f"*Path: : `{file_path}`*"))
            # 内部辅助函数：用于输出干净的左对齐无索引 HTML 表格
            def display_clean_df(df):
                styled_df = df.style.set_properties(**{'text-align': 'left'}) \
                                    .set_table_styles([{'selector': 'th', 'props': [('text-align', 'left')]}]) \
                                    .hide(axis="index")
                display(styled_df)

            # Parse Dimensions
            display(Markdown("#### Dimensions"))
            if ds.sizes:
                df_dims = pd.DataFrame(list(ds.sizes.items()), columns=["Dimension", "Size"])
                display_clean_df(df_dims)
            else:
                display(Markdown("*No dimensions found.*"))

            # Parse Variables
            display(Markdown("#### Variables"))
            var_list = []
            for var_name, var in ds.variables.items():
                # 基础信息
                var_dict = {
                    "Name": var_name,
                    "Shape": str(var.shape),
                    "Dimensions": ", ".join(var.dims) if var.dims else "scalar"
                }
                
                # 详细信息
                if detailed:
                    var_dict["Long Name"] = var.attrs.get('long_name', '')
                    var_dict["Units"] = var.attrs.get('units', '')
                
                var_list.append(var_dict)
            
            if var_list:
                df_vars = pd.DataFrame(var_list)
                # 将缺失的属性替换为空字符串而不是 NaN，表格更美观
                df_vars = df_vars.fillna('') 
                display_clean_df(df_vars)
            else:
                display(Markdown("*No variables found.*"))
                
            # Parse Global Attributes
            if detailed:
                display(Markdown("#### Global Attributes"))
                if ds.attrs:
                    # Notebook 中 HTML 会自动换行，因此不需要像终端那样截断长字符串
                    attr_list = [{"Attribute": k, "Value": str(v)} for k, v in ds.attrs.items()]
                    df_attrs = pd.DataFrame(attr_list)
                    display_clean_df(df_attrs)
                else:
                    display(Markdown("*No global attributes found.*"))

    except FileNotFoundError:
        display(Markdown(f"**Error:** The file `{file_path}` was not found."))
    except Exception as e:
        display(Markdown(f"**Error reading the NetCDF file:** {e}"))



from datetime import datetime
import calendar
import os



def nc_path(variable=None, time=None, resolution=None, info=False):
    """
    生成存储服务器中 netCDF 文件的路径。
    
    参数:
    variable (str): 变量的直观名称 (e.g., 'nitrate', 'air_temp')
    time (str): 需要的时间，格式如 'YYYY-MM-DD' 或 'YYYY-MM' (视分辨率而定)
    resolution (str): 海洋数据的时间分辨率: 'monthly', 'daily', 'hourly', 静态数据填 'static'，大气数据不需要。
    info (bool): 如果为 True, 会调用 nc_disp 打印该路径下 nc 文件的详情。
    """
    
    #  变量字典配置 
    # (内部文件分类, 可选分辨率, Grid类型/来源, 原始变量名)
    OCEAN_VARS = {
        # 生物 (biol)
        'nitrate': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'nitrate'),
        'ammonium': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'ammonium'),
        'silicon': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'silicon'),
        'diatoms': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'diatoms'),
        'flagellates': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'flagellates'),
        'microzoo': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'microzooplankton'),
        'mesozoo': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'mesozooplankton'),
        'don': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'dissolved_organic_nitrogen'),
        'pon': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'particulate_organic_nitrogen'),
        'biogenic_silicon': ('biol', ['monthly', 'daily', 'hourly'], 'T', 'biogenic_silicon'),
        
        # 化学 (chem)
        'par': ('chem', ['monthly', 'daily', 'hourly'], 'T', 'PAR'),
        'turbidity': ('chem', ['monthly', 'daily', 'hourly'], 'T', 'turbidity'),
        'dic': ('chem', ['monthly', 'daily', 'hourly'], 'T', 'dissolved_inorganic_carbon'),
        'alkalinity': ('chem', ['monthly', 'daily', 'hourly'], 'T', 'total_alkalinity'),
        'oxygen': ('chem', ['monthly', 'daily', 'hourly'], 'T', 'dissolved_oxygen'),
        'co2_flux': ('chem', ['monthly', 'daily', 'hourly'], 'T', 'CO2_flux'),
        
        # 捕食 (graz)
        'mort_flagellates': ('graz', ['monthly', 'daily'], 'T', 'MORTPHY'),
        'mort_diatoms': ('graz', ['monthly', 'daily'], 'T', 'MORTDIAT'),
        'mort_microzoo': ('graz', ['monthly', 'daily'], 'T', 'MORTMICZ'),
        
        # 基础网格 (grid)
        'ssh': ('grid', ['monthly', 'daily', 'hourly'], 'T', 'sossheig'),
        'temperature': ('grid', ['monthly', 'daily', 'hourly'], 'T', 'votemper'),
        'salinity': ('grid', ['monthly', 'daily', 'hourly'], 'T', 'vosaline'),
        'density': ('grid', ['monthly', 'daily', 'hourly'], 'T', 'sigma_theta'),
        'cell_thickness': ('grid', ['monthly', 'daily', 'hourly'], 'T', 'e3t'),
        'latitude': ('grid', ['hourly'], 'T', 'nav_lat'),
        'longitude': ('grid', ['hourly'], 'T', 'nav_lon'),
        
        # 生产力 (prod)
        'prod_diatoms': ('prod', ['monthly', 'daily'], 'T', 'PPDIAT'),
        'prod_flagellates': ('prod', ['monthly', 'daily'], 'T', 'PPPHY'),
        
        # 动力学参数 (Grid U/V/W)
        'current_u': ('grid', ['hourly'], 'U', 'vozocrtx'),
        'current_v': ('grid', ['hourly'], 'V', 'vomecrty'),
        'vertical_vel': ('grid', ['hourly'], 'W', 'vovecrtz'),
        'eddy_diff': ('grid', ['hourly'], 'W', 'vert_eddy_diff'),
        'eddy_visc': ('grid', ['hourly'], 'W', 'vert_eddy_visc'),
        'dissipation': ('grid', ['hourly'], 'W', 'dissipation'),
        
        # 静态
        'bathymetry': ('bathy', ['static'], None, 'Bathymetry')
    }

    ATMO_VARS = {
        # HRDPS
        'air_pressure': ('HRDPS', 'atmpres'),
        'cloud_fraction': ('HRDPS', 'percentcloud'),
        'precip': ('HRDPS', 'precip'),
        'humidity': ('HRDPS', 'qair'),
        'solar_rad': ('HRDPS', 'solar'),
        'air_temp': ('HRDPS', 'tair'),
        'thermal_rad': ('HRDPS', 'therm_rad'),
        'wind_u': ('HRDPS', 'u_wind'),
        'wind_v': ('HRDPS', 'v_wind'),
        
        # CaSR
        'casr_wind_u': ('CaSR', 'CaSR_v3.2_P_UUC_10m'),
        'casr_wind_v': ('CaSR', 'CaSR_v3.2_P_VVC_10m')
    }

    # 指南模式 (无输入)
    if variable is None:
        print("====== 函数使用指南 ======")
        print("调用示例: nc_path('nitrate', time='2016-03-01', resolution='monthly')")
        print("以下是所有可供选择的变量名：\n")
        
        ocean_df = pd.DataFrame([
            {'Category': 'Ocean', 'Intuitive Variable': k, 'Original Variable': v[3], 'Available Resolutions': ', '.join(v[1])}
            for k, v in OCEAN_VARS.items()
        ])
        atmo_df = pd.DataFrame([
            {'Category': 'Atmosphere', 'Intuitive Variable': k, 'Original Variable': v[1], 'Available Resolutions': 'hourly/fixed (depends on source)'}
            for k, v in ATMO_VARS.items()
        ])
        guide_df = pd.concat([ocean_df, atmo_df], ignore_index=True)
        return guide_df

    # 判断变量所属分类
    is_ocean = variable in OCEAN_VARS
    is_atmo = variable in ATMO_VARS
    
    if not is_ocean and not is_atmo:
        raise ValueError(f"未找到变量 '{variable}'，请不带参数运行本函数查看可用变量列表。")

    # 提示模式 (仅输入变量名，缺失时间/分辨率)
    if is_ocean:
        file_mid_name, valid_res, grid_type, orig_var = OCEAN_VARS[variable]
        if variable == 'bathymetry':
            # 静态数据不需要时间
            path = '/data/nsoontie/MEOPAR/NEMO-forcing/grid/bathy_meter_SalishSea2.nc'
            
        else:
            if resolution is None or time is None:
                return f"提示: '{variable}' 是海洋数据。可用的时间分辨率为 {valid_res}。请重新调用并补充 time 和 resolution 参数，例如: time='2016-03-15', resolution='{valid_res[0]}'."
            
            if resolution not in valid_res:
                raise ValueError(f"'{variable}' 不支持 '{resolution}' 分辨率，支持的选项有: {valid_res}")
            
            try:
                dt = datetime.strptime(time, "%Y-%m-%d")
            except ValueError:
                try:
                    dt = datetime.strptime(time, "%Y-%m")
                except ValueError:
                    raise ValueError("时间格式错误，请使用 'YYYY-MM-DD' 或 'YYYY-MM'")
                    
            y_str, m_str, d_str = dt.strftime("%Y"), dt.strftime("%m"), dt.strftime("%d")
            
            if resolution == 'monthly':
                last_day = calendar.monthrange(dt.year, dt.month)[1]
                path = f"/results2/SalishSea/month-avg.202111/SalishSeaCast_1m_{file_mid_name}_{grid_type}_{y_str}{m_str}01_{y_str}{m_str}{last_day}.nc"
                
            else: # daily 或 hourly
                # 生成类似 '01apr07' 的文件夹名称
                folder_name = dt.strftime("%d%b%y").lower()
                res_prefix = '1d' if resolution == 'daily' else '1h'
                date_str = f"{y_str}{m_str}{d_str}"
                
                # 针对 grid_U, grid_V, grid_W 文件命名后缀不一样
                if grid_type in ['U', 'V', 'W']:
                    path = f"/results2/SalishSea/nowcast-green.202111/{folder_name}/SalishSea_{res_prefix}_{date_str}_{date_str}_grid_{grid_type}.nc"
                else:
                    path = f"/results2/SalishSea/nowcast-green.202111/{folder_name}/SalishSea_{res_prefix}_{date_str}_{date_str}_{file_mid_name}_{grid_type}.nc"

    elif is_atmo:
        source, orig_var = ATMO_VARS[variable]
        if time is None and source != 'CaSR':
             return f"提示: '{variable}' 是大气数据 ({source})。需要提供具体的时间 (time='YYYY-MM-DD')."
        
        if source == 'CaSR':
            if 'UUC' in orig_var:
                path = "/ocean/jqiu/Atmospheric_RDPS/2008_2024_Integrated/Integrated_RDPS_P_UUC_10m_2008_2024.nc"
            else:
                path = "/ocean/jqiu/Atmospheric_RDPS/2008_2024_Integrated/Integrated_RDPS_P_VVC_10m_2008_2024.nc"
        
        elif source == 'HRDPS':
            dt = datetime.strptime(time[:10], "%Y-%m-%d")
            y_str, m_str, d_str = dt.strftime("%Y"), dt.strftime("%m"), dt.strftime("%d")
            
            # 2023-02-15 之前的路径和前缀与之后的不同
            switch_date = datetime(2023, 2, 15)
            if dt < switch_date:
                path = f"/results/forcing/atmospheric/GEM2.5/operational/ops_y{y_str}m{m_str}d{d_str}.nc"
            else:
                path = f"/results/forcing/atmospheric/continental2.5/nemo_forcing/hrdps_y{y_str}m{m_str}d{d_str}.nc"

    # Info 调用与返回
    if info:
        try:
            print(f"正在读取文件信息: {path}")
            nc_disp(path, detailed=True)
        except Exception as e:
            print(f"调用 nc_disp 失败，文件可能不存在或路径有误。报错: {e}")
            
    return path



import warnings

# 变量字典
VAR_ALIASES = {
    # Ocean 物理量
    "temperature": ["temperature", "votemper", "temp", "t", "sst", "thetao", "sea_water_temperature", "water_temp"],
    "salinity": ["salinity", "vosaline", "salt", "s", "so", "sea_water_salinity"],
    "density": ["density", "sigma_theta", "rhop", "rho", "pot_rho"],
    "ssh": ["ssh", "sossheig", "zos", "sea_surface_height", "eta"],
    "cell_thickness": ["cell_thickness", "e3t", "thkcello", "dz"],
    "current_u": ["current_u", "vozocrtx", "u", "uo", "u_velocity", "water_u"],
    "current_v": ["current_v", "vomecrty", "v", "vo", "v_velocity", "water_v"],
    "vertical_vel": ["vertical_vel", "vovecrtz", "w", "wo", "vertical_velocity"],
    "eddy_diff": ["eddy_diff", "vert_eddy_diff", "difvho", "kz"],
    "eddy_visc": ["eddy_visc", "vert_eddy_visc", "kzv"],
    "dissipation": ["dissipation", "eps"],
    "bathymetry": ["bathymetry", "Bathymetry", "deptho", "topo", "elevation"],
    
    # Ocean 生化量 (BGC)
    "nitrate": ["nitrate", "no3", "NO3"],
    "ammonium": ["ammonium", "nh4", "NH4"],
    "silicon": ["silicon", "si", "Si", "silicate"],
    "diatoms": ["diatoms", "diat", "phy2"],
    "flagellates": ["flagellates", "nano", "phy1"],
    "microzoo": ["microzoo", "microzooplankton", "zoo1"],
    "mesozoo": ["mesozoo", "mesozooplankton", "zoo2"],
    "don": ["don", "dissolved_organic_nitrogen"],
    "pon": ["pon", "particulate_organic_nitrogen"],
    "biogenic_silicon": ["biogenic_silicon", "bsi", "bSi"],
    "par": ["par", "PAR", "irradiance"],
    "turbidity": ["turbidity"],
    "dic": ["dic", "dissolved_inorganic_carbon", "DIC"],
    "alkalinity": ["alkalinity", "total_alkalinity", "talk", "ALK"],
    "oxygen": ["oxygen", "dissolved_oxygen", "o2", "O2", "dox"],
    "co2_flux": ["co2_flux", "CO2_flux", "fgco2"],
    "prod_diatoms": ["prod_diatoms", "PPDIAT", "npp_diat"],
    "prod_flagellates": ["prod_flagellates", "PPPHY", "npp_flag"],
    "mort_flagellates": ["mort_flagellates", "MORTPHY"],
    "mort_diatoms": ["mort_diatoms", "MORTDIAT"],
    "mort_microzoo": ["mort_microzoo", "MORTMICZ"],

    # Atmosphere
    "air_pressure": ["air_pressure", "atmpres", "psl", "slp", "pressure"],
    "cloud_fraction": ["cloud_fraction", "percentcloud", "clt", "cloud"],
    "precip": ["precip", "pr", "precipitation", "rain"],
    "humidity": ["humidity", "qair", "hus", "q", "rh"],
    "solar_rad": ["solar_rad", "solar", "rsds", "swrad"],
    "air_temp": ["air_temp", "tair", "tas", "t2m", "air"],
    "thermal_rad": ["thermal_rad", "therm_rad", "rlds", "lwrad"],
    "wind_u": ["wind_u", "u_wind", "uas", "u10", "uwnd"],
    "wind_v": ["wind_v", "v_wind", "vas", "v10", "vwnd"],
    "casr_wind_u": ["casr_wind_u", "CaSR_v3.2_P_UUC_10m"],
    "casr_wind_v": ["casr_wind_v", "CaSR_v3.2_P_VVC_10m"],
    
    # Grid
    "nav_lat" : ["nav_lat", "Latitude", "latitude", "lat"],
    "nav_lon" : ["nav_lon", "Longitude", "longitude", "lon"],
    "e3t" : ["e3t", "thickness"],
    "depth" : ["depth", "deptht", "depthu", "depthv", "depthw"]
}

def nc_read(filepath, vars_to_read=None, info: bool = True):
    """
    读取 NetCDF 文件。支持读取单个或多个变量，自动处理别名映射、时间轴检查与经纬度质控。
    如果未指定 vars_to_read，则自动读取整个文件并打印所有变量名。
    
    参数:
    filepath: str, NetCDF 文件路径
    vars_to_read: str 或 list, 需要读取的变量名或变量名列表。默认为 None。
    info: bool, 是否显示提示信息
    
    返回:
    xarray.DataArray (单变量) 或 xarray.Dataset (多变量/全变量)
    """
    try:
        ds = xr.open_dataset(filepath)
    except Exception as e:
        print(f"无法打开文件 {filepath}: {e}")
        return None
        
    all_vars = list(ds.variables.keys())
    
    target_vars = []
    missing_vars = []
    is_single_var = False
    
    # ==========================================
    # 变量解析与批量读取
    # ==========================================
    if not vars_to_read:
        # 如果不给变量列表，自动把所有变量读入 ds
        target_vars = all_vars
        if info:
            var_str = ", ".join(all_vars)
            display(Markdown(f"**提示**: 未指定读取的变量，已自动读取全文件。包含以下变量: `{var_str}`"))
    else:
        # 判断输入是单个字符串还是列表，统一转换为列表处理
        is_single_var = isinstance(vars_to_read, str)
        var_list = [vars_to_read] if is_single_var else vars_to_read
        
        for var in var_list:
            target_var = None
            
            if var in all_vars:
                target_var = var
            else:
                if var in VAR_ALIASES:
                    for alias in VAR_ALIASES[var]:
                        if alias in all_vars:
                            target_var = alias
                            display(Markdown(f"**提示**: 识别到类别 `{var}`，已自动抽取真实变量 `{target_var}`。"))
                            break
                if not target_var:
                    category = None
                    for cat, aliases_list in VAR_ALIASES.items():
                        if var in aliases_list:
                            category = cat
                            break
                    if category:
                        for alias in VAR_ALIASES[category]:
                            if alias in all_vars:
                                target_var = alias
                                display(Markdown(f"**提示**: 你输入了 `{var}`，但文件使用的是 `{target_var}`，已为您转换。"))
                                break

            if target_var:
                target_vars.append(target_var)
            else:
                missing_vars.append(var)

        # 如果有找不到的变量，触发报错并调用 nc_disp
        if missing_vars:
            missing_str = ", ".join(missing_vars)
            display(Markdown(f"**错误**: 在文件中找不到以下变量或其别名: `{missing_str}`。正在显示文件概览："))
            try:
                nc_disp(filepath, detailed=True)
            except NameError:
                print("找不到 nc_disp 函数，请确保它已被正确定义并导入。")
            
            # 如果连一个能读的变量都没有，直接关闭并退出
            if not target_vars:
                ds.close()
                return None

    # ==========================================
    # 时间轴信息提醒逻辑 
    # ==========================================
    time_names = ['time', 'time_counter', 'ocean_time', 't', 'date', 'MT']
    time_var = next((n for n in time_names if n in target_vars), None)
    
    if time_var:
        t_data = ds[time_var]
        t_attrs = t_data.attrs
        units = t_attrs.get('units', '未知 (可能未定义单位)')
        calendar = t_attrs.get('calendar', 'standard (标准日历)')
        
        try:
            t_start = str(t_data.values[0])[:19] 
            t_end = str(t_data.values[-1])[:19]
            range_str = f"`{t_start}` 至 `{t_end}`"
        except:
            range_str = "无法直接解析起止时间"
        if info:
            display(Markdown(f"**时间轴 (`{time_var}`)**: `{units}` | `{calendar}` | {range_str}"))
        

    # ==========================================
    # 经纬度质控逻辑
    # ==========================================
    lon_names = ['lon', 'longitude', 'nav_lon', 'x', 'x_rho', 'lon_rho']
    lat_names = ['lat', 'latitude', 'nav_lat', 'y', 'y_rho', 'lat_rho']
    
    lon_var = next((n for n in lon_names if n in target_vars), None)
    lat_var = next((n for n in lat_names if n in target_vars), None)
    
    if lon_var and lat_var:
        lon_data = ds[lon_var]
        lat_data = ds[lat_var]
        
        lon_min, lon_max = lon_data.min().item(), lon_data.max().item()
        if lon_min < 0:
             if info:
                display(Markdown(f"**经度范围**: `-180 到 180` (Min: `{lon_min:.2f}`, Max: `{lon_max:.2f}`)"))
        else:
             if info:
                display(Markdown(f"**经度范围**: `0 到 360` (Min: `{lon_min:.2f}`, Max: `{lon_max:.2f}`)"))
             
        has_origin = ((lon_data == 0) & (lat_data == 0)).any().item()
        if has_origin and info:
            display(Markdown("**警告: 坐标网格中检测到了 (0, 0) 奇点！请留意掩膜设置。**"))
    
    # ==========================================
    # 深度警告
    # ==========================================
    depth_names=['depth','deptht','depthu','depthv','depthw']
    depth_var = next((n for n in depth_names if n in target_vars), None)
    if depth_var and info:
        display(Markdown("**警告: 读取深度。请留意掩膜所需深度与读取深度是否匹配。**"))

    # ==========================================
    # 返回结果
    # ==========================================
    if not vars_to_read:
        return ds                  # 未提供列表，直接返回完整的 Dataset
    elif is_single_var:
        return ds[target_vars[0]]  # 单变量返回 DataArray
    else:
        return ds[target_vars]     # 多变量返回 Dataset 的子集
    



# 帮助文档
HELP_TEXT = """
### `nc_extract` 使用指南

**功能**: 根据坐标或格点索引，快速从 xarray 数据中提取单点、切块或剖面。

**基础输入参数**:
* `data` (xr.Dataset/DataArray): 已经读取进内存的 xarray 数据。
* `method` (str): 提取模式。可选: `'point'` (单点), `'box'` (切块), `'profile'` (剖面)。
* `lat`, `lon` (str): 纬度和经度的维度名称，默认为 `'lat'` 和 `'lon'`。
* `target_coord` (tuple): 目标地理坐标 `(目标经度, 目标纬度)`。
* `target_grid` (tuple): 目标格点索引 `(经度格点号, 纬度格点号)`。*(与coord二选一)*
* `timecounter` (int/str): 时间步长索引 (int) 或具体时间 (str，如 '2023-01-01')。

**扩展参数 (**kwargs)**:
* `box_deg` (float): method='box' 且使用坐标时，向外扩展的度数，默认 1.0。
* `box_grid` (int): method='box' 且使用格点时，向外扩展的格点数，默认 5。
* `axis` (str): method='profile' 时，沿哪个轴切剖面。可选 `'lon'` 或 `'lat'` (默认)。

**示例调用**:
> nc_extract(ds, method='point', target_coord=(110, 20), timecounter=0)

# 整体读取 nc 文件
`ds = xr.open_dataset("your_weather_data.nc")`

# 提取 2023年8月1日，北京 (经度116, 纬度40) 的单点时间序列
result = nc_extract(
    data=ds, 
    method='point', 
    target_coord=(116, 40), 
    timecounter='2023-08-01'
)
"""

def nc_extract(data=None, method=None, lat='lat', lon='lon', 
               target_coord=None, target_grid=None, timecounter=None, **kwargs):
    
    # 无输入时打印使用说明
    if data is None or method is None:
        display(Markdown(HELP_TEXT))
        return

    # 检查输入的合法性
    if target_coord is None and target_grid is None:
        raise ValueError("报错：必须提供 `target_coord` (经纬度) 或 `target_grid` (格点号) 中的一个！")

    # 处理时间维度
    if timecounter is not None:
        time_dim = kwargs.get('time_dim', 'time')
        if time_dim in data.dims:
            # 如果是整数，按索引位置提取
            if isinstance(timecounter, int):
                data = data.isel({time_dim: timecounter})
            # 如果是时间字符串或 datetime 对象，提取最接近的时间
            else:
                data = data.sel({time_dim: timecounter}, method='nearest')
        else:
            warnings.warn(f"提醒：数据中未找到名为 '{time_dim}' 的时间维度，已跳过时间筛选。")

    # 核心空间提取逻辑
    if method == 'point':
        if target_coord:
            return data.sel({lon: target_coord[0], lat: target_coord[1]}, method='nearest')
        else:
            return data.isel({lon: target_grid[0], lat: target_grid[1]})

    elif method == 'box':
        if target_coord:
            box_deg = kwargs.get('box_deg', 1.0)
            lon_slice = slice(target_coord[0] - box_deg, target_coord[0] + box_deg)
            lat_slice = slice(target_coord[1] - box_deg, target_coord[1] + box_deg)
            return data.sel({lon: lon_slice, lat: lat_slice})
        else:
            box_grid = kwargs.get('box_grid', 5)
            # 使用 max(0, ...) 防止格点索引变成负数
            lon_slice = slice(max(0, target_grid[0] - box_grid), target_grid[0] + box_grid)
            lat_slice = slice(max(0, target_grid[1] - box_grid), target_grid[1] + box_grid)
            return data.isel({lon: lon_slice, lat: lat_slice})

    elif method == 'profile':
        axis = kwargs.get('axis', 'lat')
        if target_coord:
            if axis == 'lat':
                return data.sel({lat: target_coord[1]}, method='nearest')
            elif axis == 'lon':
                return data.sel({lon: target_coord[0]}, method='nearest')
        else:
            if axis == 'lat':
                return data.isel({lat: target_grid[1]})
            elif axis == 'lon':
                return data.isel({lon: target_grid[0]})
        
    else:
        raise ValueError("报错：`method` 参数错误！必须是 'point', 'box', 或 'profile'。")
    

import numpy as np
import xarray as xr
from IPython.display import display, Markdown

# 预设的站点字典：名称 -> (经度, 纬度)
# 你可以根据自己的需求随时往里面添加站点
STATIONS = {
    "Qingdao": (120.33, 36.07),
    "Xiamen": (118.08, 24.48),
    "Sanya": (109.51, 18.25),
    "Guangzhou": (113.23, 23.16),
    "Shanghai": (121.47, 31.23),
    "Argo_Test": (-170.5, 10.5) # 用于测试负经度
}

def nc_gridnumber(ds, lon=None, lat=None, station=None, info=True):
    """
    根据经纬度或站点名称，寻找 NetCDF 中最近的格点索引。
    
    参数:
    ds: xarray.Dataset 或 xarray.DataArray
    lon: float, 目标经度
    lat: float, 目标纬度
    station: str, 站点名称（优先使用）
    info: bool, 是否打印诊断信息
    
    返回:
    dict: 包含维度和对应索引的字典，例如 {'lon': 120, 'lat': 45} 或 {'y': 10, 'x': 20}。
          可以直接配合 ds.isel(**indices) 使用。
    """
    
    # ==========================================
    # 输入校验与站点解析
    # ==========================================
    if station is not None:
        if station in STATIONS:
            input_lon, input_lat = STATIONS[station]
            input_source = f"站点 `{station}`"
        else:
            display(Markdown(f"**错误**: 未找到站点 `{station}`。可用站点字典如下："))
            return STATIONS
    elif lon is not None and lat is not None:
        input_lon, input_lat = lon, lat
        input_source = f"输入坐标 `({lon}, {lat})`"
    else:
        display(Markdown("**提示**: 未提供有效的 `(lon, lat)` 或 `station`。目前内置的站点字典如下："))
        for k, v in STATIONS.items():
            print(f" - {k}: (Lon: {v[0]}, Lat: {v[1]})")
        return STATIONS

    # ==========================================
    # 识别坐标变量
    # ==========================================
    lon_names = ['lon', 'longitude', 'nav_lon', 'x', 'x_rho', 'lon_rho']
    lat_names = ['lat', 'latitude', 'nav_lat', 'y', 'y_rho', 'lat_rho']
    
    # 在坐标和变量中寻找
    avail_vars = list(ds.coords) + list(ds.variables)
    lon_var = next((n for n in lon_names if n in avail_vars), None)
    lat_var = next((n for n in lat_names if n in avail_vars), None)
    
    if not lon_var or not lat_var:
        display(Markdown("**错误**: 无法在数据集中识别出经纬度变量。"))
        return None
        
    ds_lon = ds[lon_var]
    ds_lat = ds[lat_var]

    # ==========================================
    # 经纬度质控与格式统一 (180 vs 360)
    # ==========================================
    ds_lon_min, ds_lon_max = ds_lon.min().item(), ds_lon.max().item()
    ds_lat_min, ds_lat_max = ds_lat.min().item(), ds_lat.max().item()
    
    is_0_360 = ds_lon_max > 180 or ds_lon_min >= 0
    lon_format_str = "`0 到 360`" if is_0_360 else "`-180 到 180`"
    
    target_lon = input_lon
    
    # 如果目标经度越界，自动帮你转换
    if is_0_360 and target_lon < 0:
        target_lon += 360
    elif not is_0_360 and target_lon > 180:
        target_lon -= 360

    # ==========================================
    # 寻找最近格点 (兼容 1D 直角网格 和 2D 曲线网格)
    # ==========================================
    # 使用巧妙的求余算法，避免日界线（180/-180 或 0/360）周边的距离计算错误
    dlon = (ds_lon - target_lon + 180) % 360 - 180
    dlat = ds_lat - input_lat

    indices = {}
    if ds_lon.ndim == 1 and ds_lat.ndim == 1:
        # 1D 坐标系 (如常规的经纬度网格)
        idx_lon = int(np.nanargmin(np.abs(dlon).values))
        idx_lat = int(np.nanargmin(np.abs(dlat).values))
        
        indices = {ds_lon.dims[0]: idx_lon, ds_lat.dims[0]: idx_lat}
        found_lon = ds_lon[idx_lon].item()
        found_lat = ds_lat[idx_lat].item()
    else:
        # 2D 坐标系 (如 ROMS, NEMO 的曲线网格)
        dist_sq = dlon**2 + dlat**2
        flat_idx = np.nanargmin(dist_sq.values)
        idx_tuple = np.unravel_index(flat_idx, dist_sq.shape)
        
        # 将维度名称与算出的多维索引对应
        indices = {dim: int(idx) for dim, idx in zip(ds_lon.dims, idx_tuple)}
        found_lon = ds_lon.values[idx_tuple]
        found_lat = ds_lat.values[idx_tuple]

    # 为了应对原数据中可能用负数表示西经的情况，对找到的经度统一做显示还原
    display_found_lon = found_lon
    if display_found_lon > 180 and input_lon < 0:
        display_found_lon -= 360

    # ==========================================
    # 计算实际球面距离误差 (Haversine公式)
    # ==========================================
    def haversine(lon1, lat1, lon2, lat2):
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        a = np.sin((lat2 - lat1)/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1)/2)**2
        return 6371 * 2 * np.arcsin(np.sqrt(a)) # 地球半径约 6371 km

    dist_km = haversine(input_lon, input_lat, found_lon, found_lat)

    # ==========================================
    # 打印报告
    # ==========================================
    if info:
        report = (
            f"**目标来源**: {input_source} -> `({input_lon}, {input_lat})`\n\n"
            f"**数据经度系统**: {lon_format_str} (Min: `{ds_lon_min:.2f}`, Max: `{ds_lon_max:.2f}`)\n\n"
            f"**匹配格点坐标**: `({display_found_lon:.4f}, {found_lat:.4f})`\n\n"
            f"**球面误差距离**: `{dist_km:.2f} km`\n\n"
            f"**对应切片索引**: `{indices}`"
        )
        display(Markdown(report))

    return indices