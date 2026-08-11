"""
CaSR_to_NEMO_Forcing

This script convert CaSR forcing fields (integrated already) to the operational files required by SalishSeaCast. Please change the paths and dates inside the .py file.


Author: Junqi Qiu
"""

# ============================================================
# Import
# ============================================================

import os
import numpy as np
import pandas as pd
import xarray as xr


# ============================================================
# Paths and Dates
# ============================================================

# Path of integrated CaSR files
BASE_PATH = '/ocean/jqiu/forcings_raw/CaSR_raw/2008_2024_Integrated/'

# Path of output
OUTPUT_DIR = '/ocean/jqiu/forcings/CaSR/'

# dates of covert, both ends included
START_DATE = '2023-03-03'
END_DATE = '2023-03-31'

# Overwrite the existing files
OVERWRITE_EXISTING = True

# Examine the output 
VERIFY_OUTPUT = False


# ============================================================
# Operation variables vs CaSR file names and variable names
# ============================================================

VAR_MAPPING = {
    'solar': (
        'Integrated_RDPS_P_FB_SFC_2008_2024.nc',
        'CaSR_v3.2_P_FB_SFC'
    ),
    'precip': (
        'Integrated_RDPS_P_PR0_SFC_2008_2024.nc',
        'CaSR_v3.2_P_PR0_SFC'
    ),
    'therm_rad': (
        'Integrated_RDPS_P_FI_SFC_2008_2024.nc',
        'CaSR_v3.2_P_FI_SFC'
    ),
    'tair': (
        'Integrated_RDPS_P_TT_1.5m_2008_2024.nc',
        'CaSR_v3.2_P_TT_1.5m'
    ),
    'qair': (
        'Integrated_RDPS_P_HU_1.5m_2008_2024.nc',
        'CaSR_v3.2_P_HU_1.5m'
    ),
    'u_wind': (
        'Integrated_RDPS_P_UUC_10m_2008_2024.nc',
        'CaSR_v3.2_P_UUC_10m'
    ),
    'atmpres': (
        'Integrated_RDPS_P_PN_SFC_2008_2024.nc',
        'CaSR_v3.2_P_PN_SFC'
    ),
    'v_wind': (
        'Integrated_RDPS_P_VVC_10m_2008_2024.nc',
        'CaSR_v3.2_P_VVC_10m'
    )
}


# ============================================================
# 辅助函数
# ============================================================

def convert_units(data_array, target_var):
    """
    Convert the units of CaSR

    Variables
    ----
    data_array : xarray.DataArray
        input data

    target_var : str
        output variable name

    Return
    ----
    converted : xarray.DataArray
        Converted data

    units : str
        Converted units
    """

    if target_var == 'tair':
        # degree Celsius -> K
        converted = data_array + 273.15
        units = 'K'

    elif target_var in ('u_wind', 'v_wind'):
        # knots -> m/s
        converted = data_array * 0.514444
        units = 'm/s'

    elif target_var == 'atmpres':
        # hPa or mb -> Pa
        converted = data_array * 100.0
        units = 'Pa'

    elif target_var == 'precip':
        # m water equivalent, Quantity of precipitation (1h)   -> kg/m^2/s 
        # Note that the units for NEMO is kg/m^2/s, despite that sometime in HRDPS it denotes kg/m^2
        converted = data_array * 1000.0 / 3600
        units = 'kg/m^2/s'

    elif target_var in ('solar', 'therm_rad'):
        converted = data_array
        units = 'W/m^2'

    elif target_var == 'qair':
        converted = data_array
        units = 'kg/kg'

    else:
        converted = data_array
        units = ''

    return converted, units


def summarize_array(array):
    """
    返回数组的基本统计信息。

    对于全 NaN 数组，min 和 max 返回 NaN。
    """

    values = np.asarray(array)

    total_count = values.size
    finite_mask = np.isfinite(values)
    finite_count = int(finite_mask.sum())
    nan_count = int(np.isnan(values).sum())

    if finite_count > 0:
        minimum = float(np.nanmin(values))
        maximum = float(np.nanmax(values))
        mean = float(np.nanmean(values))
    else:
        minimum = np.nan
        maximum = np.nan
        mean = np.nan

    return {
        'total': total_count,
        'finite': finite_count,
        'nan': nan_count,
        'min': minimum,
        'max': maximum,
        'mean': mean
    }


def print_summary(variable_name, array, indent='    '):
    """
    打印变量统计信息。
    """

    stats = summarize_array(array)

    print(
        f'{indent}{variable_name}: '
        f'shape={np.shape(array)}, '
        f'finite={stats["finite"]}/{stats["total"]}, '
        f'NaN={stats["nan"]}, '
        f'min={stats["min"]:.8g}, '
        f'max={stats["max"]:.8g}, '
        f'mean={stats["mean"]:.8g}'
    )


def select_one_day(dataset, target_date):
    """
    从 Dataset 中选择指定日期。

    使用半开区间 [当天 00:00, 次日 00:00)，
    相比直接 sel(time='YYYY-MM-DD') 更明确。
    """

    day_start = pd.Timestamp(target_date)
    day_end = day_start + pd.Timedelta(days=1)

    return dataset.sel(
        time=slice(
            np.datetime64(day_start),
            np.datetime64(day_end - pd.Timedelta(nanoseconds=1))
        )
    )


def get_output_path(target_date):
    """
    根据日期构造输出文件路径。
    """

    date_str = pd.Timestamp(target_date).strftime('%Ym%md%d')
    filename = f'CaSR_y{date_str}.nc'

    return os.path.join(OUTPUT_DIR, filename)


def verify_written_file(
    output_path,
    expected_time,
    expected_shape
):
    """
    重新打开写出的 NetCDF 文件并进行基本验证。
    """

    print('  正在验证输出文件...')

    with xr.open_dataset(output_path) as ds_check:

        required_dimensions = {
            'time_counter',
            'y',
            'x'
        }

        missing_dimensions = (
            required_dimensions - set(ds_check.dims)
        )

        if missing_dimensions:
            raise ValueError(
                f'输出文件缺少维度: '
                f'{sorted(missing_dimensions)}'
            )

        written_time = ds_check['time_counter'].values

        if not np.array_equal(written_time, expected_time):
            raise ValueError(
                '输出文件中的 time_counter '
                '与参考时间不一致'
            )

        variables_to_check = list(VAR_MAPPING.keys()) + [
            'percentcloud'
        ]

        for variable_name in variables_to_check:

            if variable_name not in ds_check:
                raise KeyError(
                    f'输出文件中缺少变量 '
                    f'{variable_name}'
                )

            variable = ds_check[variable_name]

            if variable.shape != expected_shape:
                raise ValueError(
                    f'{variable_name} 输出形状为 '
                    f'{variable.shape}，'
                    f'预期为 {expected_shape}'
                )

            # 这里触发真正读取文件内容
            values = variable.values

            print_summary(
                variable_name,
                values,
                indent='    验证 '
            )

            if variable_name != 'percentcloud':
                finite_count = int(
                    np.isfinite(values).sum()
                )

                if finite_count == 0:
                    raise ValueError(
                        f'{variable_name} 写出后没有任何'
                        f'有限值，变量可能全部为 NaN'
                    )

    print('  输出文件验证通过。')


# ============================================================
# 单日转换
# ============================================================

def convert_one_day(target_date):
    """
    转换指定日期的数据，每天生成一个 NetCDF 文件。
    """

    print()
    print('=' * 78)
    print(f'开始处理日期: {target_date}')

    output_path = get_output_path(target_date)

    if os.path.exists(output_path):
        if OVERWRITE_EXISTING:
            print(f'输出文件已存在，将覆盖: {output_path}')
        else:
            print(f'输出文件已存在，跳过: {output_path}')
            return 'skipped'

    # --------------------------------------------------------
    # 使用温度文件获取参考时间、网格和经纬度
    # --------------------------------------------------------

    reference_filename = (
        'Integrated_RDPS_P_TT_1.5m_2008_2024.nc'
    )

    reference_path = os.path.join(
        BASE_PATH,
        reference_filename
    )

    if not os.path.exists(reference_path):
        raise FileNotFoundError(
            f'找不到参考文件: {reference_path}'
        )

    with xr.open_dataset(reference_path) as ds_ref_all:

        required_reference_variables = {
            'time',
            'rlat',
            'rlon',
            'lat',
            'lon'
        }

        missing_reference_variables = (
            required_reference_variables
            - set(ds_ref_all.variables)
        )

        if missing_reference_variables:
            raise KeyError(
                f'参考文件缺少变量或坐标: '
                f'{sorted(missing_reference_variables)}'
            )

        ds_ref = select_one_day(
            ds_ref_all,
            target_date
        ).load()

    time_size = ds_ref.sizes.get('time', 0)

    if time_size == 0:
        ds_ref.close()
        print(f'{target_date} 没有参考数据，跳过。')
        return 'no_data'

    y_size = ds_ref.sizes['rlat']
    x_size = ds_ref.sizes['rlon']

    reference_time = ds_ref['time'].values.copy()
    nav_lat = ds_ref['lat'].values.copy()
    nav_lon = ds_ref['lon'].values.copy()

    ds_ref.close()

    expected_shape = (
        time_size,
        y_size,
        x_size
    )

    print(f'时间步数量: {time_size}')
    print(f'网格大小: {y_size} × {x_size}')
    print(f'预期三维数组形状: {expected_shape}')
    print(f'首个时间: {reference_time[0]}')
    print(f'末个时间: {reference_time[-1]}')

    if nav_lat.shape != (y_size, x_size):
        raise ValueError(
            f'lat 形状为 {nav_lat.shape}，'
            f'预期为 {(y_size, x_size)}'
        )

    if nav_lon.shape != (y_size, x_size):
        raise ValueError(
            f'lon 形状为 {nav_lon.shape}，'
            f'预期为 {(y_size, x_size)}'
        )

    # --------------------------------------------------------
    # 创建输出 Dataset
    # --------------------------------------------------------

    ds_out = xr.Dataset(
        coords={
            'time_counter': (
                'time_counter',
                reference_time
            ),
            'y': (
                'y',
                np.arange(y_size, dtype=np.int32)
            ),
            'x': (
                'x',
                np.arange(x_size, dtype=np.int32)
            )
        }
    )

    ds_out['time_counter'].attrs = {
        'standard_name': 'time',
        'long_name': 'Time'
    }

    ds_out['y'].attrs = {
        'long_name': 'Grid y index'
    }

    ds_out['x'].attrs = {
        'long_name': 'Grid x index'
    }

    ds_out['nav_lat'] = (
        ('y', 'x'),
        nav_lat
    )

    ds_out['nav_lon'] = (
        ('y', 'x'),
        nav_lon
    )

    ds_out['nav_lat'].attrs = {
        'standard_name': 'latitude',
        'long_name': 'Latitude',
        'units': 'degrees_north'
    }

    ds_out['nav_lon'].attrs = {
        'standard_name': 'longitude',
        'long_name': 'Longitude',
        'units': 'degrees_east'
    }

    # --------------------------------------------------------
    # 逐个变量读取并转换
    # --------------------------------------------------------

    try:
        for target_var, (
            input_filename,
            source_var
        ) in VAR_MAPPING.items():

            print()
            print(f'  正在处理变量: {target_var}')

            input_path = os.path.join(
                BASE_PATH,
                input_filename
            )

            if not os.path.exists(input_path):
                raise FileNotFoundError(
                    f'找不到输入文件: {input_path}'
                )

            with xr.open_dataset(input_path) as ds_in:

                if source_var not in ds_in:
                    raise KeyError(
                        f'{input_filename} 中不存在变量 '
                        f'{source_var}'
                    )

                source_data = select_one_day(
                    ds_in[[source_var]],
                    target_date
                )[source_var].load()

            source_time_size = source_data.sizes.get(
                'time',
                0
            )

            if source_time_size == 0:
                raise ValueError(
                    f'{target_var} 在 {target_date} '
                    f'没有数据'
                )

            required_dims = {
                'time',
                'rlat',
                'rlon'
            }

            missing_dims = (
                required_dims
                - set(source_data.dims)
            )

            if missing_dims:
                raise ValueError(
                    f'{target_var} 缺少维度: '
                    f'{sorted(missing_dims)}；'
                    f'实际维度为 {source_data.dims}'
                )

            # 先将维度顺序固定下来
            source_data = source_data.transpose(
                'time',
                'rlat',
                'rlon'
            )

            if source_data.shape != expected_shape:
                raise ValueError(
                    f'{target_var} 在 {target_date} '
                    f'的形状为 {source_data.shape}，'
                    f'参考形状为 {expected_shape}'
                )

            source_time = source_data['time'].values

            if not np.array_equal(
                source_time,
                reference_time
            ):
                raise ValueError(
                    f'{target_var} 在 {target_date} '
                    f'的时间坐标与参考文件不一致'
                )

            print_summary(
                f'{target_var} 转换前',
                source_data.values
            )

            converted_data, output_units = convert_units(
                source_data,
                target_var
            )

            # ------------------------------------------------
            # 关键修复：
            #
            # 只取 NumPy 数组，然后按位置写入 ds_out。
            # 不把原始 rlat/rlon 坐标带入输出 Dataset，
            # 从而避免 xarray 按标签自动对齐后产生全 NaN。
            # ------------------------------------------------

            converted_values = np.asarray(
                converted_data.values
            )

            if converted_values.shape != expected_shape:
                raise ValueError(
                    f'{target_var} 转换后的形状为 '
                    f'{converted_values.shape}，'
                    f'预期为 {expected_shape}'
                )

            print_summary(
                f'{target_var} 转换后',
                converted_values
            )

            finite_count = int(
                np.isfinite(converted_values).sum()
            )

            if finite_count == 0:
                raise ValueError(
                    f'{target_var} 转换后没有任何有限值，'
                    f'请检查原始数据'
                )

            ds_out[target_var] = (
                ('time_counter', 'y', 'x'),
                converted_values
            )

            ds_out[target_var].attrs = {
                'units': output_units,
                'coordinates': 'nav_lat nav_lon'
            }

            # 写入 Dataset 后再次确认
            assigned_values = ds_out[target_var].values

            if assigned_values.shape != expected_shape:
                raise ValueError(
                    f'{target_var} 写入 ds_out 后形状为 '
                    f'{assigned_values.shape}，'
                    f'预期为 {expected_shape}'
                )

            assigned_finite_count = int(
                np.isfinite(assigned_values).sum()
            )

            if assigned_finite_count == 0:
                raise ValueError(
                    f'{target_var} 写入 ds_out 后全部为 '
                    f'NaN 或 Inf'
                )

            print_summary(
                f'{target_var} 写入后',
                assigned_values
            )

        # ----------------------------------------------------
        # 添加云覆盖率
        # ----------------------------------------------------

        print()
        print('  正在添加变量: percentcloud')

        percentcloud_values = np.zeros(
            expected_shape,
            dtype=np.float32
        )

        ds_out['percentcloud'] = (
            ('time_counter', 'y', 'x'),
            percentcloud_values
        )

        ds_out['percentcloud'].attrs = {
            'long_name': 'Cloud cover',
            'units': 'percent',
            'coordinates': 'nav_lat nav_lon'
        }

        # ----------------------------------------------------
        # Dataset 全局属性
        # ----------------------------------------------------

        ds_out.attrs = {
            'title': 'Converted CaSR atmospheric forcing',
            'source': (
                'Integrated RDPS/CaSR forcing files'
            ),
            'date_converted': str(target_date),
            'history': (
                'Converted with xarray; spatial dimensions '
                'stored as integer y/x indices.'
            )
        }

        # ----------------------------------------------------
        # 输出编码
        # ----------------------------------------------------

        output_encoding = {
            'time_counter': {
                'units': (
                    'seconds since 1970-01-01 00:00:00'
                ),
                'calendar': 'standard',
                'dtype': 'float64'
            },

            'y': {
                'dtype': 'int32'
            },

            'x': {
                'dtype': 'int32'
            },

            'nav_lat': {
                'zlib': True,
                'complevel': 4
            },

            'nav_lon': {
                'zlib': True,
                'complevel': 4
            }
        }

        for variable_name in VAR_MAPPING:
            output_encoding[variable_name] = {
                'zlib': True,
                'complevel': 4
            }

        output_encoding['percentcloud'] = {
            'zlib': True,
            'complevel': 4,
            'dtype': 'float32'
        }

        # ----------------------------------------------------
        # 写出前进行最终检查
        # ----------------------------------------------------

        print()
        print('  写出前变量统计:')

        for variable_name in (
            list(VAR_MAPPING.keys())
            + ['percentcloud']
        ):
            print_summary(
                variable_name,
                ds_out[variable_name].values
            )

        # ----------------------------------------------------
        # 写出文件
        # ----------------------------------------------------

        print()
        print(f'  正在写出: {output_path}')

        ds_out.to_netcdf(
            output_path,
            mode='w',
            unlimited_dims=['time_counter'],
            encoding=output_encoding
        )

    finally:
        ds_out.close()

    # --------------------------------------------------------
    # 写出后验证
    # --------------------------------------------------------

    if VERIFY_OUTPUT:
        verify_written_file(
            output_path=output_path,
            expected_time=reference_time,
            expected_shape=expected_shape
        )

    print(f'生成完毕: {output_path}')

    return 'success'


# ============================================================
# 日期范围转换
# ============================================================

def convert_date_range(start_date, end_date):
    """
    从 start_date 到 end_date 逐日转换。
    日期范围包含首日和末日。
    """

    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True
    )

    start_time = pd.Timestamp(start_date)
    end_time = pd.Timestamp(end_date)

    if end_time < start_time:
        raise ValueError(
            f'结束日期 {end_date} 不能早于'
            f'开始日期 {start_date}'
        )

    dates = pd.date_range(
        start=start_time,
        end=end_time,
        freq='D'
    )

    print('=' * 78)
    print(f'总共需要处理 {len(dates)} 天')
    print(f'日期范围: {start_date} 至 {end_date}')
    print(f'输入目录: {BASE_PATH}')
    print(f'输出目录: {OUTPUT_DIR}')
    print(f'覆盖已有文件: {OVERWRITE_EXISTING}')
    print(f'写出后验证: {VERIFY_OUTPUT}')

    results = {
        'success': [],
        'skipped': [],
        'no_data': [],
        'failed': []
    }

    for date in dates:

        target_date = date.strftime('%Y-%m-%d')

        try:
            status = convert_one_day(target_date)

            if status not in results:
                status = 'success'

            results[status].append(target_date)

        except Exception as error:
            print()
            print('!' * 78)
            print(f'处理失败: {target_date}')
            print(
                f'错误类型: '
                f'{type(error).__name__}'
            )
            print(f'错误信息: {error}')
            print('!' * 78)

            results['failed'].append(
                {
                    'date': target_date,
                    'error_type': type(error).__name__,
                    'error': str(error)
                }
            )

    # --------------------------------------------------------
    # 最终汇总
    # --------------------------------------------------------

    print()
    print('=' * 78)
    print('全部日期处理完成')
    print(f'成功: {len(results["success"])} 天')
    print(f'跳过: {len(results["skipped"])} 天')
    print(f'无数据: {len(results["no_data"])} 天')
    print(f'失败: {len(results["failed"])} 天')

    if results['success']:
        print()
        print('成功日期:')
        for date in results['success']:
            print(f'  {date}')

    if results['skipped']:
        print()
        print('跳过日期:')
        for date in results['skipped']:
            print(f'  {date}')

    if results['no_data']:
        print()
        print('无数据日期:')
        for date in results['no_data']:
            print(f'  {date}')

    if results['failed']:
        print()
        print('失败日期及原因:')

        for item in results['failed']:
            print(
                f'  {item["date"]}: '
                f'{item["error_type"]}: '
                f'{item["error"]}'
            )

    return results


# ============================================================
# 主程序
# ============================================================

if __name__ == '__main__':
    convert_date_range(
        start_date=START_DATE,
        end_date=END_DATE
    )