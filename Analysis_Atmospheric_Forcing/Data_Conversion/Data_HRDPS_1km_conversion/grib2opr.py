from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


# ============================================================
# 配置
# ============================================================

GRIB_ROOT = Path("/results/forcing/atmospheric/GEM1.0/GRIB")

START_DATE = pd.Timestamp("2023-03-16")
END_DATE = pd.Timestamp("2023-03-31")

RUN_HOUR = "00"

LEADS = range(1, 25)

# "." 表示把每天的 NetCDF 写到当前运行目录
OUTPUT_DIR = Path("/ocean/jqiu/forcings/HRDPS_1km")


# ============================================================
# HRDPS rotated grid 定义
# ============================================================

# GRIB metadata:
#
# gridType = rotated_ll
# uvRelativeToGrid = 1
#
# latitudeOfSouthernPoleInDegrees  = -33.443381
# longitudeOfSouthernPoleInDegrees = 266.463574
# angleOfRotationInDegrees         = 0.0
#
# 等价 northern pole:

ROTATED_NORTH_POLE_LAT = 33.443381
ROTATED_NORTH_POLE_LON = 86.463574


# ============================================================
# 累计量负值容忍阈值
# ============================================================

# 差分后转成最终物理单位，再判断负值。
#
# solar / therm_rad:
#     单位 W m-2
#
# precip:
#     单位 kg m-2 s-1
#
# DSWRF 的诊断已经发现约 -1.28 W m-2 的累计场数值噪声，
# 因此这里给辐射留 2 W m-2 的余量。
#
# precipitation 的阈值单独设得很小；第一次运行时建议观察
# 程序打印出来的 minimum rate 和负值个数。

NEGATIVE_RATE_TOLERANCE = {
    "solar": 16.0,
    "therm_rad": 16.0,
    "precip": 1.0e-8,
}


# ============================================================
# 输入变量定义
# ============================================================

VARIABLES = {
    "solar": {
        "grib_name": "DSWRF",
        "level": "SFC_0",
        "type": "accum",
    },

    "precip": {
        "grib_name": "APCP",
        "level": "SFC_0",
        "type": "accum",
    },

    "therm_rad": {
        "grib_name": "DLWRF",
        "level": "SFC_0",
        "type": "accum",
    },

    "tair": {
        "grib_name": "TMP",
        "level": "TGL_2",
        "type": "instant",
    },

    "qair": {
        "grib_name": "SPFH",
        "level": "TGL_2",
        "type": "instant",
    },

    "u_wind": {
        "grib_name": "UGRD",
        "level": "TGL_10",
        "type": "instant",
    },

    "v_wind": {
        "grib_name": "VGRD",
        "level": "TGL_10",
        "type": "instant",
    },

    "atmpres": {
        "grib_name": "PRMSL",
        "level": "MSL_0",
        "type": "instant",
    },
}


# ============================================================
# NetCDF 变量 metadata
# ============================================================

VARIABLE_ATTRS = {
    "solar": {
        "long_name":
            "Surface downwelling shortwave radiation flux",

        "standard_name":
            "surface_downwelling_shortwave_flux_in_air",

        "units":
            "W m-2",

        "comment":
            "Hourly mean downward shortwave radiation flux "
            "derived by differencing successive HRDPS "
            "accumulated shortwave radiation fields. Small "
            "negative differenced values attributable to "
            "numerical precision in the accumulated GRIB "
            "fields are clipped to zero. The stored forcing "
            "timestamp marks the beginning of the represented "
            "one-hour interval.",
    },

    "therm_rad": {
        "long_name":
            "Surface downwelling longwave radiation flux",

        "standard_name":
            "surface_downwelling_longwave_flux_in_air",

        "units":
            "W m-2",

        "comment":
            "Hourly mean downward longwave radiation flux "
            "derived by differencing successive HRDPS "
            "accumulated longwave radiation fields. Small "
            "negative differenced values attributable to "
            "numerical precision are clipped to zero. "
            "The stored forcing timestamp marks the beginning "
            "of the represented one-hour interval.",
    },

    "precip": {
        "long_name":
            "Surface precipitation mass flux",

        "standard_name":
            "precipitation_flux",

        "units":
            "kg m-2 s-1",

        "comment":
            "Hourly mean precipitation rate derived by "
            "differencing successive HRDPS accumulated "
            "precipitation fields. Small negative differenced "
            "values attributable to numerical precision are "
            "clipped to zero. The stored forcing timestamp "
            "marks the beginning of the represented one-hour "
            "interval.",
    },

    "tair": {
        "long_name":
            "Air temperature at 2 m above ground",

        "standard_name":
            "air_temperature",

        "units":
            "K",

        "comment":
            "HRDPS 2 m air temperature from forecast leads "
            "P001-P024 of the 00 UTC forecast cycle.",
    },

    "qair": {
        "long_name":
            "Specific humidity at 2 m above ground",

        "standard_name":
            "specific_humidity",

        "units":
            "kg kg-1",

        "comment":
            "HRDPS 2 m specific humidity from forecast leads "
            "P001-P024 of the 00 UTC forecast cycle.",
    },

    "u_wind": {
        "long_name":
            "Eastward wind at 10 m above ground",

        "standard_name":
            "eastward_wind",

        "units":
            "m s-1",

        "comment":
            "HRDPS native grid-relative 10 m wind rotated to "
            "the true geographic eastward component.",
    },

    "v_wind": {
        "long_name":
            "Northward wind at 10 m above ground",

        "standard_name":
            "northward_wind",

        "units":
            "m s-1",

        "comment":
            "HRDPS native grid-relative 10 m wind rotated to "
            "the true geographic northward component.",
    },

    "atmpres": {
        "long_name":
            "Air pressure at mean sea level",

        "standard_name":
            "air_pressure_at_mean_sea_level",

        "units":
            "Pa",

        "comment":
            "HRDPS mean sea level pressure from forecast leads "
            "P001-P024 of the 00 UTC forecast cycle.",
    },

    "percentcloud": {
        "long_name":
            "Cloud cover placeholder",

        "units":
            "percent",

        "comment":
            "Placeholder forcing field set uniformly to zero. "
            "It is not derived from HRDPS cloud forecast data.",
    },
}


# ============================================================
# 文件定位
# ============================================================

def get_file(target_date, grib_name, level_name, lead):
    """
    返回 target_date 00Z cycle 中指定 lead 的 GRIB2 文件。
    """

    folder = (
        GRIB_ROOT
        / target_date.strftime("%Y%m%d")
        / RUN_HOUR
        / f"{lead:03d}"
    )

    filename = (
        f"CMC_hrdps_west_{grib_name}_{level_name}_"
        f"rotated_latlon0.009x0.009_"
        f"{target_date:%Y%m%d}T{RUN_HOUR}Z_"
        f"P{lead:03d}-00.grib2"
    )

    path = folder / filename

    if not path.exists():
        raise FileNotFoundError(
            f"找不到 HRDPS 文件:\n{path}"
        )

    return path


# ============================================================
# GRIB 读取
# ============================================================

def read_grib(path):
    """
    读取一个 GRIB2 文件。

    返回
    ----
    data  : float32
    lat   : float64
    lon   : float64, 0-360
    attrs : GRIB variable attributes
    """

    ds = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
        },
    )

    variable_name = list(ds.data_vars)[0]

    da = ds[variable_name].squeeze(
        drop=True
    ).load()

    data = da.values.astype(
        np.float32
    )

    lat = ds["latitude"].values.astype(
        np.float64
    )

    lon = np.mod(
        ds["longitude"].values.astype(
            np.float64
        ),
        360.0,
    )

    attrs = dict(da.attrs)

    ds.close()

    return data, lat, lon, attrs


# ============================================================
# 累计量 -> 一小时平均 rate
# ============================================================

def accumulated_to_hourly_rate(
    raw,
    variable_name,
):
    """
    把从 initialization 开始累计的 HRDPS 变量转换成
    每小时平均 rate。

    输入:
        raw[0]  = P001 = 0-1 h accumulation
        raw[1]  = P002 = 0-2 h accumulation
        ...
        raw[23] = P024 = 0-24 h accumulation

    输出:
        forcing 00 = P001 / 3600
        forcing 01 = (P002 - P001) / 3600
        ...
        forcing 23 = (P024 - P023) / 3600

    forcing timestamp 使用 interval-start convention。
    """

    hourly_accum = np.empty_like(
        raw,
        dtype=np.float32,
    )

    # 第一个小时：从 initialization 到 P001
    hourly_accum[0] = raw[0]

    # 后续小时：相邻累计场差分
    hourly_accum[1:] = np.diff(
        raw,
        axis=0,
    )

    # 转成每秒 rate：
    #
    # radiation:
    #     J m-2 / s = W m-2
    #
    # precipitation:
    #     kg m-2 / s = kg m-2 s-1

    hourly_rate = (
        hourly_accum / 3600.0
    ).astype(np.float32)

    tolerance = NEGATIVE_RATE_TOLERANCE[
        variable_name
    ]

    min_rate = float(
        np.nanmin(hourly_rate)
    )

    negative_mask = (
        hourly_rate < 0.0
    )

    n_negative = int(
        np.count_nonzero(
            negative_mask
        )
    )

    print(
        f"  {variable_name}: "
        f"minimum differenced rate = "
        f"{min_rate:.8g}"
    )

    print(
        f"  {variable_name}: "
        f"negative cells before clipping = "
        f"{n_negative}"
    )

    # --------------------------------------------------------
    # 异常负值检查
    # --------------------------------------------------------

    if min_rate < -tolerance:
        raise ValueError(
            f"{variable_name}: "
            "累计量差分出现异常负值。 "
            f"minimum hourly rate = "
            f"{min_rate:.8g}; "
            f"allowed tolerance = "
            f"{tolerance:.8g}"
        )

    # --------------------------------------------------------
    # 小负值视为累计 GRIB 数值精度噪声，裁成 0
    # --------------------------------------------------------

    if n_negative > 0:

        print(
            f"  {variable_name}: "
            f"{n_negative} 个小负值被裁为 0"
        )

        hourly_rate[
            negative_mask
        ] = 0.0

    return hourly_rate


# ============================================================
# HRDPS 风旋转系数
# ============================================================

def hrdps_rotation_coefficients(
    lat,
    lon,
):
    """
    根据固定 HRDPS West rotated_ll 网格定义计算局地
    grid-relative -> true east/north 的旋转系数。

    theta 定义为:
        true east -> native grid +x
        朝 true north 为正。

    因此:

        u_east =
            u_grid*cos(theta)
            - v_grid*sin(theta)

        v_north =
            u_grid*sin(theta)
            + v_grid*cos(theta)
    """

    pole_lat = np.deg2rad(
        ROTATED_NORTH_POLE_LAT
    )

    pole_lon = np.deg2rad(
        ROTATED_NORTH_POLE_LON
    )

    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)

    dlon = (
        lon_rad
        - pole_lon
    )

    sin_theta = (
        np.cos(pole_lat)
        * np.sin(dlon)
    )

    cos_theta = (
        np.sin(pole_lat)
        * np.cos(lat_rad)
        - np.cos(pole_lat)
        * np.sin(lat_rad)
        * np.cos(dlon)
    )

    # 归一化，避免浮点误差
    norm = np.hypot(
        cos_theta,
        sin_theta,
    )

    cos_theta = (
        cos_theta / norm
    ).astype(np.float32)

    sin_theta = (
        sin_theta / norm
    ).astype(np.float32)

    return (
        cos_theta,
        sin_theta,
    )


def rotate_hrdps_wind(
    u_grid,
    v_grid,
    cos_theta,
    sin_theta,
):
    """
    HRDPS native grid-relative U/V
    -> true eastward / northward U/V.
    """

    u_east = (
        u_grid
        * cos_theta[None, :, :]
        - v_grid
        * sin_theta[None, :, :]
    )

    v_north = (
        u_grid
        * sin_theta[None, :, :]
        + v_grid
        * cos_theta[None, :, :]
    )

    return (
        u_east.astype(np.float32),
        v_north.astype(np.float32),
    )


def process_one_day(target_date):
    # ============================================================
    # forcing 时间轴
    # ============================================================

    # 保持与现有标准 forcing file 一致:
    #
    # time_counter 00:00 <- HRDPS P001
    # time_counter 01:00 <- HRDPS P002
    # ...
    # time_counter 23:00 <- HRDPS P024
    #
    # 对累计变量:
    #
    # 00:00 -> 00-01 h 平均
    # 01:00 -> 01-02 h 平均
    # ...
    #
    # 即 interval-start timestamp convention。

    target_date = pd.Timestamp(target_date).normalize()

    output_file = (
        OUTPUT_DIR
        / f"HRDPS_1km_y{target_date:%Y}m{target_date:%m}d{target_date:%d}.nc"
    )

    print()
    print("#" * 78)
    print(f"开始处理日期：{target_date:%Y-%m-%d}")
    print(f"输入 cycle：{GRIB_ROOT / target_date.strftime('%Y%m%d') / RUN_HOUR}")
    print(f"输出文件：{output_file}")
    print("#" * 78)

    times = pd.date_range(
        target_date,
        periods=24,
        freq="h",
    )


    # ============================================================
    # 读取数据
    # ============================================================

    all_data = {}

    nav_lat = None
    nav_lon = None

    reference_shape = None


    for target_var, info in VARIABLES.items():

        grib_name = info[
            "grib_name"
        ]

        level_name = info[
            "level"
        ]

        # print()
        # print("=" * 70)
        # print(
        #     f"处理 {target_var}"
        # )
        # print("=" * 70)

        lead_data = []

        for lead in LEADS:

            path = get_file(
                target_date,
                grib_name,
                level_name,
                lead,
            )

            data, lat, lon, attrs = (
                read_grib(path)
            )

            # ----------------------------------------------------
            # shape consistency
            # ----------------------------------------------------

            if reference_shape is None:

                reference_shape = (
                    data.shape
                )

            elif data.shape != reference_shape:

                raise ValueError(
                    "网格 shape 不一致:\n"
                    f"{path}\n"
                    f"{data.shape} != "
                    f"{reference_shape}"
                )

            # ----------------------------------------------------
            # lat/lon consistency
            # ----------------------------------------------------

            if nav_lat is None:

                nav_lat = lat
                nav_lon = lon

            else:

                if not np.allclose(
                    lat,
                    nav_lat,
                    rtol=0.0,
                    atol=1.0e-8,
                ):
                    raise ValueError(
                        "latitude grid 不一致:\n"
                        f"{path}"
                    )

                if not np.allclose(
                    lon,
                    nav_lon,
                    rtol=0.0,
                    atol=1.0e-8,
                ):
                    raise ValueError(
                        "longitude grid 不一致:\n"
                        f"{path}"
                    )

            # ----------------------------------------------------
            # 检查 stepType
            # ----------------------------------------------------

            step_type = attrs.get(
                "GRIB_stepType",
                "unknown",
            )

            expected_type = info[
                "type"
            ]

            if expected_type == "accum":

                if step_type != "accum":
                    raise ValueError(
                        f"{target_var}: "
                        "预期 GRIB_stepType=accum，"
                        f"实际为 {step_type}\n"
                        f"{path}"
                    )

            elif expected_type == "instant":

                if step_type != "instant":
                    raise ValueError(
                        f"{target_var}: "
                        "预期 GRIB_stepType=instant，"
                        f"实际为 {step_type}\n"
                        f"{path}"
                    )

            # ----------------------------------------------------
            # 风必须仍然是 grid-relative
            # ----------------------------------------------------

            if target_var in (
                "u_wind",
                "v_wind",
            ):

                uv_relative = attrs.get(
                    "GRIB_uvRelativeToGrid",
                    None,
                )

                if uv_relative != 1:
                    raise ValueError(
                        f"{target_var}: "
                        "GRIB_uvRelativeToGrid "
                        f"不是 1，而是 "
                        f"{uv_relative}。\n"
                        "停止处理，避免错误或重复旋转。\n"
                        f"{path}"
                    )

            lead_data.append(
                data
            )

            forcing_time = (
                times[
                    lead - 1
                ]
            )

            valid_time = (
                target_date
                + pd.Timedelta(
                    hours=lead
                )
            )

            # print(
            #     f"  forcing "
            #     f"{forcing_time:%Y-%m-%d %H:%M}"
            #     f" <- "
            #     f"P{lead:03d}"
            #     f" "
            #     f"(HRDPS valid "
            #     f"{valid_time:%Y-%m-%d %H:%M})"
            # )

        raw = np.stack(
            lead_data,
            axis=0,
        ).astype(
            np.float32
        )

        # --------------------------------------------------------
        # accumulated -> hourly rate
        # --------------------------------------------------------

        if info["type"] == "accum":

            data = (
                accumulated_to_hourly_rate(
                    raw,
                    target_var,
                )
            )

        else:

            data = raw

        all_data[
            target_var
        ] = data


    # ============================================================
    # 风旋转
    # ============================================================

    print()
    print("=" * 70)
    print(
        "旋转 HRDPS 10 m 风到 "
        "true east / true north"
    )
    print("=" * 70)

    cos_theta, sin_theta = (
        hrdps_rotation_coefficients(
            nav_lat,
            nav_lon,
        )
    )

    u_grid = all_data[
        "u_wind"
    ]

    v_grid = all_data[
        "v_wind"
    ]

    u_east, v_north = (
        rotate_hrdps_wind(
            u_grid,
            v_grid,
            cos_theta,
            sin_theta,
        )
    )


    # ============================================================
    # 风旋转 sanity check
    # ============================================================

    speed_grid = np.hypot(
        u_grid,
        v_grid,
    )

    speed_earth = np.hypot(
        u_east,
        v_north,
    )

    max_speed_error = float(
        np.nanmax(
            np.abs(
                speed_grid
                - speed_earth
            )
        )
    )

    print(
        "最大风速模旋转误差:",
        f"{max_speed_error:.6e} m s-1",
    )

    if max_speed_error > 1.0e-4:
        raise ValueError(
            "风旋转后没有保持速度模，"
            "请检查旋转公式。"
        )

    all_data[
        "u_wind"
    ] = u_east

    all_data[
        "v_wind"
    ] = v_north


    # ============================================================
    # 建立输出 Dataset
    # ============================================================

    y_size, x_size = (
        nav_lat.shape
    )

    ds_out = xr.Dataset(
        coords={
            "time_counter": (
                "time_counter",
                times,
            ),

            "y": (
                "y",
                np.arange(
                    y_size,
                    dtype=np.int32,
                ),
            ),

            "x": (
                "x",
                np.arange(
                    x_size,
                    dtype=np.int32,
                ),
            ),
        }
    )


    # ============================================================
    # 时间坐标 metadata
    # ============================================================

    ds_out[
        "time_counter"
    ].attrs = {
        "long_name":
            "Atmospheric forcing time",

        "standard_name":
            "time",

        "axis":
            "T",

        "comment":
            "Forcing timestamps follow the existing model forcing "
            "convention. Forcing hours 00-23 correspond to HRDPS "
            "forecast leads P001-P024 from the 00 UTC forecast "
            "cycle. For accumulated variables, each timestamp "
            "marks the beginning of the represented one-hour "
            "interval.",
    }


    # ============================================================
    # x/y 坐标
    # ============================================================

    ds_out[
        "y"
    ].attrs = {
        "long_name":
            "HRDPS native grid y index",
    }

    ds_out[
        "x"
    ].attrs = {
        "long_name":
            "HRDPS native grid x index",
    }


    # ============================================================
    # nav_lat / nav_lon
    # ============================================================

    ds_out[
        "nav_lat"
    ] = (
        ("y", "x"),
        nav_lat,
    )

    ds_out[
        "nav_lon"
    ] = (
        ("y", "x"),
        nav_lon,
    )


    ds_out[
        "nav_lat"
    ].attrs = {
        "long_name":
            "Latitude",

        "standard_name":
            "latitude",

        "units":
            "degrees_north",
    }


    ds_out[
        "nav_lon"
    ].attrs = {
        "long_name":
            "Longitude",

        "standard_name":
            "longitude",

        "units":
            "degrees_east",

        "comment":
            "Longitude stored using the 0 to 360 degree convention.",
    }


    # ============================================================
    # 写入 forcing 变量
    # ============================================================

    for variable_name, data in (
        all_data.items()
    ):

        ds_out[
            variable_name
        ] = (
            (
                "time_counter",
                "y",
                "x",
            ),
            data,
        )

        ds_out[
            variable_name
        ].attrs.update(
            VARIABLE_ATTRS[
                variable_name
            ]
        )

        ds_out[
            variable_name
        ].attrs[
            "coordinates"
        ] = "nav_lon nav_lat"


    # ============================================================
    # percentcloud placeholder
    # ============================================================

    ds_out[
        "percentcloud"
    ] = (
        (
            "time_counter",
            "y",
            "x",
        ),

        np.zeros(
            (
                len(times),
                y_size,
                x_size,
            ),
            dtype=np.float32,
        ),
    )

    ds_out[
        "percentcloud"
    ].attrs.update(
        VARIABLE_ATTRS[
            "percentcloud"
        ]
    )

    ds_out[
        "percentcloud"
    ].attrs[
        "coordinates"
    ] = "nav_lon nav_lat"


    # ============================================================
    # 全局 metadata
    # ============================================================

    creation_time = (
        pd.Timestamp.now(
            tz="UTC"
        ).strftime(
            "%Y-%m-%dT%H:%M:%SZ"
        )
    )

    ds_out.attrs = {

        "title":
            "HRDPS 1 km atmospheric forcing",

        "institution":
            "Environment and Climate Change Canada; "
            "processed for ocean model forcing",

        "source":
            "ECCC HRDPS West 1 km GRIB2 forecast fields",

        "Conventions":
            "CF-1.8",

        "forecast_reference_time":
            f"{target_date:%Y-%m-%d} "
            "00:00:00 UTC",

        "forecast_cycle":
            "00 UTC",

        "forecast_leads":
            "P001 through P024",

        "time_mapping":
            "Forcing timestamps 00-23 correspond respectively to "
            "HRDPS forecast leads P001-P024 from the 00 UTC "
            "initialization. HRDPS meteorological valid times are "
            "therefore one hour later than the stored forcing "
            "timestamps.",

        "processing":
            "Instantaneous variables are taken from HRDPS forecast "
            "leads P001-P024. Accumulated shortwave radiation, "
            "longwave radiation, and precipitation are converted "
            "to hourly mean rates using successive forecast-lead "
            "differences. P001 represents the accumulation from "
            "initialization through forecast hour 1. Small negative "
            "differenced values within configured numerical "
            "tolerances are clipped to zero. Native HRDPS "
            "grid-relative 10 m wind components are rotated to "
            "true geographic eastward and northward components.",

        "wind_grid":
            "HRDPS rotated latitude-longitude grid",

        "wind_rotation":
            "Native GRIB U/V fields have "
            "GRIB_uvRelativeToGrid=1 and are rotated to true "
            "geographic east/north before writing.",

        "rotated_grid_southern_pole_latitude":
            -33.443381,

        "rotated_grid_southern_pole_longitude":
            266.463574,

        "rotated_grid_angle_of_rotation":
            0.0,

        "cloud_cover_note":
            "percentcloud is a placeholder field set uniformly "
            "to zero and is not derived from HRDPS cloud data.",

        "history":
            f"{creation_time}: "
            "created from HRDPS GRIB2 files.",
    }


    # ============================================================
    # NetCDF encoding
    # ============================================================

    encoding = {

        "time_counter": {
            "units":
                "seconds since 1970-01-01 00:00:00",

            "calendar":
                "standard",

            "dtype":
                "float64",
        },

        "nav_lat": {
            "dtype":
                "float64",

            "zlib":
                True,

            "complevel":
                4,
        },

        "nav_lon": {
            "dtype":
                "float64",

            "zlib":
                True,

            "complevel":
                4,
        },
    }


    for variable_name in (
        list(VARIABLES)
        + ["percentcloud"]
    ):

        encoding[
            variable_name
        ] = {

            "dtype":
                "float32",

            "zlib":
                True,

            "complevel":
                4,

            "_FillValue":
                np.float32(
                    1.0e20
                ),
        }


    # ============================================================
    # 写出 NetCDF
    # ============================================================

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print()
    print("=" * 70)
    print(
        "写出 NetCDF"
    )
    print("=" * 70)

    ds_out.to_netcdf(
        output_file,

        engine="netcdf4",

        # 关键：
        # time_counter 是真正的 NetCDF unlimited dimension。
        unlimited_dims=[
            "time_counter",
        ],

        encoding=encoding,
    )

    ds_out.close()


    print()
    print(
        f"生成完毕："
        f"{output_file}"
    )

    print(
        "forcing 时间范围："
        f"{times[0]:%Y-%m-%d %H:%M}"
        " -> "
        f"{times[-1]:%Y-%m-%d %H:%M}"
    )

    print(
        "time_counter = unlimited dimension"
    )

# ============================================================
# 批量处理：起止日期均包含
# ============================================================

def main():
    if END_DATE < START_DATE:
        raise ValueError(
            f"END_DATE ({END_DATE:%Y-%m-%d}) 早于 "
            f"START_DATE ({START_DATE:%Y-%m-%d})"
        )

    dates = pd.date_range(START_DATE, END_DATE, freq="D")

    print()
    print("#" * 78)
    print("HRDPS 批量 forcing 转换")
    print(
        f"日期范围：{dates[0]:%Y-%m-%d} -> "
        f"{dates[-1]:%Y-%m-%d}，共 {len(dates)} 天"
    )
    print(f"每天统一使用 {RUN_HOUR} UTC cycle")
    print(f"输出目录：{OUTPUT_DIR.resolve()}")
    print("#" * 78)

    for i, target_date in enumerate(dates, start=1):
        print()
        print(
            f"[{i}/{len(dates)}] "
            f"处理 {target_date:%Y-%m-%d}"
        )
        process_one_day(target_date)

    print()
    print("#" * 78)
    print(f"全部完成：共生成 {len(dates)} 个逐日文件")
    print("#" * 78)


if __name__ == "__main__":
    main()
