import os
from datetime import datetime, timedelta

import xarray as xr


# ============================================================
# 用户配置
# ============================================================

START_DATE = "2023-12-01"
END_DATE = "2023-12-31"

INPUT_DIR = (
    "/results/forcing/atmospheric/continental2.5/"
    "nemo_forcing"
)

OUTPUT_DIR = "/ocean/jqiu/forcings/HRDPS_subsampled/"

STRIDE = 4


# ============================================================
# 函数
# ============================================================

def parse_date(date_string):
    """
    把 YYYY-MM-DD 格式的字符串转换成 datetime 对象。
    """
    try:
        return datetime.strptime(date_string, "%Y-%m-%d")
    except ValueError:
        raise ValueError(
            f"Invalid date: {date_string}\n"
            "Please use YYYY-MM-DD format, for example: 2023-03-01"
        )


def generate_dates(start_date, end_date):
    """
    生成从 start_date 到 end_date 的所有日期，
    包含首尾日期。
    """
    current_date = start_date

    while current_date <= end_date:
        yield current_date
        current_date += timedelta(days=1)


def build_encoding(ds, ds_subsampled):
    """
    根据原始 dataset 的 encoding，
    为降采样后的 dataset 创建适用于 netCDF4 的 encoding。
    """

    valid_encoding_keys = {
        "shuffle",
        "szip_coding",
        "chunksizes",
        "szip_pixels_per_block",
        "blosc_shuffle",
        "quantize_mode",
        "compression",
        "complevel",
        "fletcher32",
        "least_significant_digit",
        "contiguous",
        "significant_digits",
        "endian",
        "_FillValue",
        "dtype",
        "zlib",
    }

    encoding = {}

    for var_name in ds_subsampled.variables:

        original_encoding = ds[var_name].encoding

        # 只保留 netCDF4 后端接受的参数
        var_encoding = {
            key: value
            for key, value in original_encoding.items()
            if key in valid_encoding_keys
        }

        # datetime 变量保留 CF 时间编码
        if "units" in original_encoding:
            var_encoding["units"] = original_encoding["units"]

        if "calendar" in original_encoding:
            var_encoding["calendar"] = original_encoding["calendar"]

        # 原来的 chunksize 是针对原始数组形状的，
        # 降采样后可能不合适
        var_encoding.pop("chunksizes", None)
        var_encoding.pop("contiguous", None)

        # 避免继承其他压缩格式
        var_encoding.pop("compression", None)
        var_encoding.pop("szip_coding", None)
        var_encoding.pop("szip_pixels_per_block", None)
        var_encoding.pop("blosc_shuffle", None)

        # 使用普通 zlib 压缩
        if ds_subsampled[var_name].ndim > 0:
            var_encoding["zlib"] = True
            var_encoding["complevel"] = 4
            var_encoding["shuffle"] = True

        encoding[var_name] = var_encoding

    return encoding


def process_one_day(date):
    """
    处理一天的 HRDPS 文件。
    """

    date_tag = date.strftime("y%Ym%md%d")

    input_filename = f"hrdps_{date_tag}.nc"
    output_filename = f"hrdps_subsampled_{date_tag}.nc"

    input_path = os.path.join(
        INPUT_DIR,
        input_filename,
    )

    output_path = os.path.join(
        OUTPUT_DIR,
        output_filename,
    )

    print()
    print("=" * 70)
    print(f"Processing date: {date.strftime('%Y-%m-%d')}")
    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")
    print("=" * 70)

    # --------------------------------------------------------
    # 检查输入文件是否存在
    # --------------------------------------------------------

    if not os.path.exists(input_path):
        print("WARNING: Input file does not exist.")
        print("Skipping this date.")
        return "missing"

    # --------------------------------------------------------
    # 如果输出已经存在，就跳过
    # --------------------------------------------------------

    if os.path.exists(output_path):
        print("Output file already exists.")
        print("Skipping this date.")
        return "exists"

    # --------------------------------------------------------
    # 打开并处理文件
    # --------------------------------------------------------

    with xr.open_dataset(input_path) as ds:

        print()
        print("Original dimensions:")
        print(dict(ds.sizes))

        # 检查 x 和 y 维度
        for dim in ["x", "y"]:
            if dim not in ds.dims:
                raise ValueError(
                    f"Dimension '{dim}' was not found in:\n"
                    f"{input_path}\n"
                    f"Available dimensions: {list(ds.dims)}"
                )

        # x、y 方向每 STRIDE 个点取 1 个
        ds_subsampled = ds.isel(
            x=slice(None, None, STRIDE),
            y=slice(None, None, STRIDE),
        )

        print()
        print("Subsampled dimensions:")
        print(dict(ds_subsampled.sizes))

        # 创建 encoding
        encoding = build_encoding(
            ds,
            ds_subsampled,
        )

        # time_counter 如果存在，则保持 unlimited
        if "time_counter" in ds_subsampled.dims:
            unlimited_dims = ["time_counter"]
        else:
            unlimited_dims = None

        print()
        print("Writing output file...")

        ds_subsampled.to_netcdf(
            output_path,
            mode="w",
            format="NETCDF4",
            engine="netcdf4",
            encoding=encoding,
            unlimited_dims=unlimited_dims,
        )

    print(f"Finished: {output_filename}")

    return "success"


# ============================================================
# 主程序
# ============================================================

def main():

    # 确保输出目录存在
    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    # 解析日期
    start_date = parse_date(START_DATE)
    end_date = parse_date(END_DATE)

    if end_date < start_date:
        raise ValueError(
            "END_DATE cannot be earlier than START_DATE."
        )

    number_of_days = (
        end_date - start_date
    ).days + 1

    print("=" * 70)
    print("HRDPS daily subsampling")
    print("=" * 70)

    print(f"Start date: {START_DATE}")
    print(f"End date:   {END_DATE}")
    print(f"Total days: {number_of_days}")
    print(f"Stride:     {STRIDE}")
    print(f"Input dir:  {INPUT_DIR}")
    print(f"Output dir: {OUTPUT_DIR}")

    # --------------------------------------------------------
    # 逐日处理
    # --------------------------------------------------------

    success_count = 0
    missing_count = 0
    existing_count = 0

    for date in generate_dates(
        start_date,
        end_date,
    ):

        status = process_one_day(date)

        if status == "success":
            success_count += 1

        elif status == "missing":
            missing_count += 1

        elif status == "exists":
            existing_count += 1

    # --------------------------------------------------------
    # 最终总结
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print("All processing completed.")
    print("=" * 70)

    print(f"Requested days:       {number_of_days}")
    print(f"Newly processed:      {success_count}")
    print(f"Missing input files:  {missing_count}")
    print(f"Existing outputs:     {existing_count}")
    print(f"Output folder:        {OUTPUT_DIR}")


if __name__ == "__main__":
    main()