import csv
import pandas as pd
from datetime import datetime


def _safe_float(s):
    s = str(s).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _safe_int(s):
    s = str(s).strip()
    if not s:
        return None
    try:
        return int(s)
    except ValueError:
        return None


def _parse_station_line(line):
    """
    FORMAT (A10,5X,A20,5X,A10)
    """
    line = line.rstrip("\n")

    return {
        "station_type": line[0:10].strip(),
        "station_name": line[15:35].strip(),
        "station_id": line[40:50].strip(),
    }


def _parse_admin_line(line):
    """
    FORMAT (2F10.4,F8.1,I4,2I2,I6,F8.1,E12.3,2X,A2,I4,2I3,I4) NOT Reliable!! Recheck the format for time!!
    """
    line = line.rstrip("\n")

    latitude = _safe_float(line[0:10])
    longitude = _safe_float(line[10:20])
    water_depth = _safe_float(line[20:28])

    year = _safe_int(line[29:33])
    month = _safe_int(line[33:35])
    day = _safe_int(line[35:37])

    time_raw = line[37:42].strip()

    record_length_min = _safe_float(line[42:50])
    sampling_frequency_hz = _safe_float(line[50:62])

    quality_code = line[64:66].strip()

    n_additional_params = _safe_int(line[66:70]) or 0
    n_wave_heights = _safe_int(line[70:73]) or 0
    n_wave_periods = _safe_int(line[73:76]) or 0
    n_spectral_estimates = _safe_int(line[76:80]) or 0

    # 处理 HHMM，例如 0, 30, 130, 1230
    time_hhmm = time_raw.zfill(4)

    hour = None
    minute = None
    datetime_utc = None

    if len(time_hhmm) >= 4:
        try:
            hour = int(time_hhmm[:-2])
            minute = int(time_hhmm[-2:])

            if year is not None and month is not None and day is not None:
                datetime_utc = datetime(
                    year, month, day, hour, minute
                ).isoformat() + "Z"
        except ValueError:
            pass

    return {
        "latitude": latitude,
        "longitude": longitude,
        "water_depth_m": water_depth,
        "year": year,
        "month": month,
        "day": day,
        "time_hhmm": time_hhmm,
        "hour": hour,
        "minute": minute,
        "datetime_utc": datetime_utc,
        "record_length_min": record_length_min,
        "sampling_frequency_hz": sampling_frequency_hz,
        "quality_code": quality_code,
        "n_additional_params": n_additional_params,
        "n_wave_heights": n_wave_heights,
        "n_wave_periods": n_wave_periods,
        "n_spectral_estimates": n_spectral_estimates,
    }


def _parse_additional_params(lines, n_params):
    """
    FORMAT (5(E12.5,A4))

    每行最多 5 组：
    value: E12.5, 12 chars
    code: A4, 4 chars
    """
    params = {}
    count = 0

    for line in lines:
        line = line.rstrip("\n")

        for j in range(5):
            if count >= n_params:
                break

            start = j * 16
            value_str = line[start:start + 12]
            code = line[start + 12:start + 16].strip()

            if code:
                key = code
                if key in params:
                    suffix = 2
                    while f"{key}_{suffix}" in params:
                        suffix += 1
                    key = f"{key}_{suffix}"

                params[key] = _safe_float(value_str)

            count += 1

    return params


def _parse_wave_height_period_lines(lines, n_heights, n_periods):
    """
    FORMAT (8(F6.2,A4))

    每行最多 8 组：
    value: F6.2, 6 chars
    code: A4, 4 chars

    前 n_heights 个是 wave heights，
    后 n_periods 个是 wave periods。
    """
    params = {}
    total = n_heights + n_periods
    count = 0

    for line in lines:
        line = line.rstrip("\n")

        for j in range(8):
            if count >= total:
                break

            start = j * 10
            value_str = line[start:start + 6]
            code = line[start + 6:start + 10].strip()

            if code:
                key = code
                if key in params:
                    suffix = 2
                    while f"{key}_{suffix}" in params:
                        suffix += 1
                    key = f"{key}_{suffix}"

                params[key] = _safe_float(value_str)

            count += 1

    return params


def _parse_spectral_lines(lines, n_spectral_estimates):
    """
    FORMAT (6E12.4)

    每个 spectral estimate 是三个数：
    frequency, bandwidth, density

    每行 6 个数，所以每行最多 2 个 spectral estimates。
    """
    numbers = []

    for line in lines:
        line = line.rstrip("\n")

        for j in range(6):
            start = j * 12
            value_str = line[start:start + 12]
            value = _safe_float(value_str)

            if value is not None:
                numbers.append(value)

    needed = n_spectral_estimates * 3
    numbers = numbers[:needed]

    spectra = []

    for k in range(0, len(numbers), 3):
        if k + 2 < len(numbers):
            spectra.append({
                "spectral_index": k // 3 + 1,
                "frequency_hz": numbers[k],
                "bandwidth_hz": numbers[k + 1],
                "spectral_density_m2_per_hz": numbers[k + 2],
            })

    return spectra


def _lines_needed(n_items, items_per_line):
    if n_items <= 0:
        return 0
    return (n_items + items_per_line - 1) // items_per_line


def read_formatb(path, encoding="utf-8"):
    """
    读取 MEDS ASCII FormatB 非方向波谱数据，返回 pandas.DataFrame。

    Parameters
    ----------
    path : str
        FormatB 原始文本文件路径。
    encoding : str, default "utf-8"
        文件编码。如果 utf-8 读不出来，可以试试 "latin1"。

    Returns
    -------
    pandas.DataFrame
        长格式 DataFrame。
        每一行对应一个 spectral estimate，即一个 frequency/bandwidth/density 三元组。
    """

    rows = []

    with open(path, "r", encoding=encoding, errors="replace") as f:
        # 保留固定宽度位置，但去掉换行符
        lines = [line.rstrip("\n") for line in f if line.strip()]

    i = 0
    record_number = 0

    while i < len(lines):
        record_number += 1

        if i + 1 >= len(lines):
            raise ValueError(
                f"文件在 record {record_number} 附近提前结束：缺少 station/admin line。"
            )

        station_line = lines[i]
        i += 1
        station = _parse_station_line(station_line)

        admin_line = lines[i]
        i += 1
        admin = _parse_admin_line(admin_line)

        n_add = admin["n_additional_params"]
        n_h = admin["n_wave_heights"]
        n_p = admin["n_wave_periods"]
        n_spec = admin["n_spectral_estimates"]

        n_add_lines = _lines_needed(n_add, 5)
        additional_lines = lines[i:i + n_add_lines]
        i += n_add_lines
        additional_params = _parse_additional_params(additional_lines, n_add)

        n_hp_lines = _lines_needed(n_h + n_p, 8)
        hp_lines = lines[i:i + n_hp_lines]
        i += n_hp_lines
        wave_params = _parse_wave_height_period_lines(
            hp_lines, n_h, n_p
        )

        n_spec_lines = _lines_needed(n_spec, 2)
        spectral_lines = lines[i:i + n_spec_lines]
        i += n_spec_lines
        spectra = _parse_spectral_lines(spectral_lines, n_spec)

        base = {
            "record_number": record_number,
            **station,
            **admin,
            **additional_params,
            **wave_params,
        }

        for spec in spectra:
            rows.append({
                **base,
                **spec,
            })

    df = pd.DataFrame(rows)

    # 如果有 datetime_utc，顺手转成 pandas datetime
    if "datetime_utc" in df.columns:
        df["datetime_utc"] = pd.to_datetime(
            df["datetime_utc"],
            utc=True,
            errors="coerce"
        )

    return df


def formatb_to_csv(input_path, output_path, encoding="utf-8"):
    """
    读取 FormatB 文件并保存为 CSV。

    Parameters
    ----------
    input_path : str
        FormatB 原始文本文件路径。
    output_path : str
        输出 CSV 文件路径。
    encoding : str, default "utf-8"
        输入文件编码。

    Returns
    -------
    pandas.DataFrame
        转换后的 DataFrame。
    """
    df = read_formatb(input_path, encoding=encoding)
    df.to_csv(output_path, index=False)
    return df