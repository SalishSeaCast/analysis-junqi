#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""为 NEMO 生成大气强迫插值权重；修改下面配置后在 VS Code 点 ▶ 运行。

环境：Python >= 3.10，numpy、scipy、netCDF4（在 VS Code 选中的解释器中安装）。
本脚本替代 get_weight_nemo 的权重生成流程，不移植 map.F90 的场读取/时间处理。

方法：保留结构网格单元拓扑，在 geographic lon/lat 平面进行完整四边形双线性
反解；外域点采用球面大圆距离最近源点，权重 [1,0,0,0]，四个索引均合法。
使用所有单元包围盒的空间搜索，不将有限个近邻搜索失败误认为外域。
仅支持区域网格（展开后经度跨度 <180 度），包括跨日期变更线的区域网格。
不支持全球周期接缝、极点、折叠/凹/退化单元；这些情况明确报错。
源坐标必须是地理经纬度，不能使用 HRDPS 的旋转坐标 rlon/rlat。

输出沿用所给 Fortran 的 NetCDF 布局：
  src01..04: int32 (y,x), wgt01..04: float64 (y,x)
  索引从 1 开始，src = (i_python+1) + nx*j_python。
  Python field[y,x].ravel(order='C')[src-1] 等于 Fortran field(i,j) 的取值。
  角点顺序 (00,01,10,11)，与旧代码 i 外循环、j 内循环一致。
  源二维坐标 nav_lon/nav_lat 的维度为 (lat,lon)，一维则为 lon/lat。
  额外输出 extrapolated(y,x)、nearest_distance_km(y,x) 和目标经纬度。
  最近距离对所有目标计算；不是双线性插值的误差估计。

注意：不做陆海掩膜、不读取气象场值、不旋转风矢量、不做面积守恒重映射。
不同气象变量只有网格坐标、维度、裁剪和排列完全一致才可共用权重。
缺测气象值需要在上游处理（0*NaN 仍是 NaN）。
所有目标（包括陆点）均计算；坐标缺测或坏网格报错，不冒充合法外推。
先完整校验，再写临时文件并原子替换；异常不会留下看似成功的半成品。
原始代码头注明 CeCILL；如对外分发，请保留并核实上游许可要求。
"""
from pathlib import Path
from datetime import datetime, timezone
import os
import tempfile
import numpy as np
from scipy.spatial import cKDTree
from netCDF4 import Dataset

# ====================== 只需修改这一区域 ======================
# 相对路径相对于脚本所在目录，不依赖 VS Code 的工作目录。
BASE_DIR = Path(__file__).resolve().parent
# SOURCE_FILE = BASE_DIR / "hrdps_1km_sample.nc"  # 任一同网格的 forcing 文件
# TARGET_FILE = BASE_DIR / "bathy_meter.nc"       # 或含 glamt/gphit 的 mesh 文件
SOURCE_FILE = '/ocean/jqiu/forcings/HRDPS_1km/HRDPS_1km_y2023m03d01.nc'
TARGET_FILE = '/home/jqiu/analysis-junqi/Analysis_Atmospheric_Forcing/Analysis_weights/HRDPS_1km_Weights/bathy_meter.nc'
OUTPUT_FILE = BASE_DIR / "met_gem_weight.nc"

SOURCE_LON = "nav_lon"   # 一维文件可改为 "lon" / "lat"
SOURCE_LAT = "nav_lat"
TARGET_LON = "nav_lon"   # mesh 文件可改为 "glamt" / "gphit"
TARGET_LAT = "nav_lat"

# 必须填写 NetCDF 中真实的水平维度名，按 (y, x) 顺序；并非变量名。
# 明确指定维度，避免方形网格转置后仍能“正常运行”的隐蔽错误。
SOURCE_DIMS = ("y", "x")
TARGET_DIMS = ("y", "x")
# 坐标若另含 time_counter、t 等维度，明确选择某一帧；默认只接受单元素额外维度。
SOURCE_SLICES = {}        # 例如 {"time_counter": 0}
TARGET_SLICES = {}        # 例如 {"t": 0}

OVERWRITE = True
BATCH_SIZE = 10000       # 降低此值可减少目标搜索的临时内存
# None 表示不限制外推距离；设数字则超出时拒绝写文件（单位 km）。
MAX_EXTRAPOLATION_KM = None
# =============================================================
EARTH_RADIUS_KM = 6371.0088
XY_TOL = 1.0e-9          # 经纬度反解残差容差，单位 degree
UV_TOL = 1.0e-7          # 单元局部坐标的边界舍入容差
JAC_REL_TOL = 1.0e-12


def read_coordinate(ds, name, horizontal_dims, selections):
    if name not in ds.variables:
        raise ValueError(f"{ds.filepath()}: 没有变量 {name!r}")
    v = ds.variables[name]
    if len(set(horizontal_dims)) != 2:
        raise ValueError("水平维度必须为两个不同名称，顺序 (y, x)")
    key, kept = [], []
    for dim, size in zip(v.dimensions, v.shape):
        if dim in horizontal_dims:
            key.append(slice(None))
            kept.append(dim)
        else:
            index = selections.get(dim, 0 if size == 1 else None)
            if not isinstance(index, (int, np.integer)) or not 0 <= index < size:
                raise ValueError(f"{name}: 请在 SLICES 中选择维度 {dim!r} (长度 {size})")
            key.append(int(index))
    data = np.asarray(np.ma.filled(np.ma.asarray(v[tuple(key)], dtype=float), np.nan))
    if not np.isfinite(data).all():
        raise ValueError(f"{name}: 坐标含缺测/NaN/Inf，不能作为外推处理")
    standard = str(getattr(v, "standard_name", "")).lower()
    units = str(getattr(v, "units", "")).lower()
    if standard in {"grid_longitude", "grid_latitude"} or "radian" in units:
        raise ValueError(f"{name}: 必须提供以度为单位的地理经纬度")
    return data, kept


def read_grid(path, lon_name, lat_name, dims, selections):
    with Dataset(path) as ds:
        lon, ld = read_coordinate(ds, lon_name, dims, selections)
        lat, ad = read_coordinate(ds, lat_name, dims, selections)
    if lon.ndim == lat.ndim == 1:
        if ld != [dims[1]] or ad != [dims[0]]:
            raise ValueError("一维 lon 必须沿 x 维，lat 必须沿 y 维；请检查 DIMS")
        original = (lon.copy(), lat.copy())
        lon, lat = np.meshgrid(lon, lat)
    elif lon.ndim == lat.ndim == 2 and set(ld) == set(ad) == set(dims):
        lon = lon.transpose([ld.index(d) for d in dims])
        lat = lat.transpose([ad.index(d) for d in dims])
        original = (lon.copy(), lat.copy())
    else:
        raise ValueError("需要 lon(x)/lat(y)，或水平维度相同的二维经纬度")
    if lon.size == 0 or np.any(np.abs(lat) >= 90) or np.any(np.abs(lon) > 720):
        raise ValueError("坐标为空、超出经纬度范围或位于极点；请核对坐标变量")
    return np.ascontiguousarray(lon), np.ascontiguousarray(lat), original


def cross(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def xyz(lon, lat):
    lon, lat = np.deg2rad(lon), np.deg2rad(lat)
    c = np.cos(lat)
    return np.column_stack((c * np.cos(lon), c * np.sin(lon), np.sin(lat)))


class RegionalGrid:
    """一次构建搜索索引；所有 cell 保留原始行列编号。"""
    def __init__(self, lon, lat):
        if lon.shape != lat.shape or lon.ndim != 2 or min(lon.shape) < 2:
            raise ValueError("源网格必须至少为 2×2，且经纬度 shape 相同")
        if not np.isfinite(lon).all() or not np.isfinite(lat).all():
            raise ValueError("源坐标不完整")
        self.ny, self.nx = lon.shape
        if lon.size > np.iinfo(np.int32).max:
            raise ValueError("源点数超过旧 NEMO int32 索引范围")
        # 先连续展开，再将目标映射到同一经度分支。
        unwrapped = np.rad2deg(np.unwrap(np.unwrap(np.deg2rad(lon), axis=1), axis=0))
        if np.ptp(unwrapped) >= 180:
            raise ValueError("本脚本仅处理经度跨度 <180° 的区域网格，不支持全球周期网格")
        self.reference = (unwrapped.min() + unwrapped.max()) / 2
        self.lon, self.lat = unwrapped, lat
        p = np.stack((unwrapped, lat), axis=-1)
        # perimeter: 00 -> 10 -> 11 -> 01
        self.quads = np.stack((p[:-1, :-1], p[:-1, 1:], p[1:, 1:], p[1:, :-1]), axis=2).reshape(-1, 4, 2)
        edge = np.roll(self.quads, -1, axis=1) - self.quads
        turns = cross(edge, np.roll(edge, -1, axis=1))
        scale2 = np.max(np.sum(edge * edge, axis=2), axis=1)
        good_pos = np.all(turns > JAC_REL_TOL * scale2[:, None], axis=1)
        good_neg = np.all(turns < -JAC_REL_TOL * scale2[:, None], axis=1)
        good = good_pos | good_neg
        if not good.all():
            k = int(np.flatnonzero(~good)[0])
            j, i = divmod(k, self.nx - 1)
            raise ValueError(f"存在凹/自交/退化单元，共 {(~good).sum()} 个；首个左下角 Python (j,i)=({j},{i})")
        if not (good_pos.all() or good_neg.all()):
            raise ValueError("网格单元朝向不一致，疑似折叠或坐标维度错误")
        self.orientation = 1 if good_pos.all() else -1
        self.lower, self.upper = self.quads.min(axis=1), self.quads.max(axis=1)
        centers = (self.lower + self.upper) / 2
        # Chebyshev 半径取所有 bbox 的最大半宽，保证不漏任何可能包含目标的单元。
        self.radius = float(np.max((self.upper - self.lower) / 2) + 2 * XY_TOL)
        self.cell_tree = cKDTree(centers)
        self.node_tree = cKDTree(xyz(lon.ravel(), lat.ravel()))

    def invert(self, cells, points):
        q = self.quads[cells]
        p0 = q[:, 0]
        a, b = q[:, 1] - p0, q[:, 3] - p0
        c = q[:, 2] - q[:, 1] - q[:, 3] + p0
        r = points - p0
        det = cross(a, b)
        u, v = cross(r, b) / det, cross(a, r) / det
        for _ in range(30):
            residual = p0 + u[:, None]*a + v[:, None]*b + (u*v)[:, None]*c - points
            if np.max(np.abs(residual), initial=0) <= XY_TOL:
                break
            du, dv = a + v[:, None]*c, b + u[:, None]*c
            jac = cross(du, dv)
            scale = np.linalg.norm(du, axis=1) * np.linalg.norm(dv, axis=1)
            if np.any(np.abs(jac) <= JAC_REL_TOL * scale):
                raise ValueError("域内双线性反解遇到奇异 Jacobian；未将其冒充外推")
            u -= cross(residual, dv) / jac
            v -= cross(du, residual) / jac
        residual = p0 + u[:, None]*a + v[:, None]*b + (u*v)[:, None]*c - points
        ok = (np.isfinite(u) & np.isfinite(v) &
              (u >= -UV_TOL) & (u <= 1+UV_TOL) & (v >= -UV_TOL) & (v <= 1+UV_TOL) &
              (np.max(np.abs(residual), axis=1) <= XY_TOL))
        if not ok.all():
            raise ValueError("已判定在单元内，但双线性反解不收敛；停止而不是错误外推")
        u, v = np.clip(u, 0, 1), np.clip(v, 0, 1)
        w = np.column_stack(((1-u)*(1-v), (1-u)*v, u*(1-v), u*v))
        # clip 只允许边界舍入；确认 clip 后仍能复原坐标。
        if np.max(np.abs(np.einsum('nk,nkd->nd', w, q[:, [0, 3, 1, 2]]) - points)) > 3*XY_TOL:
            raise ValueError("权重不能复原目标坐标")
        return w

    def weights(self, lon, lat, batch_size=BATCH_SIZE):
        if lon.shape != lat.shape or not np.isfinite(lon).all() or not np.isfinite(lat).all():
            raise ValueError("目标经纬度不完整或 shape 不一致")
        if batch_size < 1 or lon.size == 0:
            raise ValueError("BATCH_SIZE 必须为正，目标不能为空")
        x = self.reference + (lon.ravel() - self.reference + 180) % 360 - 180
        points = np.column_stack((x, lat.ravel()))
        n = len(points)
        indices = np.empty((n, 4), dtype=np.int32)
        weights = np.zeros((n, 4), dtype=np.float64)
        extrapolated = np.ones(n, dtype=np.int8)
        distance = np.empty(n)
        for start in range(0, n, batch_size):
            stop = min(n, start + batch_size)
            p = points[start:stop]
            chord, nearest = self.node_tree.query(xyz(lon.ravel()[start:stop], lat.ravel()[start:stop]))
            distance[start:stop] = 2*EARTH_RADIUS_KM*np.arcsin(np.clip(chord/2, 0, 1))
            indices[start:stop] = nearest[:, None] + 1
            weights[start:stop, 0] = 1
            candidates = self.cell_tree.query_ball_point(p, self.radius, p=np.inf, return_sorted=True)
            counts = np.fromiter((len(c) for c in candidates), dtype=np.int64, count=len(p))
            if counts.sum():
                owners = np.repeat(np.arange(len(p)), counts)
                cells = np.concatenate(candidates).astype(np.int64)
                inside_bbox = np.all((p[owners] >= self.lower[cells]-XY_TOL) &
                                     (p[owners] <= self.upper[cells]+XY_TOL), axis=1)
                owners, cells = owners[inside_bbox], cells[inside_bbox]
                q = self.quads[cells]
                edges = np.roll(q, -1, axis=1) - q
                signed = self.orientation * cross(edges, p[owners, None] - q)
                inside = np.all(signed >= -XY_TOL*np.linalg.norm(edges, axis=2), axis=1)
                owners, cells = owners[inside], cells[inside]
                if len(owners):
                    # 共用边界可属于多个 cell；稳定选取原网格顺序中的第一个。
                    local, first = np.unique(owners, return_index=True)
                    chosen = cells[first]
                    w = self.invert(chosen, p[local])
                    jj, ii = np.divmod(chosen, self.nx-1)
                    base = jj*self.nx + ii + 1
                    indices[start+local] = base[:, None] + np.array([0, self.nx, 1, self.nx+1])
                    weights[start+local] = w
                    extrapolated[start+local] = 0
            print(f"已处理 {stop:,}/{n:,} 个目标点", flush=True)
        validate(indices, weights, self.nx*self.ny)
        return indices, weights, extrapolated, distance


def validate(indices, weights, source_size):
    if indices.shape != weights.shape or indices.ndim != 2 or indices.shape[1] != 4:
        raise ValueError("权重数组 shape 错误")
    if not np.issubdtype(indices.dtype, np.integer) or np.any((indices < 1) | (indices > source_size)):
        raise ValueError("输出存在非法源索引")
    if not np.isfinite(weights).all() or np.any((weights < 0) | (weights > 1)):
        raise ValueError("输出权重不是有限的 [0,1] 数值")
    if not np.allclose(weights.sum(axis=1), 1, rtol=0, atol=1e-12):
        raise ValueError("存在未计算的目标点或权重和不等于 1")


def write_weights(path, source_original, target_lon, target_lat, result, source_path, target_path,
                  overwrite=True):
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"输出已存在：{path}；请更改 OUTPUT_FILE 或 OVERWRITE")
    slon, slat = source_original
    ny, nx = (len(slat), len(slon)) if slon.ndim == 1 else slon.shape
    indices, weights, extrapolated, distance = result
    validate(indices, weights, nx*ny)
    shape = target_lon.shape
    if shape != target_lat.shape or indices.shape[0] != target_lon.size:
        raise ValueError("目标 shape 与权重不一致")
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.stem + ".", suffix=".tmp.nc", dir=path.parent)
    os.close(fd)
    try:
        # NetCDF3 64-bit offset，无 HDF5/压缩依赖，保留旧程序变量布局。
        with Dataset(tmp, "w", format="NETCDF3_64BIT_OFFSET") as ds:
            for name, size in (("x", shape[1]), ("y", shape[0]), ("lon", nx), ("lat", ny), ("numwgt", 4)):
                ds.createDimension(name, size)
            coords = (("lon", slon, ("lon",)), ("lat", slat, ("lat",))) if slon.ndim == 1 else (
                ("nav_lon", slon, ("lat", "lon")), ("nav_lat", slat, ("lat", "lon")))
            for name, data, dims in coords:
                var = ds.createVariable(name, "f8", dims, fill_value=False)
                var[:] = data
                var.units = "degrees_east" if "lon" in name else "degrees_north"
            for k in range(4):
                iv = ds.createVariable(f"src{k+1:02d}", "i4", ("y", "x"), fill_value=False)
                iv[:] = indices[:, k].reshape(shape)
                iv.index_base = np.int32(1)
                wv = ds.createVariable(f"wgt{k+1:02d}", "f8", ("y", "x"), fill_value=False)
                wv[:] = weights[:, k].reshape(shape)
            ev = ds.createVariable("extrapolated", "i1", ("y", "x"), fill_value=False)
            ev[:] = extrapolated.reshape(shape)
            ev.flag_values = np.array([0, 1], dtype=np.int8)
            ev.flag_meanings = "bilinear nearest_neighbor_extrapolation"
            dv = ds.createVariable("nearest_distance_km", "f8", ("y", "x"), fill_value=False)
            dv[:] = distance.reshape(shape)
            dv.units = "km"
            for name, data, units in (("target_lon", target_lon, "degrees_east"), ("target_lat", target_lat, "degrees_north")):
                var = ds.createVariable(name, "f8", ("y", "x"), fill_value=False)
                var[:] = data
                var.units = units
            ds.source_file = str(source_path)
            ds.target_file = str(target_path)
            ds.history = datetime.now(timezone.utc).isoformat() + " generated by get_weight_nemo.py"
            ds.method = "Regional geographic-lon/lat quadrilateral bilinear; spherical nearest outside domain"
            ds.source_index_formula = "1 + i_zero_based + nx * j_zero_based; field[y,x] C-order"
            ds.corner_order = "00,01,10,11"
            ds.extrapolated_count = np.int32(extrapolated.sum())
            ds.earth_radius_km = EARTH_RADIUS_KM
        # 重新打开临时文件，核对实际写入的索引和权重，而非仅内存数组。
        with Dataset(tmp) as ds:
            ri = np.column_stack([ds[f"src{k:02d}"][:].ravel() for k in range(1, 5)])
            rw = np.column_stack([ds[f"wgt{k:02d}"][:].ravel() for k in range(1, 5)])
            validate(ri, rw, nx*ny)
            if not np.array_equal(ri, indices) or not np.array_equal(rw, weights):
                raise IOError("写入后校验不一致")
        if overwrite:
            os.replace(tmp, path)
        else:
            os.link(tmp, path)  # 原子拒绝覆盖，即使其他进程刚创建了同名文件
            os.unlink(tmp)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def main():
    paths = [Path(p).expanduser().resolve() for p in (SOURCE_FILE, TARGET_FILE, OUTPUT_FILE)]
    source_path, target_path, output_path = paths
    if output_path in (source_path, target_path):
        raise ValueError("输出不能覆盖输入文件")
    if output_path.exists() and not OVERWRITE:
        raise FileExistsError(output_path)
    if MAX_EXTRAPOLATION_KM is not None and (
            not np.isfinite(MAX_EXTRAPOLATION_KM) or MAX_EXTRAPOLATION_KM < 0):
        raise ValueError("MAX_EXTRAPOLATION_KM 应为 None 或非负有限数值")
    slon, slat, original = read_grid(source_path, SOURCE_LON, SOURCE_LAT, SOURCE_DIMS, SOURCE_SLICES)
    tlon, tlat, _ = read_grid(target_path, TARGET_LON, TARGET_LAT, TARGET_DIMS, TARGET_SLICES)
    print(f"源网格 (ny,nx)={slon.shape}；目标网格={tlon.shape}", flush=True)
    grid = RegionalGrid(slon, slat)
    result = grid.weights(tlon, tlat, BATCH_SIZE)
    _, _, extra, distance = result
    n_extra = int(extra.sum())
    print(f"双线性点：{len(extra)-n_extra:,}；最近邻外推点：{n_extra:,} ({n_extra/len(extra):.2%})")
    if n_extra:
        d = distance[extra.astype(bool)]
        print(f"外推最近点距离 km：中位数={np.median(d):.3f}，最大={d.max():.3f}")
        if MAX_EXTRAPOLATION_KM is not None and d.max() > MAX_EXTRAPOLATION_KM:
            raise ValueError("外推距离超过配置上限，未写出文件")
    write_weights(output_path, original, tlon, tlat, result, source_path, target_path, OVERWRITE)
    print(f"完成并通过权重/索引校验：{output_path}")


if __name__ == "__main__":
    main()
