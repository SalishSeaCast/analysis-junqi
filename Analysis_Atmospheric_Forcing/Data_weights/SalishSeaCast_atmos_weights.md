# SalishSeaCast / NEMO 大气 forcing weights 文件说明

基于以下材料：

- Salish Sea MEOPAR 文档：**Atmospheric Forcing**
- `get_weight_nemo.F90`
- `map.F90`

以说明：

1. weights 文件是怎么构建出来的；
2. `src01..src04` 和 `wgt01..wgt04` 分别是什么；
3. 如何在 Python/xarray 里手动把大气 forcing 映射到 SalishSeaCast/NEMO 网格。

---

## 1. 这个 weights 文件是干什么的？

NEMO 的 CORE bulk atmospheric forcing 接口可以使用 **Interpolation On the Fly (IOF)**，也就是在运行时把原始大气 forcing 网格上的变量插值到 NEMO 模型网格。

SalishSeaCast 的文档明确说，`get_weight_nemo` 会生成一份 weighting factors，使一个 atmospheric forcing grid 上的变量能够插值到由 bathymetry dataset 定义的 model grid。

所以 weights 文件本身**不包含插值后的大气变量**。它只保存：

- 对于每一个 NEMO 目标网格点，要去原始大气网格拿哪 4 个点；
- 这 4 个点分别乘什么权重。

对于当前这类文件，例如：

```text
Dimensions:  (lat: 105, lon: 70, y: 898, x: 398)

nav_lat  (lat, lon)
nav_lon  (lat, lon)
src01    (y, x)
wgt01    (y, x)
...
src04    (y, x)
wgt04    (y, x)
```

可以理解为：

```text
大气 source grid:   (lat, lon) = (105, 70)
                              |
                              |  src01..04 + wgt01..04
                              v
NEMO target grid:   (y, x)   = (898, 398)
```

---

## 2. 映射过去的到底是什么网格？

### 源码能确认：目标网格来自 `bathy_meter.nc`

`get_weight_nemo.F90` 首先打开：

```fortran
NF90_OPEN('bathy_meter.nc', ...)
```

然后读取：

```fortran
nav_lon
nav_lat
```

并把它们存到：

```fortran
glamt(jpi,jpj)
gphit(jpi,jpj)
```

随后 `get_weight()` 把这两个数组作为目标坐标 `xdest, ydest` 传给插值权重计算程序。

因此可以确定：

> **weights 的目标 `(y,x)` 网格就是 `bathy_meter.nc` 中 `nav_lon/nav_lat` 所定义的 NEMO ocean/model grid。**

SalishSeaCast 文档的措辞也是 “the model grid (as defined in the bathymetry dataset)”。

### 是不是 T grid？

源码里目标经纬度变量叫：

```fortran
glamt
gphit
```

按照 NEMO 的命名习惯，这表示 **T-point longitude/latitude**。而且 bathymetry 文件的 `nav_lon/nav_lat` 通常也是 NEMO T 点的位置。

所以对这份代码来说，可以把目标理解为：

> **SalishSeaCast 的 NEMO T 网格 / tracer grid，也就是 bathymetry 的 `(y,x)` 网格。**



---

## 3. weights 文件是怎么构建的？

### 3.1 读取 NEMO 目标网格

`get_weight_nemo.F90` 从：

```text
bathy_meter.nc
```

读取 NEMO 网格的：

```text
nav_lon
nav_lat
```

并获得：

```fortran
jpi = NEMO x dimension
jpj = NEMO y dimension
```

对应 NetCDF 输出中的：

```text
x = jpi
y = jpj
```

---

### 3.2 读取 atmospheric source grid

程序根据 `namelist` 找到一个 atmospheric forcing 文件，然后调用：

```fortran
CALL get_atmo_grid(sf(1))
```

文档特别说明：构建 weights 时，从 atmospheric forcing 文件真正需要的是**网格点位置**，而不是具体的大气变量值。

得到 source grid 的大小：

```fortran
nxgr
nygr
```

并在输出 weights 文件中写成：

```text
lon = nxgr
lat = nygr
```

对于二维 source grid，还会把原始大气网格坐标写为：

```text
nav_lon(lat, lon)
nav_lat(lat, lon)
```

因此 weights 文件里的 `nav_lon/nav_lat` 是 **source atmospheric grid** 的经纬度，不是 NEMO `(y,x)` 网格坐标。

---

### 3.3 明确使用 4 点双线性插值

`get_weight_nemo.F90` 写得非常直接：

```fortran
refwgt%numwgt = 4 ! only consider bilinear case
```

然后分配：

```fortran
refwgt%data_jpij(4,jpi,jpj)
refwgt%data_wgt (4,jpi,jpj)
```

并调用：

```fortran
CALL get_weight(refwgt)
```

所以每一个 NEMO `(y,x)` 点对应：

- 4 个 source index；
- 4 个 bilinear interpolation weights。

---

## 4. `src01..src04` 到底是什么？

`map.F90` 给出了明确公式。

对于 source grid 上二维位置 `(i1,j1)`，写入 weights 文件的线性索引是：

```fortran
indexij = i1 + nx * (j1 - 1)
```

因此：

\[
\mathrm{src}=i+N_x(j-1)
\]

其中：

- `i = 1 .. nx`
- `j = 1 .. ny`
- `src` 是 **Fortran 1-based linear index**。

这意味着在 Python 里使用时必须：

```python
idx = src - 1
```

才能变成 0-based index。

### 对 `(lat, lon)` 的 xarray 数组意味着什么？

如果你的大气变量在 Python 中是：

```python
field.dims == ("lat", "lon")
field.shape == (105, 70)
```

那么 `lon` 是最后一维，也是变化最快的一维。

源码公式中的 `i` 对应 source grid 的 x/lon 方向，所以 Python 中可直接：

```python
src_flat = field.values.ravel(order="C")
```

然后：

```python
src_flat[src01 - 1]
```

就能拿到对应的 source value。

换句话说：

```text
Fortran:
    src = i + nx*(j-1)

Python 对 (lat, lon):
    flat_index = ilat*nlon + ilon
```

两者只差 Fortran 从 1 开始、Python 从 0 开始。

---

## 5. `wgt01..wgt04` 是怎么得到的？

对于规则网格，`map.F90` 先找到包围目标点的 source-grid cell，然后计算目标点在这个 cell 中的相对位置：

```fortran
a = (x - xgr(k0)) / dx
b = (y - ygr(l0)) / dy
```

再定义：

```fortran
xx(0) = 1-a
xx(1) = a

yy(0) = 1-b
yy(1) = b
```

随后四个权重为：

```fortran
weight = xx(i) * yy(j)
```

也就是标准 bilinear interpolation：

\[
w_1=(1-a)(1-b)
\]

\[
w_2=(1-a)b
\]

\[
w_3=a(1-b)
\]

\[
w_4=ab
\]

最后源码还会把 4 个 weights 再归一化，使其和为 1（只要总权重大于 0）：

```fortran
w = sum(weight)
weight = weight / w
```

因此正常点应满足：

```python
wgt01 + wgt02 + wgt03 + wgt04 ~= 1.0
```

---

## 6. `src01..04` 和 `wgt01..04` 在文件里是怎么写出来的？

在内存中：

```fortran
data_jpij(4,jpi,jpj)
data_wgt (4,jpi,jpj)
```

输出到 NetCDF 时，程序逐个 `jk=1..4` 写成：

```text
src01(y,x)
wgt01(y,x)

src02(y,x)
wgt02(y,x)

src03(y,x)
wgt03(y,x)

src04(y,x)
wgt04(y,x)
```

所以对于任意目标格点 `(y,x)`：

\[
F_{NEMO}(y,x)
=
\sum_{k=1}^{4}
wgt_k(y,x)\,
F_{atmos}[src_k(y,x)]
\]

注意 `src_k` 是 1-based linear index，因此 Python 中要减 1。

---

## 7. Python 中如何手动映射一个二维 atmospheric field

假设：

```python
field
```

是某一个时刻的大气变量，维度为：

```text
(lat, lon) = (105, 70)
```

weights dataset 为：

```python
weights
```

其中：

```text
src01..src04 (y,x)
wgt01..wgt04 (y,x)
```

最直接的写法：

```python
import numpy as np
import xarray as xr


def remap_one_field(field, weights):
    """
    将一个二维 atmospheric field 从 source (lat, lon)
    映射到 SalishSeaCast/NEMO (y, x) 网格。
    """

    # 确保输入顺序正确：lat, lon
    field = field.transpose("lat", "lon")

    # lon/x 方向变化最快，对应源码中的
    # src = i + nx*(j-1)
    src_flat = field.values.ravel(order="C")

    out = np.zeros(
        (weights.sizes["y"], weights.sizes["x"]),
        dtype=np.result_type(field.dtype, np.float64),
    )

    for k in range(1, 5):
        idx = weights[f"src0{k}"].values - 1   # Fortran -> Python
        wgt = weights[f"wgt0{k}"].values
        out += src_flat[idx] * wgt

    return xr.DataArray(
        out,
        dims=("y", "x"),
        name=field.name,
        attrs=field.attrs,
    )
```

使用：

```python
mapped = remap_one_field(casr_var, weights)
```

输入：

```text
(105, 70)
```

输出：

```text
(898, 398)
```

---

## 8. 带 time 维的大气文件怎么映射？

如果大气变量是：

```text
(time, lat, lon)
```

可以一次处理整个时间维：

```python
import numpy as np
import xarray as xr


def remap_atmos_to_nemo(field, weights):
    """
    field: DataArray, dims (..., lat, lon)
    weights: Dataset, src01..04 / wgt01..04 dims (y, x)
    """

    # 把 lat/lon 放到最后两维
    other_dims = [d for d in field.dims if d not in ("lat", "lon")]
    field = field.transpose(*other_dims, "lat", "lon")

    nlat = field.sizes["lat"]
    nlon = field.sizes["lon"]

    # (..., lat, lon) -> (..., source_point)
    src_flat = field.values.reshape(*field.shape[:-2], nlat * nlon)

    target_shape = (
        *field.shape[:-2],
        weights.sizes["y"],
        weights.sizes["x"],
    )

    out = np.zeros(
        target_shape,
        dtype=np.result_type(field.dtype, np.float64),
    )

    for k in range(1, 5):
        idx = weights[f"src0{k}"].values - 1
        wgt = weights[f"wgt0{k}"].values

        # advanced indexing:
        # (..., source_point) -> (..., y, x)
        out += src_flat[..., idx] * wgt

    return xr.DataArray(
        out,
        dims=(*other_dims, "y", "x"),
        coords={d: field.coords[d] for d in other_dims if d in field.coords},
        name=field.name,
        attrs=field.attrs,
    )
```

例如：

```python
tair_nemo = remap_atmos_to_nemo(casr["tair"], weights)
```

如果输入是：

```text
(time, lat, lon)
```

输出就是：

```text
(time, y, x)
```

---

## 9. 使用之前建议做的 sanity checks

### 9.1 source index 范围

对于 source grid：

```python
nsource = weights.sizes["lat"] * weights.sizes["lon"]
```

检查：

```python
for k in range(1, 5):
    src = weights[f"src0{k}"]
    print(k, src.min().item(), src.max().item())
```

正常的 1-based index 应满足：

```text
1 <= src <= nsource
```

---

### 9.2 权重和

```python
wsum = sum(weights[f"wgt0{k}"] for k in range(1, 5))

print(wsum.min().item())
print(wsum.max().item())
```

正常有效区域应非常接近：

```text
1.0
```

---

### 9.3 source grid 大小是否与 forcing 一致

```python
assert field.sizes["lat"] == weights.sizes["lat"]
assert field.sizes["lon"] == weights.sizes["lon"]
```

这一步很重要：weights 是针对特定 atmospheric grid geometry 预计算的，不能随便拿给另一套不同 source grid 使用。

---

### 9.4 最好检查 source 经纬度

如果 atmospheric file 自己也有二维经纬度，可以进一步确认它和 weights 中：

```text
weights.nav_lat
weights.nav_lon
```

一致。

这样可以避免“shape 一样但实际网格位置不同”的情况。

---

## 10. 一个容易混淆的地方：weights 文件里没有 target nav_lon/nav_lat

这份 weights 文件中的：

```text
nav_lat(lat, lon)
nav_lon(lat, lon)
```

是 **source atmospheric grid**。

而目标 NEMO grid 的经纬度并没有作为 `(y,x)` 的 `nav_lat/nav_lon` 一并写进这个 weights 文件。

目标网格是在构建 weights 时从：

```text
bathy_meter.nc
```

读取的。

因此如果映射完成以后希望给输出 `DataArray` 加上真正的 SalishSeaCast 经纬度，应该另外从相同/对应的 NEMO bathymetry 或 mesh/grid 文件读取：

```text
nav_lon(y,x)
nav_lat(y,x)
```

然后作为坐标附上。

例如：

```python
mapped = mapped.assign_coords(
    nav_lon=(("y", "x"), bathy["nav_lon"].values),
    nav_lat=(("y", "x"), bathy["nav_lat"].values),
)
```

---

## 11. 整个流程可以概括成

```text
             atmospheric forcing file
              source grid (lat, lon)
                       |
                       | get_atmo_grid()
                       |
                       v
              source lon/lat geometry
                       |
                       |
       bathy_meter.nc  |       get_weight_nemo
       NEMO nav_lon ---+------------------+
       NEMO nav_lat ---+                  |
                                          v
                               bilinear interpolation
                                     map.F90
                                          |
                                          v
                               met_gem_weight.nc
                         src01..04 + wgt01..04
                                          |
                                          | NEMO IOF
                                          v
                              atmospheric variables
                                  on NEMO (y,x)
```

手动在 Python 里做的事情，其实就是重复最后一步：

```text
source atmospheric variable
        -> flatten source grid
        -> src01..04 取四个值
        -> 乘 wgt01..04
        -> 相加
        -> NEMO (y,x)
```

---

## 12. 总结

1. 这是 **4-point bilinear interpolation weights**。
2. `src01..src04` 是 atmospheric source grid 上的 **Fortran 1-based linear indices**。
3. 索引公式由源码明确给出：

   ```fortran
   src = i + nx*(j-1)
   ```

4. 对 Python 中 `(lat, lon)` 排列的数组，可使用：

   ```python
   field.values.ravel(order="C")
   ```

   并用：

   ```python
   src - 1
   ```

   进行索引。

5. `wgt01..wgt04` 是四个对应 source points 的 bilinear weights，正常情况下和为 1。
6. 映射目标是 **`bathy_meter.nc` 的 `nav_lon/nav_lat` 所定义的 NEMO model grid**；从源码命名 `glamt/gphit` 看，是 NEMO T-point/tracer grid。
7. weights 文件自身的 `nav_lon(lat,lon)` / `nav_lat(lat,lon)` 描述的是 **source atmospheric grid**，不是目标 NEMO 网格。

