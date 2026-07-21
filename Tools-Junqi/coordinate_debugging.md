# 使用 `pcolormesh` 绘制曲线网格时的伪影问题

## 1. 问题现象

在使用二维经纬度坐标绘制模型结果时，图上可能出现以下异常：

- 横跨整个绘图区的长条；
- 三角形或楔形色块；
- 从海域边缘延伸到远处的巨大网格单元；
- 明明已经把异常数据设成 `NaN`，伪影却仍然存在。

本例中的根本原因是：经纬度数组中存在类似

```python
lon = -1
lat = -1
```

的填充值或坏坐标。

`(-1, -1)` 从全球经纬度范围来看并不是非法值，因此仅使用下面这种全球范围检查无法识别它：

```python
(lon < -180) | (lon > 180) | (lat < -90) | (lat > 90)
```

但是对于 Salish Sea 模型网格，`(-1, -1)` 显然不属于合理的模型区域。

---

## 2. 为什么把数据设成 `NaN` 仍然不够

假设某个网格点的坐标为：

```python
lon[j, i] = -1
lat[j, i] = -1
```

即使执行：

```python
field[j, i] = np.nan
```

也不能完全解决问题。

原因是 `matplotlib.pyplot.pcolormesh()` 不仅使用数据值，还会使用经纬度数组构造网格单元的几何形状。

也就是说：

- 数据值决定网格单元的颜色；
- 经纬度决定网格单元画在哪里、形状多大。

如果某个坐标点突然从 Salish Sea 跳到 `(-1, -1)`，`pcolormesh` 可能会把这个坏点与周围的正常点连接，构造出一个异常巨大的四边形。

即使对应的数据已经是 `NaN`，坏坐标仍可能影响相邻网格单元的边界，因此产生条带、三角形或跨区域伪影。

核心结论是：

> 只 mask 数据，不足以阻止坏坐标参与 `pcolormesh` 的网格几何计算。

---

## 3. 推荐的处理方法

处理这类问题时，建议分为三步。

### 3.1 使用模型区域范围识别坏坐标

不要只检查全球经纬度范围，而应根据模型区域设置一个较宽松但合理的范围。

例如，Salish Sea 模型可以使用：

```python
valid_coordinates = (
    np.isfinite(lon)
    & np.isfinite(lat)
    & (lon >= -130)
    & (lon <= -120)
    & (lat >= 45)
    & (lat <= 53)
)
```

这样可以识别：

- `NaN` 和无穷值；
- `(-1, -1)`；
- `(0, 0)`；
- 其他虽属于全球合法经纬度、但明显不属于模型网格的填充值。

模型有效范围应略大于实际绘图范围，以避免误删边界附近的正常网格点。

### 3.2 Mask 坏坐标及其邻近数据

一个坏坐标可能影响与它相连的多个网格单元。

因此，除了遮掉坏点本身，通常还应遮掉它周围一圈数据：

```python
from scipy.ndimage import binary_dilation

bad_neighbourhood = binary_dilation(
    bad_coordinates,
    structure=np.ones((3, 3), dtype=bool),
    iterations=1,
)
```

然后：

```python
field_masked = np.ma.masked_where(
    bad_neighbourhood,
    field,
)
```

这可以避免与坏坐标相邻的网格单元继续产生不规则几何形状。

### 3.3 不要把坏坐标直接传给 `pcolormesh`

即使数据已经 masked，也应避免让 `pcolormesh` 直接看到异常坐标。

一种实用方法是：

- 数据保持 masked；
- 仅为绘图坐标创建一个副本；
- 用最近的有效坐标临时替换坏坐标。

这些替换后的坐标不会产生实际颜色，因为对应数据仍然处于 masked 状态；它们只用于保持网格几何稳定。

---

## 4. 简洁示例代码

下面的示例展示了一个可复用的处理流程。

```python
import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import binary_dilation, distance_transform_edt


def fill_bad_coordinates_nearest(lon, lat, bad_mask):
    # 用最近的有效网格点填补坏坐标。
    # 对应数据仍保持 masked；填补坐标只是为了稳定网格几何。

    if not np.any(bad_mask):
        return lon.copy(), lat.copy()

    if np.all(bad_mask):
        raise ValueError("All coordinates are invalid.")

    nearest_indices = distance_transform_edt(
        bad_mask,
        return_distances=False,
        return_indices=True,
    )

    lon_plot = lon[tuple(nearest_indices)]
    lat_plot = lat[tuple(nearest_indices)]

    return lon_plot, lat_plot


# lon、lat 和 field 均为形状相同的二维数组。

valid_coordinates = (
    np.isfinite(lon)
    & np.isfinite(lat)
    & (lon >= -130)
    & (lon <= -120)
    & (lat >= 45)
    & (lat <= 53)
)

bad_coordinates = ~valid_coordinates


# 坏坐标可能影响相邻网格单元，因此扩大一圈 mask。
bad_neighbourhood = binary_dilation(
    bad_coordinates,
    structure=np.ones((3, 3), dtype=bool),
    iterations=1,
)


# Mask 数据。
field_masked = np.ma.masked_where(
    bad_neighbourhood | ~np.isfinite(field),
    field,
)


# 创建不会破坏网格几何的绘图坐标。
lon_plot, lat_plot = fill_bad_coordinates_nearest(
    lon,
    lat,
    bad_coordinates,
)


plt.figure(figsize=(8, 6))

plt.pcolormesh(
    lon_plot,
    lat_plot,
    field_masked,
    shading="auto",
    cmap="viridis",
)

plt.colorbar(label="Field value")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.show()
```

---

## 5. 建议加入的诊断检查

绘图前可以快速检查常见填充值：

```python
minus_one_points = (
    np.isclose(lon, -1)
    | np.isclose(lat, -1)
)

print(
    "Points with lon=-1 or lat=-1:",
    np.count_nonzero(minus_one_points),
)
```

也可以输出坏点的模型索引：

```python
j_bad, i_bad = np.where(minus_one_points)

print(
    list(zip(j_bad[:20], i_bad[:20]))
)
```

另外，建议始终检查坐标范围：

```python
print("Longitude range:", np.nanmin(lon), np.nanmax(lon))
print("Latitude range:", np.nanmin(lat), np.nanmax(lat))
```

如果结果中出现与模型区域完全无关的值，例如：

```text
Longitude range: -125.4  -1.0
Latitude range:   -1.0  50.9
```

基本可以确定坐标数组中存在填充值。

---

## 6. 如何避免类似问题

今后使用二维曲线网格和 `pcolormesh` 时，可以遵循以下原则：

1. 不要只按全球经纬度范围判断坐标是否合法。
2. 根据具体模型区域设置合理的有效坐标范围。
3. 数据 mask 和坐标清理必须同时进行。
4. 坏坐标附近的网格单元也可能受到影响，应适当扩大 mask。
5. 不要将 `NaN`、无穷值或明显错误的填充坐标直接传给 `pcolormesh`。
6. 绘图前打印坐标范围和异常点数量。
7. 优先使用 `shading="auto"`，减少数组形状与网格解释方式不一致带来的问题。

---

## 7. 总结

本次伪影不是速度差值计算造成的，也不是两个模拟网格不一致造成的。

真正的问题是：

```text
经纬度数组中存在 (-1, -1) 等模型填充值
        ↓
坏坐标仍被传入 pcolormesh
        ↓
pcolormesh 使用坏坐标构造异常巨大的网格单元
        ↓
图中出现条带、三角形或跨区域伪影
```

正确的处理方式是：

```text
识别模型区域之外的坏坐标
        ↓
mask 坏点及其邻近数据
        ↓
为绘图坐标临时填补坏点
        ↓
使用清理后的坐标调用 pcolormesh
```
