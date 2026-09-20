# `map.F90` 工作原理、HRDPS 1 km 权重文件故障分析与修改建议

## 1. 文档目的

本文针对旧版 Fortran `map.F90` 中的大气强迫网格读取、规则/非规则网格插值以及权重文件生成逻辑进行代码审阅，重点解释：

1. 原始程序如何读取 atmospheric forcing grid；
2. `grid_type=1` 与 `grid_type=2` 时分别走什么插值路径；
3. `weight_interpolation_irreg()` 如何搜索 HRDPS 二维网格并生成四点双线性权重；
4. 为什么同一套代码在旧网格上可能长期正常，而 HRDPS West 1 km 网格会触发严重故障；
5. 当前已经可以从源码确认的 bug；
6. 建议的最低限度修复和更稳健的重构方案；
7. 新权重文件生成后必须执行的验收检查。

本分析基于当前提供的 `map.F90` 原始源码以及对两份权重文件的对比检查：一份已知可用的较粗 HRDPS/continental 权重文件，以及一份明显异常的 HRDPS 1 km 权重文件。

---

## 2. `map.F90` 的总体职责

`map.F90` 定义 `MODULE map`，主要负责以下几类工作：

- 保存 forcing 文件信息 (`FLD_N`, `FLD`)；
- 保存插值权重 (`WGT`)；
- 从 NetCDF forcing 文件识别并读取 atmospheric grid；
- 对规则网格执行双线性插值；
- 对二维非规则/曲线网格执行双线性插值；
- 预先生成 source-grid index 与 weight，使运行时可以重复使用权重而无需每个时间步重新搜索。

其中 `WGT` 的核心成员为：

```fortran
TYPE :: WGT
   INTEGER :: numwgt
   INTEGER, DIMENSION(:,:,:), POINTER :: data_jpij
   double precision, DIMENSION(:,:,:), POINTER :: data_wgt
END TYPE WGT
```

对 bilinear interpolation，`numwgt=4`。因此目标 NEMO 网格上的每一个点都会保存：

- 4 个 source grid flattened index；
- 对应的 4 个 interpolation weights。

其数学目标是：

```text
F_target = w1*F(src1) + w2*F(src2) + w3*F(src3) + w4*F(src4)
```

正常情况下应满足：

```text
0 <= wi <= 1
w1 + w2 + w3 + w4 ≈ 1
```

---

## 3. Atmospheric grid 的识别流程

### 3.1 `get_atmo_grid()`

`get_atmo_grid(sd)` 打开某个 forcing NetCDF 文件，寻找 longitude、latitude 和 forcing variable。

代码会尝试识别以下坐标名：

```text
LON / lon / nav_lon
LAT / lat / nav_lat
```

随后根据 longitude 和 latitude 的维数判断网格类型：

```fortran
IF (grid_dim(1)==grid_dim(2)) THEN
   grid_type=grid_dim(1)
ENDIF
```

因此大致有两类：

### `grid_type == 1`

经纬度分别是一维坐标：

```text
lon(x)
lat(y)
```

这属于规则/rectilinear grid。

### `grid_type == 2`

经纬度均为二维：

```text
nav_lon(x,y)
nav_lat(x,y)
```

这属于二维曲线/旋转网格。

HRDPS rotated latitude-longitude 文件在转换成 geographic longitude/latitude 后属于这一类，因此会进入 irregular-grid 路径。

---

## 4. 插值路径的选择

`get_weight()` 的逻辑非常直接：

```fortran
IF ( grid_type == 2 ) THEN
  CALL weight_interpolation_irreg(...)
  RETURN
ENDIF
```

否则进入：

```fortran
CALL weight_interpolation_reg(...)
```

因此：

```text
1-D lon/lat forcing grid
        ↓
weight_interpolation_reg()

2-D nav_lon/nav_lat forcing grid
        ↓
weight_interpolation_irreg()
```

重要的是：已知正常的权重文件本身也可以来自二维 `nav_lon/nav_lat` 网格，因此不能简单认为 “`grid_type=2` 一定坏”。真正的问题在于 `weight_interpolation_irreg()` 的搜索算法存在多个潜伏缺陷，而 HRDPS West 1 km 网格更容易触发这些缺陷。

---

# 5. `weight_interpolation_irreg()` 的原始工作原理

## 5.1 对每一个 NEMO target point 单独搜索

外层循环：

```fortran
do j2=1,ny2
  do i2=1,nx2
```

因此程序按 `(i2,j2)` 顺序扫描整个目标 NEMO 网格。

目标点经纬度为：

```fortran
x=xdest(i2,j2)
y=ydest(i2,j2)
```

如果 longitude 为负数，则：

```fortran
if (x<0) x=x+360
```

把 `[-180,180]` 转成 `[0,360]` 表示。

---

## 5.2 先做一个矩形 bounding-box 检查

程序计算：

```fortran
grid_xmin=minval(xgr)
grid_xmax=maxval(xgr)
grid_ymin=minval(ygr)
grid_ymax=maxval(ygr)
```

然后检查目标点是否满足：

```text
xmin <= x <= xmax
ymin <= y <= ymax
```

如果不满足，直接 `return`。

注意：对于 rotated/curvilinear grid，这只是整个网格的轴对齐矩形包围盒，并不等价于真正的 source domain 边界。

---

## 5.3 从 source grid 中心开始搜索

每个 target point 都从：

```fortran
k0=nx/2
l0=ny/2
```

开始。

`k0,l0` 是当前猜测的 source grid point。

程序取当前 source point 的两个局部方向：

```fortran
u = grid(k0+1,l0) - grid(k0,l0)
v = grid(k0,l0+1) - grid(k0,l0)
```

在经纬度二维空间中写成：

```text
u = (ux,uy)
v = (vx,vy)
```

目标点相对当前位置为：

```text
d = target - grid(k0,l0)
```

然后解局部线性系统：

```text
d ≈ a*u + b*v
```

源码中：

```fortran
d=ux*vy-uy*vx

a=(dx*vy-dy*vx)/d
b=(ux*dy-uy*dx)/d
```

这相当于利用当前网格的局部线性方向估计目标点距离当前 source point 大约有多少个 `k` 和 `l` 网格间隔。

然后执行：

```fortran
k0 = k0 + nint(a)
l0 = l0 + nint(b)
```

最多重复 10 次，期望逐步逼近目标附近的 source grid point。

可以把它理解为一种简化的 Newton-like / local-linear iterative search。

---

## 5.4 在 `k0,l0` 周围检查四个 sector

搜索结束后，程序以 `grid(k0,l0)` 为中心，依次测试四个方向组合：

```text
sector 1: (+k,+l)
sector 2: (-k,+l)
sector 3: (-k,-l)
sector 4: (+k,-l)
```

对于每个 sector，再计算一次 `(a,b)`。

如果认为目标点属于该 sector，则设置：

```fortran
step_k = ±1
step_l = ±1
```

并进入 interpolation。

---

## 5.5 生成四个 source index

找到 sector 后：

```fortran
i1=step_k*i+k0
j1=step_l*j+l0
indexij(count,i2,j2)=i1+nx*(j1-1)
```

这里采用 Fortran column-major 风格的一维 index：

```text
index = i + nx*(j-1)
```

四个 corner 对应：

```text
(k0,          l0)
(k0+step_k,   l0)
(k0,          l0+step_l)
(k0+step_k,   l0+step_l)
```

---

## 5.6 生成双线性权重

程序使用：

```fortran
xx(0)=1.-a
xx(1)=a

yy(0)=1.-b
yy(1)=b
```

四个权重为：

```text
(1-a)(1-b)
(1-a)b
a(1-b)
ab
```

随后把四个权重的总和记为 `w`，若 `w>0` 则除以 `w` 做归一化。

---

# 6. 已确认的关键 bug

以下问题均可以直接从当前源码确认，不依赖猜测。

---

## Bug 1：`k0` 的边界检查被错误复制成了 `l0` 检查

### 原始代码

```fortran
k0 = k0 + nint (a)
l0 = l0 + nint (b)

if(l0<1) l0=1
if(l0>ny) l0=ny-1

if(l0<1) l0=1
if(l0>ny) l0=ny-1
```

很明显，第二组检查极可能原本应为：

```fortran
if(k0<1)  k0=1
if(k0>=nx) k0=nx-1
```

但实际代码把 `l0` 的两行复制了两遍。

### 后果

`l0` 被限制，而 `k0` 可以变成：

```text
0
负数
nx
大于 nx
```

之后程序会直接访问：

```fortran
xgr(k0+1,l0)
xgr(k0-1,l0)
```

以及最终：

```fortran
index = i1 + nx*(j1-1)
```

因此错误的 `k0` 可以直接产生非法 source index。

### 与坏权重文件的对应

异常 HRDPS 1 km weights 中出现了负 source indices，例如成对的负值。由于相邻两个 corner 的 `i` 相差 1，这与错误的 `k0` 经 flatten 后形成相邻负 index 的特征一致。

### 严重程度

**Critical**。

这是当前源码中最明确、最直接的索引错误之一。

---

## Bug 2：10 次搜索结束后又更新一次 `k0,l0`，但完全不做边界检查

迭代循环结束后程序再次计算：

```fortran
a=(dx*vy-dy*vx)/d
b=(ux*dy-uy*dx)/d

k0 = k0 + nint (a)
l0 = l0 + nint (b)
```

之后立刻访问：

```fortran
xgr(k0+1,l0)
xgr(k0-1,l0)
xgr(k0,l0+1)
xgr(k0,l0-1)
```

此处既没有 `k0` check，也没有 `l0` check。

即使修复 Bug 1，如果这里不修，仍然可能越界。

### 严重程度

**Critical**。

---

## Bug 3：sector 判定只检查 `a >= 0`、`b >= 0`，没有检查上界

原始代码：

```fortran
if (a>-eps .and. b>-eps) goto 1000
```

这并不能判断 target 是否真正位于当前四边形 cell 内。

对局部平行四边形坐标，正确条件至少应近似为：

```text
0 <= a <= 1
0 <= b <= 1
```

考虑浮点误差：

```fortran
a >= -eps .and. a <= 1.d0+eps .and. &
b >= -eps .and. b <= 1.d0+eps
```

### 为什么当前判断危险

例如：

```text
a = 3
b = 2
```

当前代码仍然认为这个 sector 合法。

四个权重变成：

```text
(1-3)(1-2) =  2
(1-3)(2)   = -4
3(1-2)     = -3
3*2        =  6
```

它们的总和仍然恰好是：

```text
2 - 4 - 3 + 6 = 1
```

因此后面的：

```fortran
w=sum(weight)
weight=weight/w
```

完全无法识别这种错误。

### 与坏权重文件的对应

异常权重文件中存在：

- 负权重；
- 大于 1 的权重；
- 极端权重约十几。

这种现象与 `a>1` 或 `b>1` 被错误接受完全一致。

### 严重程度

**Critical**。

---

## Bug 4：一个 target point 找不到 sector 时，直接 `RETURN` 整个 subroutine

原始代码：

```fortran
if (a>-eps .and. b>-eps) then
   goto 1000
else
   write(*,*) 'issues with interpolation in sectors'
   return
endif
```

这里的 `return` 不是“跳过当前目标点”，而是：

> 立即退出整个 `weight_interpolation_irreg()`。

由于 target grid 的循环在 subroutine 内部：

```fortran
do j2=1,ny2
  do i2=1,nx2
```

任何一个 target point 失败都会导致后续全部 target points 不再计算。

### 与坏权重文件的对应

异常 HRDPS 1 km weights 中只有目标网格第一行开头很少一部分点被实际处理，后面的绝大部分权重为零/未完成状态。

这与：

```text
前 N 个 target 被计算
第 N+1 个 target sector search 失败
→ RETURN
→ 后面的 target 全部未计算
```

完全吻合。

### 严重程度

**Critical**。

这是“单点失败扩大成整个权重文件报废”的直接原因。

---

## Bug 5：生成 `indexij` 前没有任何 source index bounds check

代码：

```fortran
i1=step_k*i+k0
j1=step_l*j+l0
indexij(count,i2,j2)=i1+nx*(j1-1)
```

没有检查：

```text
1 <= i1 <= nx
1 <= j1 <= ny
```

如果 `k0/l0` 已经错误或 target 靠近网格边缘，代码会把非法 `(i1,j1)` 直接 flatten。

### 严重程度

**Critical**。

---

## Bug 6：除法前没有检查 determinant `d` 是否接近 0

多个位置都直接执行：

```fortran
d=ux*vy-uy*vx

a=(dx*vy-dy*vx)/d
b=(ux*dy-uy*dx)/d
```

但没有：

```fortran
if (abs(d) < tolerance) ...
```

如果局部两个网格方向在 geographic lon/lat 空间近乎共线，可能产生：

```text
极大 a/b
Inf
NaN
```

对于正常 HRDPS 网格这未必是本次主因，但属于必须修复的数值稳定性问题。

### 严重程度

**High**。

---

## Bug 7：bounding-box 检查不等于 curvilinear-domain 检查

代码只使用：

```text
[min(lon), max(lon)] × [min(lat), max(lat)]
```

判断 target 是否在 atmospheric grid 中。

对于 rotated grid，其真实覆盖区域通常是倾斜的四边形/曲线区域。某个点可以位于 global bounding box 内，却在实际网格外。

因此 bounding box 最多只能作为快速排除，不能作为真正的 inside-domain 判定。

### 严重程度

**Medium / design limitation**。

---

## Bug 8：`weight_interpolation_irreg()` 的所谓 bilinear inverse mapping 实际使用局部平行四边形近似

算法用：

```text
P ≈ P00 + a*(P10-P00) + b*(P01-P00)
```

然后使用标准：

```text
(1-a)(1-b), a(1-b), (1-a)b, ab
```

对于严格平行四边形是成立的。

但是一般 curvilinear quadrilateral 的第四个点 `P11` 不一定满足：

```text
P11 = P00 + (P10-P00) + (P01-P00)
```

因此这不是严格的一般四边形 bilinear inverse mapping。

在局部网格非常规则、曲率较小时可以是合理近似，但对于更大的 rotated domain 不应默认它永远可靠。

### 严重程度

**Medium / algorithmic approximation**。

---

# 7. 为什么旧网格可能一直正常，而 HRDPS West 1 km 出问题

这次故障不要求 forcing physical variables 有任何不一致。

权重生成实际上只关心：

```text
source nav_lon/nav_lat
source nx,ny
NEMO glamt/gphit
```

而不关心：

```text
solar
precip
tair
qair
wind
pressure
```

因此 “forcing fields 都由同一套 Python 代码一致生成” 并不能保护旧 weights generator。

## 7.1 搜索算法本身依赖局部线性近似

每个 target point 都从 source grid 中心开始，通过：

```fortran
k0 += nint(a)
l0 += nint(b)
```

做跳跃搜索。

对一个较小、旋转较弱、目标区域相对接近 source-grid 中心的网格，局部线性估计可能长期表现良好，因此 `k0` 从未越界，copy/paste bug 一直没有暴露。

HRDPS West 1 km 网格则具有：

- 更大的 source dimensions；
- 更大的地理覆盖范围；
- geographic lon/lat 空间中更明显的整体倾斜；
- 从 source-grid 中心到 Salish Sea target 的 index 距离可能较大。

因此一次 `nint(a)` / `nint(b)` jump 更容易过冲，使潜伏的 bounds-check bug 第一次真正被触发。

## 7.2 错误不是一步发生，而是多个缺陷叠加

最合理的故障链为：

```text
HRDPS 1 km rotated grid
    ↓
中心起步的 local-linear search 对某些 target 过冲
    ↓
k0 更新错误/越界
    ↓
BUG：k0 没有 bounds check
    ↓
BUG：最后一次 k0/l0 更新也没有 bounds check
    ↓
搜索到错误 source neighborhood
    ↓
BUG：sector test 只要求 a,b >= 0
    ↓
错误 cell 仍可能被接受
    ↓
产生负 weight、weight > 1、极端 weight
    ↓
产生非法 flattened src index
    ↓
某个 target 错得太远，四个 sector 全部失败
    ↓
BUG：RETURN 整个 subroutine
    ↓
后续整个 NEMO target grid 未计算
    ↓
上层程序未验证 completeness
    ↓
半成品 weights 被正常写成 NetCDF
```

因此 “以前别人能跑” 只能说明以前输入没有触发这些条件，不能证明 irregular-grid 算法本身健壮。

---

# 8. `map_interpolation_irreg()` 也有同类问题

完整文件中不仅 `weight_interpolation_irreg()` 使用该搜索方法，直接 interpolation 的：

```fortran
map_interpolation_irreg()
```

也存在相同的：

- 重复 `l0` bounds check、漏掉 `k0`；
- 最后一次 `k0/l0` 更新没有 bounds check；
- sector 判定缺少 `a<=1,b<=1`；
- sector failure 后直接 `return`；
- determinant 没有保护。

因此这不是单独的 “weights writer bug”，而是整个旧版 irregular-grid interpolation implementation 的通用风险。

如果生产模式曾直接调用 `map_interpolation_irreg()` 而不是预计算 weights，也应视为需要修复。

---

# 9. 最低限度修改建议

以下改动属于“继续沿用当前算法前提下的最低安全修复”。

## 9.1 每次更新 `k0/l0` 后都做严格 bounds clamp

至少需要保证后续可安全访问 `k0±1, l0±1`。

建议定义：

```fortran
k0 = max(2, min(k0, nx-1))
l0 = max(2, min(l0, ny-1))
```

如果实际 sector search 必须允许位于边界第一/最后一个点，则应单独设计边界逻辑；不能只靠任意 clamp。

至少应替换原来的重复代码：

```fortran
if(l0<1) l0=1
if(l0>ny) l0=ny-1

if(l0<1) l0=1
if(l0>ny) l0=ny-1
```

---

## 9.2 最后一轮 `k0/l0` 更新后必须再次检查

在：

```fortran
k0 = k0 + nint(a)
l0 = l0 + nint(b)
```

之后、任何 `xgr(k0±1,l0±1)` 访问之前，必须立即做 bounds check。

---

## 9.3 sector 条件改成 `[0,1] × [0,1]`

原代码：

```fortran
if (a>-eps .and. b>-eps) goto 1000
```

建议：

```fortran
if ( a >= -eps .and. a <= 1.d0 + eps .and. &
     b >= -eps .and. b <= 1.d0 + eps ) goto 1000
```

并可在接受后把极小的 round-off 裁回区间：

```fortran
a = max(0.d0, min(1.d0, a))
b = max(0.d0, min(1.d0, b))
```

只有在已经确认误差只是 `eps` 级别时才应该 clamp；不能用 clamp 掩盖真正的大幅越界。

---

## 9.4 检查 determinant

建议：

```fortran
if (abs(d) < d_tol) then
   ! current sector invalid
endif
```

`d_tol` 应结合 coordinate magnitude/grid spacing 设置，不能简单依赖机器最小数。

---

## 9.5 生成 flattened index 前检查 `i1/j1`

```fortran
if (i1 < 1 .or. i1 > nx .or. &
    j1 < 1 .or. j1 > ny) then
    ! mark current target invalid / fail loudly
endif
```

绝不应让非法 `(i,j)` 被转换成 flattened index。

---

## 9.6 单点失败不要 `RETURN` 整张网格

当前：

```fortran
write(*,*) 'issues with interpolation in sectors'
return
```

建议改为：

```fortran
indexij(:,i2,j2) = 0
weight(:,i2,j2)  = mask
n_failed = n_failed + 1
cycle
```

但仅仅继续跑并不够。

例程结束时应输出失败数，并由调用者决定：

```text
如果 expected target domain 应全部被 source 覆盖，n_failed > 0 就 STOP / error stop。
```

这样既能获得完整 diagnostics，也不会静默输出半成品。

---

## 9.7 输出权重前增加全局 sanity checks

对每一个 target point：

### index

```text
1 <= src <= nx*ny
```

### finite

```text
all weights are finite
```

### weight range

对于真正 interpolation：

```text
wi >= -small_tolerance
wi <= 1 + small_tolerance
```

### sum

```text
abs(sum(w)-1) < tolerance
```

### completeness

```text
processed_points == nx2*ny2
```

如果 target domain 理论上全部位于 source domain 内，则任何 failure 都应该使权重生成任务失败，而不是产生一个看似正常的 `.nc` 文件。

---

# 10. 建议的代码级修复框架

下面不是完整可编译替换版，而是推荐控制逻辑。

```fortran
! initialize outputs explicitly
indexij = 0
weight  = mask
n_failed = 0

DO j2 = 1, ny2
  DO i2 = 1, nx2

    ! target coordinate
    x = xdest(i2,j2)
    y = ydest(i2,j2)

    ! establish initial guess
    k0 = nx/2
    l0 = ny/2

    found = .FALSE.

    DO attempt = 1, max_attempt

      ! k0/l0 must allow neighbor access
      IF (k0 < 2 .OR. k0 > nx-1 .OR. &
          l0 < 2 .OR. l0 > ny-1) EXIT

      ... compute local a,b ...

      IF (ABS(d) < d_tol) EXIT

      k_new = k0 + NINT(a)
      l_new = l0 + NINT(b)

      k_new = MAX(2, MIN(k_new, nx-1))
      l_new = MAX(2, MIN(l_new, ny-1))

      IF (k_new == k0 .AND. l_new == l0) EXIT

      k0 = k_new
      l0 = l_new
    END DO

    ! Test all four neighboring cells.
    ! Accept only if both coordinates are within [0,1] (with eps).

    IF (.NOT. found) THEN
       n_failed = n_failed + 1
       CYCLE
    ENDIF

    ! Validate corners before flattening
    ...

    ! Compute weights and validate them
    ...

  END DO
END DO

IF (n_failed > 0) THEN
   WRITE(*,*) 'Failed targets:', n_failed
   ERROR STOP 'Weight generation incomplete'
ENDIF
```

---

# 11. 更推荐的长期方案：不要继续依赖这套手写 cell-search

即使把上述明显 bug 全修掉，当前 irregular-grid 算法仍然具有结构性局限：

1. 每个 target 都从 source grid 中心重新开始搜索，效率低；
2. 使用 geographic lon/lat 上的局部平行四边形近似；
3. 只通过局部 `u/v` 线性化寻找 cell；
4. 没有严格的 curvilinear quadrilateral containment test；
5. 对大 rotated grid、边缘区域和强曲率不够稳健。

对于现代 workflow，更可靠的方案是使用成熟 regridding library 生成 weights，例如：

- ESMF / ESMPy；
- xESMF；
- SCRIP-compatible remapping；
- 或基于投影坐标/空间索引 (`KDTree`) 的显式 cell search，再做严格四边形 inverse mapping。

如果必须保持现有 NEMO `src01...src04 / wgt01...wgt04` 文件格式，也可以：

1. 用现代工具/算法计算 source neighbor 和 weight；
2. 再按旧 NetCDF schema 写出兼容文件。

这样可以保留下游 Fortran 使用方式，而不必继续依赖脆弱的旧搜索算法。

---

# 12. 对 HRDPS 1 km 权重文件的故障结论

根据源码与实际异常权重文件的对应关系，目前可以较有把握地作如下判断：

## 可以确认

- `weight_interpolation_irreg()` 存在明确的 `k0` bounds-check copy/paste bug；
- 最后一轮 `k0/l0` 更新无边界检查；
- sector membership 条件错误地缺少 `a<=1,b<=1`；
- flattened source index 写入前没有范围检查；
- single target failure 会 `RETURN` 整个 routine；
- determinant 没有数值保护；
- 上层旧程序没有足够的权重完整性验证。

## 与实际坏文件高度一致的表现

- 非法负 source index；
- 负权重；
- 权重大于 1，甚至达到十几；
- 只有 target grid 开头极少数点被处理；
- 后面几乎整张 target grid 未生成有效权重。

## 因此最合理的解释

HRDPS West 1 km forcing grid 本身不必有任何错误；更可能是它的尺寸和 rotated-grid 几何第一次明显触发了旧 `weight_interpolation_irreg()` 中长期潜伏的 search/index bugs。

“其它 forcing/其它 grid 以前都能跑”并不能证明这段代码是正确的，只能说明以前的输入没有进入这些失败路径。

---

# 13. 推荐实施顺序

建议按以下顺序处理：

1. **不要再使用当前 `weights-HRDPS-1km_202108.nc`。** 该文件已有明显结构性异常。
2. 修复 `k0/l0` bounds checks。
3. 修复 sector 的 `[0,1]` 判定。
4. 添加 `d`、`i1/j1`、finite-value 检查。
5. 删除单点失败导致全局 `RETURN` 的行为，改成统计 failure 并最终 fail loudly。
6. 初始化 `indexij` 和 `weight`，避免半成品数组含未定义内容。
7. 重新生成 HRDPS 1 km weights。
8. 对新文件执行自动验收测试。
9. 再用实际 forcing 做一两个字段的离线 interpolation，与 xESMF/ESMF 或独立 Python 方法交叉验证。
10. 若未来仍长期使用 HRDPS 1 km，建议用现代 regridding 工具重新实现 weights generator，而不是继续维护该 old Fortran cell-search。

---

# 14. 新 weights 文件最低验收标准

建议生成完成后自动检查：

```text
Ntarget = jpi*jpj

所有需要覆盖的 target：
    4 个 src index 全部合法
    4 个 weight 全部 finite
    min(weight) >= -1e-10（或合理浮点阈值）
    max(weight) <= 1+1e-10
    |sum(weight)-1| < 1e-10（double precision 情况）
```

同时检查：

```text
min(src) >= 1
max(src) <= nxgr*nygr
```

如果文件约定使用 0-based source index，则相应改为：

```text
0 <= src < nxgr*nygr
```

必须首先确认下游读取程序使用的是哪一种 index convention。当前 Fortran 生成公式：

```fortran
index = i + nx*(j-1)
```

是 **1-based flattened index**。

此外建议输出统计：

```text
number of target points
number processed
number failed
min/max src index
min/max weight
max |sum(weight)-1|
number of negative weights
number of weights > 1
```

只要出现：

```text
failed > 0
src out of bounds
large negative weight
weight >> 1
```

就不应发布该 weights 文件。

---

## 15. 最终结论

HRDPS 1 km 的异常权重文件可以由 `map.F90` 中的源码缺陷得到一致解释。最重要的不是 forcing physical fields 的一致性，而是旧 irregular-grid search 对 source-grid geometry 的处理不够安全。

其中最具决定性的三个问题是：

1. **`k0` 边界检查因 copy/paste 错误完全缺失；**
2. **sector 判断没有要求 `a,b <= 1`，错误 cell 可被当成合法 interpolation cell；**
3. **一个 target point 失败即 `RETURN` 整个 routine，使后续整张目标网格停止计算。**

这三者叠加后，可以自然产生本次观察到的负 source index、极端权重以及只生成少数 target weights 的现象。

因此，对 HRDPS 1 km 的优先行动应当是修复/替换 weights generator，并重新生成和严格验证权重，而不是首先怀疑新 Python forcing 生成脚本改变了网格。
