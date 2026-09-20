# HRDPS 1 km → NEMO 权重文件故障：最终诊断报告与修改建议

## 1. 结论摘要

本次问题的核心并不是 HRDPS 1 km forcing 物理场本身有问题，而是旧版 `map.F90` 中针对二维曲线/旋转网格的 `weight_interpolation_irreg()` 存在多个长期潜伏的搜索、边界检查和错误处理缺陷。

结合源码、异常权重文件的表现，以及对 HRDPS 南侧边界与 NEMO 南侧边界的几何检查，目前最合理、最完整的故障链为：

```text
NEMO 南侧部分 target point 落在 HRDPS 1 km 实际覆盖范围之外
        ↓
旧程序只做矩形 bounding-box 检查，不能正确识别真实 curvilinear domain 边界
        ↓
这些点仍被送入 irregular-grid cell search
        ↓
sector 判定只要求 a>=0、b>=0，没有要求 a<=1、b<=1
        ↓
本应判定为“域外”的点被当成合法 cell
        ↓
插值实际上退化成不受控 extrapolation
        ↓
出现负权重、权重大于 1，甚至极端权重
        ↓
搜索进一步偏离，叠加 k0 边界检查 bug / 最后一次更新无 bounds check
        ↓
可能产生非法 source index
        ↓
某一个 target point 最终四个 sector 全部失败
        ↓
代码直接 RETURN 整个 weight_interpolation_irreg()
        ↓
后续全部 NEMO target points 不再计算
        ↓
上层没有检查权重生成是否完整
        ↓
半成品数组仍被写成一个“成功生成”的 NetCDF 权重文件
        ↓
模型运行时使用坏权重，最终出问题
```

因此，本次故障不是单一 bug，而是：

1. **真实网格边界差异是触发条件；**
2. **不受控 extrapolation 是早期异常；**
3. **索引与搜索缺陷进一步放大问题；**
4. **单点失败后全局 `RETURN` 是灾难性控制流错误；**
5. **缺乏最终完整性检查使半成品文件成功落盘。**

---

## 2. 几何背景：NEMO 南侧确实部分超出 HRDPS 1 km

当前边界检查显示：

- NEMO 南侧边界的一部分位于 HRDPS 1 km 实际覆盖范围之外；
- 之后随着经度变化，NEMO 南侧边界重新进入 HRDPS 覆盖范围；
- 超出部分主要位于陆地区域，因此在最终海洋模型积分中很可能不属于真正需要高质量 atmospheric forcing 的海洋点。

这意味着：

**“NEMO 网格有部分点不在 HRDPS 域内”本身并不一定应该导致整个系统失败。**

真正的问题是旧权重生成程序没有显式区分：

```text
正常 interpolation
有限、允许的 extrapolation
真正的搜索失败
```

而是把三者混在同一套逻辑中。

---

## 3. 当前源码中可以确认的关键问题

## 3.1 `k0` bounds check 存在明显 copy/paste bug

当前代码更新：

```fortran
k0 = k0 + nint(a)
l0 = l0 + nint(b)
```

之后写成：

```fortran
if(l0<1) l0=1
if(l0>ny) l0=ny-1

if(l0<1) l0=1
if(l0>ny) l0=ny-1
```

第二组显然应当是对 `k0` 的检查。

结果是：

- `l0` 被限制；
- `k0` 完全可能变成 `0`、负数、`nx` 或大于 `nx`；
- 后续仍会直接访问 `xgr(k0±1,l0)`；
- 还可能把非法 `(i,j)` flatten 成错误 source index。

这可以直接解释异常权重文件中的负 source index 或越界 index。

---

## 3.2 搜索结束后的最后一次 `k0/l0` 更新没有任何 bounds check

迭代搜索结束后，程序还会再次执行：

```fortran
k0 = k0 + nint(a)
l0 = l0 + nint(b)
```

随后直接访问邻点：

```fortran
xgr(k0+1,l0)
xgr(k0-1,l0)
xgr(k0,l0+1)
xgr(k0,l0-1)
```

这一段完全没有边界保护。

因此，即使修复前面的 copy/paste bug，这里仍然可以再次越界。

---

## 3.3 sector membership 判定错误：允许无限制 extrapolation

当前 sector 判断本质上是：

```fortran
if (a > -eps .and. b > -eps) then
    ! accept
endif
```

但双线性 interpolation 的局部坐标正常应满足：

```text
0 <= a <= 1
0 <= b <= 1
```

当前逻辑没有检查 `a<=1`、`b<=1`。

例如：

```text
a = 3
b = 2
```

仍会被接受。

对应四个权重：

```text
(1-a)(1-b) =  2
(1-a)b     = -4
a(1-b)     = -3
ab         =  6
```

总和仍然恰好是 1。

所以仅仅检查：

```text
sum(weight) ≈ 1
```

完全无法发现这种错误。

这正是为什么异常权重文件可以同时出现：

- 负权重；
- 大于 1 的权重；
- 极端大权重；
- 但权重和仍然接近 1。

---

## 3.4 bounding-box 检查不能判断真正的 HRDPS curvilinear domain

程序只利用：

```text
[min(lon), max(lon)] × [min(lat), max(lat)]
```

判断 target 是否“在 source domain 内”。

对于 rotated / curvilinear grid，这只是一个轴对齐矩形包围盒。

因此完全可能出现：

```text
target 在 bounding box 内
但 target 实际已经在 HRDPS 网格边界之外
```

当前 NEMO 南侧红色区域正是需要重点考虑的这种情况。

---

## 3.5 单个 target point 失败时直接 `RETURN` 整个 routine

当前代码在四个 sector 都失败后执行：

```fortran
write(*,*) 'issues with interpolation in sectors'
return
```

这里的 `RETURN` 不是：

```text
跳过当前 target point
```

而是：

```text
立即退出整个 weight_interpolation_irreg()
```

由于所有 target point 都在该 subroutine 内部双循环中处理：

```fortran
do j2 = 1, ny2
    do i2 = 1, nx2
```

所以故障行为是：

```text
前面 N 个点被处理
第 N+1 个点失败
        ↓
RETURN
        ↓
剩余全部 target point 永远不再计算
```

这可以直接解释异常权重文件中：

```text
开头少量点有数据
之后绝大部分区域为 0 / 空状态 / 未完成状态
```

---

## 3.6 source index 写入前没有 bounds validation

当前逻辑直接执行：

```fortran
i1 = ...
j1 = ...

indexij(count,i2,j2) = i1 + nx*(j1-1)
```

却没有先检查：

```fortran
1 <= i1 <= nx
1 <= j1 <= ny
```

因此错误的二维 index 可以被直接转换成一个 flattened integer 并写入文件。

---

## 3.7 determinant `d` 没有数值保护

局部坐标计算中：

```fortran
d = ux*vy - uy*vx

a = (...) / d
b = (...) / d
```

没有检查：

```fortran
abs(d) < tolerance
```

如果局部两个网格方向在 geographic lon/lat 空间中过于接近平行，则可能得到：

```text
极大 a/b
Inf
NaN
```

虽然这未必是本次事故的主要触发因素，但属于必须修复的数值稳定性问题。

---

## 4. 对本次坏权重文件的最终解释

当前最合理的过程不是：

```text
程序一开始就完全找不到 HRDPS 网格
```

而更可能是：

```text
NEMO 南侧第一个或前几个域外点
        ↓
搜索还能找到一个“看起来像”的 source neighborhood
        ↓
因为 a,b 上界没有检查
        ↓
程序把 extrapolation 当 interpolation 接受
        ↓
出现负权重、>1 权重
        ↓
继续扫描南侧 target row
        ↓
某个点偏离太大 / 搜索失败 / 索引状态恶化
        ↓
四个 sector 全失败
        ↓
RETURN 整个 routine
        ↓
剩余 target points 不再处理
```

如果 target grid 是：

```fortran
do j2=1,ny2
    do i2=1,nx2
```

而南边界对应 `j2=1`，那么这还会自然形成一种非常有规律的异常文件：

```text
南侧第一行前面少量点：
    有值，但可能是错误 extrapolation

某一点：
    search 失败

从该点开始：
    同一行后续 + 所有其它行都不再计算
```

这与当前观察到的坏权重文件行为高度一致。

---

# 5. 推荐的修改目标

不建议简单把程序改成：

```text
只允许严格 interpolation
任何 target 出 source domain 就失败
```

因为在当前 NEMO–HRDPS 配置下，NEMO 南侧确实存在少量主要位于陆地的 target points 超出 HRDPS 范围。

更适合这个模型的目标是：

```text
1. 域内点：正常 interpolation
2. 边界外少量点：允许受控、有限的 extrapolation
3. 真正无法解析的点：记录 failure
4. 所有 target 都扫描完以后统一报告
5. 只要存在 unresolved failure，就绝不发布正式 weights 文件
```

也就是：

```text
合法外推
正常报错
不乱输出
```

---

# 6. 建议的三状态处理逻辑

建议明确把每一个 target point 分类成：

```text
INTERPOLATION
CONTROLLED_EXTRAPOLATION
FAILURE
```

---

## 6.1 正常 interpolation

严格条件：

```fortran
a >= -eps .and. a <= 1.d0 + eps .and. &
b >= -eps .and. b <= 1.d0 + eps
```

对于只是浮点误差导致的极小越界，可以在确认之后裁剪：

```fortran
a = max(0.d0, min(1.d0, a))
b = max(0.d0, min(1.d0, b))
```

此时：

```text
4 个 weight 应基本位于 [0,1]
sum(weight) ≈ 1
```

---

## 6.2 Controlled extrapolation

如果明确允许边界外推，可以定义一个最大允许距离：

```fortran
real(kind=8), parameter :: aext = 0.25d0
```

例如接受：

```fortran
a >= -aext .and. a <= 1.d0 + aext .and. &
b >= -aext .and. b <= 1.d0 + aext
```

这意味着最多允许超出当前 cell 大约 `0.25` 个局部网格间距。

具体阈值不应凭经验永久写死，建议结合：

- HRDPS 1 km 实际网格尺度；
- NEMO 红色域外区域距离；
- 海陆 mask；
- 一次诊断运行的统计结果；

再决定最终值。

建议初期先打印所有 extrapolation 的：

```text
i2, j2
lon, lat
a, b
source cell
```

观察实际需要的 extrapolation 幅度。

---

# 7. 两种可选的 extrapolation 策略

## 7.1 方案 A：真正的线性 extrapolation

保留原始 `a,b`：

```fortran
w00 = (1.d0-a)*(1.d0-b)
w01 = (1.d0-a)*b
w10 = a*(1.d0-b)
w11 = a*b
```

优点：

- 数学上保持局部一阶变化趋势；
- 对非常小的域外距离比较自然。

缺点：

- 会产生负权重和 >1 权重；
- 外推越远越容易放大极值；
- 必须严格限制最大 extrapolation 距离。

如果使用该方案，则：

```text
“出现负权重”不能再自动视为错误
```

必须结合 extrapolation 状态判断。

---

## 7.2 方案 B：nearest-edge / constant extrapolation

如果 target 稍微超出 domain，则将局部坐标限制到边界：

```fortran
a_eff = max(0.d0, min(1.d0, a))
b_eff = max(0.d0, min(1.d0, b))
```

然后使用：

```fortran
w00 = (1.d0-a_eff)*(1.d0-b_eff)
...
```

这相当于：

```text
HRDPS 域外点使用最近 source cell 边缘上的值
```

优点：

- 权重仍保持 `[0,1]`；
- 不放大 atmospheric field；
- 对主要位于陆地、只是为了填满 NEMO rectangular grid 的点更安全；
- 后续 sanity checks 更简单。

缺点：

- 不是真正意义上的线性 extrapolation；
- 边界外场变成 constant extension。

### 对当前模型的建议

如果确认 NEMO 超出 HRDPS 的区域：

- 很小；
- 主要位于陆地；
- 不承担实际海洋动力计算；

那么 **nearest-edge / constant extrapolation 更保守，也更适合生产环境**。

如果后来发现有真实海洋 wet point 位于 HRDPS 外，则应重新考虑 source forcing domain，而不应依赖较大的 extrapolation 来掩盖覆盖不足。

---

# 8. 必须修复的 bounds handling

无论采用哪一种 extrapolation，以下修复都必须完成。

## 8.1 每次更新 `k0/l0` 后都保证可安全访问邻点

因为后面需要：

```text
k0-1
k0+1
l0-1
l0+1
```

建议保证：

```fortran
k0 = max(2, min(k0, nx-1))
l0 = max(2, min(l0, ny-1))
```

如果业务逻辑确实需要使用最外一圈 source point，则应单独设计 boundary-cell 逻辑，而不是让 `k0±1` 自由越界。

---

## 8.2 最后一次更新后必须再次 bounds check

任何：

```fortran
k0 = k0 + nint(a)
l0 = l0 + nint(b)
```

之后，在访问 source grid 之前都必须重新检查。

---

## 8.3 flatten index 前严格检查二维索引

必须在：

```fortran
indexij = i1 + nx*(j1-1)
```

之前执行：

```fortran
if (i1 < 1 .or. i1 > nx .or. &
    j1 < 1 .or. j1 > ny) then
    ! current target is FAILURE
endif
```

绝不能允许非法 `(i1,j1)` 被写进权重文件。

---

# 9. determinant 数值保护

任何计算：

```fortran
a = ... / d
b = ... / d
```

之前都应该做：

```fortran
if (abs(d) < d_tol) then
    ! current sector invalid
endif
```

`d_tol` 应根据：

- 经纬度量级；
- 网格尺度；
- double precision；

合理设定，而不是直接使用机器最小数。

---

# 10. 最大控制流修复：禁止单点失败导致全局 RETURN

这是本次最需要修复的行为之一。

当前：

```fortran
write(*,*) 'issues with interpolation in sectors'
return
```

建议改成：

```fortran
if (.not. found) then
    n_failed = n_failed + 1

    failed_i(n_failed) = i2
    failed_j(n_failed) = j2

    indexij(:,i2,j2) = 0
    weight(:,i2,j2)  = 0.d0

    cycle
endif
```

这里的 `0` 只是内部失败占位符。

**它绝不能被视为一个合法最终权重。**

程序应该继续扫描所有 target points，以获得完整 diagnostics。

---

# 11. 全部 target 扫描结束后统一决定成功或失败

建议统计：

```fortran
n_interp = 0
n_extrap = 0
n_failed = 0
```

最后打印：

```fortran
write(*,*) 'Interpolation points : ', n_interp
write(*,*) 'Extrapolation points : ', n_extrap
write(*,*) 'Failed points        : ', n_failed
```

如果：

```fortran
n_failed > 0
```

则：

```fortran
error stop 'Weight generation failed: unresolved target points'
```

正式权重文件不得发布。

---

# 12. 输出文件必须采用“验证通过后发布”的模式

当前最大的工程问题之一是：

```text
routine 没有算完
但程序仍然成功留下一个 weights.nc
```

推荐流程：

```text
计算全部 weights
        ↓
运行全局 sanity checks
        ↓
检查 n_failed
        ↓
全部通过
        ↓
写正式 weights 文件
```

更稳健的做法是：

```text
先写：
weights-HRDPS-1km.nc.tmp

全部验证通过后：

rename →
weights-HRDPS-1km.nc
```

这样即使程序：

- 中途 `ERROR STOP`；
- 被系统杀死；
- 出现 I/O failure；

也不会留下一个看起来正式可用的半成品文件。

---

# 13. 推荐的最终 sanity checks

正式发布权重文件前，至少检查：

## 13.1 source indices

对于所有有效 target：

```text
1 <= src <= nx*ny
```

当前 Fortran flatten 公式：

```fortran
index = i + nx*(j-1)
```

是 1-based index。

---

## 13.2 finite

所有权重必须：

```text
finite
not NaN
not Inf
```

---

## 13.3 权重和

无论 interpolation 还是 linear extrapolation：

```text
abs(sum(w)-1) < tolerance
```

都必须成立。

---

## 13.4 interpolation 权重范围

对于严格 interpolation：

```text
weight >= -small_tol
weight <= 1 + small_tol
```

---

## 13.5 extrapolation 权重限制

如果采用线性 extrapolation：

不能简单要求 `[0,1]`。

但应该限制：

```text
max(abs(weight)) < wgt_limit
```

同时限制：

```text
a,b 不得超出设定的 aext
```

例如：

```text
aext <= 0.25 或经诊断确认后的其它值
```

---

## 13.6 completeness

必须保证：

```text
n_interp + n_extrap + n_failed == nx2 * ny2
```

并且正式文件要求：

```text
n_failed == 0
```

---

# 14. 如果能够使用 NEMO land/sea mask，建议进一步优化

这是当前模型最合理的长期处理方式。

如果 weights generator 可以获得 target mask，则可以区分：

```text
wet ocean point
land point
```

推荐规则：

```text
wet point:
    优先要求严格 interpolation
    最多允许非常小的 controlled extrapolation
    如果明显超出 HRDPS domain → FAILURE

land point:
    可以允许更宽松的 nearest-edge extrapolation
    或使用明确的 harmless filler
```

这样可以避免：

```text
为了填充根本不会用于海洋计算的陆地点
反而让整个权重生成失败
```

同时又能保证真正的 ocean forcing 不被“大尺度外推”掩盖。

---

# 15. 推荐的代码控制框架

下面是建议的总体结构，而非可直接替换的完整源码。

```fortran
indexij = 0
weight  = 0.d0

n_interp = 0
n_extrap = 0
n_failed = 0

do j2 = 1, ny2
    do i2 = 1, nx2

        x = xdest(i2,j2)
        y = ydest(i2,j2)

        ! ---------------------------------
        ! 1. establish initial source guess
        ! ---------------------------------

        k0 = nx/2
        l0 = ny/2

        found = .false.

        ! ---------------------------------
        ! 2. iterative search
        ! ---------------------------------

        do iter = 1, max_iter

            if (k0 < 2 .or. k0 > nx-1 .or. &
                l0 < 2 .or. l0 > ny-1) exit

            ! compute local u/v and determinant d

            if (abs(d) < d_tol) exit

            ! compute a,b

            k_new = k0 + nint(a)
            l_new = l0 + nint(b)

            k_new = max(2, min(k_new, nx-1))
            l_new = max(2, min(l_new, ny-1))

            if (k_new == k0 .and. l_new == l0) exit

            k0 = k_new
            l0 = l_new

        enddo

        ! ---------------------------------
        ! 3. test neighboring sectors
        ! ---------------------------------

        ! First try strict interpolation:
        !
        ! 0 <= a <= 1
        ! 0 <= b <= 1
        !
        ! If not, try controlled extrapolation:
        !
        ! -aext <= a <= 1+aext
        ! -aext <= b <= 1+aext

        if (strict_cell_found) then

            n_interp = n_interp + 1

        else if (controlled_extrap_cell_found) then

            n_extrap = n_extrap + 1

            ! either:
            !   A) retain a,b for true linear extrapolation
            !
            ! or:
            !   B) clamp a,b to [0,1] for nearest-edge extension

        else

            n_failed = n_failed + 1
            cycle

        endif

        ! ---------------------------------
        ! 4. validate source corner indices
        ! ---------------------------------

        if (i1 < 1 .or. i1 > nx .or. &
            j1 < 1 .or. j1 > ny) then

            n_failed = n_failed + 1
            cycle

        endif

        ! ---------------------------------
        ! 5. calculate weights
        ! ---------------------------------

        ! ...

        ! ---------------------------------
        ! 6. validate weights
        ! ---------------------------------

        if (.not. finite(weights)) then
            n_failed = n_failed + 1
            cycle
        endif

        if (abs(sum(weights)-1.d0) > sum_tol) then
            n_failed = n_failed + 1
            cycle
        endif

    enddo
enddo

write(*,*) 'Interpolation points : ', n_interp
write(*,*) 'Extrapolation points : ', n_extrap
write(*,*) 'Failed points        : ', n_failed

if (n_failed > 0) then
    error stop 'Weight generation incomplete'
endif
```

---

# 16. 建议增加的 diagnostics

每次生成 weights 时建议自动报告：

```text
Ntarget
Ninterpolation
Nextrapolation
Nfailed

min/max source index
min/max weight
max(abs(weight))
max(abs(sum(weight)-1))

max extrapolation in a
max extrapolation in b

number of negative weights
number of weights > 1
```

如果使用 land/sea mask，还应分别报告：

```text
wet extrapolated points
land extrapolated points
wet failed points
land failed points
```

其中：

```text
wet failed > 0
```

应该直接视为严重错误。

---

# 17. 推荐的调试验证步骤

修改代码后，建议不要直接进入生产模型。

先按以下顺序验证：

## Step 1：定位旧代码第一个失败点

在旧版 `RETURN` 前打印：

```fortran
write(*,*) 'FAILED target:', &
           i2, j2, &
           xdest(i2,j2), ydest(i2,j2), &
           'k0/l0=', k0, l0, &
           'a/b=', a, b
```

将该点画到 HRDPS/NEMO 边界图上。

如果它位于南侧 red/outside 区域附近，则可以对当前故障链做进一步闭环确认。

---

## Step 2：新代码仅运行 weight generation

检查：

```text
n_failed
n_extrap
max extrapolation
weight range
source index range
```

先不要跑 NEMO。

---

## Step 3：把 extrapolated points 全部画出来

特别检查：

```text
是否都集中在 HRDPS 南侧边缘
是否主要位于 NEMO land points
是否存在 interior ocean extrapolation
```

如果 extrapolated point 大量出现在 domain 内部，则说明搜索算法仍然有问题。

---

## Step 4：独立验证少量点

挑选：

```text
正常 interior point
HRDPS 边缘 point
extrapolated land point
```

用 Python / xESMF / 独立几何方法交叉检查。

---

## Step 5：检查新 NetCDF 权重文件

确认：

```text
没有负 source index
没有 out-of-range source index
没有 NaN/Inf
没有未处理区域
没有大片无理由的 0 weight
```

---

## Step 6：再运行 NEMO

先做较短测试运行，并检查 atmospheric forcing 应用后的：

```text
min/max
NaN
局部异常跳变
边界异常
```

确认无异常之后再进入正式积分。

---

# 18. 对“允许 extrapolation”这件事的最终建议

对于当前 HRDPS 1 km → NEMO 配置：

**允许有限、明确、可统计的边界 extrapolation 是合理的。**

但必须满足：

```text
1. extrapolation 有明确距离上限
2. 每个 extrapolated point 都被记录
3. source indices 永远合法
4. 权重永远 finite
5. 不允许算法搜索失败后偷偷继续
6. unresolved failure 必须导致整个生成任务失败
7. 正式权重文件只在全局验证通过后发布
```

对于当前主要位于陆地区域的超界点，优先建议：

```text
nearest-edge / constant extrapolation
```

而不是无限制线性外推。

原因是：

- 目标不是在 HRDPS 域外重建高质量 atmospheric field；
- 只是为 rectangular NEMO grid 的少量边缘陆地点提供稳定、有限的 forcing；
- nearest-edge 不会产生极端负权重或放大 forcing 极值；
- 工程风险更低。

---

# 19. 最终推荐实施顺序

建议实际修改时按以下顺序执行：

1. 修复 `k0` copy/paste bounds bug；
2. 修复最后一次 `k0/l0` 更新后的 bounds check；
3. 添加 determinant `d` 保护；
4. 添加 `i1/j1` source index validation；
5. 删除单点失败导致整个 routine `RETURN` 的行为；
6. 增加 `INTERPOLATION / EXTRAPOLATION / FAILURE` 三状态；
7. 明确 controlled extrapolation 最大范围；
8. 对当前模型优先实现 nearest-edge extrapolation；
9. 添加 `n_interp / n_extrap / n_failed` 统计；
10. 添加完整的 weight/index sanity checks；
11. `n_failed > 0` 时 `ERROR STOP`；
12. 正式 `.nc` 只在全部 validation 通过后生成或 rename；
13. 重新生成 HRDPS 1 km weights；
14. 绘制 extrapolated target points 做空间检查；
15. 用独立方法交叉验证典型点；
16. 最后再进入 NEMO 测试运行。

---

# 20. 最终结论

本次 HRDPS 1 km 权重故障可以由旧版 `map.F90` 的源码缺陷和 NEMO/HRDPS 边界几何关系一致解释。

最关键的事实是：

```text
NEMO 南侧存在少量 HRDPS 域外点
```

这本来只是一个普通的边界处理问题。

但旧代码：

```text
无法严格识别真实 curvilinear domain
允许无限制 extrapolation
缺少完整 bounds validation
单点失败后直接 RETURN 整个 routine
最终又不检查 completeness
```

使一个局部、可控的问题演变成了整张坏权重文件。

因此推荐的修复原则不是简单地“禁止一切外推”，而是：

> **允许明确、有限、可追踪的边界 extrapolation；真正失败时完整报告并终止；任何情况下都不得输出未经全局验证的半成品权重文件。**

对于当前模型，更具体地说：

> **域内正常插值；少量域外陆地点采用受控 nearest-edge extrapolation；任何 wet-point 覆盖不足、非法 source index、数值异常或无法找到合理 source cell 的情况都必须显式失败；正式 weights 文件只在所有检查通过后发布。**

这应当作为后续修改 `map.F90` 和重新生成 HRDPS 1 km weights 的最终设计原则。
