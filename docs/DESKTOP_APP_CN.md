# CartMesh2D 桌面前端

## 1. 分层

| 层 | 内容 | 代码 |
|---|---|---|
| 输入 | .xy / .dxf / .svg / .csv / .txt，内置 11 个样例 | `desktop/src/core/geometry.js`、`svg.js`、`samples.js` |
| 方法 | 纯 Cut-cell（稳定）、贴体边界层 + Cut-cell（Beta） | `capabilities.js` |
| 尺寸场 | 远场倍体长、壁面单元、每级带宽、曲率、间隙、尾迹 | `sizing.js`、`job.js` |
| 结果 | 分层级着色的 CM2D 预览、三门质量面板、层级直方图 | `renderer/viewport.js`、`app.js` |

主进程只做编排：几何转换、拉起 CLI、读报告。**所有参数校验和 argv 构造都在
`core/` 里，可以在 Node 下单测**（`npm test`，36 项）。渲染进程 `contextIsolation`
打开、CSP 收紧到 `default-src 'self'`，因此 `.cm2d` 在主进程解析后再传过去，格式
读取器只有一份。

## 2. 输入格式

- **.xy** 原生，按米，直通；
- **.dxf** 走 `cartmesh2d_dxf_cli`（单位换算和实体诊断都在那边，fail-closed）；
- **.svg** 进程内解析 `path`/`polygon`/`polyline`/`circle`/`ellipse`/`rect`，
  三次与二次 Bézier 按弦高容差自适应细分，椭圆弧走中心参数化。**不支持
  `transform`**——会警告而不是静默错位。y 轴翻转，使网格与图看起来一致；
- **.csv / .txt** 每行一对坐标，分隔符任意，空行分环，首行表头自动跳过。

无物理单位的格式（svg/csv/txt）会把最大跨度归一到 1 m。尺寸场的每个量都以体长为
单位，所以这个归一不影响任何 sizing 默认值；只有 hybrid 的首层厚度是绝对量。

## 3. 层级预算：远场和壁面抢同一个深度

树的 level 是相对**计算域**量的，所以

```
wallLevel = ceil(log2((1 + 2·farFieldSpans) · wallCellsPerSpan))
```

体长约掉了。界面把这条算式实时算出来并和实测安全上限 11 比较，超了就在生成前
拒绝，并同时给出两条出路（远场降到多少、或壁面降到多少）。`desktop/tests/sizing.test.js`
把这个预测钉在几组真实 `--size-field-only` 读数上。

曲率和间隙要求的深度取决于几何，算不出来，所以有 `预检尺寸场` 按钮：它调
`cartmesh2d_cli --size-field-only`，解析尺寸场后立刻退出，把壁面 / 曲率 / 间隙
三个层级要求和距离带、分段带、尾迹箱的数量报回来。被拒绝的请求**也会**报这些数字，
因为那正是修请求需要的信息。

## 4. 质量面板：三门分开报

界面永远请求 OpenFOAM case 目录，不是因为总要导出，而是因为 `cartmesh2d_cli` 把
solver 划分、solver 质量和 Q1 合同全部放在 `if (openFoamCase)` 里——不要 case 就
等于不要质量报告。

- **拓扑不变量**：两条 CLI 在审计失败时都返回非零，所以走到有结果就是通过；
- **Solver 质量**：这是真正拦住 OpenFOAM 导出的硬门，用户该看的就是它；
- **Q1 合同**：比 solver 门和 `checkMesh` 都严，是诊断。界面明确写出它的分类型
  计数是在 solver 凸划分之后统计的，所以归到 `cartesian` 的项可能含划分碎片。

生成器写网格、尺寸场和全部质量报告都在尝试 OpenFOAM 导出**之前**，所以导出阶段
失败时界面仍然显示网格和三门，状态改成"网格已生成，后续步骤失败"。

## 5. 样例的实测设置

11 个样例各自带一组跑通过的 `(远场, 壁面单元, α)`。**这些组合在 α 和壁面分辨率上
都不单调**——靠近层级上限时，决定 solver 门的是壁面顶点有没有擦过格线，不是参数调得
好不好。一个实测例子：48 顶点的 CSV 轮廓在远场 8× 和 9× 上过不了门，10× 和 11× 上
过，壁面分辨率一个字没动。

因此失败提示直接说这件事：改远场 ±1 倍体长会移动整个格点阵列，通常比调壁面分辨率
有效。

## 6. 已知边界

- `narrow_gap` 样例的 OpenFOAM 导出失败（`embedded boundary edge has no unique
  input-loop identity`），网格本身通过全部门。既有缺陷，不在前端范围内。
- 普通 SVG 也能撞上 `CURRENT_STATE_CN.md` 第 4 节 A 的切点擦格点缺陷：一个
  右端点恰好落在 y = -0.6 的手绘 blob 在测过的 15 组设置上全部 `source global
  topology audit failed`。这说明那个缺陷不只出现在压力网格上。
- 贴体边界层路径没有尺寸场（`--size-field` 只在 `cartmesh2d_cli` 上），实测上限
  level 8，界面直接把输入范围卡在 8。

## 7. 验证

```sh
cd desktop
npm test                                  # 36 项，core/ 的纯函数
sh scripts/build-macos.sh                 # 会校验打包出的三个 CLI 不依赖 mesasdk
node_modules/.bin/electron . --smoke=circle --out=/tmp/x --shot=/tmp/x/circle.png
```

`--smoke=<样例 id>` 走真实渲染进程：选样例、生成、把状态 / 计数 / 三门 / 直方图
打成 JSON，可选截图。`--method=hybrid` 换路径。打包脚本用 `/usr/bin/clang++` 而不是
默认编译器——mesasdk 的 g++ 在 PATH 上，它产出的二进制加载
`@rpath/lib/libstdc++.6.dylib`，进了 .app 就起不来。
