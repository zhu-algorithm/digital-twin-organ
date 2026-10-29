# stomach3d — 三维解剖学数字胃

把这个目录连同 `stomach3d.html` 一起放进 `docs/`，页面就能直接用：没有构建步骤，
没有后端，全部相对路径。`assetBase` 是 `./stomach3d/assets/`，import map 指向
`./stomach3d/vendor/three/three.module.js`，与已上线的心脏/肝脏/肾脏/小肠页一致。

```
docs/
  stomach3d.html        ← 页面
  stomach3d/
    js/                 场景模块 (main, motility, vessels, flow, gland, molecular, materials, hud)
    vendor/three/       固定版本的 three.js + 场景用到的 jsm 模块
    assets/             页面读取的全部生成产物
    ACCEPTANCE.md       逐项实测验收，含"未达成"一节
    GEOMETRY_PROVENANCE.md  来源、许可、受试者筛查、逐项 contour_anchored
```

## 资产与生成它们的脚本

| 文件 | 由谁生成 | 是什么 |
|---|---|---|
| `stomach.glb` | `build_stomach_fields.py` | 真实 CT 分割表面 + 烘进顶点的场（156,492 面，5.09 MB） |
| `stomach_geometry.npz` | `build_stomach_geometry.py` | 同样的表面，外加打包掩膜、分区标签、中心线，供下游工具用 |
| `stomach_landmarks.json` | `build_stomach_geometry.py` | 贲门、幽门、胃底顶点、中心线、大小弯迹线、角切迹、腹腔干起始、脾门 |
| `geometry_report.json` | `build_stomach_geometry.py` | 受试者筛查、面数、分区体积、实测尺寸对参考范围 |
| `stomach_vessels.json` | `build_stomach_vessels.py` | 五棵动脉树 + 胃静脉树 |
| `vessel_report.json` | `build_stomach_vessels.py` | Murray 残差、四阶段净空求解、逐对间隙、终末口径审计 |
| `stomach_flow.json` / `flow_report.json` | `build_stomach_flow.py` | 烘进路径的 v(s)=Q(s)/A(s) 与逐采样点逆流审计 |
| `stomach_motion.json` / `motion_report.json` | `build_stomach_motion.py` | 胃窦蠕动波运动学、回推、排空曲线参数 |
| `stomach_micro.json` | `build_stomach_motion.py` | 真实尺寸胃腺柱与黏液屏障 + 声明的放大倍数 |
| `stomach_model.json` | `build_stomach_model.py` | 生理与药理参数集，每个数都带来源与 `reference_located` 标记 |
| `stomach_expression.json` | `build_stomach_expression.py` (+ `repair_...py`) | 131 细胞类型 × 61 ADME 基因 |
| `stomach_proton_pump_expression.json` | 同上 | ATP4A/ATP4B，**面板外附加基因，单独报告** |
| `stomach_census_provenance.json` | 同上 | Census 溯源、区块覆盖、修复记录 |
| `stomach_fields_report.json` | `build_stomach_fields.py` | 顶点场布局、glb 预算、实测 vs 构造清单 |

## 页面内实测

页面启动时自动在**实时顶点缓冲**上跑一次运动学核验（`api.measureMotion()`），
HUD 显示的"浏览器实测"数字就是它量出来的，不是从生成文件里读的声明值。
`tools/benchmark_stomach3d.py` 在无头 Chromium 里重复这套核验，并加测
着色器链接状态、逐单步粒子逆流、回推计数、每帧渲染成本与 DOM 溢出。
