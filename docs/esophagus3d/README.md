# 三维解剖学数字食管 — esophagus3d

直接解进 `docs/` 即可：`esophagus3d.html` 与 `esophagus3d/` 并列放在 `docs/` 下，
全部路径为相对路径，`assetBase = './esophagus3d/assets/'`，
importmap 指向 `./esophagus3d/vendor/three/three.module.js`。无需构建步骤、无后端。

```
docs/
  esophagus3d.html
  esophagus3d/
    js/      main · peristalsis · vessels · flow · wall · molecular · hud
    vendor/  three.js r1xx (MIT, LICENSE 随包)
    assets/  glb + 全部 json/csv 资产
    ACCEPTANCE.md
    GEOMETRY_PROVENANCE.md
```

五个视图：器官 / 血管 / 剖面 / 微观 / 分子。
动画：吞咽（原发性蠕动）、扩张触发（继发性蠕动）、门脉高压（栅栏区血流反向 → 静脉曲张）。

**读这一页之前请先读 `ACCEPTANCE.md` 第 0 节和第 14 节。** 受试者不是其它四页的 s0236；
中心线长度超出已发表范围；膈食管裂孔处的第三处狭窄未测到；GPU 帧率未测到。
