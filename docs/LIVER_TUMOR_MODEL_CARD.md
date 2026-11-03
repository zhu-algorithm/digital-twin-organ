# 肝脏肿瘤虚拟器官模型卡

## 用途与边界

这是一个可重复运行的机制研究原型，用于连接肝脏双重入流、组织氧供、正常肝细胞和肿瘤细胞群的动态。它没有使用患者影像进行个体化，也没有经过临床结局校准，不能用于诊断、预后或治疗决策。

## 状态与机制

模型用常微分方程追踪八个非负状态：正常肝细胞、肿瘤活性边缘、肿瘤核心、坏死肿瘤、功能血管，以及正常组织、肿瘤边缘和核心的氧供。门静脉和肝动脉分别贡献默认总入流的 75% 和 25%。肿瘤负荷压迫有效灌注；低氧促进核心形成、坏死与血管生成信号；正常细胞与肿瘤细胞竞争有限组织空间。

页面提供三个可比较的机制情景：自然进展、75% 血管生成抑制、总入流降至 65%。这些是反事实扰动，不代表具体药物或剂量。

## 肝肿瘤–TDLN 闭环

`liver_tdln_coupling.py` 将肝肿瘤活性负荷映射为引流抗原，驱动 TDLN 中的成熟 DC、Tpex 与效应 T 细胞；由 TDLN 产生并回流肿瘤的效应细胞再作用于肿瘤活性边缘。抗 PD-1 的血浆／TDLN 暴露和受体占有率沿用 TDLN 模型。细胞量尺度、免疫杀伤增益、坏死及缺氧抗原增益均为显式桥接参数，默认值是可审计假设，不是临床标定值。

闭环反事实检查分别验证受体占有率、肿瘤内效应细胞时间积分和最终活肿瘤负荷；模型并不假定每个晚期时间点的效应细胞数都随治疗单调增加。

肝来源抗原还通过一个非负、时间变化的输入进入六区 TDLN 网络：`afferent_lymph → subcapsular_sinus → paracortex`。抗原在被膜下窦被摄取并产生成熟 DC，旁皮质内发生 Tpex 扩增和分化，效应细胞经髓质进入输出淋巴。运输层在关闭反应时满足“区室总量变化 = 累计输入”的质量平衡测试；把引流增益设为零会同时消除肝来源的窦内抗原与旁皮质 DC 信号。

## 可重复性

- 实现：`src/tdln_pd1/liver_tumor.py`
- 肝肿瘤–TDLN 闭环：`src/tdln_pd1/liver_tdln_coupling.py`
- 六区空间桥接：`src/tdln_pd1/liver_tdln_spatial.py`
- 纵向桥接校准：`src/tdln_pd1/liver_bridge_assimilation.py`
- 确定性数据导出：`scripts/export_liver_tumor_dashboard.py`
- 固定网页数据：`docs/data/liver_tumor_dashboard.json`
- 自动测试：`tests/test_liver_tumor.py` 与 `tests/test_liver_tumor_dashboard.py`

## 已知限制

- 使用归一化细胞量和氧供，而非影像体积、绝对细胞数或 mmHg。
- 没有空间网格、肿瘤克隆异质性、免疫浸润或药代动力学。
- 页面上的细胞点是状态量的可视编码，不是一细胞一测序记录；正常肝脏页面中的单细胞层仍是本仓库的实测数据入口。
- 参数用于产生稳定、可审计的参考轨迹；需要以纵向影像、灌注和病理数据进行外部标定后，才可讨论预测性能。

## 纵向校准输入

校准器接受长表 CSV：`day,target,value,sd`。系统目标包括 `liver_tumour_fraction`、`tumour_antigen`、`effector_tumour` 和 `pd1_occupancy`；空间目标写作 `species@zone`，例如 `antigen@subcapsular_sinus` 或 `mature_dc@paracortex`。每一行必须带正的测量标准差，时间–目标组合必须唯一。可拟合参数以命名空间区分，例如 `bridge.antigen_drainage_gain` 和 `coupling.immune_killing_gain`。

```bash
python scripts/calibrate_liver_tdln_bridge.py \
  --observations observations.csv \
  --fit bridge.antigen_drainage_gain,coupling.immune_killing_gain \
  --dose-days 0,21,42 --output artifacts/liver_bridge
```

输出包括完整参数、加权残差、局部可辨识性诊断和绑定输入／输出 SHA-256 的审计清单。缺乏足够观测或雅可比秩不足会设置 `identifiability_warning`；优化器成功不等于生物学、外部或临床验证。

仓库中的 `examples/liver_tdln_bridge_observations.csv` 仅用于 CI 参数恢复。它由模型在引流增益 0.8 下合成，每行均标记 `synthetic_demo_not_biological_evidence`；它不是患者数据、真实队列或机制有效性的证据。

## 单细胞患者特征映射

`liver_bridge_features.py` 将标签无关的单细胞汇总特征有界映射到三个机制：DC 激活影响抗原引流和 DC 产出，CD8 细胞毒性／比例提高免疫杀伤，耗竭特征降低杀伤。优先使用谱系限制的 `*_cd8_t`、`*_dc` 特征，再回退到全局签名；真正缺失的特征保持群体基线，并记录在审计的 `fallback_features` 中。所有乘数限制在 0.25–4 倍，防止任意表达尺度生成负值或无界参数。映射系数仍是待校准假设。

```bash
python scripts/simulate_patient_liver_tdln.py \
  --features patient_features.csv --patient-id P001 \
  --dose-days 0,21,42 --paired-counterfactual \
  --output artifacts/P001_liver_tdln
```

命令输出系统层与六区空间轨迹，以及记录所用特征、回退项、映射系数、输入文件哈希和输出哈希的患者级审计。任何疑似治疗结局字段都会在模拟前被拒绝，防止应答标签泄漏进机制特征。

`--paired-counterfactual` 使用完全相同的患者特征、初始状态、桥接参数和数值设置，同时运行自然进展与给定抗 PD-1 方案。附加 JSON 报告最终肿瘤负荷的绝对变化／相对降低、肿瘤效应细胞时间积分、峰值 PD-1 占有率以及输出淋巴效应细胞端点。它是模型内部的机制处理效应，不是 RECIST 分类、临床应答概率或因果治疗效应估计。

队列训练前可运行 `scripts/augment_liver_tdln_counterfactual.py`。它逐患者生成 11 个 `liver_cf_*` 机制特征，并保留原始单细胞特征。审计文件绑定输入、输出、剂量方案和映射系数，验证时会重新模拟整个队列并比较确定性输出哈希。CI 的 PD-1 训练管线实际使用这张增强表，因此这些反事实不是仅供展示的旁路结果。

