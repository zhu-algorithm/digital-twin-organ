# 肝脏肿瘤虚拟器官模型卡

## 用途与边界

这是一个可重复运行的机制研究原型，用于连接肝脏双重入流、组织氧供、正常肝细胞和肿瘤细胞群的动态。它没有使用患者影像进行个体化，也没有经过临床结局校准，不能用于诊断、预后或治疗决策。

## 状态与机制

模型用常微分方程追踪八个非负状态：正常肝细胞、肿瘤活性边缘、肿瘤核心、坏死肿瘤、功能血管，以及正常组织、肿瘤边缘和核心的氧供。门静脉和肝动脉分别贡献默认总入流的 75% 和 25%，可通过 `portal_flow_multiplier` 与 `arterial_flow_multiplier` 独立扰动；保留的 `flow_multiplier` 同比缩放两路血流并与旧接口数值等价。氧输送按两路入流各自的相对氧含量加权，默认状态归一为 1。当只给出门静脉扰动时，模型以显式 `arterial_buffer_fraction=0.25` 部分增加肝动脉流量，并受 `arterial_buffer_max_multiplier=2.5` 限制；总灌注仍不能完全恢复。若两路血流均由患者清单提供，则认为测得动脉值已经包含任何缓冲，模型不会重复补偿。默认补偿比例受麻醉猫实验约 25.5% 的结果启发（[Lautt et al., 1985](https://pubmed.ncbi.nlm.nih.gov/2579585/)），兔实验报告约 18.7%（[Mathie et al., 2003](https://pmc.ncbi.nlm.nih.gov/articles/PMC305370/)）；这些动物急性实验不是人肝肿瘤参数标定。肿瘤负荷压迫有效灌注；低氧促进核心形成、坏死与血管生成信号；正常细胞与肿瘤细胞竞争有限组织空间。

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

校准器接受长表 CSV：`day,target,value,sd`。系统目标包括 `liver_tumour_fraction`、`tumour_antigen`、`effector_tumour`、`pd1_occupancy`，以及 `liver_portal_flow`、`liver_arterial_flow`、`liver_arterial_buffer_flow`、`liver_oxygenated_inflow` 和 `liver_perfusion`；空间目标写作 `species@zone`，例如 `antigen@subcapsular_sinus` 或 `mature_dc@paracortex`。每一行必须带正的测量标准差，时间–目标组合必须唯一。可拟合参数以命名空间区分，例如 `bridge.antigen_drainage_gain`、`coupling.immune_killing_gain`、`liver.arterial_buffer_fraction`、`liver.arterial_buffer_max_multiplier` 和两路氧含量。

```bash
python scripts/calibrate_liver_tdln_bridge.py \
  --observations observations.csv \
  --fit bridge.antigen_drainage_gain,coupling.immune_killing_gain \
  --dose-days 0,21,42 --output artifacts/liver_bridge
```

HABR 参数只能在“孤立的门静脉扰动”实验语义下拟合：命令必须提供 `--portal-flow-multiplier`，并且不能同时提供 `--arterial-flow-multiplier`，否则实测动脉流已包含缓冲反应，模型会拒绝重复归因。例如：

```bash
python scripts/calibrate_liver_tdln_bridge.py \
  --observations portal_restriction_observations.csv \
  --fit liver.arterial_buffer_fraction \
  --portal-flow-multiplier 0.4 --output artifacts/habr_calibration
```

这使肝参数与桥接参数一起进入拟合文件、残差、局部协方差、不确定性抽样和队列重放，但不解决实验可辨识性：若缓冲上限在观测范围内从未激活，`arterial_buffer_max_multiplier` 仍无法由数据估计；诊断中的秩、条件数和 `identifiability_warning` 必须审查。当前恢复测试使用模型生成的合成真值，只验证软件实现，不是人类 HABR 标定。

输出包括完整参数、加权残差、局部可辨识性诊断和绑定输入／输出 SHA-256 的审计清单。缺乏足够观测或雅可比秩不足会设置 `identifiability_warning`；优化器成功不等于生物学、外部或临床验证。

仓库中的 `examples/liver_tdln_bridge_observations.csv` 仅用于 CI 参数恢复。它由模型在引流增益 0.8 下合成，每行均标记 `synthetic_demo_not_biological_evidence`；它不是患者数据、真实队列或机制有效性的证据。

## 单细胞患者特征映射

`liver_bridge_features.py` 将标签无关的单细胞汇总特征有界映射到三个机制：DC 激活影响抗原引流和 DC 产出，CD8 细胞毒性／比例提高免疫杀伤，耗竭特征降低杀伤。优先使用谱系限制的 `*_cd8_t`、`*_dc` 特征，再回退到全局签名；真正缺失的特征保持群体基线，并记录在审计的 `fallback_features` 中。所有乘数限制在 0.25–4 倍，防止任意表达尺度生成负值或无界参数。映射系数仍是待校准假设。

映射系数现在可用独立的机制测量长表进行加权鲁棒校准，而不接触 PD-1 应答标签。输入按患者提供单细胞特征，以及 `antigen_drainage_multiplier`、`dc_yield_multiplier` 或 `immune_killing_multiplier` 的正值与标准差；每位患者的重复特征必须一致。校准输出包含局部可辨识性、协方差、残差和输入／输出哈希。随后把拟合文件通过 `--feature-bridge-mapping` 传入队列反事实增强，其文件哈希、完整映射和 label-free 契约均进入可重放审计。

```bash
python scripts/calibrate_liver_feature_bridge.py \
  --observations liver_feature_bridge_observations.csv \
  --fit killing_cytotoxic_coefficient,killing_exhaustion_coefficient \
  --output artifacts/liver_feature_bridge

python scripts/augment_liver_tdln_counterfactual.py \
  --features patient_features.csv --treatment-manifest treatment_manifest.csv \
  --feature-bridge-mapping artifacts/liver_feature_bridge/fitted_liver_feature_bridge_mapping.json \
  --output artifacts/liver_counterfactual_features.csv
```

仓库示例观测是由映射方程生成的合成恢复夹具，只证明校准、审计和下游加载可执行；它不是单细胞信号与人体免疫杀伤之间的实测换算。

```bash
python scripts/simulate_patient_liver_tdln.py \
  --features patient_features.csv --patient-id P001 \
  --dose-days 0,21,42 --paired-counterfactual \
  --portal-flow-multiplier 0.8 --arterial-flow-multiplier 0.6 \
  --output artifacts/P001_liver_tdln
```

两路血流倍率也进入患者审计设置。验证器不仅重新核对输入和输出哈希，还会从绑定的患者特征行、给药方案、两路血流设置和配对反事实标志重新运行系统层与六区空间层，并比较规范化 CSV/JSON 哈希；修改审计中的血流设置而不重新生成轨迹会失败。

命令输出系统层与六区空间轨迹，以及记录所用特征、回退项、映射系数、输入文件哈希和输出哈希的患者级审计。任何疑似治疗结局字段都会在模拟前被拒绝，防止应答标签泄漏进机制特征。

`--paired-counterfactual` 使用完全相同的患者特征、初始状态、桥接参数和数值设置，同时运行自然进展与给定抗 PD-1 方案。附加 JSON 报告最终肿瘤负荷的绝对变化／相对降低、肿瘤效应细胞时间积分、峰值 PD-1 占有率以及输出淋巴效应细胞端点。它是模型内部的机制处理效应，不是 RECIST 分类、临床应答概率或因果治疗效应估计。

队列训练前可运行 `scripts/augment_liver_tdln_counterfactual.py`。它逐患者生成 11 个 `liver_cf_*` 机制特征，并保留原始单细胞特征。`--bridge-parameters` 可读取 `calibrate_liver_tdln_bridge.py` 产生的拟合参数；文件必须明确标记 `research_only=true`，内容及 SHA-256 会写入增强审计。验证时会重新模拟整个队列并比较确定性输出哈希。CI 按“桥接校准 → 队列反事实 → 预测训练”的顺序实际使用这张增强表，因此这些反事实不是仅供展示的旁路结果。

同时提供 `--bridge-diagnostics` 与至少 2 个 `--uncertainty-samples` 时，增强器按诊断文件中的局部对数参数协方差进行相关抽样。涉及 HABR 时采用联合拒绝抽样，将整体相关分布截断在同一套可执行生理边界内（缓冲比例 0–1、动脉扩张上限至少为 1），而不是逐参数裁剪；若协方差在合理尝试次数内无法产生足量合法样本则失败并要求检查边界接触与可辨识性。每个点估计特征增加 `_q025`、`_q50`、`_q975` 列；种子、样本数和诊断文件哈希都进入可重放审计。区间仅传播局部桥接参数的不确定性，不包含模型结构误差、患者抽样误差、域偏移或临床结局误差。

真实队列应以 `--treatment-manifest` 取代一个全队列共用的 `--dose-days`。CSV 必须逐患者给出 `regimen_id,agent_class,dose_days,model_dose_amount,dose_scale,response_definition`；患者集合必须与特征表完全一致，当前机制只接受 `agent_class=anti_pd1` 和明确的 `arbitrary_model_exposure` 尺度，同一训练队列必须共享一种应答定义。每位患者可有不同的给药日与模型暴露量。该清单及解析后的治疗契约进入增强审计，防止在标签治疗背景与模拟方案不一致时静默训练。

队列还可提供 `--perfusion-manifest`，逐患者声明 `portal_flow_multiplier` 与 `arterial_flow_multiplier`，以及全队列一致的 `measurement_definition` 和 `normalization_reference`。倍率必须为有限正数且患者集合与特征表严格一致；它们同时应用于同一患者的自然史和 anti-PD-1 两臂，也应用于每个桥接参数不确定性抽样。原始清单哈希及解析后的无量纲灌注契约进入增强审计并在验证时重放。仓库示例仅是合成软件夹具，不是灌注测量或生理标定证据。

后续 `augment-twin` 会识别 `liver-counterfactual-audit.json`，核对它与输入特征表的哈希绑定，并把完整载荷嵌入混合增强审计。锁定预测模型因此保留从治疗清单、桥接参数、协方差诊断到最终特征的嵌套血缘，而不是只记录一个失去上游语义的中间 CSV 哈希。

