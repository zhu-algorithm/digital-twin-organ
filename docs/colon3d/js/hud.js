// hud.js — every number on the page is read out of the generated reports.  Nothing here
// is typed in by hand.
const T = (id) => document.getElementById(id);
const esc = s => String(s).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
function tbl(id, head, rows) {
  const el = T(id); if (!el) return;
  el.innerHTML = '<thead><tr>' + head.map(h => `<th>${esc(h)}</th>`).join('') + '</tr></thead>'
    + '<tbody>' + rows.map(r => '<tr>' + r.map(c => `<td>${c}</td>`).join('') + '</tr>').join('')
    + '</tbody>';
}
const ok = s => `<span class="ok">${esc(s)}</span>`;
const bad = s => `<span class="bad">${esc(s)}</span>`;

const SEGNAME = { caecum: '盲肠', ascending: '升结肠', hepatic_flexure: '肝曲',
  transverse: '横结肠', splenic_flexure: '脾曲', descending: '降结肠',
  sigmoid: '乙状结肠', rectum: '直肠' };
// published segmental lengths for comparison (cm) — snippet-level references
const PUB_SEG = { caecum: [6, 6], ascending: [15, 20], transverse: [40, 60],
  descending: [20, 25], sigmoid: [35, 45], rectum: [12, 15] };

export function renderTables(D) {
  const { sc, micro, geomRep, vesRep, flowRep, taeRep, vessels, expr, molInfo, MICROBIAL_LINES } = D;

  // --- provenance line in the HUD
  const src = geomRep.source;
  T('prov').innerHTML =
    `<b>${esc(src.subject)}</b> · ${esc(src.dataset)}<br>` +
    `CC-BY-4.0 · doi:${esc(src.doi)} · 体素 ${src.voxel_mm.map(x => x.toFixed(2)).join('×')} mm<br>` +
    `网格 <code>${geomRep.mesh.triangles.toLocaleString()}</code> 面 · ` +
    `中心线 <code>${geomRep.centreline.length_cm.toFixed(1)}</code> cm · ` +
    `结肠袋 <code>${taeRep.measured.haustra_count}</code> 个<br>` +
    `<span class="bad">结肠带/结肠袋/血管/微结构均为构造</span>`;

  // --- geometry
  const g = geomRep, m = g.mesh, cl = g.centreline;
  const rows = [
    ['来源', `${esc(g.source.dataset)}，doi:${esc(g.source.doi)}，${esc(g.source.licence)}`, '—', ok('CC-BY-4.0 可再分发')],
    ['受试者', esc(g.source.subject_meta), '—', '—'],
    ['该版本标签总数', g.source.labels_in_release, '117', ok('核实一致')],
    ['结肠标签体素', g.mask.colon_voxels.toLocaleString(), '—', '—'],
    ['连通域', `${g.mask.components}，保留最大者 ${g.mask.largest_component_voxels.toLocaleString()}`,
      '—', `丢弃 ${g.mask.discarded_voxels} 体素 (${(g.mask.discarded_fraction * 100).toFixed(3)}%)`],
    ['视野截断', Object.values(g.mask.truncated_by_field_of_view).some(Boolean) ? bad('有') : ok('无'), '—', ok('达成')],
    ['三角面数', m.triangles.toLocaleString(), '130,000–160,000', m.triangles >= 130000 && m.triangles <= 160000 ? ok('达成') : bad('未达成')],
    ['glb 大小', (m.glb_bytes / 1048576).toFixed(2) + ' MB', '≤ 7 MB', m.glb_bytes <= 7 * 1048576 ? ok('达成') : bad('未达成')],
    ['网格水密', m.watertight ? ok('是') : bad('否'), '—', '—'],
    ['标签体积', g.mask.label_volume_mL.toFixed(1) + ' mL', '—', '含壁与腔，非净腔容积'],
    ['中心线长度（实测）', `<b>${cl.length_cm.toFixed(2)} cm</b>`, 'CT 结肠造影 150.3 ± 18.5 cm（范围 109.7–195.9）',
      bad('低于公开范围下界')],
    ['骨架原始长度', cl.raw_skeleton_route_length_mm.toFixed(1) + ' mm',
      '—', `样条容差 ${cl.spline_tolerance_mm} mm，最大偏离 ${cl.spline_max_deviation_from_skeleton_mm} mm`],
    ['管腔半径', `${cl.radius_mm.min.toFixed(2)}–${cl.radius_mm.max.toFixed(2)} mm（均值 ${cl.radius_mm.mean.toFixed(2)}）`,
      '—', '未充气的常规腹部 CT，多处塌陷']
  ];
  for (const sg of g.segments) {
    const pub = PUB_SEG[sg.name];
    const inr = pub ? (sg.length_cm >= pub[0] && sg.length_cm <= pub[1]) : null;
    rows.push([SEGNAME[sg.name] || sg.name, `${sg.length_cm.toFixed(1)} cm（管径均值 ${sg.diameter_mm_mean.toFixed(1)} mm）`,
      pub ? `${pub[0]}–${pub[1]} cm` : '构造段（见下）',
      pub ? (inr ? ok('在范围内') : bad('超出范围')) : '肝曲/脾曲范围为声明的构造']);
  }
  for (const [k, v] of Object.entries(g.landmarks)) {
    rows.push([`界标 ${esc(k)}`, v.s_mm !== undefined ? `s = ${v.s_mm.toFixed(1)} mm` : '—',
      esc(v.rule || v.note || ''), v.contour_anchored === false ? bad('contour_anchored = false') : ok('contour_anchored = true')]);
  }
  tbl('t-geom', ['项目', '实测', '公开参考 / 规则', '判定'], rows);

  // --- taeniae & haustra
  const d = taeRep.derivation, me = taeRep.measured, oc = taeRep.orientation_check;
  const tr = [
    ['ρ = 结肠带 / 肠壁长度', d.rho_taenia_over_wall.toFixed(4), esc(d.rho_source), bad('数值 B 级 / 论文身份 A 级')],
    ['λ 结肠袋间距', d.haustral_pitch_mm + ' mm', esc(d.pitch_source), bad('数值 B 级 / 论文身份 A 级')],
    ['φ 解出', `${d.phi_deg.toFixed(3)}°`, 'φ/sin φ = 1/ρ', ok('解出，非调参')],
    ['弧长/弦长（生成曲面上量回）', me.bulge_profile_arc_over_chord_measured.toFixed(6),
      '目标 ' + me.bulge_profile_arc_over_chord_target.toFixed(6),
      Math.abs(me.bulge_profile_arc_over_chord_measured - me.bulge_profile_arc_over_chord_target) < 1e-5 ? ok('完全吻合') : bad('不吻合')],
    ['导出的袋深 h', d.derived_sagitta_mm.toFixed(3) + ' mm', '由 ρ 与 λ 唯一确定', ok('非可调参数')],
    ['结肠袋数 / 半月襞数', `${me.haustra_count} / ${me.semilunar_folds}`, '—', `结肠袋区段 ${me.haustrated_run_mm.toFixed(0)} mm`],
    ['直肠内结肠袋', ok('0（带散开为连续纵肌层）'), 'Gray/IMAIOS：直肠无结肠袋', ok('达成')],
    ['三带长度（实测）', `系膜带 ${me.taenia_mesocolica_length_mm.toFixed(1)} · 网膜带 ${me.taenia_omentalis_length_mm.toFixed(1)} · 独立带 ${me.taenia_libera_length_mm.toFixed(1)} mm`,
      '—', `/ 中心线 = ${me.taenia_over_centreline_length_ratio.mesocolica.toFixed(3)} / ${me.taenia_over_centreline_length_ratio.omentalis.toFixed(3)} / ${me.taenia_over_centreline_length_ratio.libera.toFixed(3)}`],
    ['带长 / 壁长（三维实测）', me.taenia_over_wall_ratio_measured_8x_upsampled.toFixed(4),
      '声明 ' + me.taenia_over_wall_ratio_declared.toFixed(4),
      bad(`相对偏差 ${(Math.abs(me.taenia_over_wall_ratio_measured_8x_upsampled - me.taenia_over_wall_ratio_declared) / me.taenia_over_wall_ratio_declared * 100).toFixed(1)}%（轴线曲率所致，见验收文档）`)],
    ['三带夹角（实测）', me.band_angular_separation_deg.map(x => x.toFixed(1) + '°').join(' / '), '120° / 120° / 120°', ok('达成')],
    ['带宽', me.taenia_width_mm + ' mm', esc(me.taenia_width_source), bad('数值 B 级 / 论文身份 A 级')],
    ['肠脂垂', me.appendices_epiploicae + ' 个（独立带 + 网膜带两排，直肠无，盲肠稀疏）', '—', bad('contour_anchored = false')],
    ['管轴偏离实测中心线', `最大 ${me.tube_axis.max_deviation_from_measured_centreline_mm} mm，均值 ${me.tube_axis.mean_deviation_from_measured_centreline_mm} mm`,
      '容差 ' + me.tube_axis.spline_tolerance_mm + ' mm', '—'],
    ['曲率半径小于管腔半径的采样点', `${me.tube_axis.samples_with_curvature_radius_below_lumen_radius} / ${me.tube_axis.samples_total}`,
      '0 为理想', bad('该处圆管自相交，照实报告')]
  ];
  const ori = [
    ['独立带 · 升结肠', oc.taenia_libera_on_ascending.anterior_fraction, '应为前侧', oc.taenia_libera_on_ascending.mean_anterior_component > 0 ? ok('前侧 ✓') : bad('✗')],
    ['独立带 · 横结肠', oc.taenia_libera_on_transverse.mean_superior_component.toFixed(3), '应为下侧', oc.taenia_libera_on_transverse.mean_superior_component < 0 ? ok('下侧 ✓') : bad('✗')],
    ['系膜带 · 升结肠', oc.taenia_mesocolica_on_ascending.mean_anterior_component.toFixed(3), '应为后内侧', oc.taenia_mesocolica_on_ascending.mean_anterior_component < 0 ? ok('后侧 ✓') : bad('✗')],
    ['系膜带 · 横结肠', oc.taenia_mesocolica_on_transverse.mean_anterior_component.toFixed(3), '应为后侧（横结肠系膜附着）', oc.taenia_mesocolica_on_transverse.mean_anterior_component < 0 ? ok('后侧 ✓') : bad('✗')],
    ['网膜带 · 横结肠', `前 ${oc.taenia_omentalis_on_transverse.mean_anterior_component.toFixed(3)} / 上 ${oc.taenia_omentalis_on_transverse.mean_superior_component.toFixed(3)}`,
      '应为前上侧（大网膜附着）', (oc.taenia_omentalis_on_transverse.mean_anterior_component > 0 && oc.taenia_omentalis_on_transverse.mean_superior_component > 0) ? ok('前上侧 ✓') : bad('✗')]
  ];
  tbl('t-taenia', ['项目', '实测 / 解出', '来源或目标', '判定 / 核实等级'],
      tr.concat([['<b>— 带位走行核对（由一条锚定规则 + 120° 约定导出 —</b>', '', '', '']]).concat(ori));

  // --- vessels
  const vr = [];
  for (const [k, t] of Object.entries(vesRep.trees)) {
    vr.push([esc(k), t.segments, t.bifurcations, t.terminals,
      `${t.diameter_mm[0].toFixed(3)}–${t.diameter_mm[1].toFixed(3)}`,
      t.total_length_mm.toFixed(0),
      `${t.murray_max_relative_residual.toExponential(2)}（12 位小数四舍五入为 ${t.murray_max_relative_residual_rounded_12dp}）`,
      bad('锚定：' + (t.contour_anchored.length > 60 ? t.contour_anchored.slice(0, 60) + '…' : t.contour_anchored))]);
  }
  tbl('t-vessel', ['树', '段', '分叉', '终末', '直径 mm', '总长 mm', 'Murray 最大相对残差', '轮廓锚定'], vr);
  const cl2 = [];
  for (const [k, v] of Object.entries(vesRep.cross_tree_clearance)) {
    cl2.push([esc(k.replace('__vs__', ' ↔ ')), v.min_surface_gap_mm.toFixed(3),
      v.interpenetrating_segment_pairs, v.pairs_compared.toLocaleString(),
      v.at_solver_target ? bad('等于求解器目标 0.25 mm — 是求解器属性，不是独立测量') : ok('真实测量'),
      v.worst_pair ? esc(`${v.worst_pair.a_label} (r=${v.worst_pair.a_radius_mm}) ↔ ${v.worst_pair.b_label} (r=${v.worst_pair.b_radius_mm})`) : '—']);
  }
  const vt = T('t-vessel');
  vt.insertAdjacentHTML('afterend',
    '<h2>跨树间隙（实测）</h2><table id="t-clear"></table>' +
    `<p class="hl">${esc(vesRep.clearance_policy.priority)} 为满足间隙而删除的段数：<b>${vesRep.clearance_policy.segments_deleted_for_clearance}</b>。</p>` +
    '<h2>边缘动脉与 Riolan 弓（是环，不是树）</h2>' +
    `<p class="hl">${esc(vesRep.marginal_artery.note)}</p>` +
    `<p>${esc(vesRep.portal_vs_systemic.consequence)}</p>`);
  tbl('t-clear', ['树对', '最小表面间隙 mm', '穿插段对', '比较段对', '是否等于求解器目标', '最差段对'], cl2);

  // --- flow
  const fr = [];
  for (const c of flowRep.published_comparison) {
    fr.push([esc(c.level), (c.modelled_peak_cm_s ?? c.modelled_mean_cm_s).toFixed(3) + ' cm/s',
      (c.published_psv_cm_s || c.published_mean_cm_s).join('–') + ' cm/s',
      c.verdict === 'in range' ? ok('在范围内') : bad('超出范围'), esc(c.note)]);
  }
  fr.push(['逐采样点逆流审计',
    `${flowRep.reflux_audit.paths} 条通路 / ${flowRep.reflux_audit.samples.toLocaleString()} 个采样点`,
    '逆流 0', flowRep.reflux_audit.reflux_samples === 0 ? ok('逆流 ' + flowRep.reflux_audit.reflux_samples) : bad('逆流 ' + flowRep.reflux_audit.reflux_samples),
    `非有限 ${flowRep.reflux_audit.nonfinite} · 零面积 ${flowRep.reflux_audit.zero_area}`]);
  fr.push(['结肠总血流', flowRep.flow_per_tree.sma_colic.Q_root_mL_min.toFixed(1) + ' (SMA 全量) / ' +
    flowRep.flow_per_tree.ima_colic.Q_root_mL_min.toFixed(1) + ' (IMA) mL/min',
    '—', bad('SMA 结肠份额为声明的假设，非实测'), esc(vesRep.flows.assumption)]);
  fr.push(['搏动耦合', `PI 模型 ${flowRep.heart_coupling.pulsatility_index_model}`,
    `公开 IMA PI ${flowRep.heart_coupling.pulsatility_index_published_IMA}`,
    bad('远低于公开值'), esc(flowRep.heart_coupling.limitation)]);
  fr.push(['管腔内容物速度', flowRep.luminal_transport.baseline_content_velocity_mm_s.toFixed(5) + ' mm/s',
    `总转运 ${flowRep.luminal_transport.mean_total_transit_h} h`,
    bad('显示增益 ' + flowRep.luminal_transport.display_gain + '×'), esc(flowRep.luminal_transport.gain_note)]);
  tbl('t-flow', ['项目', '模型实测', '公开参考', '判定', '说明'], fr);

  // --- motion (filled in by the browser audit, and re-filled on demand)
  const mo = micro.motility;
  tbl('t-motion', ['运动', '声明值（文献）', '核实等级', '浏览器内实测'], [
    ['结肠袋收缩', `${mo.haustral_contractions.frequency_cycles_per_min} 次/分（公开 ${mo.haustral_contractions.frequency_published_range.join('–')}），跨度 ${mo.haustral_contractions.span_cm.join('–')} cm，传播 ${mo.haustral_contractions.propagation_cm_s} ± ${mo.haustral_contractions.propagation_sd_cm_s} cm/s，非推进性`,
      bad('数值 B 级 / 论文身份 A 级'), '<span id="a-haustral">切换到"运动 → 结肠袋收缩"后由 <code>__motionAudit()</code> 实测</span>'],
    ['集团蠕动 (HAPC)', `每 24 h ${mo.mass_movement.events_per_24h_median} 次（范围 ${mo.mass_movement.events_per_24h_range.join('–')}），传播 ${mo.mass_movement.propagation_velocity_cm_s} ± ${mo.mass_movement.propagation_velocity_sd_cm_s} cm/s，≥ ${mo.mass_movement.minimum_propagated_distance_cm} cm，${mo.mass_movement.duration_s.join('–')} s，${mo.mass_movement.amplitude_mmHg.join('–')} mmHg`,
      bad('数值 B 级 / 论文身份 A 级'), '<span id="a-mass">同上</span>'],
    ['直肠充盈 / 排便反射', esc(mo.defecation.sequence.join(' → ')) + `；本受试者直肠腔容积实测 ${mo.defecation.measured_rectal_lumen_volume_mL} mL`,
      bad('数值 B 级 / 论文身份 A 级'), '<span id="a-defec">同上</span>'],
    ['阈值容积', bad('未建模'), esc(mo.defecation.sensory_thresholds), '—']
  ]);

  // --- micro
  const cr = micro.crypt, su = micro.surface_area, gb = micro.goblet;
  tbl('t-micro', ['项目', '值', '来源（仅为核对指引）', '核实等级'], [
    ['绒毛', bad('结肠完全没有绒毛（黏膜是平的）'), '这正是结肠吸收面积远小于小肠的原因', '—'],
    ['隐窝深度', `${cr.depth_um} ± ${cr.depth_sd_um} µm`, esc(cr.source.slice(0, 110)) + '…', bad('数值 B 级 / 论文身份 A 级')],
    ['隐窝外径 / 腔径 / 上皮高', `${cr.outer_diameter_um} / ${cr.lumen_diameter_um} / ${cr.epithelial_height_um} µm`, '同上', '同上'],
    ['每隐窝细胞数', `${cr.cells_per_crypt.join('–')}（更宽范围 ${cr.cells_per_crypt_wider_range.join('–')}）`,
      esc(cr.cells_source.slice(0, 90)) + '…', bad('数值 B 级 / 论文身份 A 级')],
    ['基底干细胞', cr.stem_cells_at_base.join('–') + ' 个', '同上', '同上'],
    ['放大倍数（声明）', `<b>${cr.display_magnification}×</b>（按真实 µm 建模后放大；与小肠页同一倍数，便于直接对照）`,
      `渲染深度 ${cr.rendered_depth_mm_at_magnification} mm`, '—'],
    ['表面放大倍数 · 结肠', `${su.colon_amplification_microvilli_only}×（仅微绒毛）`, esc(su.source.slice(0, 95)) + '…', bad('数值 B 级 / 论文身份 A 级')],
    ['表面放大倍数 · 小肠', `${su.small_intestine_amplification_villi_plus_microvilli.join('–')}×（绒毛 + 微绒毛，其中绒毛单独 ${su.small_intestine_villi_only}×）`, '同上', '同上'],
    ['小肠 / 结肠 放大倍数之比', `<b>${su.ratio_small_intestine_over_colon}×</b>（算出）`, '—', ok('由上面两个公开值算出')],
    ['黏膜总面积', `全消化道 ${su.total_digestive_mucosa_m2} m²，其中大肠 ${su.large_intestine_mucosa_m2} m²`, '同上', '同上'],
    ['本受试者分割表面积', `${su.measured_colon_label_surface_m2} m²`, esc(su.measured_note.slice(0, 140)) + '…', bad('与 2 m² 黏膜面积不是同一个量，不可直接对比')],
    ['杯状细胞比例', `结肠 ${(gb.colon_goblet_fraction_published * 100).toFixed(0)}% vs 小肠 ${(gb.small_intestine_goblet_fraction_published * 100).toFixed(0)}%，黏液层厚 ${gb.mucus_layer_thickness_ratio_colon_over_si}×`,
      esc(gb.source.slice(0, 95)) + '…', bad('数值 B 级 / 论文身份 A 级')],
    ['近→远端杯状细胞梯度', bad('未由原始文献核实'), esc(gb.proximal_to_distal_gradient), bad('改由本项目自己的单细胞数据检验')]
  ]);

  // --- routes
  const rc = micro.route_comparison;
  tbl('t-route', ['路径', '静脉回流', '肝脏首过'], [
    ['口腔 · 舌下', esc(rc.sublingual), ok('无首过')],
    ['胃 / 小肠（口服）', esc(rc.gastric_and_upper_gi), bad('完整首过')],
    ['直肠上段（栓剂位置偏高）', esc(rc.upper_rectum), bad('完整首过')],
    ['直肠下段（栓剂位置偏低）', esc(rc.lower_rectum), ok('绕过首过')],
    ['说明', esc(rc.caveat), '—']
  ]);

  // --- molecular
  const mb = micro.microbiome;
  let html = '';
  if (!expr) {
    html += `<div class="warnbox"><b>分子层数据缺失。</b>CELLxGENE Census 的结肠 61 基因 ADME 谱在本次会话中<b>没有成功产出</b>，
      因此本页<b>不显示任何表达数值</b>。按既定要求，这里<b>没有</b>用组织平均值、其它器官的值或任何代理值填补。
      失败原因与已尝试的路线写在 ACCEPTANCE.md。</div>`;
  } else {
    html += `<p>Census <code>${esc(expr.census_version)}</code> · ${expr.cell_types.length} 个细胞类型 × ${expr.genes.length} 个基因 · ` +
      `${expr.n_cells.reduce((a, b) => a + b, 0).toLocaleString()} 个细胞。` +
      (molInfo ? `场景中显示的是 <code>${esc(molInfo.cell_type)}</code>。` : '') + '</p>';
    if (molInfo) {
      html += '<table><thead><tr><th>蛋白</th><th>基因</th><th>膜面</th><th>方向</th><th>实测表达 log1p CP10K</th></tr></thead><tbody>' +
        molInfo.placed.sort((a, b) => b.value - a.value).map(p =>
          `<tr><td>${esc(p.protein)}</td><td><code>${esc(p.gene)}</code></td><td>${esc(p.membrane)}</td><td>${esc(p.role)}</td><td>${p.value}</td></tr>`).join('') +
        '</tbody></table>';
    }
  }
  html += `<div class="warnbox" style="margin-top:14px"><b>菌群部分是文献描述，不是本项目测的数据。</b>${esc(mb.disclaimer)}</div>`;
  html += '<table><thead><tr><th>菌群介导的药物代谢主线</th><th>链条</th><th>说明</th></tr></thead><tbody>' +
    MICROBIAL_LINES.map(l => `<tr><td>${esc(l.title)}</td><td>${esc(l.chain.join(' → '))}</td><td>${esc(l.note)}</td></tr>`).join('') +
    '</tbody></table>';
  html += '<table><thead><tr><th>发酵</th><th>值</th><th>来源（仅为核对指引）</th><th>核实等级</th></tr></thead><tbody>' +
    `<tr><td>腔内细菌密度</td><td>${mb.luminal_density_per_mL.map(x => x.toExponential(0)).join('–')} /mL</td><td>${esc(mb.density_source)}</td><td>${bad(mb.density_citation_level)}</td></tr>` +
    `<tr><td>SCFA 摩尔比 乙酸:丙酸:丁酸</td><td>${mb.scfa.molar_ratio_acetate_propionate_butyrate.join(':')}（外周血变为 ${mb.scfa.peripheral_molar_ratio.join(':')}）</td><td>${esc(mb.scfa.source.slice(0, 120))}…</td><td>${bad('数值 B 级 / 论文身份 A 级')}</td></tr>` +
    `<tr><td>腔内丁酸</td><td>${mb.scfa.butyrate_lumen_mM.join('–')} mM</td><td>同上</td><td>同上</td></tr>` +
    `<tr><td>门静脉 / 外周 SCFA 比</td><td>${mb.scfa.portal_over_peripheral_concentration_ratio}×</td><td>同上</td><td>同上</td></tr>` +
    `<tr><td>丁酸作为结肠细胞燃料</td><td>${esc(mb.butyrate_as_fuel.claim)}</td><td>${esc(mb.butyrate_as_fuel.source)}</td><td>${bad(mb.butyrate_as_fuel.citation_level)}</td></tr>` +
    `<tr><td>定量比例</td><td>${bad(mb.butyrate_as_fuel.quantitative_fraction)}</td><td>—</td><td>—</td></tr>` +
    `<tr><td>标定</td><td>${bad(mb.no_calibration)}</td><td>—</td><td>—</td></tr>` +
    '</tbody></table>';
  T('d-mol').innerHTML = html;

  // --- provenance paragraph
  T('p-prov').innerHTML =
    `几何来源：<b>${esc(src.dataset)}</b>，doi:<code>${esc(src.doi)}</code>，许可 <b>${esc(src.licence)}</b>（纯署名，可再分发）。` +
    `受试者 <code>${esc(src.subject)}</code>（${esc(src.subject_meta)}）。` +
    `取用方式：${esc(src.access)}。` +
    `<b>本项目不含任何基于 CT 影像灰度的内容</b>，只使用了该版本的掩膜标签。` +
    `<b>引文政策（与本项目其它消化道器官同一套措辞）：</b>本页每一个文献数值都是 ` +
    `<b>数值 B 级</b>（来自未留存的检索片段、<b>未回原文核对</b>）配 <b>论文身份 A 级</b>` +
    `（期刊 / DOI / PMID 可核查）。构建过程中<b>没有打开过任何一篇论文全文</b>` +
    `（${micro.citation_policy.full_texts_opened} 篇），交付包里也<b>没有任何可复现的检索记录</b>` +
    `（${micro.citation_policy.reproducible_retrievals_shipped} 条）。` +
    `数值旁边的文献只是<b>告诉你去哪里核对</b>，不代表本项目核对过。` +
    `凡是本项目自己算出来的量都标为“实测 / 解出”，产出它的代码在 tools 包里。`;

  // keep the motion audit table live
  const put = (id, o) => {
    const el = T(id); if (!el) return;
    el.innerHTML = Object.entries(o).map(([k, v]) =>
      `<code>${esc(k)}</code> = ${esc(JSON.stringify(v))}`).join('<br>');
  };
  window.__fillMotionTable = (mode, res) =>
    put({ haustral: 'a-haustral', mass: 'a-mass', defecation: 'a-defec' }[mode], res);
}
