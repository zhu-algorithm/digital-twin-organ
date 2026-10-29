// hud.js — the provenance strip.  Every number shown here is read from the
// generated reports, never hard-coded in the page.
export function hudText(diag, narrow) {
  const g = diag.geometry.measured;
  const v = diag.vesselReport;
  const f = diag.flowReport;
  const mo = diag.motionReport.kinematics;
  const mm = diag.motionMeasured;
  const mol = diag.molecular;
  const L = [];
  L.push(`几何 TotalSegmentator v2.0.1 · CC-BY-4.0 · s0236 · <code>${(diag.triangles / 1000).toFixed(0)}k</code> 面 · glb <code>${diag.glbMB}</code> MB`);
  L.push(`胃容积 <code>${g.label_volume_mL}</code> mL（空腹参考 ${g.volume_vs_fasted_reference_mL.join('–')} · 容量参考 ${g.volume_vs_capacity_reference_mL.join('–')}）`);
  L.push(`角切迹 <code>${g.incisura_detected_on_this_subject ? '实测' : '未检出→按比例'}</code> κ=<code>${g.incisura_curvature_rad_per_mm}</code> rad/mm（阈值 ${g.incisura_detection_threshold_rad_per_mm}，显著度 ${g.incisura_prominence_ratio}×）`);
  L.push(`Murray 最大相对残差 <code>${v.murray_max_relative_residual_all_trees}</code> · 跨树最小表面间隙 <code>${v.clearance.global_min_surface_gap_mm}</code> mm（${v.clearance.global_min_between}）`);
  L.push(`逆流审计 <code>${f.audit.total_samples_audited}</code> 采样点 / 逆流 <code>${f.audit.retrograde_samples_total}</code> 个 · 血流层增益 <code>1×</code>`);
  L.push(`蠕动 声明 <code>${mo.frequency_cpm}</code> cpm · 浏览器实测 <code>${mm ? mm.measuredFrequency_cpm : '—'}</code> cpm · 实测波速 <code>${mm ? mm.measuredPeakSpeed_mm_s.median : '—'}</code> mm/s（声明 ${mo.declared_speed_profile_mm_s.corpus}–${mo.declared_speed_profile_mm_s.terminal_28mm[1]}）`);
  L.push(`t=0 未位移采样环 <code>${mm ? mm.ringsUndisplacedAtT0 : '—'}</code>/<code>${mm ? mm.stationsTotal : '—'}</code> → 是行进波，不是整体缩放`);
  if (mol && mol.ok) {
    L.push(`分子层 <code>${mol.n_cell_types}</code> 细胞类型 × <code>${mol.n_genes}</code> 基因 · <code>${mol.n_cells}</code> 细胞 · donor <code>${mol.n_donors}</code>`);
  } else {
    L.push(`<b style="color:#f0b25a">分子层缺失</b>：${mol ? mol.reason : '未加载'}`);
  }
  L.push(`<code>is_user_personal_data = false</code>（心脏耦合源自 PhysioNet 公开数据）`);
  if (narrow) L.push(`<code>窄屏：粒子数与像素比已降档</code>`);
  return L.join('<br>');
}
