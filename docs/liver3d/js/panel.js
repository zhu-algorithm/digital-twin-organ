/**
 * panel.js — the per-cell readout, and the resolution label on every row.
 *
 * The whole point of this panel is that the three kinds of number in it do NOT
 * have the same resolution, and the reader has to be able to see which is which
 * at a glance. The heart page established the wording; it is reused unchanged:
 *
 *   本细胞实测   this cell's own measurement. Only the 61 gene values qualify:
 *                they are the raw counts recorded for this one soma_joinid in
 *                the CELLxGENE Census, with its own raw_sum.
 *   本细胞类型   a property of the cell TYPE this cell belongs to, not of this
 *                cell. All antibody-IHC protein evidence is here: HPA reports a
 *                staining level per cell type, never per cell.
 *   整块组织     a property of a liver homogenate — every cell type mixed. The
 *                shotgun iBAQ values and the absolute pmol/mg values are here.
 *
 * No row ever crosses a tier. There is no expression-to-protein coefficient and
 * no abundance-to-activity coefficient anywhere in this file, because this
 * project has not calibrated one; where a tier has no value for a gene the panel
 * says 无数据 rather than borrowing from the tier above or below.
 */

const TIER = {
  cell: { label: '本细胞实测', cls: 'tier-cell',
          hint: '这个 soma_joinid 自己的测序读数' },
  celltype: { label: '本细胞类型', cls: 'tier-type',
              hint: '该细胞类型的属性，不是这个细胞的属性' },
  tissue: { label: '整块组织', cls: 'tier-tissue',
            hint: '肝匀浆，所有细胞类型混在一起，无细胞分辨率' },
  othercell: { label: '他人单细胞', cls: 'tier-other',
               hint: '是单个肝细胞的蛋白测量，但那是另一队列的 792 个细胞，'
                   + '不是你点的这个 Census 细胞' },
  none: { label: '无数据', cls: 'tier-none', hint: '该层没有这个基因的数据' }
};

const CLASS_ZH = {
  phase_I: 'I 相代谢', phase_II: 'II 相代谢',
  uptake_transporter: '摄取转运体', efflux_transporter: '外排转运体',
  nuclear_receptor: '核受体', housekeeping_reference: '管家参照'
};

export class CellPanel {
  /**
   * @param {HTMLElement} host    container the panel is appended to
   * @param {object} cellData     parsed liver_cells.json
   * @param {object} protein      parsed liver_protein_evidence.json
   * @param {object} expression   parsed liver_expression.json (cell-type means)
   */
  constructor (host, cellData, protein, expression) {
    this.genes = cellData.genes;
    this.geneClass = cellData.gene_class;
    this.protein = protein || { genes: {}, tiers: {} };
    this.expression = expression;
    this.sampling = cellData.sampling;
    this.seed = cellData.seed;
    this.normalization = cellData.normalization;
    this.el = document.createElement('div');
    this.el.id = 'cellpanel';
    this.el.style.display = 'none';
    host.appendChild(this.el);
    this.sortBy = 'gene';
    this._current = null;
  }

  hide () { this.el.style.display = 'none'; this._current = null; }

  /** CP10K-log1p of one cell from its own counts and its own raw_sum. */
  static cp10k (count, rawSum) {
    if (!rawSum) return 0;
    return Math.log1p(count / rawSum * 1e4);
  }

  /**
   * @param {object} rec  {cls, cell, zone01} from CellTier.intersect
   */
  show (rec) {
    if (!rec || !rec.cell) {
      this.el.innerHTML = '<div class="cp-head">没有绑定到真实细胞</div>'
        + '<div class="cp-note">这个渲染实例所对应的 Census 细胞类型在抽样文件中没有可用细胞，'
        + '因此它没有任何逐细胞数据。它不会显示任何基因值——不会用细胞类型均值代替。</div>';
      this.el.style.display = 'block';
      return;
    }
    this._current = rec;
    const c = rec.cell;
    const counts = new Int32Array(this.genes.length);
    for (const [i, v] of c.counts) counts[i] = v;
    const nnz = c.counts.length;

    const pg = this.protein.genes || {};
    const rows = this.genes.map((g, i) => {
      const raw = counts[i];
      const norm = CellPanel.cp10k(raw, c.raw_sum);
      const p = pg[g] || {};
      const ihc = p.ihc || {};
      const ihcKey = Object.keys(ihc).find(k => matchCellType(k, rec.cls)) ||
                     Object.keys(ihc)[0];
      return {
        gene: g, cls: this.geneClass[i], raw, norm,
        ihcCellType: ihcKey || null,
        ihc: ihcKey ? ihc[ihcKey] : null,
        ibaq: p.ms_detected ? p.ibaq_ppm_of_total : null,
        abs: (p.absolute && p.absolute.length) ? p.absolute[0] : null,
        sc: p.sc_protein || null
      };
    });
    rows.sort((a, b) => this.sortBy === 'expr' ? b.norm - a.norm
                                              : a.gene.localeCompare(b.gene));

    const tierTag = (t) => `<span class="tier ${TIER[t].cls}" title="${TIER[t].hint}">${TIER[t].label}</span>`;
    const esc = s => String(s === null || s === undefined ? '' : s)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

    const body = rows.map(r => `<tr>
      <td class="g">${r.gene}<span class="gc">${CLASS_ZH[r.cls] || r.cls}</span></td>
      <td class="n">${r.raw}</td>
      <td class="n">${r.norm.toFixed(3)}</td>
      <td class="ihc">${r.ihc
        ? `${esc(r.ihc.level || '—')}<span class="sub">${esc(r.ihcCellType)}${
            r.ihc.intensity ? ' · ' + esc(r.ihc.intensity) : ''}${
            r.ihc.quantity ? ' · ' + esc(r.ihc.quantity) : ''}</span>`
        : '<span class="na">无</span>'}</td>
      <td class="n">${r.ibaq !== null ? r.ibaq.toFixed(2) : '<span class="na">未检出</span>'}</td>
      <td class="n">${r.abs ? `${esc(r.abs.mean)}<span class="sub">${esc(r.abs.unit)}</span>`
                            : '<span class="na">无</span>'}</td>
      <td class="n">${r.sc && r.sc.n_cells_quantified
        ? `${r.sc.pct.toFixed(0)}%<span class="sub">${r.sc.n_cells_quantified}/${r.sc.n_cells_total}</span>`
        : '<span class="na">未定量</span>'}</td>
    </tr>`).join('');

    const ident = this.protein.liver_experiment_identification || {};
    const sc = (this.protein.tiers && this.protein.tiers.single_cell) || {};

    this.el.innerHTML = `
      <div class="cp-head">
        <b>${esc(c.cell_type)}</b>
        <button id="cp-close" title="关闭">×</button>
      </div>
      <div class="cp-id">
        <div><span>soma_joinid</span><code>${c.soma_joinid}</code>${tierTag('cell')}</div>
        <div><span>donor_id</span><code>${esc(c.donor_id)}</code>${tierTag('cell')}</div>
        <div><span>dataset_id</span><code>${esc(c.dataset_id).slice(0, 18)}…</code>${tierTag('cell')}</div>
        <div><span>assay</span><code>${esc(c.assay)}</code>${tierTag('cell')}</div>
        <div><span>tissue</span><code>${esc(c.tissue)}</code>${tierTag('cell')}</div>
        ${c.disease ? `<div><span>disease</span><code>${esc(c.disease)}</code>${tierTag('cell')}</div>` : ''}
        ${c.sex ? `<div><span>sex</span><code>${esc(c.sex)}</code>${tierTag('cell')}</div>` : ''}
        <div><span>raw_sum</span><code>${c.raw_sum}</code>${tierTag('cell')}</div>
        <div><span>61 基因中非零</span><code>${nnz}</code>${tierTag('cell')}</div>
        ${rec.zone01 !== undefined
          ? `<div><span>渲染位置</span><code>zone ${rec.zone01.toFixed(2)}</code>
             <span class="tier tier-none" title="这是它在本页渲染中的位置，不是 Census 记录的位置">渲染坐标</span></div>`
          : ''}
      </div>
      <div class="cp-banner">
        这一行行数字<b>不是"某类细胞的平均值"</b>。counts 与 CP10K-log1p 两列是
        <b>这一个</b> <code>soma_joinid</code> 自己的测序读数，
        CP10K-log1p = <code>log1p(count / raw_sum × 1e4)</code> 由它自己的 counts 与 raw_sum 当场算出。
        右边三列<b>分辨率更粗</b>，逐列标注。
      </div>
      <div class="cp-sort">
        排序 <button class="${this.sortBy === 'gene' ? 'on' : ''}" data-sort="gene">基因名</button>
        <button class="${this.sortBy === 'expr' ? 'on' : ''}" data-sort="expr">本细胞表达量</button>
      </div>
      <table class="cp-tab">
        <thead><tr>
          <th>基因</th>
          <th>raw counts<br>${tierTag('cell')}</th>
          <th>CP10K-log1p<br>${tierTag('cell')}</th>
          <th>HPA IHC 染色等级<br>${tierTag('celltype')}</th>
          <th>iBAQ ppm<br>${tierTag('tissue')}</th>
          <th>绝对定量<br>${tierTag('tissue')}</th>
          <th>单细胞蛋白检出率<br>${tierTag('othercell')}</th>
        </tr></thead>
        <tbody>${body}</tbody>
      </table>
      <div class="cp-note">
        <b>每一列是什么、不是什么。</b><br>
        · <b>raw counts / CP10K-log1p</b> — 本细胞实测。CELLxGENE Census
          <code>${esc(this.protein.census_version || '')}</code>，抽样方式：${esc(this.sampling)}，
          随机种子 <code>${this.seed}</code>。<br>
        · <b>HPA IHC 染色等级</b> — 本细胞类型级。抗体免疫组化给的是
          <b>有序类别</b>（high / medium / low / not detected）与染色强度、阳性细胞比例、
          亚细胞定位，<b>不是丰度</b>，没有每细胞拷贝数。它描述的是"肝细胞"或"胆管细胞"
          这一整类，不是你点的这一个细胞。<br>
        · <b>iBAQ ppm</b> — 整块组织级。29 种健康人体组织的 label-free 鸟枪蛋白质组
          （${esc((this.protein.tiers.tissue_homogenate || {}).source || '')}）。
          iBAQ 在同一样品内正比于摩尔量，可用于蛋白之间的相对比较，
          <b>完全没有细胞分辨率</b>：肝匀浆把肝细胞和所有非实质细胞混在一起。
          该数据集未附"实验编号→组织"对照表，肝脏实验是
          <code>${esc(ident.experiment_assigned || '')}</code>，
          由${esc(ident.panel_genes_peaking_in_assigned_experiment || '')} 个肝细胞特异（非分泌型）标志物峰值一致判定，
          与次高组织相差 <b>${esc(ident.separation_fold || '')}×</b>——这是<b>推断</b>，不是作者给的标签。<br>
        · <b>绝对定量</b> — 整块组织级，但<b>是绝对量</b>（pmol/mg 蛋白）。
          来自本项目 <code>protein_abundance_anchors.csv</code> 中
          <code>measure = absolute_protein_abundance</code> 的 36 行（覆盖 61 个基因中的 23 个），
          基质为人肝细胞、人肝微粒体或肝总膜的<b>混合制备物</b>。
          分母是"每毫克蛋白"，<b>不是"每个细胞"</b>——本项目没有"每细胞蛋白量"换算系数，
          所以它无法变成"这个肝细胞里有多少分子"。<br>
        · <b>单细胞蛋白检出率</b> — 他人单细胞级。人肝单细胞蛋白质组<b>确实已经存在</b>，
          任务简报里"基本不存在"的说法经 PRIDE 检索后<b>被证伪</b>：
          ${esc((sc.projects || []).map(p => p.accession).join(' / '))}，
          ${esc(((sc.projects || [])[0] || {}).reference || '')}。
          本页从中读出的是实测覆盖度：<b>${esc(sc.n_single_cells)}</b> 个激光显微切割的单个肝细胞，
          每个细胞中位定量 <b>${esc((sc.proteins_per_single_cell || {}).median)}</b> 个蛋白，
          61 个 ADME 基因中 <b>${esc(sc.adme_panel_genes_quantified)}</b> 个在至少一个单细胞里被定量到。
          <b>但这一列不是你点的这个细胞的蛋白量</b>：那是另一个队列的另一批细胞。
          本项目里<b>没有任何一个细胞同时测了转录组和蛋白质组</b>。
          该数据集附了 18 个按个体分的分区打分文件（共 792 行，与 792 个样品列数目一致），
          但<b>投稿没有给出打分行与 mzML 样品列的对应键</b>，
          按行序硬配会是一个<b>被发明出来的映射</b>，所以本页<b>没有</b>由它导出任何门脉／中央分区蛋白值。<br>
        <b>本项目没有任何"表达量→蛋白丰度"或"丰度→活性"的标定系数</b>，
        因此这张表<b>不会</b>把左边的基因读数换算成右边的蛋白量，也不会反过来。
      </div>`;

    const close = this.el.querySelector('#cp-close');
    if (close) close.onclick = () => this.hide();
    for (const b of this.el.querySelectorAll('[data-sort]')) {
      b.onclick = () => { this.sortBy = b.dataset.sort; this.show(this._current); };
    }
    this.el.style.display = 'block';
    this.el.scrollTop = 0;
  }
}

/** Does an HPA liver cell-type label correspond to the rendered class? */
function matchCellType (hpaLabel, cls) {
  const l = (hpaLabel || '').toLowerCase();
  if (cls && cls.startsWith('hepatocyte')) return l.includes('hepatocyte');
  if (cls === 'cholangiocyte') return l.includes('bile duct');
  return false;
}

export const PANEL_CSS = `
#cellpanel{position:absolute;right:10px;top:10px;bottom:10px;width:min(96vw,430px);
  background:rgba(8,11,18,.94);border:1px solid rgba(120,140,170,.3);border-radius:10px;
  z-index:30;overflow:auto;padding:0 0 12px;font-size:12px;backdrop-filter:blur(8px)}
#cellpanel .cp-head{position:sticky;top:0;background:rgba(10,14,22,.97);padding:10px 12px;
  border-bottom:1px solid rgba(120,140,170,.26);font-size:13.5px;display:flex;
  justify-content:space-between;align-items:center;z-index:2}
#cellpanel .cp-head button{background:none;border:none;color:#9fb4cf;font-size:19px;
  cursor:pointer;line-height:1;padding:0 4px;min-height:0}
#cellpanel .cp-id{padding:8px 12px;display:grid;grid-template-columns:1fr;gap:3px}
#cellpanel .cp-id div{display:flex;gap:6px;align-items:baseline;flex-wrap:wrap}
#cellpanel .cp-id span:first-child{color:#8fa0b8;min-width:88px;font-size:11px}
#cellpanel code{font-family:ui-monospace,Menlo,monospace;font-size:11px;color:#dce6f3}
#cellpanel .cp-banner{margin:6px 12px 8px;padding:8px 10px;border-radius:7px;
  background:#152033;border:1px solid #2d4668;color:#cfe0f5;
  font-size:11.5px;line-height:1.6}
#cellpanel .cp-sort{padding:0 12px 6px;color:#8fa0b8;font-size:11px}
#cellpanel .cp-sort button{font-size:11px;padding:3px 7px;min-height:0}
#cellpanel table.cp-tab{width:100%;border-collapse:collapse;font-size:11px}
#cellpanel .cp-tab th{position:sticky;top:38px;background:rgba(10,14,22,.97);
  text-align:left;padding:5px 6px;border-bottom:1px solid rgba(120,140,170,.3);
  font-weight:600;color:#9fb4cf;font-size:10px;line-height:1.5;vertical-align:bottom}
#cellpanel .cp-tab td{padding:4px 6px;border-bottom:1px solid rgba(120,140,170,.12);
  vertical-align:top;font-variant-numeric:tabular-nums}
#cellpanel td.g{color:#e9eef6;font-weight:600;white-space:nowrap}
#cellpanel .gc{display:block;color:#7f8fa6;font-weight:400;font-size:9.5px}
#cellpanel td.n{text-align:right;color:#dce6f3;white-space:nowrap}
#cellpanel td.ihc{color:#dce6f3}
#cellpanel .sub{display:block;color:#7f8fa6;font-size:9.5px}
#cellpanel .na{color:#6b7787}
#cellpanel .tier{display:inline-block;padding:1px 5px;border-radius:4px;font-size:9px;
  white-space:nowrap;font-weight:600;letter-spacing:.02em}
#cellpanel .tier-cell{background:#13402c;color:#7fe0ab;border:1px solid #2e7a53}
#cellpanel .tier-type{background:#3d3312;color:#e5c45f;border:1px solid #7d6a24}
#cellpanel .tier-tissue{background:#3b1f2c;color:#e59aad;border:1px solid #7d3b53}
#cellpanel .tier-other{background:#1d2f45;color:#86c0e8;border:1px solid #3a6a93}
#cellpanel .tier-none{background:#232a36;color:#8d9aad;border:1px solid #3c4659}
#cellpanel .cp-note{margin:10px 12px 0;padding-top:9px;border-top:1px solid rgba(120,140,170,.24);
  color:#a9b7c9;font-size:10.5px;line-height:1.75}
@media(max-width:620px){
  #cellpanel{left:10px;right:10px;width:auto;top:auto;height:62vh;bottom:10px}
  #cellpanel .cp-tab th{top:36px}
}
`;
