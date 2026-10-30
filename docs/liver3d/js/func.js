/**
 * func.js — the dynamic liver-function channels, and the grade of each one.
 *
 * Every channel is driven by the SAME geometric zone coordinate as the rest of
 * the scene: zone01 = d(nearest portal terminal) / (d_portal + d_hepatic-vein
 * terminal), computed from this subject's own two vascular trees. What differs
 * between channels is where the NUMBERS on that axis come from, and that is the
 * only thing the reader needs to be able to check. So each channel declares:
 *
 *   evidence  'measured_this_project' — a quantity this project computed from
 *                  its own data, reported as computed
 *             'published'            — a number taken from the literature, with
 *                  the citation and a value grade
 *             'schematic'            — a shape chosen to be legible, carrying no
 *                  measurement at all
 *
 * Citation grading follows the convention the stomach page established:
 *   value grade A  the number was read back from the source text
 *   value grade B  the number is reported from memory of the field and has NOT
 *                  been re-read from the source in this session
 *   paper grade A  the paper identity (DOI/PMID) is checkable
 *   paper grade B  the paper identity is not independently checkable here
 *
 * HARD LINE, stated in every place the clearance channel appears: this project
 * has NO expression-to-CLint calibration. The substrate-decay channel is driven
 * by a USER-SET extraction ratio. It is an input, not a result, and it is not
 * derived from any expression value on this page.
 */

export const FUNCTION_CHANNELS = [
  {
    id: 'oxygen',
    mode: 0,
    zh: '氧梯度',
    label: 'pO2 periportal -> pericentral',
    evidence: 'published',
    summary_zh: '门静脉周围约 60–65 mmHg → 中央静脉周围约 30–35 mmHg，沿窦状隙单调下降。',
    numbers: { periportal_mmHg: [60, 65], pericentral_mmHg: [30, 35] },
    shape: 'linear in the zone coordinate — the PROFILE is not measured, only the '
         + 'two endpoints are taken from the literature',
    citation: 'Jungermann K, Kietzmann T. Oxygen: modulator of metabolic zonation and '
            + 'disease of the liver. Hepatology 2000;31:255-260. '
            + 'DOI 10.1002/hep.510310201',
    value_grade: 'B — the 60-65 / 30-35 mmHg endpoints are reported from the field\'s '
               + 'standard values and were NOT re-read from the source text in this '
               + 'session',
    paper_grade: 'A — DOI is checkable',
    measured_by_this_project: false
  },
  {
    id: 'zonation',
    mode: 1,
    zh: '代谢分区',
    label: 'zone 1 vs zone 3 pathway assignment',
    evidence: 'published',
    summary_zh: '门静脉周围（zone 1）：糖异生、尿素合成、β 氧化、胆汁酸摄取。'
              + '中央静脉周围（zone 3）：糖酵解、脂肪生成、谷氨酰胺合成、异源物代谢（CYP）。',
    zone1: ['gluconeogenesis', 'urea synthesis', 'beta-oxidation', 'bile acid uptake'],
    zone3: ['glycolysis', 'lipogenesis', 'glutamine synthesis', 'xenobiotic metabolism'],
    shape: 'categorical labels on the two ends of the zone axis; no quantity is plotted',
    citation: 'Jungermann K, Kietzmann T. Zonation of parenchymal and nonparenchymal '
            + 'metabolism in liver. Annu Rev Nutr 1996;16:179-203. '
            + 'DOI 10.1146/annurev.nu.16.070196.001143; '
            + 'Kietzmann T. Metabolic zonation of the liver: The oxygen gradient '
            + 'revisited. Redox Biol 2017;11:622-630. DOI 10.1016/j.redox.2017.01.012',
    value_grade: 'not applicable — this channel carries no numeric value',
    paper_grade: 'A — DOIs are checkable',
    measured_by_this_project: false,
    caveat_zh: '这是<b>通路归属</b>的文献定性，不是本项目的测量。'
             + '本项目的 61 基因面板只覆盖 ADME 基因，不含 CPS1/GLUL 等尿素与谷氨酰胺合成酶，'
             + '所以糖代谢与尿素两条通路在本页<b>没有本项目自己的数据支撑</b>。'
  },
  {
    id: 'cyp',
    mode: 2,
    zh: 'CYP 分区（逐基因实测）',
    label: 'per-gene measured central/periportal ratio',
    evidence: 'measured_this_project',
    summary_zh: '逐基因显示本项目实测的中央/门脉周围比值，颜色带的方向跟着比值走。'
              + '<b>不做"CYP 一律中央高"的全局断言</b>。',
    shape: 'the ramp direction is the sign of log(ratio); the amplitude is scaled by '
         + '|log(ratio)| so a flat gene renders flat and a reversed gene ramps the '
         + 'other way',
    citation: 'this project: CELLxGENE Census liver, disease == normal, three '
            + 'zone-annotated hepatocyte populations',
    value_grade: 'computed in this session from the Census pull; see liver_zonation.json',
    paper_grade: 'not applicable — own measurement',
    measured_by_this_project: true
  },
  {
    id: 'clearance',
    mode: 3,
    zh: '沿窦状隙的底物衰减',
    label: 'substrate concentration along the sinusoid',
    evidence: 'schematic',
    summary_zh: '底物浓度从门静脉端到中央静脉端的衰减。'
              + '<b>这条曲线不是由本项目的表达数据算出来的。</b>',
    shape: 'C(zone) = exp(-zone * (-ln(1 - E))) with E the extraction ratio set by the '
         + 'slider. E is an INPUT the reader chooses; it is not fitted and not derived',
    citation: 'none — no parameter in this channel comes from a publication either',
    value_grade: 'not applicable — schematic',
    paper_grade: 'not applicable — schematic',
    measured_by_this_project: false,
    hard_statement_zh: '本项目<b>没有</b>表达量 → CLint 的标定系数，也没有做 IVIVE 外推。'
                     + '把这条衰减曲线读成"由 CYP 表达量推出的清除率"是<b>错误的</b>。'
                     + '滑杆上的抽提比 E 是你自己设的输入值。',
    hard_statement_en: 'This project has NO expression-to-CLint calibration and has not '
                     + 'done IVIVE. Reading this curve as a clearance predicted from the '
                     + 'CYP expression on this page is WRONG. The extraction ratio E is '
                     + 'an input the reader sets.'
  },
  {
    id: 'bile',
    mode: 4,
    zh: '胆汁逆流',
    label: 'canalicular bile, counter-current',
    evidence: 'measured_this_project',
    summary_zh: '胆汁在毛细胆管内<b>反向</b>流动：实质 → 肝门，与同一肝板里的血流相反。'
              + '真实流速约 0.07 cm/s，约为门静脉血的三百分之一。',
    shape: 'the ramp runs from the pericentral end to the periportal end, i.e. opposite '
         + 'to the blood; the particle speed is the derived canalicular velocity times a '
         + 'DECLARED display gain',
    display_gain: 120,
    citation: 'derived in this project from bile output per lobule over the canalicular '
            + 'cross-section; the bile flow rate and liver mass inputs are flagged in '
            + 'liver_model.json as textbook-magnitude, not this project\'s measurements',
    value_grade: 'computed in this session',
    paper_grade: 'not applicable — own derivation from flagged inputs',
    measured_by_this_project: true,
    caveat_zh: '胆汁比血流慢约 300 倍，不加增益在屏幕上看不出在动，'
             + '所以粒子速度乘了一个<b>声明过的 120× 显示增益</b>。'
  }
];

export const CHANNEL_BY_ID = Object.fromEntries(FUNCTION_CHANNELS.map(c => [c.id, c]));
