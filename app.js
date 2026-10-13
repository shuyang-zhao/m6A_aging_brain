// Front end of the Aging Brain Gene Browser. Reads data/meta.json (units, m6A maps, reader binding sources and
// references), data/genes.json (search index: key, human symbol, mouse symbol) and data/genes/<key>.json (one gene,
// rows of [unit, layer, % change, CI low, CI high, P, q, samples]), all written by pipeline/62_website.py. Text of the
// page is in English and Chinese: static text in index.html (.i-en, .i-zh), the text written here in TEXT
const CONFIG = {
  zenodo: "", // Zenodo record of the full tables (Additional files 3 and 4); the download stays inactive while empty
  code: "", // repository of the analysis code
  paper: "", // the article
};
const DEFAULT_GENE = "YTHDF3";
const VIEWS = ["gene", "landscape", "cohorts", "definitions", "sources"];
const CONTEXTS = ["human_protein_P", "human_protein_A", "human_cortex_R", "human_striatum_R", "human_hippocampus_R",
  "human_cerebellum_R", "human_neurons_R", "mouse_protein_P", "mouse_protein_A", "mouse_mrna_R", "mouse_neurons_R"];
const DEFAULT_CONTEXT = "mouse_protein_P";
const LAYER_KEYS = ["P", "A", "R", "S"];
const COLOR = { P: "var(--protein)", A: "var(--adjusted)", R: "var(--bulk)", S: "var(--sn)" };
const AGE_RANGE = { human: [0, 105], mouse: [0, 36] };
const M6A = "m<sup>6</sup>A";

const TEXT = {
  en: {
    site: "Aging Brain Gene Browser",
    layer: { P: "Protein", A: "mRNA-adjusted protein", R: "Bulk mRNA", S: "Single-nucleus mRNA" },
    species: { human: "Human", mouse: "Mouse" },
    axis: { human: "Change per decade of age (%)", mouse: "Change, late vs early age group (%)" },
    scale: { human: "% per decade", mouse: "% late vs early" },
    ageTerm: { "per decade": "per decade", "late vs early": "late vs early" },
    ageUnit: { years: "years", months: "months" },
    region: {},
    kind: { list: "gene list of the original study", peaks: "CLIP peaks (nf-core/clipseq, PureCLIP)",
      mapped_peaks: "CLIP peaks of the other species" },
    mappedFrom: (other) => `, mapped from ${other}`,
    viaOrthologs: ", mapped through one-to-one orthologs",
    search: "Search a gene: human or mouse symbol, or Ensembl id",
    filter: "Filter by name",
    allLayers: "All layers",
    allRegions: "All regions",
    mouseOnly: "mouse only",
    noGene: (q) => `No gene named <b>${q}</b> was quantified in any unit. Search by the official human or mouse
      symbol, or by the Ensembl gene id.`,
    notQuantified: (s) => s === "human" ? "Not quantified in any human unit."
      : "Not quantified in any mouse unit, or without a one-to-one mouse ortholog.",
    unitLayers: "unit layers quantified",
    sig: "q &lt; 0.05 (down, up)",
    load: { mouse: `${M6A} sites`, human: `${M6A} peaks` },
    loadTitle: { mouse: "single-nucleotide m⁶A sites on the canonical transcript, four mouse brain maps",
      human: "m⁶A-seq peaks on the canonical transcript, five human brain regions" },
    binding: `${M6A} reader binding (primary source of each reader)`,
    notBound: `Not bound by any of the 14 ${M6A} readers.`,
    reader: "reader",
    marked: `${M6A} within 50 nt of binding`,
    unmarked: `bound, no ${M6A} nearby`,
    chip: (r, marked, source) => `${r}: ${marked ? "bound with m⁶A within 50 nt" : "bound without m⁶A nearby"}. Source: ${source}`,
    plot: "Age effect in every unit",
    filled: "q &lt; 0.05", open: "q ≥ 0.05", ci: "95% CI",
    tipEffect: "Age effect", tipQ: "q (BH, within unit)", tipSamples: "Samples", tipAge: "Age term", to: "to",
    allValues: "All values", nLayers: (n) => `${n} unit layers`, csv: "Download CSV",
    geneCols: ["Species", "Layer", "Unit", "Region", "Samples", "Change (%)", "95% CI (%)", "P", "q"],
    cohortCols: ["Unit", "Species", "Layers", "Region", "Donors", "Age", "Genes", "Model covariates", "Source"],
    summary: (n, d, layers, adjusted) => `<b>${n}</b> units · <b>${d}</b> datasets · <b>${layers}</b> unit layers of
      protein or RNA${adjusted ? `, and <b>${adjusted}</b> of mRNA-adjusted protein` : ""}`,
    hidden: (k) => `${k} hidden factor${k > 1 ? "s" : ""}`,
    maps: `${M6A} maps`, mapCols: ["Species", "Map", "Sites or peaks"],
    readers: `${M6A} reader binding: primary sources`, readerCols: ["Reader", "Mouse", "Human"],
    refGroups: ["Aging datasets", `${M6A} maps`, `${M6A} reader binding`],
    zenodo: "Open the Zenodo record",
    codeAt: (url) => `The analysis code is at <a href="${url}">${url}</a>.`,
    foot: ["Article", "Code", "Data"],
    theme: { dark: "Switch to light mode", light: "Switch to dark mode" },
    ctx: { human_protein_P: "Protein", human_protein_A: "mRNA-adjusted protein", human_cortex_R: "Cortex mRNA",
      human_striatum_R: "Striatum mRNA", human_hippocampus_R: "Hippocampus mRNA", human_cerebellum_R: "Cerebellum mRNA",
      human_neurons_R: "Neuronal mRNA (single-nucleus)", mouse_protein_P: "Protein", mouse_protein_A: "mRNA-adjusted protein",
      mouse_mrna_R: "Bulk mRNA", mouse_neurons_R: "Neuronal mRNA (single-nucleus)" },
    ovTitle: "Genes with q &lt; 0.05 in each context", down: "down", up: "up",
    nCohorts: (n) => `${n} cohort${n > 1 ? "s" : ""}`, nGenes: (n) => `${n} genes`,
    volcano: "Combined age effect of every gene", color: { dir: "Direction", m6a: `${M6A} sites` },
    ns: "q ≥ 0.05", sigDown: "q &lt; 0.05, down", sigUp: "q &lt; 0.05, up", m6aClass: { mouse: `${M6A} sites`, human: `${M6A} peaks` },
    xAxis: { human: "Combined age effect (% per decade)", mouse: "Combined age effect (% late vs early)" },
    yAxis: "−log₁₀ P", qLine: "q = 0.05", beyond: "▲ beyond the axis (values in the tooltip and table)",
    rank: "Genes ranked", dir: { "": "All", down: "Down", up: "Up" }, minLabel: "|change| ≥",
    sort: { q: "Sort by q", effect: "Sort by effect size" }, more: "Show 50 more",
    shownOf: (a, b) => `${a} of ${b} genes`, csvContext: "Download this context (CSV)", tsvAll: "All contexts (TSV)",
    rankCols: ["#", "Gene", "Change (%)", "95% CI (%)", "P", "q", "Same direction", "I²", `${M6A}`],
    sameTitle: "cohorts whose mean has the sign of the combined estimate, of the cohorts measuring the gene",
    tipCombined: "Combined effect", tipSame: "Same direction", tipCohorts: (a, k) => `${a} of ${k} cohorts`,
    contextCols: ["Context", "Cohorts", "Genes"],
  },
  zh: {
    site: "脑衰老基因浏览器",
    layer: { P: "蛋白", A: "mRNA 校正蛋白", R: "Bulk mRNA", S: "单核 mRNA" },
    species: { human: "人", mouse: "小鼠" },
    axis: { human: "每十年的变化 (%)", mouse: "晚期组相对早期组的变化 (%)" },
    scale: { human: "每十年 %", mouse: "晚期对早期 %" },
    ageTerm: { "per decade": "每十年", "late vs early": "晚期对早期" },
    ageUnit: { years: "岁", months: "月" },
    region: {
      Cerebellum: "小脑", Striatum: "纹状体", Hippocampus: "海马", Amygdala: "杏仁核",
      "Thalamus and hypothalamus": "丘脑与下丘脑", "Midbrain, brainstem and spinal cord": "中脑、脑干与脊髓",
      "Olfactory bulb": "嗅球", Cerebrum: "大脑", "White matter and ventricles": "白质与脑室", Cortex: "皮层",
      "Whole brain": "全脑",
    },
    kind: { list: "原研究的基因列表", peaks: "CLIP 峰（nf-core/clipseq、PureCLIP）", mapped_peaks: "另一物种的 CLIP 峰" },
    mappedFrom: (other) => `，由${other === "human" ? "人" : "小鼠"}映射`,
    viaOrthologs: "，经一对一直系同源映射",
    search: "搜索基因：人或小鼠基因符号，或 Ensembl 编号",
    filter: "按名称筛选",
    allLayers: "全部分子层",
    allRegions: "全部脑区",
    mouseOnly: "仅小鼠",
    noGene: (q) => `没有在任何单元中定量到名为 <b>${q}</b> 的基因。请使用人或小鼠的官方基因符号，或 Ensembl 基因编号。`,
    notQuantified: (s) => s === "human" ? "未在任何人类单元中定量。" : "未在任何小鼠单元中定量，或没有一对一的小鼠直系同源基因。",
    unitLayers: "已定量的单元层",
    sig: "q &lt; 0.05（下降、上升）",
    load: { mouse: `${M6A} 位点`, human: `${M6A} 峰` },
    loadTitle: { mouse: "规范转录本上的单核苷酸 m⁶A 位点，四张小鼠脑图谱", human: "规范转录本上的 m⁶A-seq 峰，五个人脑区" },
    binding: `${M6A} 阅读蛋白结合（各阅读蛋白的主来源）`,
    notBound: `未被 14 个 ${M6A} 阅读蛋白中的任何一个结合。`,
    reader: "阅读蛋白",
    marked: `结合区 50 nt 内有 ${M6A}`,
    unmarked: `有结合，附近无 ${M6A}`,
    chip: (r, marked, source) => `${r}：${marked ? "结合区 50 nt 内有 m⁶A" : "有结合，附近无 m⁶A"}。来源：${source}`,
    plot: "各单元中的年龄效应",
    filled: "q &lt; 0.05", open: "q ≥ 0.05", ci: "95% 置信区间",
    tipEffect: "年龄效应", tipQ: "q（单元内 BH）", tipSamples: "样本数", tipAge: "年龄项", to: "至",
    allValues: "全部数值", nLayers: (n) => `${n} 个单元层`, csv: "下载 CSV",
    geneCols: ["物种", "分子层", "单元", "脑区", "样本数", "变化 (%)", "95% 置信区间 (%)", "P", "q"],
    cohortCols: ["单元", "物种", "分子层", "脑区", "供体数", "年龄", "基因数", "模型协变量", "来源"],
    summary: (n, d, layers, adjusted) => `<b>${n}</b> 个单元 · <b>${d}</b> 个数据集 · <b>${layers}</b> 个蛋白或 RNA 单元层${adjusted ? `，另有 <b>${adjusted}</b> 个 mRNA 校正蛋白层` : ""}`,
    hidden: (k) => `${k} 个隐藏因子`,
    maps: `${M6A} 图谱`, mapCols: ["物种", "图谱", "位点或峰数"],
    readers: `${M6A} 阅读蛋白结合：主来源`, readerCols: ["阅读蛋白", "小鼠", "人"],
    refGroups: ["衰老数据集", `${M6A} 图谱`, `${M6A} 阅读蛋白结合`],
    zenodo: "打开 Zenodo 记录",
    codeAt: (url) => `分析代码见 <a href="${url}">${url}</a>。`,
    foot: ["论文", "代码", "数据"],
    theme: { dark: "切换到日间模式", light: "切换到夜间模式" },
    ctx: { human_protein_P: "蛋白", human_protein_A: "mRNA 校正蛋白", human_cortex_R: "皮层 mRNA", human_striatum_R: "纹状体 mRNA",
      human_hippocampus_R: "海马 mRNA", human_cerebellum_R: "小脑 mRNA", human_neurons_R: "神经元 mRNA（单核）",
      mouse_protein_P: "蛋白", mouse_protein_A: "mRNA 校正蛋白", mouse_mrna_R: "Bulk mRNA", mouse_neurons_R: "神经元 mRNA（单核）" },
    ovTitle: "各情境中 q &lt; 0.05 的基因数", down: "下降", up: "上升",
    nCohorts: (n) => `${n} 个队列`, nGenes: (n) => `${n} 个基因`,
    volcano: "每个基因的合并年龄效应", color: { dir: "方向", m6a: `${M6A} 位点` },
    ns: "q ≥ 0.05", sigDown: "q &lt; 0.05，下降", sigUp: "q &lt; 0.05，上升", m6aClass: { mouse: `${M6A} 位点`, human: `${M6A} 峰` },
    xAxis: { human: "合并年龄效应（每十年 %）", mouse: "合并年龄效应（晚期对早期 %）" },
    yAxis: "−log₁₀ P", qLine: "q = 0.05", beyond: "▲ 超出坐标轴（确切数值见悬停提示与表格）",
    rank: "基因排名", dir: { "": "全部", down: "下降", up: "上升" }, minLabel: "|变化| ≥",
    sort: { q: "按 q 排序", effect: "按效应大小排序" }, more: "再显示 50 个",
    shownOf: (a, b) => `显示 ${a} / ${b} 个基因`, csvContext: "下载本情境 (CSV)", tsvAll: "全部情境 (TSV)",
    rankCols: ["#", "基因", "变化 (%)", "95% 置信区间 (%)", "P", "q", "同方向队列", "I²", `${M6A}`],
    sameTitle: "均值符号与合并估计相同的队列数 / 测到该基因的队列数",
    tipCombined: "合并效应", tipSame: "同方向队列", tipCohorts: (a, k) => `${k} 个队列中的 ${a} 个`,
    contextCols: ["情境", "队列", "基因数"],
  },
};

const ICONS = {
  dark: '<svg viewBox="0 0 20 20" aria-hidden="true"><circle cx="10" cy="10" r="3.6"/><path d="M10 1.8v2.2M10 16v2.2M1.8 10H4M16 10h2.2M4.2 4.2l1.6 1.6M14.2 14.2l1.6 1.6M4.2 15.8l1.6-1.6M14.2 5.8l1.6-1.6"/></svg>',
  light: '<svg viewBox="0 0 20 20" aria-hidden="true"><path d="M16.5 12.4A7 7 0 0 1 7.6 3.5a7 7 0 1 0 8.9 8.9z"/></svg>',
};

let meta, index, lang, shown = null, current = null, refNo = {};
const cache = {};
const $ = (s, el = document) => el.querySelector(s);
const T = () => TEXT[lang];
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
const rich = (s) => esc(s).replace(/m6A/g, M6A).replace(/\bN6-/g, "N<sup>6</sup>-"); // m6A and N6 as superscripts
const plain = (s) => String(s).replace(/m6A/g, "m⁶A"); // the same in attributes, which take no markup
const minus = (s) => String(s).replace("-", "−");
const pct = (v) => (v > 0 ? "+" : "") + minus(v.toFixed(1));
const num = (v) => v.toLocaleString("en-US");
const region = (r) => T().region[r] || r;
const store = (key, value) => { try { localStorage.setItem(key, value); } catch (e) { /* private window */ } };
const getJSON = (path) => fetch(path).then((r) => r.json());
const fmtP = (p) => {
  if (p >= 0.001) return p.toPrecision(2);
  const [m, e] = p.toExponential(1).split("e");
  return `${m} × 10<sup>${minus(e)}</sup>`;
};
const ref = (key) => `<sup class="ref"><a href="#sources/ref-${refNo[key]}" title="${esc(plain(meta.refs[key].text))}">${refNo[key]}</a></sup>`;

async function init() {
  [meta, index] = await Promise.all([getJSON("data/meta.json"), getJSON("data/genes.json")]);
  const order = [...meta.units.map((u) => u.ref), ...["mouse", "human"].flatMap((s) => meta.maps[s].map((m) => m.ref))];
  order.push(...["mouse", "human"].flatMap((s) => meta.readers.map((r) => meta.sources[s][r].ref)));
  [...new Set(order)].forEach((key, i) => (refNo[key] = i + 1));
  initTheme();
  initSearch();
  initCohorts();
  initLandscape();
  $("#lang").addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) setLang(b.dataset.v); });
  addEventListener("hashchange", route);
  addEventListener("resize", () => { clearTimeout(init.timer); init.timer = setTimeout(() => { drawPlots(); drawVolcano(); }, 120); });
  setLang(document.documentElement.dataset.lang === "zh" ? "zh" : "en");
  route();
}

function setLang(next) {
  lang = next;
  const root = document.documentElement;
  root.dataset.lang = next;
  root.lang = next === "zh" ? "zh-CN" : "en";
  store("lang", next);
  $("#lang").querySelectorAll("button").forEach((b) => b.classList.toggle("on", b.dataset.v === next));
  $("#query").placeholder = T().search;
  $("#f-text").placeholder = T().filter;
  if (CONFIG.zenodo) Object.assign($("#zenodo"), { href: CONFIG.zenodo, textContent: T().zenodo }).removeAttribute("aria-disabled");
  $("#code-link").innerHTML = CONFIG.code ? T().codeAt(CONFIG.code) : "";
  $("#foot-links").innerHTML = [CONFIG.paper, CONFIG.code, CONFIG.zenodo].map((url, i) => [T().foot[i], url])
    .filter(([, url]) => url).map(([name, url]) => `<a href="${url}">${name}</a>`).join("");
  paintTheme();
  fillFilters();
  renderCohorts();
  renderSources();
  renderContextTable();
  renderLandscape();
  if (shown) renderGene();
  document.title = shown && !$("#view-gene").hidden ? `${shown.name} · ${T().site}` : T().site;
}

function initTheme() {
  const system = matchMedia("(prefers-color-scheme: dark)");
  $("#theme").addEventListener("click", () => {
    const next = isDark() ? "light" : "dark";
    document.documentElement.dataset.theme = next;
    store("theme", next);
    paintTheme();
  });
  system.addEventListener("change", paintTheme);
  isDark.system = system;
}

function isDark() {
  const chosen = document.documentElement.dataset.theme;
  return chosen ? chosen === "dark" : isDark.system.matches;
}

function paintTheme() {
  const mode = isDark() ? "dark" : "light";
  Object.assign($("#theme"), { innerHTML: ICONS[mode], title: T().theme[mode] });
  $("#theme").setAttribute("aria-label", T().theme[mode]);
  drawVolcano();
}

function route() {
  const [view, arg] = decodeURIComponent(location.hash.slice(1)).split(/\/(.*)/s);
  const v = VIEWS.includes(view) ? view : "gene";
  for (const name of VIEWS) $(`#view-${name}`).hidden = name !== v;
  document.querySelectorAll(".tabs a").forEach((a) => a.classList.toggle("on", a.dataset.view === v));
  document.title = T().site;
  if (v === "gene") {
    if (arg || !current) showGene(arg || DEFAULT_GENE);
    else if (shown) { document.title = `${shown.name} · ${T().site}`; drawPlots(); }
  }
  if (v === "landscape") {
    const id = CONTEXTS.includes(arg) ? arg : land.id;
    if (id !== land.id || !land.rows) openContext(id); else drawVolcano();
  }
  const same = v === route.view;
  route.view = v;
  document.querySelectorAll(".hl").forEach((el) => el.classList.remove("hl"));
  const target = v !== "gene" && arg && document.getElementById(arg);
  if (target) { target.classList.add("hl"); target.scrollIntoView({ block: v === "sources" ? "center" : "start" }); }
  else if (!same || v === "gene") scrollTo(0, 0);
}

// ---------- search ----------
function lookup(q) {
  q = q.trim().toUpperCase();
  return index.find((g) => g[1].toUpperCase() === q) || index.find((g) => g[2].toUpperCase() === q) ||
    index.find((g) => g[0] === q);
}

function matches(q) {
  q = q.trim().toUpperCase();
  if (!q) return [];
  const score = (g) => {
    const names = [g[1].toUpperCase(), g[2].toUpperCase()];
    if (names.includes(q) || g[0] === q) return 0;
    if (names.some((n) => n.startsWith(q)) || g[0].startsWith(q)) return 1;
    return q.length > 2 && names.some((n) => n.includes(q)) ? 2 : 9;
  };
  return index.map((g) => [score(g), g]).filter(([s]) => s < 9)
    .sort((a, b) => a[0] - b[0] || (a[1][1] || a[1][2]).length - (b[1][1] || b[1][2]).length)
    .slice(0, 8).map(([, g]) => g);
}

function initSearch() {
  const input = $("#query"), list = $("#suggest");
  let items = [], active = 0;
  const go = (g) => { list.hidden = true; location.hash = `gene/${g[1] || g[2]}`; input.blur(); };
  const paint = () => {
    list.innerHTML = items.map((g, i) => `<li role="option" class="${i === active ? "on" : ""}" data-i="${i}">
      <b>${esc(g[1] || g[2])}</b><span class="alt">${g[1] && g[2] ? esc(g[2]) + " · " : ""}${g[1] ? "" : T().mouseOnly + " · "}${g[0]}</span></li>`).join("");
    list.hidden = !items.length;
  };
  input.addEventListener("input", () => { items = matches(input.value); active = 0; paint(); });
  input.addEventListener("keydown", (e) => {
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      active = (active + (e.key === "ArrowDown" ? 1 : items.length - 1)) % Math.max(items.length, 1);
      paint(); e.preventDefault();
    } else if (e.key === "Enter") {
      const g = items[active] || lookup(input.value);
      if (g) go(g); else { list.hidden = true; location.hash = `gene/${input.value.trim()}`; }
    } else if (e.key === "Escape") list.hidden = true;
  });
  list.addEventListener("mousedown", (e) => { const li = e.target.closest("li"); if (li) go(items[li.dataset.i]); });
  input.addEventListener("blur", () => setTimeout(() => (list.hidden = true), 100));
}

// ---------- gene ----------
async function showGene(q) {
  const hit = lookup(q);
  current = q;
  if (!hit) {
    shown = null;
    $("#gene").innerHTML = `<div class="card message">${T().noGene(esc(q))}</div>`;
    return;
  }
  cache[hit[0]] ??= await getJSON(`data/genes/${hit[0]}.json`);
  if (current !== q) return;
  const gene = cache[hit[0]];
  const rows = gene.rows.map(([u, layer, pct, lo, hi, p, q, n]) => ({ unit: meta.units[u], u, layer, pct, lo, hi, p, q, n }))
    .sort((a, b) => (a.unit.species > b.unit.species ? 1 : a.unit.species < b.unit.species ? -1 : 0) ||
      LAYER_KEYS.indexOf(a.layer) - LAYER_KEYS.indexOf(b.layer) || a.u - b.u);
  shown = { gene, rows, name: hit[1] || hit[2] };
  $("#query").value = shown.name;
  document.title = `${shown.name} · ${T().site}`;
  renderGene();
}

function renderGene() {
  const { gene, rows, name } = shown, box = $("#gene");
  const open = $("#gene details")?.open;
  box.innerHTML = geneHead(gene, name) + `<div class="facts">${["human", "mouse"].map((s) => fact(gene, rows, s)).join("")}</div>` +
    plotCard(rows) + tableCard(rows);
  if (open) $("#gene details").open = true;
  drawPlots();
  $("#csv").addEventListener("click", (e) => { e.preventDefault(); download(rows, name); });
  bindTips(box, (el) => tipHTML(rows[el.dataset.i]));
}

function geneHead(gene, name) {
  const ids = [["human", "Homo_sapiens"], ["mouse", "Mus_musculus"]].filter(([s]) => gene[s])
    .map(([s, path]) => `<a href="https://www.ensembl.org/${path}/Gene/Summary?g=${gene[s].id}">${gene[s].id}</a>`);
  const alias = gene.human && gene.mouse ? `<span class="alias">${esc(gene.mouse.symbol)}</span>` : "";
  return `<div class="gene-head"><h2>${esc(name)}</h2>${alias}<span class="ids">${ids.join("")}</span></div>`;
}

function fact(gene, rows, s) {
  const g = gene[s], tx = T();
  if (!g) return `<div class="card pad"><div class="fact-title"><h3>${tx.species[s]}</h3></div>
    <p class="muted">${tx.notQuantified(s)}</p></div>`;
  const all = meta.units.filter((u) => u.species === s).reduce((n, u) => n + Object.keys(u.layers).length, 0);
  const mine = rows.filter((r) => r.unit.species === s), sig = mine.filter((r) => r.q < 0.05);
  const down = sig.filter((r) => r.pct < 0).length;
  const chips = meta.readers.filter((r) => r in g.readers).map((r) => {
    const src = meta.sources[s][r];
    const source = tx.kind[src.kind] + (src.mapped ? tx.viaOrthologs : "");
    return `<span class="chip ${g.readers[r] ? "marked" : ""}" title="${esc(tx.chip(r, g.readers[r], source))}">${r}${ref(src.ref)}</span>`;
  }).join("");
  const mapRefs = [...new Set(meta.maps[s].map((m) => m.ref))].map(ref).join("");
  return `<div class="card pad">
    <div class="fact-title"><h3>${tx.species[s]} <span class="muted">${esc(g.symbol)}</span></h3><span class="muted">${tx.scale[s]}</span></div>
    <div class="stats">
      <div class="stat"><div class="v">${mine.length}<span class="of"> / ${all}</span></div><div class="k">${tx.unitLayers}</div></div>
      <div class="stat"><div class="v">${sig.length ? `↓${down}&thinsp; ↑${sig.length - down}` : "0"}</div><div class="k">${tx.sig}</div></div>
      <div class="stat" title="${esc(tx.loadTitle[s])}"><div class="v">${g.m6a ?? "–"}</div><div class="k">${tx.load[s]}${mapRefs}</div></div>
    </div>
    <div class="row-label">${tx.binding}</div>
    <div class="chips">${chips || `<span class="muted">${tx.notBound}</span>`}</div>
    <div class="legend-inline"><span><span class="chip marked">${tx.reader}</span> ${tx.marked}</span><span><span class="chip">${tx.reader}</span> ${tx.unmarked}</span></div>
  </div>`;
}

function plotCard(rows) {
  const tx = T(), species = ["human", "mouse"].filter((s) => rows.some((r) => r.unit.species === s));
  const legend = LAYER_KEYS.filter((k) => rows.some((r) => r.layer === k))
    .map((k) => `<span><span class="dot" style="background:${COLOR[k]}"></span>${tx.layer[k]}</span>`).join("");
  return `<div class="card plot-card">
    <div class="plot-head"><h3>${tx.plot}</h3>
      <div class="legend">${legend}<span><span class="dot fill"></span>${tx.filled}</span><span><span class="dot open"></span>${tx.open}</span><span>— ${tx.ci}</span></div></div>
    <div class="panels ${species.length === 1 ? "one" : ""}">${species.map((s) => `<div class="panel"><h4>${tx.species[s]}<span>${tx.scale[s]}</span></h4>
      <div class="forest-box" data-species="${s}"></div></div>`).join("")}</div></div>`;
}

function drawPlots() {
  if (!shown) return;
  document.querySelectorAll(".forest-box").forEach((box) => {
    const rows = shown.rows.map((r, i) => ({ ...r, i })).filter((r) => r.unit.species === box.dataset.species);
    box.innerHTML = forest(rows, box.dataset.species, box.clientWidth);
  });
}

function niceTicks(lim) {
  const raw = lim / 2, p = 10 ** Math.floor(Math.log10(raw)), f = raw / p;
  const step = (f <= 1 ? 1 : f <= 2 ? 2 : f <= 2.5 ? 2.5 : f <= 5 ? 5 : 10) * p;
  const n = Math.ceil(lim / step - 1e-9);
  return Array.from({ length: 2 * n + 1 }, (_, i) => +((i - n) * step).toPrecision(6));
}

// forest plot of one species: rows grouped by layer, symmetric axis around zero set by the bulk of the estimates;
// confidence intervals beyond the axis (small units) end in an arrow, and estimates beyond it become a triangle at the
// edge (values in the tooltip and the table)
function forest(rows, species, width) {
  if (!width) return "";
  const L = width < 480 ? 116 : 172, R = 12, top = 26, rowH = 17, headH = 26, bottom = 46;
  const groups = LAYER_KEYS.map((k) => [k, rows.filter((r) => r.layer === k)]).filter(([, r]) => r.length);
  const quantile = (v, q) => v.sort((a, b) => a - b)[Math.floor((v.length - 1) * q)];
  const points = quantile(rows.map((r) => Math.abs(r.pct)), 0.9);
  const ends = quantile(rows.flatMap((r) => [Math.abs(r.lo), Math.abs(r.hi)]), 0.75);
  const ticks = niceTicks(Math.max(0.5, points * 1.4, Math.min(ends, points * 4)));
  const lim = ticks.at(-1), span = width - L - R;
  const x = (v) => (L + ((Math.max(-lim, Math.min(lim, v)) + lim) / (2 * lim)) * span).toFixed(1);
  const height = top + groups.reduce((s, [, r]) => s + headH + r.length * rowH, 0) + bottom;
  const maxChars = Math.floor((L - 14) / 6.1);
  const out = [];
  for (const t of ticks) {
    out.push(`<line class="${t === 0 ? "zero" : "grid"}" x1="${x(t)}" x2="${x(t)}" y1="${top - 4}" y2="${height - bottom + 4}"/>`);
    for (const ty of [top - 10, height - bottom + 18])
      out.push(`<text class="tick" x="${x(t)}" y="${ty}" text-anchor="middle">${t > 0 ? "+" : ""}${minus(t)}</text>`);
  }
  let y = top;
  for (const [k, rs] of groups) {
    const c = COLOR[k];
    out.push(`<text class="group" x="0" y="${y + 17}">${T().layer[k]}<tspan class="tick" dx="6">${rs.length}</tspan></text>`);
    y += headH;
    for (const r of rs) {
      const cy = y + rowH / 2, label = r.unit.unit.length > maxChars ? r.unit.unit.slice(0, maxChars - 1) + "…" : r.unit.unit;
      const fill = r.q < 0.05 ? c : "var(--surface)";
      const arrows = [[r.lo < -lim, -1], [r.hi > lim, 1]].filter(([beyond]) => beyond)
        .map(([, d]) => `<path d="M${x(d * lim)} ${cy}l${-d * 5} -3.5v7z" style="fill:${c}"/>`).join("");
      const mark = Math.abs(r.pct) > lim
        ? `<path d="M${x(r.pct)} ${cy}l${r.pct > 0 ? -9 : 9} -5v10z" style="fill:${fill};stroke:${c}" stroke-width="1.4" stroke-linejoin="round"/>`
        : `<circle cx="${x(r.pct)}" cy="${cy}" r="4.2" style="fill:${fill};stroke:${c}" stroke-width="1.6"/>`;
      out.push(`<g class="row" data-i="${r.i}"><rect class="hit" x="0" y="${y}" width="${width}" height="${rowH}" rx="3"/>
        <text x="${L - 10}" y="${cy + 4}" text-anchor="end">${esc(label)}</text>
        <line x1="${x(r.lo)}" x2="${x(r.hi)}" y1="${cy}" y2="${cy}" style="stroke:${c}" stroke-width="1.6" stroke-linecap="round"/>${arrows}${mark}</g>`);
      y += rowH;
    }
  }
  out.push(`<text class="axis-title" x="${L + span / 2}" y="${height - 6}" text-anchor="middle">${T().axis[species]}</text>`);
  return `<svg class="forest" viewBox="0 0 ${width} ${height}" height="${height}" role="img"
    aria-label="${esc(T().plot)}: ${T().species[species]}">${out.join("")}</svg>`;
}

function tipHTML(r) {
  const u = r.unit, tx = T();
  return `<b>${esc(u.unit)}</b><div>${tx.layer[r.layer]} · ${esc(region(u.region))}</div>
    <div class="line" style="margin-top:6px"><span>${tx.tipEffect}</span><b>${pct(r.pct)}%</b></div>
    <div class="line"><span>${tx.ci}</span><span>${pct(r.lo)} ${tx.to} ${pct(r.hi)}%</span></div>
    <div class="line"><span>P</span><span>${fmtP(r.p)}</span></div>
    <div class="line"><span>${tx.tipQ}</span><span>${fmtP(r.q)}</span></div>
    <div class="line"><span>${tx.tipSamples}</span><span>${r.n}</span></div>
    <div class="line"><span>${tx.tipAge}</span><span>${tx.ageTerm[u.age_term]}</span></div>
    <div class="src">[${refNo[u.ref]}] ${rich(meta.refs[u.ref].text)}</div>`;
}

function bindTips(root, html) {
  const tip = $("#tip");
  root.onmousemove = (e) => {
    const el = e.target.closest("[data-i]");
    if (!el || !root.contains(el)) { tip.hidden = true; return; }
    tip.innerHTML = html(el);
    tip.hidden = false;
    const w = tip.offsetWidth, h = tip.offsetHeight;
    tip.style.left = `${Math.min(e.clientX + 16, innerWidth - w - 8)}px`;
    tip.style.top = `${e.clientY + 16 + h > innerHeight ? e.clientY - h - 12 : e.clientY + 16}px`;
  };
  root.onmouseleave = () => (tip.hidden = true);
}

function tableCard(rows) {
  const tx = T();
  const body = rows.map((r) => `<tr>
    <td>${tx.species[r.unit.species]}</td>
    <td class="nowrap"><span class="layer-key ${r.layer}"></span>${tx.layer[r.layer]}</td>
    <td>${esc(r.unit.unit)}${ref(r.unit.ref)}</td>
    <td>${esc(region(r.unit.region))}</td>
    <td class="num">${r.n}</td>
    <td class="num nowrap ${r.q < 0.05 ? "sig" : ""}">${pct(r.pct)}</td>
    <td class="num nowrap muted">${pct(r.lo)} ${tx.to} ${pct(r.hi)}</td>
    <td class="num nowrap">${fmtP(r.p)}</td>
    <td class="num nowrap">${fmtP(r.q)}</td></tr>`).join("");
  const head = tx.geneCols.map((c, i) => `<th class="${i >= 4 ? "num" : ""}">${c}</th>`).join("");
  return `<details class="card table-card"><summary class="table-head"><h3>${tx.allValues} <span class="muted">${tx.nLayers(rows.length)}</span></h3>
    <button id="csv" class="button ghost">${tx.csv}</button></summary>
    <div class="scroll"><table class="data"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div></details>`;
}

function download(rows, name) {
  const head = ["species", "layer", "unit_id", "unit", "dataset", "region", "samples", "effect_pct", "ci_low_pct",
    "ci_high_pct", "p", "q", "age_term", "reference", "doi"];
  const cell = (v) => (/[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : v);
  const lines = rows.map((r) => [r.unit.species, TEXT.en.layer[r.layer], r.unit.unit_id, r.unit.unit, r.unit.dataset,
    r.unit.region, r.n, r.pct, r.lo, r.hi, r.p, r.q, r.unit.age_term, meta.refs[r.unit.ref].text,
    meta.refs[r.unit.ref].doi].map(cell).join(","));
  const a = Object.assign(document.createElement("a"), {
    href: URL.createObjectURL(new Blob([[head.join(","), ...lines].join("\n")], { type: "text/csv" })),
    download: `${name}_age_effects.csv`,
  });
  a.click();
  URL.revokeObjectURL(a.href);
}

// ---------- cohorts ----------
const layersOf = (u) => LAYER_KEYS.filter((k) => k in u.layers);
const COLUMNS = [ // sort keys of the cohort table, in the order of cohortCols; true: numeric
  [(u) => u.unit], [(u) => u.species], [(u) => layersOf(u).join("")], [(u) => u.region], [(u) => u.donors, true],
  [(u) => u.age_min, true], [(u) => Math.max(...Object.values(u.layers)), true], [(u) => u.covariates],
];
let sortBy = null, sortDir = 1;

function initCohorts() {
  $("#f-species").addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (!b) return;
    $("#f-species").querySelectorAll("button").forEach((x) => x.classList.toggle("on", x === b));
    renderCohorts();
  });
  ["#f-layer", "#f-region", "#f-text"].forEach((s) => $(s).addEventListener("input", renderCohorts));
  $("#cohort-table").addEventListener("click", (e) => {
    const th = e.target.closest("th.sortable");
    if (!th) return;
    const c = +th.dataset.c;
    sortDir = sortBy === c ? -sortDir : 1;
    sortBy = c;
    renderCohorts();
  });
}

function fillFilters() {
  const layer = $("#f-layer").value, chosen = $("#f-region").value;
  const regions = [...new Set(meta.units.map((u) => u.region))].sort((a, b) => region(a).localeCompare(region(b), lang));
  $("#f-layer").innerHTML = `<option value="">${T().allLayers}</option>` +
    LAYER_KEYS.map((k) => `<option value="${k}">${T().layer[k]}</option>`).join("");
  $("#f-region").innerHTML = `<option value="">${T().allRegions}</option>` +
    regions.map((r) => `<option value="${esc(r)}">${esc(region(r))}</option>`).join("");
  $("#f-layer").value = layer;
  $("#f-region").value = chosen;
}

function ageBar(u) {
  const [lo, hi] = AGE_RANGE[u.species], at = (v) => ((v - lo) / (hi - lo)) * 100;
  return `<div class="age-bar"><i style="left:${at(u.age_min)}%;width:${Math.max(at(u.age_max) - at(u.age_min), 2)}%"></i></div>`;
}

function renderCohorts() {
  const tx = T(), species = $("#f-species .on").dataset.v, layer = $("#f-layer").value, chosen = $("#f-region").value;
  const text = $("#f-text").value.trim().toLowerCase();
  let units = meta.units.filter((u) => (!species || u.species === species) && (!layer || layer in u.layers) &&
    (!chosen || u.region === chosen) && (!text || `${u.unit} ${u.dataset} ${u.unit_id}`.toLowerCase().includes(text)));
  if (sortBy !== null) {
    const [key] = COLUMNS[sortBy];
    units = [...units].sort((a, b) => (key(a) > key(b) ? 1 : key(a) < key(b) ? -1 : 0) * sortDir);
  }
  const layers = units.flatMap((u) => layersOf(u).filter((k) => !layer || k === layer));
  const adjusted = layers.filter((k) => k === "A").length;
  $("#cohort-summary").innerHTML = tx.summary(units.length, new Set(units.map((u) => u.ref)).size, layers.length - adjusted, adjusted);
  const head = tx.cohortCols.map((name, i) => i < COLUMNS.length
    ? `<th class="sortable ${COLUMNS[i][1] ? "num" : ""} ${sortBy === i ? "sorted" : ""}" data-c="${i}">${name}${sortBy === i ? (sortDir > 0 ? " ↑" : " ↓") : ""}</th>`
    : `<th>${name}</th>`).join("");
  const body = units.map((u) => `<tr>
    <td class="unit"><b>${esc(u.unit)}</b><div class="muted">${esc(u.unit_id)}</div></td>
    <td>${tx.species[u.species]}</td>
    <td>${layersOf(u).map((k) => `<div class="layer-tag"><span class="layer-key ${k}"></span>${tx.layer[k]}</div>`).join("")}</td>
    <td>${esc(region(u.region))}</td>
    <td class="num">${u.donors}</td>
    <td class="nowrap">${+u.age_min.toFixed(1)}–${+u.age_max.toFixed(1)} ${tx.ageUnit[u.age_unit]}<div class="muted">${tx.ageTerm[u.age_term]}</div>${ageBar(u)}</td>
    <td class="num">${layersOf(u).map((k) => num(u.layers[k])).join("<br>")}</td>
    <td class="small">${esc(u.covariates.replaceAll(",", ", "))}${u.hidden_factors ? `; ${tx.hidden(u.hidden_factors)}` : ""}</td>
    <td class="nowrap">${esc(u.dataset)}${ref(u.ref)}${u.doi ? `<div><a href="https://doi.org/${u.doi}">doi</a></div>` : ""}</td></tr>`).join("");
  $("#cohort-table").innerHTML = `<thead><tr>${head}</tr></thead><tbody>${body}</tbody>`;
}

// ---------- sources ----------
function renderSources() {
  const tx = T();
  const item = (key) => {
    const r = meta.refs[key];
    return `<li id="ref-${refNo[key]}"><span class="n">${refNo[key]}</span><div>${rich(r.text)}<div class="meta">
      ${r.doi ? `<a href="https://doi.org/${r.doi}">doi:${esc(r.doi)}</a>` : ""}${r.accession ? `<span>${esc(r.accession)}</span>` : ""}</div></div></li>`;
  };
  const keys = Object.keys(refNo);
  const nUnits = new Set(meta.units.map((u) => u.ref)).size;
  const nMaps = new Set(["mouse", "human"].flatMap((s) => meta.maps[s].map((m) => m.ref))).size;
  const group = (title, from, to) => `<h4>${title}</h4><ol>${keys.slice(from, to).map(item).join("")}</ol>`;
  const maps = ["mouse", "human"].flatMap((s) => meta.maps[s].map((m) => `<tr><td>${tx.species[s]}</td>
    <td>${rich(m.source)}${ref(m.ref)}</td><td class="num">${num(m.sites)}</td></tr>`)).join("");
  const readers = meta.readers.map((r) => `<tr><td><b>${r}</b></td>${["mouse", "human"].map((s) => {
    const src = meta.sources[s][r];
    const mapped = src.mapped && src.kind === "list" ? tx.mappedFrom(s === "mouse" ? "human" : "mouse") : "";
    return `<td>${tx.kind[src.kind]}${mapped}${ref(src.ref)}<div class="muted">${esc(src.source)}</div></td>`;
  }).join("")}</tr>`).join("");
  const table = (title, cols, body) => `<div class="card table-card"><div class="table-head"><h3>${title}</h3></div>
    <div class="scroll"><table class="data"><thead><tr>${cols.map((c, i) => `<th class="${i === 2 && cols === tx.mapCols ? "num" : ""}">${c}</th>`).join("")}</tr></thead>
    <tbody>${body}</tbody></table></div></div>`;
  $("#source-tables").innerHTML = `<div class="grid-2">${table(tx.maps, tx.mapCols, maps)}${table(tx.readers, tx.readerCols, readers)}</div>`;
  $("#references").innerHTML = `<div class="refs">${group(tx.refGroups[0], 0, nUnits)}${group(tx.refGroups[1], nUnits, nUnits + nMaps)}
    ${group(tx.refGroups[2], nUnits + nMaps)}</div>`;
}

// ---------- landscape ----------
const land = { id: DEFAULT_CONTEXT, rows: null, color: "dir", dir: "", min: 0, sort: "q", n: 50 };
const landCache = {};
const quantileOf = (values, q) => [...values].sort((a, b) => a - b)[Math.floor((values.length - 1) * q)];
const contextName = (id) => `${T().species[id.split("_")[0]]} · ${T().ctx[id]}`;

function initLandscape() {
  const go = (id) => (location.hash = `landscape/${id}`);
  $("#l-context").addEventListener("change", (e) => go(e.target.value));
  $("#overview").addEventListener("click", (e) => { const row = e.target.closest(".ov-row"); if (row) go(row.dataset.id); });
  const segment = (sel, key, after) => $(sel).addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (!b) return;
    land[key] = b.dataset.v;
    land.n = 50;
    after();
  });
  segment("#l-color", "color", () => { paintSegments(); drawVolcano(); });
  segment("#l-dir", "dir", () => { paintSegments(); renderRank(); });
  $("#l-min").addEventListener("input", (e) => { land.min = Math.abs(+e.target.value || 0); land.n = 50; renderRank(); });
  $("#l-sort").addEventListener("change", (e) => { land.sort = e.target.value; land.n = 50; renderRank(); });
  $("#l-more").addEventListener("click", () => { land.n += 50; renderRank(); });
  $("#l-csv").addEventListener("click", downloadContext);
  const canvasBox = $("#volcano");
  canvasBox.addEventListener("mousemove", (e) => hoverVolcano(e));
  canvasBox.addEventListener("mouseleave", () => { $("#tip").hidden = true; $("#v-ring")?.setAttribute("r", 0); });
  canvasBox.addEventListener("click", (e) => { const r = hoverVolcano(e); if (r) location.hash = `gene/${r.symbol}`; });
}

async function openContext(id) {
  land.id = id;
  land.n = 50;
  landCache[id] ??= (await getJSON(`data/landscape/${id}.json`))
    .map(([symbol, pct, lo, hi, p, q, k, same, i2, m6a]) => ({ symbol, pct, lo, hi, p, q, k, same, i2, m6a }));
  if (land.id !== id) return;
  land.rows = landCache[id];
  renderLandscape();
}

function paintSegments() {
  const tx = T();
  $("#l-color").innerHTML = ["dir", "m6a"].map((v) => `<button data-v="${v}" class="${land.color === v ? "on" : ""}">${tx.color[v]}</button>`).join("");
  $("#l-dir").innerHTML = ["", "down", "up"].map((v) => `<button data-v="${v}" class="${land.dir === v ? "on" : ""}">${tx.dir[v]}</button>`).join("");
}

function renderLandscape() {
  if (!meta) return;
  const tx = T();
  $("#ov-title").innerHTML = tx.ovTitle;
  $("#ov-legend").innerHTML = `<span><span class="dot down"></span>${tx.down}</span><span><span class="dot up"></span>${tx.up}</span>`;
  $("#volcano-title").innerHTML = tx.volcano;
  $("#rank-title").innerHTML = tx.rank;
  $("#l-min-label").innerHTML = tx.minLabel;
  $("#l-sort").innerHTML = ["q", "effect"].map((v) => `<option value="${v}">${tx.sort[v]}</option>`).join("");
  $("#l-sort").value = land.sort;
  $("#l-csv").textContent = tx.csvContext;
  $("#l-tsv").textContent = tx.tsvAll;
  $("#l-more").textContent = tx.more;
  $("#l-context").innerHTML = CONTEXTS.map((id) => `<option value="${id}">${contextName(id)}</option>`).join("");
  $("#l-context").value = land.id;
  paintSegments();
  renderOverview();
  drawVolcano();
  renderRank();
}

function renderOverview() {
  const tx = T(), entries = CONTEXTS.map((id) => meta.landscape.find((c) => c.id === id));
  const max = Math.max(1, ...entries.flatMap((c) => [c.down, c.up]));
  const side = (n, cls) => `<div class="ov-side ${cls}"><i style="width:${(n / max) * 82}%"></i><span>${num(n)}</span></div>`;
  $("#overview").innerHTML = ["human", "mouse"].map((s) => `<div class="ov-group">${tx.species[s]}</div>` +
    entries.filter((c) => c.species === s).map((c) => `<div class="ov-row ${c.id === land.id ? "on" : ""}" data-id="${c.id}">
      <div class="ov-label">${tx.ctx[c.id]}<span class="muted">${tx.nCohorts(c.cohorts.length)} · ${tx.nGenes(num(c.genes))}</span></div>
      <div class="ov-bars">${side(c.down, "left")}<div class="axis"></div>${side(c.up, "right")}</div></div>`).join("")).join("");
}

function classOf(r) {
  if (land.color === "dir") return r.q < 0.05 ? (r.pct < 0 ? "down" : "up") : "ns";
  return !r.m6a ? "ns" : r.m6a === 1 ? "m1" : r.m6a <= 3 ? "m2" : "m3";
}

// volcano: points on a canvas, axes and labels in an SVG on top; the axes cover the bulk of the genes and the few
// beyond are drawn on the edge (values in the tooltip and the table)
function drawVolcano() {
  const box = $("#volcano");
  if (!land.rows || $("#view-landscape").hidden || !box.clientWidth) return;
  const tx = T(), rows = land.rows, species = land.id.split("_")[0];
  const W = box.clientWidth, H = Math.round(Math.min(500, Math.max(320, W * 0.48))), M = { l: 52, r: 18, t: 18, b: 46 };
  const ys = rows.map((r) => -Math.log10(Math.max(r.p, 1e-300)));
  const xt = niceTicks(Math.max(1, quantileOf(rows.map((r) => Math.abs(r.pct)), 0.995) * 1.1)), xlim = xt.at(-1);
  const ytop = Math.max(3, quantileOf(ys, 0.999) * 1.08), ystep = niceStep(ytop / 5);
  const yt = Array.from({ length: Math.ceil(ytop / ystep) + 1 }, (_, i) => +(i * ystep).toPrecision(6)), ylim = yt.at(-1);
  const sx = (v) => M.l + ((Math.max(-xlim, Math.min(xlim, v)) + xlim) / (2 * xlim)) * (W - M.l - M.r);
  const sy = (v) => M.t + (1 - Math.min(v, ylim) / ylim) * (H - M.t - M.b);
  box.style.height = `${H}px`;
  const canvas = box.querySelector("canvas") || box.appendChild(document.createElement("canvas"));
  const dpr = devicePixelRatio || 1;
  Object.assign(canvas, { width: W * dpr, height: H * dpr });
  Object.assign(canvas.style, { width: `${W}px`, height: `${H}px` });
  const g = canvas.getContext("2d");
  g.scale(dpr, dpr);
  const css = getComputedStyle(document.documentElement), color = (name) => css.getPropertyValue(name).trim();
  const palette = { ns: color("--ns"), down: color("--down"), up: color("--up"), m1: color("--m6a-1"), m2: color("--m6a-2"), m3: color("--m6a-3") };
  // a point beyond an axis is a triangle on the edge, pointing out
  const pts = rows.map((r, i) => ({ x: sx(r.pct), y: sy(ys[i]), r, c: classOf(r),
    out: ys[i] > ylim ? "up" : Math.abs(r.pct) > xlim ? (r.pct > 0 ? "right" : "left") : null }));
  const mark = (p, size) => {
    const t = { up: [[0, -1.6], [-1.3, 0.8], [1.3, 0.8]], right: [[1.6, 0], [-0.8, -1.3], [-0.8, 1.3]],
      left: [[-1.6, 0], [0.8, -1.3], [0.8, 1.3]] }[p.out];
    if (!t) { g.moveTo(p.x + size, p.y); g.arc(p.x, p.y, size, 0, 2 * Math.PI); return; }
    g.moveTo(p.x + t[0][0] * size * 1.5, p.y + t[0][1] * size * 1.5);
    for (const [dx, dy] of t.slice(1)) g.lineTo(p.x + dx * size * 1.5, p.y + dy * size * 1.5);
    g.closePath();
  };
  for (const layer of ["ns", "m1", "m2", "m3", "down", "up"]) {
    g.fillStyle = palette[layer];
    g.globalAlpha = layer === "ns" ? 0.55 : 0.85;
    g.beginPath();
    for (const p of pts) if (p.c === layer) mark(p, layer === "ns" ? 1.9 : 2.3);
    g.fill();
  }
  land.grid = new Map();
  for (const p of pts) {
    const key = `${Math.floor(p.x / 10)},${Math.floor(p.y / 10)}`;
    land.grid.set(key, [...(land.grid.get(key) || []), p]);
  }
  const hits = rows.filter((r) => r.q < 0.05);
  const svg = [];
  for (const t of xt) svg.push(`<line class="${t === 0 ? "zero" : "grid"}" x1="${sx(t)}" x2="${sx(t)}" y1="${M.t}" y2="${H - M.b}"/>`,
    `<text class="tick" x="${sx(t)}" y="${H - M.b + 16}" text-anchor="middle">${t > 0 ? "+" : ""}${minus(t)}</text>`);
  for (const t of yt) svg.push(`<line class="grid" x1="${M.l}" x2="${W - M.r}" y1="${sy(t)}" y2="${sy(t)}"/>`,
    `<text class="tick" x="${M.l - 8}" y="${sy(t) + 4}" text-anchor="end">${t}</text>`);
  if (hits.length) {
    const yq = sy(-Math.log10(Math.max(...hits.map((r) => r.p))));
    svg.push(`<line class="threshold" x1="${M.l}" x2="${W - M.r}" y1="${yq}" y2="${yq}"/>`,
      `<text class="tick" x="${M.l + 6}" y="${yq - 5}">${tx.qLine}</text>`);
  }
  // labels of the five most significant genes of each direction, skipped where they would overlap a placed label
  const placed = [];
  for (const dir of [-1, 1]) {
    for (const r of hits.filter((h) => Math.sign(h.pct) === dir).slice(0, 5)) {
      const x = sx(r.pct), y = Math.max(sy(-Math.log10(Math.max(r.p, 1e-300))) - 4, M.t + 10), w = r.symbol.length * 6.6;
      const right = x > W - 90, x0 = right ? x - 6 - w : x + 6;
      if (placed.some((b) => x0 < b.x + b.w + 4 && b.x < x0 + w + 4 && Math.abs(b.y - y) < 13)) continue;
      placed.push({ x: x0, y, w });
      svg.push(`<text class="label" x="${right ? x - 6 : x + 6}" y="${y}" text-anchor="${right ? "end" : "start"}">${esc(r.symbol)}</text>`);
    }
  }
  svg.push(`<text x="${M.l + (W - M.l - M.r) / 2}" y="${H - 8}" text-anchor="middle">${tx.xAxis[species]}</text>`,
    `<text transform="translate(14 ${M.t + (H - M.t - M.b) / 2}) rotate(-90)" text-anchor="middle">${tx.yAxis}</text>`,
    `<circle id="v-ring" class="ring" r="0"/>`);
  const old = box.querySelector("svg");
  if (old) old.remove();
  box.insertAdjacentHTML("beforeend", `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">${svg.join("")}</svg>`);
  const legend = land.color === "dir"
    ? [["down", tx.sigDown], ["up", tx.sigUp], ["ns", tx.ns]]
    : [["ns", "0"], ["m1", "1"], ["m2", "2–3"], ["m3", "≥ 4"]];
  $("#volcano-legend").innerHTML = (land.color === "m6a" ? `<em class="legend-text">${tx.m6aClass[species]}:</em>` : "") +
    legend.map(([c, label]) => `<span><span class="dot ${c}"></span>${label}</span>`).join("") + `<em class="legend-text muted">${tx.beyond}</em>`;
}

function niceStep(raw) {
  const p = 10 ** Math.floor(Math.log10(raw)), f = raw / p;
  return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 2.5 ? 2.5 : f <= 5 ? 5 : 10) * p;
}

function hoverVolcano(e) {
  const box = $("#volcano"), rect = box.getBoundingClientRect(), x = e.clientX - rect.left, y = e.clientY - rect.top;
  let best = null, dist = 64;
  for (let i = Math.floor(x / 10) - 1; i <= Math.floor(x / 10) + 1; i++)
    for (let j = Math.floor(y / 10) - 1; j <= Math.floor(y / 10) + 1; j++)
      for (const p of land.grid?.get(`${i},${j}`) || []) {
        const d = (p.x - x) ** 2 + (p.y - y) ** 2;
        if (d < dist) { dist = d; best = p; }
      }
  const tip = $("#tip"), ring = $("#v-ring");
  box.style.cursor = best ? "pointer" : "default";
  if (!best) { tip.hidden = true; ring?.setAttribute("r", 0); return null; }
  const r = best.r, tx = T(), s = land.id.split("_")[0];
  ring.setAttribute("cx", best.x); ring.setAttribute("cy", best.y); ring.setAttribute("r", 5);
  tip.innerHTML = `<b>${esc(r.symbol)}</b><div>${contextName(land.id)}</div>
    <div class="line" style="margin-top:6px"><span>${tx.tipCombined}</span><b>${pct(r.pct)}%</b></div>
    <div class="line"><span>${tx.ci}</span><span>${pct(r.lo)} ${tx.to} ${pct(r.hi)}%</span></div>
    <div class="line"><span>P</span><span>${fmtP(r.p)}</span></div>
    <div class="line"><span>q</span><span>${fmtP(r.q)}</span></div>
    <div class="line"><span>${tx.tipSame}</span><span>${tx.tipCohorts(r.same, r.k)}</span></div>
    <div class="line"><span>I²</span><span>${r.i2 ?? "–"}${r.i2 == null ? "" : "%"}</span></div>
    <div class="line"><span>${tx.m6aClass[s]}</span><span>${r.m6a ?? "–"}</span></div>`;
  tip.hidden = false;
  tip.style.left = `${Math.min(e.clientX + 16, innerWidth - tip.offsetWidth - 8)}px`;
  tip.style.top = `${e.clientY + 16 + tip.offsetHeight > innerHeight ? e.clientY - tip.offsetHeight - 12 : e.clientY + 16}px`;
  return r;
}

function rankedRows() {
  const rows = land.rows.filter((r) => (!land.dir || (land.dir === "down" ? r.pct < 0 : r.pct > 0)) && Math.abs(r.pct) >= land.min);
  return land.sort === "effect" ? [...rows].sort((a, b) => Math.abs(b.pct) - Math.abs(a.pct)) : rows; // stored by P
}

function renderRank() {
  if (!land.rows) return;
  const tx = T(), rows = rankedRows(), shownRows = rows.slice(0, land.n);
  const head = tx.rankCols.map((c, i) => `<th class="${i >= 2 ? "num" : ""}" ${i === 6 ? `title="${tx.sameTitle}"` : ""}>${c}</th>`).join("");
  const body = shownRows.map((r, i) => `<tr>
    <td class="muted">${i + 1}</td>
    <td><a class="gene" href="#gene/${encodeURIComponent(r.symbol)}">${esc(r.symbol)}</a></td>
    <td class="num nowrap ${r.q < 0.05 ? "sig" : ""}">${pct(r.pct)}</td>
    <td class="num nowrap muted">${pct(r.lo)} ${tx.to} ${pct(r.hi)}</td>
    <td class="num nowrap">${fmtP(r.p)}</td>
    <td class="num nowrap">${fmtP(r.q)}</td>
    <td class="num">${r.same} / ${r.k}</td>
    <td class="num">${r.i2 == null ? "–" : `${r.i2}%`}</td>
    <td class="num">${r.m6a ?? "–"}</td></tr>`).join("");
  $("#rank-table").innerHTML = `<thead><tr>${head}</tr></thead><tbody>${body}</tbody>`;
  $("#rank-count").textContent = tx.shownOf(num(shownRows.length), num(rows.length));
  $("#l-more").hidden = shownRows.length >= rows.length;
}

function downloadContext() {
  const head = ["gene", "effect_pct", "ci_low_pct", "ci_high_pct", "p", "q", "n_cohorts", "n_same_direction", "i2", "m6a_sites"];
  const lines = land.rows.map((r) => [r.symbol, r.pct, r.lo, r.hi, r.p, r.q, r.k, r.same, r.i2 ?? "", r.m6a ?? ""].join(","));
  const a = Object.assign(document.createElement("a"), {
    href: URL.createObjectURL(new Blob([[head.join(","), ...lines].join("\n")], { type: "text/csv" })),
    download: `${land.id}_landscape.csv`,
  });
  a.click();
  URL.revokeObjectURL(a.href);
}

function renderContextTable() {
  const tx = T();
  $("#context-table").innerHTML = `<thead><tr>${tx.contextCols.map((c, i) => `<th class="${i === 2 ? "num" : ""}">${c}</th>`).join("")}</tr></thead>
    <tbody>${CONTEXTS.map((id) => meta.landscape.find((c) => c.id === id)).map((c) => `<tr><td class="nowrap">${contextName(c.id)}</td>
    <td>${c.cohorts.map(esc).join(", ")}</td><td class="num">${num(c.genes)}</td></tr>`).join("")}</tbody>`;
}

init();
