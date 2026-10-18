// Proteome clock of the Clock page, computed in the browser (nothing is uploaded). Weights: data/clock_<species>.json
// (columns gene_id, symbol, weight; pipeline/62_website.py, from script 52). Score of a sample: each protein is
// standardised across the uploaded samples (mean 0, SD 1 over the quantified values, missing = 0), the score is the
// weighted sum over the genes with a weight, then standardised across the samples. The same calculation as
// donor_scores of pipeline/52_proteome_clock.py; also loaded by dev/verification/clock_check.js under node.

// delimited text to rows of cells; quoted cells may hold the delimiter and doubled quotes
function parseDelimited(text, delim) {
  const rows = [];
  let row = [], cell = "", quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c !== '"') cell += c;
      else if (text[i + 1] === '"') { cell += '"'; i++; } else quoted = false;
    } else if (c === '"' && cell === "") quoted = true;
    else if (c === delim) { row.push(cell); cell = ""; }
    else if (c === "\n" || c === "\r") {
      if (c === "\r" && text[i + 1] === "\n") i++;
      row.push(cell); cell = "";
      if (row.length > 1 || row[0] !== "") rows.push(row);
      row = [];
    } else cell += c;
  }
  row.push(cell);
  if (row.length > 1 || row[0] !== "") rows.push(row);
  return rows;
}

// first column: gene symbol or Ensembl gene id (version and the ids after a ";" dropped), other columns: samples; a row
// whose first cell is "age" holds the ages. Empty and non-numeric cells (NA, NaN, Inf) are missing.
function parseTable(text) {
  text = text.replace(/^﻿/, "");
  const head = text.slice(0, text.search(/[\r\n]|$/));
  const rows = parseDelimited(text, head.includes("\t") ? "\t" : ",");
  const samples = rows[0].slice(1).map((s) => s.trim());
  const value = (x) => { const v = x.trim() === "" ? NaN : Number(x); return Number.isFinite(v) ? v : NaN; }; // NA, NaN, Inf: missing
  let age = null;
  const ids = [], matrix = [];
  for (const r of rows.slice(1)) {
    const id = r[0].trim();
    const numbers = Float64Array.from(samples, (_, j) => value(r[j + 1] ?? ""));
    if (id.toLowerCase() === "age") age = Array.from(numbers);
    else { ids.push(id); matrix.push(numbers); }
  }
  if (age?.every(Number.isNaN)) age = null;
  return { samples, ids, matrix, age };
}

const ENSEMBL = /^ENS(MUS)?G\d{6,}/i;
const normalise = (id, byId) => { const first = id.split(";")[0].trim(); return (byId ? first.replace(/\.\d+$/, "") : first).toUpperCase(); };

// weight of each uploaded row (NaN where the gene has none); ids are Ensembl gene ids when most rows look like them,
// else symbols (case-insensitive). A gene present in several rows takes the first.
function matchWeights(ids, clock) {
  const byId = ids.filter((id) => ENSEMBL.test(id)).length > ids.length / 2;
  const weight = new Map(clock.weight.map((w, i) => [byId ? clock.gene_id[i].toUpperCase() : clock.symbol[i].toUpperCase(), w]));
  const seen = new Set();
  let duplicates = 0;
  const w = ids.map((id) => {
    const key = normalise(id, byId);
    if (!weight.has(key)) return NaN;
    if (seen.has(key)) { duplicates++; return NaN; }
    seen.add(key);
    return weight.get(key);
  });
  return { weights: w, byId, matched: seen.size, duplicates };
}

// score of every sample and the number of weighted proteins quantified in it
function clockScores(matrix, weights) {
  const n = matrix[0].length, raw = new Float64Array(n), used = new Int32Array(n);
  matrix.forEach((x, g) => {
    if (Number.isNaN(weights[g])) return;
    let k = 0, sum = 0;
    for (const v of x) if (!Number.isNaN(v)) { k++; sum += v; }
    if (k < 2) return;
    const mean = sum / k;
    let ss = 0;
    for (const v of x) if (!Number.isNaN(v)) ss += (v - mean) ** 2;
    const sd = Math.sqrt(ss / (k - 1));
    if (!(sd > 0)) return;
    for (let j = 0; j < n; j++) if (!Number.isNaN(x[j])) { raw[j] += (weights[g] * (x[j] - mean)) / sd; used[j]++; }
  });
  const mean = raw.reduce((a, b) => a + b, 0) / n;
  const sd = Math.sqrt(raw.reduce((a, b) => a + (b - mean) ** 2, 0) / (n - 1));
  return { scores: Array.from(raw, (v) => (v - mean) / sd), used: Array.from(used) };
}

function ranks(x) {
  const order = x.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]), r = new Array(x.length);
  for (let i = 0; i < order.length;) {
    let j = i;
    while (j + 1 < order.length && order[j + 1][0] === order[i][0]) j++;
    for (let k = i; k <= j; k++) r[order[k][1]] = (i + j) / 2 + 1; // ties share the mean rank
    i = j + 1;
  }
  return r;
}

// Spearman rho over the samples with both values
function spearman(x, y) {
  const keep = x.map((_, i) => i).filter((i) => !Number.isNaN(x[i]) && !Number.isNaN(y[i]));
  const a = ranks(keep.map((i) => x[i])), b = ranks(keep.map((i) => y[i])), n = keep.length;
  const ma = a.reduce((s, v) => s + v, 0) / n, mb = b.reduce((s, v) => s + v, 0) / n;
  let ab = 0, aa = 0, bb = 0;
  for (let i = 0; i < n; i++) { ab += (a[i] - ma) * (b[i] - mb); aa += (a[i] - ma) ** 2; bb += (b[i] - mb) ** 2; }
  return { rho: ab / Math.sqrt(aa * bb), n };
}

if (typeof module !== "undefined") module.exports = { parseTable, matchWeights, clockScores, spearman };
