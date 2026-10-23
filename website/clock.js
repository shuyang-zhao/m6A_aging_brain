// m6Brain clock of the Clock page, computed in the browser (nothing is uploaded). Weights: data/clock_<species>.json
// (columns gene_id, symbol, weight and the intercept; pipeline/62_website.py, from script 55). Predicted age of a
// sample: each protein is standardised across the uploaded samples (mean 0, SD 1 over the quantified values, missing
// = 0); the age is the intercept plus the weighted sum over the clock proteins. Months in mouse, years in human.
// The same calculation as predict_age of the m6brain package; also loaded by dev/verification/clock_check.js under node.

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

// predicted age of every sample and the number of clock proteins quantified in it
function clockAges(matrix, weights, intercept) {
  const n = matrix[0].length, predicted = new Float64Array(n).fill(intercept), used = new Int32Array(n);
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
    for (let j = 0; j < n; j++) if (!Number.isNaN(x[j])) { predicted[j] += (weights[g] * (x[j] - mean)) / sd; used[j]++; }
  });
  return { predicted: Array.from(predicted), used: Array.from(used) };
}

// Pearson r over the samples with both values
function pearson(x, y) {
  const keep = x.map((_, i) => i).filter((i) => !Number.isNaN(x[i]) && !Number.isNaN(y[i])), n = keep.length;
  const mx = keep.reduce((s, i) => s + x[i], 0) / n, my = keep.reduce((s, i) => s + y[i], 0) / n;
  let xy = 0, xx = 0, yy = 0;
  for (const i of keep) { xy += (x[i] - mx) * (y[i] - my); xx += (x[i] - mx) ** 2; yy += (y[i] - my) ** 2; }
  return { r: xy / Math.sqrt(xx * yy), n };
}

// age acceleration: residual of the predicted age on the age in a least-squares line over the samples with an age
function acceleration(predicted, age) {
  const keep = age.map((_, i) => i).filter((i) => !Number.isNaN(age[i])), n = keep.length;
  const ma = keep.reduce((s, i) => s + age[i], 0) / n, mp = keep.reduce((s, i) => s + predicted[i], 0) / n;
  let ap = 0, aa = 0;
  for (const i of keep) { ap += (age[i] - ma) * (predicted[i] - mp); aa += (age[i] - ma) ** 2; }
  return predicted.map((p, i) => p - (mp + (ap / aa) * (age[i] - ma))); // NaN where the age is missing
}

if (typeof module !== "undefined") module.exports = { parseTable, matchWeights, clockAges, pearson, acceleration };
