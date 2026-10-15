# Additional files 3 and 4 of the article as Excel workbooks
# supplementary_data_1.xlsx (Additional file 3): README, the manifest of the raw files of data_sources.tsv with size and
# MD5, one sheet per result table of scripts 10-38 and 50 (the per-gene age effects of script 10 are in workbook 2), the
# supplementary tables (ST1-ST7) and the source data (SD)
# supplementary_data_2.xlsx (Additional file 4): README and the per-gene age effects of every unit in the three layers
# (script 10)
# the README of each workbook comes from docs/supplementary_data_columns.tsv (one row per sheet, then one per column;
# sheet "*" holds the columns shared by all sheets, "* source data" those of the source data); a column without a
# description, or a description without its column, stops the script
# reads data_sources.tsv, the raw files it lists, output/10-50, source_data/ and the column table
# writes output/61_supplementary_data/
import hashlib
import os
import time
from fnmatch import fnmatch
from glob import glob
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed
from openpyxl import Workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

import utils

OUT = utils.OUT / "61_supplementary_data"
RAW = utils.DATA / "raw"
WORKERS = 16
CHUNK = 100_000  # rows converted and written at a time
BOLD = Font(bold=True)
WRAP = Alignment(wrap_text=True, vertical="top")
# tables that are not sheets: permuted estimates (too large), the reader CLIP targets (input of 38 and 50), the model
# inputs and GPU runs of 50
NOT_SHEETS = [
    "12_permutation_fdr/null_*",
    "40_target_sets/*",
    "50_transcript_ageing_model/inputs",
    "50_transcript_ageing_model/runs",
]

described = utils.read(utils.ROOT / "docs" / "supplementary_data_columns.tsv", keep_default_na=False, dtype=str)
entries = described[described.column == ""]  # one row per sheet, in workbook order
columns_of = {sheet: rows.set_index("column")[["meaning", "unit"]] for sheet, rows in described.groupby("sheet")}


def readme_rows(book, columns):
    # a row per sheet and, below it, a row per column; every column of every table must be described, and every
    # description must belong to a column
    rows = []
    for entry in entries[entries.workbook == book].itertuples():
        script, source = f"pipeline/{entry.file.split('/')[0]}.py", f"output/{entry.file}"
        if entry.file == "data_sources.tsv":
            script, source = "pipeline/61_supplementary_data.py", entry.file
        if entry.file.startswith("source_data/"):
            script, source = "pipeline/60_source_data.py", entry.file
        note = " Values rounded to 4 significant digits." if entry.file.startswith("source_data/") else ""
        rows.append((entry.sheet, None, entry.meaning + note, None, script, source))
        if entry.in_workbook == "no":
            continue
        own = columns_of[entry.sheet].drop(index="")
        shared = [columns_of["*"]] + ([columns_of["* source data"]] if entry.file.startswith("source_data/") else [])
        known = pd.concat([own, *shared[::-1]])
        known = known[~known.index.duplicated()]  # the sheet's own meaning first, then source data, then shared
        missing = [column for column in columns[entry.sheet] if column not in known.index]
        extra = [column for column in own.index if column not in columns[entry.sheet]]
        assert not missing and not extra, f"{entry.sheet}: not described {missing}, described but absent {extra}"
        rows += [(entry.sheet, column, *known.loc[column], None, None) for column in columns[entry.sheet]]
    return rows


def md5(path):
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "md5").hexdigest()


def files_of(path):
    # the files a path of data_sources.tsv stands for, under data/raw (or, for a file of the repository, the repository)
    for base in (RAW, utils.ROOT):
        target = base / path
        if target.is_file():
            return [target]
        if target.is_dir():
            return sorted(p for p in target.rglob("*") if p.is_file())
        found = sorted(Path(p) for p in glob(str(target)) if os.path.isfile(p))
        if found:
            return found
    return []


def read_table(file):
    path = utils.ROOT / file if file.startswith("source_data/") else utils.OUT / file
    return utils.read(path, keep_default_na=False, na_values=[""], low_memory=False)  # only empty fields are missing


def add_sheet(wb, name, table):
    # header row frozen, bold and filtered; numbers stay numbers, missing values stay empty; columns sized from the head
    ws = wb.create_sheet(name)
    ws.freeze_panes = "A2"
    header = []
    for i, column in enumerate(table, 1):
        lengths = table[column].head(1000).dropna().astype(str).str.len()
        text = min(max(len(column), lengths.max() if len(lengths) else 0), 60)
        numeric = pd.api.types.is_numeric_dtype(table[column])
        ws.column_dimensions[get_column_letter(i)].width = max(len(column), 10) + 2 if numeric else text + 2
        header.append(WriteOnlyCell(ws, value=column))
        header[-1].font = BOLD
    ws.append(header)
    for start in range(0, len(table), CHUNK):
        part = table.iloc[start : start + CHUNK]
        for row in part.astype(object).where(part.notna(), None).itertuples(index=False, name=None):
            ws.append(row)
    ws.auto_filter.ref = f"A1:{get_column_letter(table.shape[1])}{len(table) + 1}"


def add_readme(wb, rows):
    ws = wb.create_sheet("README")
    ws.freeze_panes = "A2"
    for letter, width in zip("ABCDEF", (32, 30, 110, 44, 40, 70), strict=True):
        ws.column_dimensions[letter].width = width

    def cells(values, bold):
        out = [WriteOnlyCell(ws, value=value) for value in values]
        for cell in out:
            cell.alignment = WRAP
            cell.font = BOLD if bold else Font()
        return out

    ws.append(cells(["Sheet", "Column", "Meaning", "Measurement unit", "Source script", "Output file"], True))
    for row in rows:
        ws.append(cells(row, row[1] is None))
    ws.auto_filter.ref = f"A1:F{len(rows) + 1}"


def workbook(book, load, columns):
    # README first, then a sheet for every entry with a table, each loaded when it is written
    path = OUT / f"supplementary_data_{book}.xlsx"
    wb = Workbook(write_only=True)
    add_readme(wb, readme_rows(book, columns))
    for entry in entries[(entries.workbook == book) & (entries.in_workbook == "yes")].itertuples():
        add_sheet(wb, entry.sheet, load(entry.sheet))
    wb.save(path)
    print(f"{path.name}: {path.stat().st_size / 1e6:.0f} MB, {time.time() - start:.0f} s", flush=True)


start = time.time()
OUT.mkdir(parents=True, exist_ok=True)
sources = utils.read(utils.ROOT / "data_sources.tsv", keep_default_na=False)
assert entries.sheet.is_unique and entries.sheet.str.len().max() <= 31, "sheet names: unique, 31 characters at most"

# every table of the scripts 10-50 and of source_data/ has a sheet or an entry
folders = [p.stem for p in (utils.ROOT / "pipeline").glob("[1-5]?_*.py") if (utils.OUT / p.stem).is_dir()]
found = {f"{folder}/{p.name}" for folder in folders for p in (utils.OUT / folder).iterdir()}
found |= {f"source_data/{p.name}" for p in (utils.ROOT / "source_data").glob("*.tsv")}
unlisted = [file for file in sorted(found - set(entries.file)) if not any(fnmatch(file, p) for p in NOT_SHEETS)]
assert not unlisted, f"tables without a sheet or an entry: {unlisted}"

# Supplementary Data 1
book1 = entries[(entries.workbook == "1") & (entries.in_workbook == "yes")]
tables = {e.sheet: read_table(e.file) for e in book1.itertuples() if e.file != "data_sources.tsv"}

# raw data manifest: a row per file of every source, MD5 of the distinct files in parallel
listing = {row.id: files_of(row.path) for row in sources.itertuples()}
missing = [source for source, files in listing.items() if not files]
assert not missing, f"sources without files: {missing}"
distinct = sorted({file for files in listing.values() for file in files}, key=lambda f: -f.stat().st_size)
checksum = dict(zip(distinct, Parallel(n_jobs=WORKERS, batch_size=1)(delayed(md5)(f) for f in distinct), strict=True))
rows = []
for row in sources.itertuples():
    head = [row.id, row.type, row.citation or None, row.doi or None, row.accession or None, row.url or None]
    for f in listing[row.id]:
        path = os.path.relpath(f, RAW if f.is_relative_to(RAW) else utils.ROOT)  # data/raw, or the repository
        rows.append([*head, path, f.stat().st_size, checksum[f]])
columns = "source id,type,citation,DOI,accession,URL,file path,size in bytes,MD5".split(",")
tables["Raw data manifest"] = pd.DataFrame(rows, columns=columns)
workbook("1", tables.get, {sheet: list(table) for sheet, table in tables.items()})

# Supplementary Data 2: one table at a time (the mRNA layer has about 890,000 rows)
files = dict(entries[entries.workbook == "2"][["sheet", "file"]].itertuples(index=False))
header = {sheet: list(utils.read(utils.OUT / file, nrows=0)) for sheet, file in files.items()}
workbook("2", lambda sheet: read_table(files[sheet]), header)
