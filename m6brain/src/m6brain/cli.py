"""Command line: m6brain age|score FILE --species mouse|human [--out result.csv]"""

import argparse
import sys

from . import __version__, predict_age, score, weights_version


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="m6brain",
        description="m6Brain clock of brain proteome age (months in mouse, years in human) and its label-free score.",
    )
    parser.add_argument(
        "--version", action="version", version=f"m6brain {__version__} (weights {weights_version()['date']})"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    for name, text in [
        ("age", "predict the age of the samples of a proteome table (genes x samples, log abundances)"),
        ("score", "label-free score of the samples (weights predicted from the transcript, no age fitted)"),
    ]:
        cmd = sub.add_parser(name, help=text)
        cmd.add_argument(
            "file", help="CSV/TSV: first column gene symbol or Ensembl id, one column per sample; optional row 'age'"
        )
        cmd.add_argument("--species", required=True, choices=["mouse", "human"])
        cmd.add_argument("--ids", default="auto", choices=["auto", "symbol", "ensembl"])
        cmd.add_argument("--out", help="output CSV (default: stdout)")
        if name == "age":
            cmd.add_argument("--controls", help="comma-separated sample names on which the age acceleration is fitted")
    args = parser.parse_args(argv)

    if args.command == "age":
        controls = args.controls.split(",") if args.controls else None
        result = predict_age(args.file, args.species, controls=controls, ids=args.ids)
    else:
        result = score(args.file, args.species, ids=args.ids)
    result.to_csv(args.out or sys.stdout, index=False)
    a = result.attrs
    total = f"{a['clock_genes']} clock proteins" if args.command == "age" else f"{a['genes_input']} rows"
    print(f"{a['genes_matched']} of {total} matched ({a['ids']} ids)", file=sys.stderr)
    if "pearson_r" in a:
        r, p = a["pearson_r"], a["pearson_p"]
        print(f"Pearson r with age = {r:.3f} (P = {p:.2g}, n = {len(result)})", file=sys.stderr)
