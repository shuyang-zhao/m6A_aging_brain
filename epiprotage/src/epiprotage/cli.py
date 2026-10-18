"""Command line: epiprotage score FILE --species mouse|human [--out scores.csv]"""

import argparse
import sys

from . import __version__, score, weights_version


def main(argv=None):
    parser = argparse.ArgumentParser(prog="epiprotage", description="Zero-shot clock of brain proteome aging.")
    parser.add_argument(
        "--version", action="version", version=f"epiprotage {__version__} (weights {weights_version()['date']})"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    cmd = sub.add_parser("score", help="score the samples of a proteome table (genes x samples, log abundances)")
    cmd.add_argument(
        "file", help="CSV/TSV: first column gene symbol or Ensembl id, one column per sample; a row 'age' is optional"
    )
    cmd.add_argument("--species", required=True, choices=["mouse", "human"])
    cmd.add_argument("--ids", default="auto", choices=["auto", "symbol", "ensembl"])
    cmd.add_argument("--out", help="output CSV (default: stdout)")
    args = parser.parse_args(argv)

    result = score(args.file, args.species, ids=args.ids)
    result.to_csv(args.out or sys.stdout, index=False)
    print(
        f"{result.attrs['genes_matched']} of {result.attrs['genes_input']} rows matched ({result.attrs['ids']} ids)",
        file=sys.stderr,
    )
    if "spearman_rho" in result.attrs:
        rho, p = result.attrs["spearman_rho"], result.attrs["spearman_p"]
        print(f"Spearman rho with age = {rho:.3f} (P = {p:.2g}, n = {len(result)})", file=sys.stderr)
