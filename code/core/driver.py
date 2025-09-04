#!/usr/bin/env python3

import argparse
import csv
import os
from datetime import datetime


OUT_TEX_DIR = os.path.join("out", "tex")
OUT_CSV_DIR = os.path.join("out", "csv")
OUT_FIG_DIR = os.path.join("out", "figs")
RUN_LOG = os.path.join("out", "run.log")


def ensure_dirs():
    os.makedirs(OUT_TEX_DIR, exist_ok=True)
    os.makedirs(OUT_CSV_DIR, exist_ok=True)
    os.makedirs(OUT_FIG_DIR, exist_ok=True)


def read_table_csv(path, required_columns):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for col in required_columns:
            if col not in reader.fieldnames:
                raise ValueError(f"Missing column '{col}' in {path}. Found {reader.fieldnames}")
        rows = list(reader)
    return rows


def write_tex_table_leptons(rows, out_path):
    # Expected columns: Ratio, Predicted, Sigma, PDG, DeviationPPM
    header = (
        "\\begin{center}\n"
        "\\begin{tabular}{lcccc}\n"
        "\\toprule\n"
        "Ratio & Predicted & $\\pm$ & PDG & Deviation (ppm) \\\\ \n"
        "\\midrule\n"
    )
    body_lines = []
    for r in rows:
        body_lines.append(
            f"{r['Ratio']} & {r['Predicted']} & $\\pm$ {r['Sigma']} & {r['PDG']} & {r['DeviationPPM']} \\\\"  # noqa: E501
        )
    footer = (
        "\\bottomrule\n"
        "\\end{tabular}\n"
        "\\end{center}\n"
    )
    with open(out_path, "w") as f:
        f.write(header)
        f.write("\n".join(body_lines) + "\n")
        f.write(footer)


def write_tex_table_quarks(rows, out_path, mz_variant=False):
    # Expected columns: Ratio, Predicted, Sigma, PDGAtScale, Deviation
    title = "Ratio @ $M_Z$" if mz_variant else "Ratio"
    header = (
        "\\begin{center}\n"
        "\\begin{tabular}{lcccc}\n"
        "\\toprule\n"
        f"{title} & Predicted & $\\pm$ & PDG (scheme) & Deviation \\\\ \n"
        "\\midrule\n"
    )
    body_lines = []
    for r in rows:
        body_lines.append(
            f"{r['Ratio']} & {r['Predicted']} & $\\pm$ {r['Sigma']} & {r['PDGAtScale']} & {r['Deviation']} \\\\"  # noqa: E501
        )
    footer = (
        "\\bottomrule\n"
        "\\end{tabular}\n"
        "\\end{center}\n"
    )
    with open(out_path, "w") as f:
        f.write(header)
        f.write("\n".join(body_lines) + "\n")
        f.write(footer)


def write_run_log(note):
    ensure_dirs()
    with open(RUN_LOG, "a") as f:
        f.write(f"[{datetime.utcnow().isoformat()}Z] {note}\n")


def build_tables(args):
    ensure_dirs()

    # Leptons
    leptons_csv = args.leptons or os.path.join("data", "leptons.csv")
    if os.path.exists(leptons_csv):
        lep_rows = read_table_csv(
            leptons_csv, ["Ratio", "Predicted", "Sigma", "PDG", "DeviationPPM"]
        )
        write_tex_table_leptons(lep_rows, os.path.join(OUT_TEX_DIR, "leptons_table.tex"))
        write_run_log(f"Wrote {OUT_TEX_DIR}/leptons_table.tex from {leptons_csv}")
    else:
        raise FileNotFoundError(f"Leptons CSV not found: {leptons_csv}")

    # Quarks at mu*
    q_mu_csv = args.quarks_mu or os.path.join("data", "quarks_muStar.csv")
    if os.path.exists(q_mu_csv):
        qmu_rows = read_table_csv(
            q_mu_csv, ["Ratio", "Predicted", "Sigma", "PDGAtScale", "Deviation"]
        )
        write_tex_table_quarks(qmu_rows, os.path.join(OUT_TEX_DIR, "quarks_muStar.tex"))
        write_run_log(f"Wrote {OUT_TEX_DIR}/quarks_muStar.tex from {q_mu_csv}")
    elif args.both_scales:
        raise FileNotFoundError(f"Quarks@mu* CSV not found: {q_mu_csv}")

    # Quarks at MZ
    q_mz_csv = args.quarks_mz or os.path.join("data", "quarks_MZ.csv")
    if os.path.exists(q_mz_csv):
        qmz_rows = read_table_csv(
            q_mz_csv, ["Ratio", "Predicted", "Sigma", "PDGAtScale", "Deviation"]
        )
        write_tex_table_quarks(
            qmz_rows, os.path.join(OUT_TEX_DIR, "quarks_MZ.tex"), mz_variant=True
        )
        write_run_log(f"Wrote {OUT_TEX_DIR}/quarks_MZ.tex from {q_mz_csv}")
    elif args.both_scales:
        raise FileNotFoundError(f"Quarks@MZ CSV not found: {q_mz_csv}")


def main():
    parser = argparse.ArgumentParser(description="Build tables/figures for the manuscript")
    parser.add_argument("--tables", action="store_true", help="Write TeX tables to out/tex/")
    parser.add_argument("--figs", action="store_true", help="(Reserved) Build figures to out/figs/")
    parser.add_argument("--both-scales", action="store_true", help="Emit quark tables at mu* and MZ")
    parser.add_argument("--leptons", help="Path to leptons.csv")
    parser.add_argument("--quarks-mu", help="Path to quarks_muStar.csv")
    parser.add_argument("--quarks-mz", help="Path to quarks_MZ.csv")
    args = parser.parse_args()

    if args.tables:
        build_tables(args)
    else:
        print("No action selected. Use --tables to generate TeX tables.")


if __name__ == "__main__":
    main()


