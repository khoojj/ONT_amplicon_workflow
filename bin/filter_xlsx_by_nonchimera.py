#!/usr/bin/env python3
import argparse
import pandas as pd


def read_nonchimera_ids(fasta_path: str) -> set[str]:
    keep = set()
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            h = line[1:].strip()

            # Remove anything after first whitespace (if any)
            h = h.split()[0]

            # Remove anything after ';' (e.g. ;size=2481;)
            h = h.split(";")[0]

            # Keep only the part before the first underscore:
            # b06c0_f5bf... -> b06c0
            cid = h.split("_", 1)[0]

            if cid:
                keep.add(cid)
    return keep


def main():
    ap = argparse.ArgumentParser(
        description="Filter NanoCLUST XLSX summary to keep only consensus IDs present in a non-chimeric FASTA."
    )
    ap.add_argument("--xlsx", required=True, help="Input .xlsx (NanoCLUST summary)")
    ap.add_argument("--nonchimera_fa", required=True, help="FASTA produced by vsearch --nonchimeras")
    ap.add_argument("--out", required=True, help="Output filtered .xlsx")
    ap.add_argument("--id_col", default="consensus_id", help="Column containing consensus IDs (default: consensus_id)")
    args = ap.parse_args()

    keep = read_nonchimera_ids(args.nonchimera_fa)

    df = pd.read_excel(args.xlsx, engine="openpyxl")
    if args.id_col not in df.columns:
        raise SystemExit(f"Expected column '{args.id_col}' in {args.xlsx}, found: {list(df.columns)}")

    df_f = df[df[args.id_col].astype(str).isin(keep)].copy()
    df_f.to_excel(args.out, index=False, engine="openpyxl")


if __name__ == "__main__":
    main()
