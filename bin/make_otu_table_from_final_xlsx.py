#!/usr/bin/env python3
import argparse
import os
import pandas as pd

def safe_taxid(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    # sometimes "123,456" -> keep first
    if "," in s:
        s = s.split(",", 1)[0]
    # sometimes floats from excel
    s = s.replace(".0", "")
    return int(s) if s.isdigit() else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", nargs="+", required=True, help="Input *.final.xlsx files")
    ap.add_argument("--out", required=True, help="Output OTU table .xlsx")
    ap.add_argument("--add_taxonomy", action="store_true", help="Add Taxonomy column using ete3")
    ap.add_argument("--barcode_from", choices=["filename_prefix"], default="filename_prefix",
                    help="How to extract barcode name (default: filename_prefix -> barcodeXX from basename)")
    args = ap.parse_args()

    otu = None

    for f in args.xlsx:
        base = os.path.basename(f)
        # expects barcode06...final.xlsx -> "barcode06"
        barcode = base.split(".")[0].split("_")[0]

        df = pd.read_excel(f, engine="openpyxl")

        required = {"taxid", "reads_in_cluster"}
        miss = required - set(df.columns)
        if miss:
            raise SystemExit(f"{f}: missing columns {sorted(miss)}")

        df["taxid_clean"] = df["taxid"].apply(safe_taxid)
        df = df.dropna(subset=["taxid_clean"]).copy()
        df["taxid_clean"] = df["taxid_clean"].astype("int64")

        t = df.groupby("taxid_clean", as_index=False)["reads_in_cluster"].sum()
        t.columns = ["taxid", barcode]

        otu = t if otu is None else pd.merge(otu, t, on="taxid", how="outer")

    otu.fillna(0, inplace=True)

    if args.add_taxonomy:
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()

        relevant = {"kingdom","phylum","class","order","family","genus","species"}
        ranks_order = ["kingdom","phylum","class","order","family","genus","species"]

        cache = {}
        def get_tax(taxid: int) -> str:
            if taxid in cache:
                return cache[taxid]
            try:
                lineage = ncbi.get_lineage(int(taxid))
                names = ncbi.get_taxid_translator(lineage)
                rank_map = ncbi.get_rank(lineage)
                lineage_dict = {rank_map[t]: names.get(t, "Unknown") for t in lineage if rank_map.get(t) in relevant}
                tax = "; ".join([lineage_dict.get(r, "Unknown") for r in ranks_order])
            except Exception:
                tax = "Unknown"
            cache[taxid] = tax
            return tax

        otu["Taxonomy"] = otu["taxid"].apply(get_tax)

    otu.to_excel(args.out, index=False, engine="openpyxl")

if __name__ == "__main__":
    main()
