#!/usr/bin/env python3
import argparse
import logging
import os

import pandas as pd
import taxopy


def setup_logging(log_file: str):
    logging.basicConfig(
        filename=log_file,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s\t%(levelname)s\t%(message)s",
    )
    # Also print to stdout (Nextflow captures this too)
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger("").addHandler(console)


def _normalize_rank(rank: str | None) -> str | None:
    """
    Treat anything below species / undefined as 'species' for the purpose of 'move up one rank'.
    """
    if rank is None:
        return "species"
    r = str(rank).strip().lower()
    below_species = {"subspecies", "strain", "isolate", "no rank", "forma", "varietas", "serovar"}
    if r in below_species:
        return "species"
    return r


def _target_rank_move_up_one(cur_rank: str | None) -> str | None:
    """
    Move up exactly one rank according to ladder:
      species -> genus -> family -> order
    If already at order (or higher), returns None (no change).
    """
    cur = _normalize_rank(cur_rank)
    move_up = {
        "species": "genus",
        "genus": "family",
        "family": "order",
        "order": None,
    }
    return move_up.get(cur, None)


def _find_rank_name_in_lineage(taxon: taxopy.Taxon, wanted_rank: str) -> str | None:
    for rank, name in taxon.ranked_name_lineage:
        if rank == wanted_rank:
            return name
    return None


def update_sciname_and_taxid_low_ident(
    df: pd.DataFrame,
    taxdb: taxopy.TaxDb,
    threshold: float,
    id_col: str,
):
    # Force numeric per_ident
    df["_per_ident_num"] = pd.to_numeric(df["per_ident"], errors="coerce")

    changed = 0
    for idx, row in df.iterrows():
        per_ident = row["_per_ident_num"]
        if pd.isna(per_ident) or per_ident >= threshold:
            continue  # only touch < threshold

        consensus_id = str(row.get(id_col, f"row{idx+1}"))
        old_taxid = row.get("taxid", None)
        old_sciname = row.get("sciname", None)

        if pd.isna(old_taxid) or not str(old_taxid).isdigit():
            logging.warning(f"SKIP\t{consensus_id}\tper_ident={per_ident}\tinvalid_taxid={old_taxid}")
            continue

        try:
            taxon = taxopy.Taxon(int(old_taxid), taxdb)

            old_rank = getattr(taxon, "rank", None)
            target_rank = _target_rank_move_up_one(old_rank)

            if not target_rank:
                logging.info(
                    f"NOCHANGE\t{consensus_id}\tper_ident={per_ident}\told_rank={old_rank}\told_taxid={old_taxid}\t(reason=no_higher_rank_in_ladder)"
                )
                continue

            target_name = _find_rank_name_in_lineage(taxon, target_rank)
            if not target_name:
                logging.info(
                    f"NOCHANGE\t{consensus_id}\tper_ident={per_ident}\told_rank={old_rank}\ttarget_rank={target_rank}\told_taxid={old_taxid}\t(reason=target_rank_missing_in_lineage)"
                )
                continue

            target_taxid = taxopy.taxid_from_name(target_name, taxdb)
            target_taxid = target_taxid[0] if target_taxid else None

            # Apply change
            df.at[idx, "sciname"] = target_name
            df.at[idx, "taxid"] = target_taxid
            changed += 1

            logging.info(
                "CHANGE\t%s\tper_ident=%.3f\told_rank=%s\told_sciname=%s\told_taxid=%s\tnew_rank=%s\tnew_sciname=%s\tnew_taxid=%s",
                consensus_id,
                float(per_ident),
                old_rank,
                old_sciname,
                old_taxid,
                target_rank,
                target_name,
                target_taxid,
            )

        except Exception as e:
            logging.error(f"ERROR\t{consensus_id}\tper_ident={per_ident}\ttaxid={old_taxid}\t{e}")

    df.drop(columns=["_per_ident_num"], inplace=True)
    return changed


def main():
    ap = argparse.ArgumentParser(
        description=(
            "For per_ident below threshold: move taxonomy up ONE rank using ladder "
            "species->genus->family->order, updating sciname/taxid (taxopy)."
        )
    )
    ap.add_argument("--xlsx", required=True, help="Input XLSX (filtered)")
    ap.add_argument("--out", required=True, help="Output XLSX")
    ap.add_argument("--log", required=True, help="Log file path")
    ap.add_argument("--threshold", type=float, default=97.0, help="Only modify rows where per_ident < threshold (default 97.0)")
    ap.add_argument("--id_col", default="consensus_id", help="Consensus ID column name (default consensus_id)")
    ap.add_argument("--backup", action="store_true", help="Create _backup.xlsx next to input before writing output")
    args = ap.parse_args()

    setup_logging(args.log)
    logging.info(f"START\tinput={args.xlsx}\tout={args.out}\tthreshold={args.threshold}\tid_col={args.id_col}")

    df = pd.read_excel(args.xlsx, engine="openpyxl")

    required = {"sciname", "taxid", "per_ident", args.id_col}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns {sorted(missing)} in {args.xlsx}")

    if args.backup:
        backup_file = args.xlsx.replace(".xlsx", "_backup.xlsx")
        os.rename(args.xlsx, backup_file)
        logging.info(f"BACKUP\t{backup_file}")

    taxdb = taxopy.TaxDb()
    changed = update_sciname_and_taxid_low_ident(df, taxdb, args.threshold, args.id_col)

    df.to_excel(args.out, index=False, engine="openpyxl")
    logging.info(f"DONE\tchanged_rows={changed}")


if __name__ == "__main__":
    main()
