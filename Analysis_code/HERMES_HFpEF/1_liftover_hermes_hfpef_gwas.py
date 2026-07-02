#!/usr/bin/env python3
"""Lift HERMES HFpEF GWAS summary statistics from hg19 to hg38.

Input GWAS columns:
  rsID, chr, pos_b37, A1, A2, A1_beta, se, pval

Output keeps the original columns and adds:
  CHR37, POS37, CHR38, POS38, variant_hg38, rsid

Fast default:
  The input already contains rsID, so the script keeps that rsID by default and
  skips the extra 1000G rsID scan. Use --add-1000g-rsid only when hg38 rsID
  replacement is required.
"""

import argparse
import gzip
import os
from collections import defaultdict

import pandas as pd


DEFAULT_PAPER_DIR = os.environ.get(
    "PAPER_SMIXCAN_DIR",
    "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan",
)
DEFAULT_INPUT_FILE = os.path.join(
    DEFAULT_PAPER_DIR,
    "Heart",
    "HERMES_HFpEF",
    "FORMAT-METAL_Pheno4_EUR.tsv.gz",
)
DEFAULT_RSID_FILE = os.path.join(
    DEFAULT_PAPER_DIR,
    "Heart",
    "Data",
    "1000g_b38_snpIDs.txt",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Lift HERMES HFpEF hg19 GWAS to hg38."
    )
    parser.add_argument("--input-file", default=DEFAULT_INPUT_FILE)
    parser.add_argument("--rsid-file", default=DEFAULT_RSID_FILE)
    parser.add_argument("--output-file", default=None)
    parser.add_argument("--chunksize", type=int, default=1_000_000)
    parser.add_argument(
        "--add-1000g-rsid",
        action="store_true",
        help="Also scan 1000G hg38 rsID annotation and replace input rsID when available.",
    )
    return parser.parse_args()


def strip_compression_suffix(path):
    base = os.path.basename(path)
    for suffix in (".tsv.gz", ".tsv.bgz", ".bgz", ".gz", ".tsv", ".txt"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def default_output_path(input_file):
    output_dir = os.path.dirname(input_file)
    base = strip_compression_suffix(input_file)
    return os.path.join(output_dir, f"{base}_hg38_rsid.tsv.gz")


def make_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def normalize_hermes_chunk(chunk):
    """Standardize HERMES columns before liftover.

    HERMES reports A1_beta as the effect for A1. Downstream scripts keep A1 as
    the GWAS effect allele and A2 as the baseline allele.
    """
    required = ["rsID", "chr", "pos_b37", "A1", "A2", "A1_beta", "se", "pval"]
    missing = [col for col in required if col not in chunk.columns]
    if missing:
        raise ValueError(f"HERMES HFpEF GWAS is missing columns: {', '.join(missing)}")

    chunk = chunk.copy()
    chunk["CHR37"] = chunk["chr"].astype(str).str.replace("^chr", "", regex=True)
    chunk["POS37"] = pd.to_numeric(chunk["pos_b37"], errors="coerce").astype("Int64")
    chunk["A1"] = chunk["A1"].astype(str).str.upper()
    chunk["A2"] = chunk["A2"].astype(str).str.upper()
    chunk["A1_beta"] = pd.to_numeric(chunk["A1_beta"], errors="coerce")
    chunk["se"] = pd.to_numeric(chunk["se"], errors="coerce")
    chunk["pval"] = pd.to_numeric(chunk["pval"], errors="coerce")
    chunk["rsid"] = chunk["rsID"].astype(str)
    return chunk


def liftover_coordinate(lo, chrom, pos37):
    lifted = lo.convert_coordinate(f"chr{chrom}", int(pos37) - 1)
    if not lifted:
        return None, None
    chrom38 = str(lifted[0][0]).replace("chr", "")
    pos38 = int(lifted[0][1]) + 1
    return chrom38, pos38


def add_liftover_columns(chunk, lo, liftover_cache):
    # Lift each unique hg19 position once, then merge the hg38 coordinates back
    # to the full chunk. This avoids repeated pyliftover calls for duplicate
    # variants at the same position.
    positions = chunk[["CHR37", "POS37"]].drop_duplicates()

    map_items = []
    for row in positions.itertuples(index=False):
        key = (str(row.CHR37), int(row.POS37))
        if key not in liftover_cache:
            chrom38, pos38 = liftover_coordinate(lo, key[0], key[1])
            liftover_cache[key] = (chrom38, pos38)
        chrom38, pos38 = liftover_cache[key]
        map_items.append((key[0], key[1], chrom38, pos38))

    map_df = pd.DataFrame(map_items, columns=["CHR37", "POS37", "CHR38", "POS38"])
    map_df["CHR37"] = map_df["CHR37"].astype(str)
    map_df["POS37"] = map_df["POS37"].astype("int64")
    chunk["CHR37"] = chunk["CHR37"].astype(str)
    chunk["POS37"] = chunk["POS37"].astype("int64")
    chunk = chunk.merge(map_df, on=["CHR37", "POS37"], how="left")
    chunk = chunk.dropna(subset=["CHR38", "POS38"]).copy()
    chunk["POS38"] = chunk["POS38"].astype("int64")
    chunk["variant_hg38"] = (
        chunk["CHR38"].astype(str)
        + ":"
        + chunk["POS38"].astype(str)
        + ":"
        + chunk["A2"].astype(str)
        + ":"
        + chunk["A1"].astype(str)
    )
    return chunk


def write_lifted_gwas(input_file, output_file, chunksize, add_1000g_rsid):
    try:
        from pyliftover import LiftOver
    except ImportError as exc:
        raise SystemExit(
            "Missing dependency: pyliftover. Install it with:\n"
            "  python3 -m pip install pyliftover\n"
        ) from exc

    lo = LiftOver("hg19", "hg38")
    liftover_cache = {}
    needed_positions = defaultdict(set)
    total_rows = 0
    lifted_rows = 0
    wrote_header = False

    with gzip.open(output_file, "wt") as handle:
        for chunk in pd.read_csv(
            input_file,
            sep="\t",
            compression="gzip",
            chunksize=chunksize,
            low_memory=False,
        ):
            total_rows += len(chunk)
            chunk = normalize_hermes_chunk(chunk)
            chunk = chunk.dropna(subset=["POS37", "A1_beta", "se", "pval"]).copy()
            chunk = add_liftover_columns(chunk, lo, liftover_cache)

            if add_1000g_rsid:
                for chrom, positions in chunk.groupby("CHR38")["POS38"]:
                    needed_positions[str(chrom)].update(positions.astype(int).tolist())

            chunk.to_csv(
                handle,
                sep="\t",
                index=False,
                header=not wrote_header,
                na_rep="",
            )
            wrote_header = True
            lifted_rows += len(chunk)
            print(
                f"Lifted {lifted_rows:,} rows so far "
                f"({len(liftover_cache):,} unique positions cached)..."
            )

    return needed_positions, total_rows, lifted_rows


def build_rsid_lookup(rsid_file, needed_positions, chunksize):
    """Read only rsID annotation rows needed by lifted GWAS positions."""
    rsid_lookup = {}
    names = ["CHR38", "POS38", "rsid_from_1000g"]
    use_chromosomes = set(needed_positions)

    for chunk in pd.read_csv(
        rsid_file,
        sep="\t",
        names=names,
        dtype={"CHR38": str, "POS38": "int64", "rsid_from_1000g": str},
        chunksize=chunksize,
    ):
        chunk = chunk[chunk["CHR38"].isin(use_chromosomes)]
        if chunk.empty:
            continue

        keep_parts = []
        for chrom, sub in chunk.groupby("CHR38"):
            positions = needed_positions.get(str(chrom))
            if positions:
                keep_parts.append(sub[sub["POS38"].isin(positions)])

        if not keep_parts:
            continue

        matched = pd.concat(keep_parts, ignore_index=True)
        for row in matched.itertuples(index=False):
            key = (str(row.CHR38), int(row.POS38))
            old = rsid_lookup.get(key)
            if old is None:
                rsid_lookup[key] = row.rsid_from_1000g
            elif row.rsid_from_1000g not in old.split(";"):
                rsid_lookup[key] = f"{old};{row.rsid_from_1000g}"

        print(f"Collected {len(rsid_lookup):,} rsID positions so far...")

    return rsid_lookup


def replace_with_1000g_rsid(input_file, output_file, rsid_lookup, chunksize):
    temp_file = f"{output_file}.tmp_1000g_rsid.tsv.gz"
    with gzip.open(temp_file, "wt") as handle:
        wrote_header = False
        written_rows = 0
        for chunk in pd.read_csv(input_file, sep="\t", compression="gzip", chunksize=chunksize):
            hg38_rsid = [
                rsid_lookup.get((str(chrom), int(pos)), "")
                for chrom, pos in zip(chunk["CHR38"], chunk["POS38"])
            ]
            chunk["rsid_hg38"] = hg38_rsid
            chunk["rsid"] = [
                hg38 if isinstance(hg38, str) and hg38 else original
                for original, hg38 in zip(chunk["rsid"], chunk["rsid_hg38"])
            ]
            chunk.to_csv(
                handle,
                sep="\t",
                index=False,
                header=not wrote_header,
                na_rep="",
            )
            wrote_header = True
            written_rows += len(chunk)
            print(f"Annotated {written_rows:,} rows with 1000G rsID so far...")

    os.replace(temp_file, output_file)


def main():
    args = parse_args()
    output_file = args.output_file or default_output_path(args.input_file)
    make_parent_dir(output_file)

    needed_positions, total_rows, lifted_rows = write_lifted_gwas(
        args.input_file,
        output_file,
        args.chunksize,
        args.add_1000g_rsid,
    )
    print(f"Rows read: {total_rows:,}")
    print(f"Rows lifted: {lifted_rows:,}")

    if args.add_1000g_rsid:
        rsid_lookup = build_rsid_lookup(args.rsid_file, needed_positions, args.chunksize)
        replace_with_1000g_rsid(output_file, output_file, rsid_lookup, args.chunksize)

    print(f"Wrote: {output_file}")


if __name__ == "__main__":
    main()
