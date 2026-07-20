import argparse
import gzip
import os
import tempfile
from collections import defaultdict

import pandas as pd


DEFAULT_PAPER_DIR = os.environ.get(
    "PAPER_SMIXCAN_DIR",
    "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan",
)
DEFAULT_INPUT_FILE = os.path.join(
    DEFAULT_PAPER_DIR,
    "Heart",
    "Data",
    "I21_MI_.gwas.imputed_v3.both_sexes.tsv.bgz",
)
DEFAULT_RSID_FILE = os.path.join(
    DEFAULT_PAPER_DIR,
    "Heart",
    "Data",
    "1000g_b38_snpIDs.txt",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Lift MI GWAS variants from GRCh37/hg19 to GRCh38/hg38 and add rsID annotation."
    )
    parser.add_argument("--input-file", default=DEFAULT_INPUT_FILE)
    parser.add_argument("--rsid-file", default=DEFAULT_RSID_FILE)
    parser.add_argument(
        "--output-file",
        default=None,
        help="Default: input filename with _hg38_rsid.tsv.gz suffix in the same directory.",
    )
    parser.add_argument(
        "--coord-file",
        default=None,
        help="Default: input filename with _hg38_rsid_chr_pos.txt suffix in the same directory.",
    )
    parser.add_argument(
        "--coord-only",
        action="store_true",
        help="Write only rsid/CHR38/POS38 coordinate output, not the full annotated GWAS file.",
    )
    parser.add_argument("--chunksize", type=int, default=250_000)
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep the temporary lifted file used before rsID annotation.",
    )
    return parser.parse_args()


def default_output_path(input_file):
    output_dir = os.path.dirname(input_file)
    base = os.path.basename(input_file)
    for suffix in (".tsv.bgz", ".tsv.gz", ".bgz", ".gz", ".tsv", ".txt"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return os.path.join(output_dir, f"{base}_hg38_rsid.tsv.gz")


def default_coord_path(input_file):
    output_dir = os.path.dirname(input_file)
    base = os.path.basename(input_file)
    for suffix in (".tsv.bgz", ".tsv.gz", ".bgz", ".gz", ".tsv", ".txt"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return os.path.join(output_dir, f"{base}_hg38_rsid_chr_pos.txt")


def make_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def parse_variant_column(df):
    variant_parts = df["variant"].astype(str).str.split(":", expand=True)
    if variant_parts.shape[1] < 4:
        raise ValueError("Expected variant column in chr:pos:ref:alt format.")

    df["CHR37"] = variant_parts[0].str.replace("^chr", "", regex=True)
    df["POS37"] = pd.to_numeric(variant_parts[1], errors="coerce").astype("Int64")
    df["REF"] = variant_parts[2]
    df["ALT"] = variant_parts[3]
    return df


def liftover_coordinate(lo, chrom, pos):
    if pd.isna(pos):
        return pd.NA, pd.NA

    # GWAS coordinates are 1-based. pyliftover maps 0-based chain coordinates.
    lifted = lo.convert_coordinate(f"chr{chrom}", int(pos) - 1)
    if not lifted:
        return pd.NA, pd.NA
    chrom38 = str(lifted[0][0]).replace("chr", "")
    pos38 = int(lifted[0][1]) + 1
    return chrom38, pos38


def lift_gwas(input_file, temp_file, chunksize):
    try:
        from pyliftover import LiftOver
    except ImportError as exc:
        raise SystemExit(
            "Missing dependency: pyliftover. Install it with:\n"
            "  python3 -m pip install pyliftover\n"
        ) from exc

    lo = LiftOver("hg19", "hg38")
    needed_positions = defaultdict(set)
    total_rows = 0
    lifted_rows = 0
    wrote_header = False

    for chunk in pd.read_csv(
        input_file,
        sep="\t",
        compression="gzip",
        chunksize=chunksize,
        low_memory=False,
    ):
        total_rows += len(chunk)
        chunk = parse_variant_column(chunk)
        chunk = chunk.dropna(subset=["POS37"]).copy()
        lifted = [
            liftover_coordinate(lo, chrom, pos)
            for chrom, pos in zip(chunk["CHR37"], chunk["POS37"])
        ]
        chunk["CHR38"] = [item[0] for item in lifted]
        chunk["POS38"] = [item[1] for item in lifted]
        chunk = chunk.dropna(subset=["CHR38", "POS38"]).copy()
        chunk["POS38"] = chunk["POS38"].astype("int64")
        chunk["variant_hg38"] = (
            chunk["CHR38"].astype(str)
            + ":"
            + chunk["POS38"].astype(str)
            + ":"
            + chunk["REF"].astype(str)
            + ":"
            + chunk["ALT"].astype(str)
        )

        for chrom, positions in chunk.groupby("CHR38")["POS38"]:
            needed_positions[str(chrom)].update(positions.astype(int).tolist())

        chunk.to_csv(
            temp_file,
            sep="\t",
            index=False,
            mode="a",
            header=not wrote_header,
            na_rep="",
        )
        wrote_header = True
        lifted_rows += len(chunk)
        print(f"Lifted {lifted_rows:,} rows so far...")

    return needed_positions, total_rows, lifted_rows


def build_rsid_lookup(rsid_file, needed_positions, chunksize):
    rsid_lookup = {}
    names = ["CHR38", "POS38", "rsid"]
    use_chromosomes = set(needed_positions)

    for chunk in pd.read_csv(
        rsid_file,
        sep="\t",
        names=names,
        dtype={"CHR38": str, "POS38": "int64", "rsid": str},
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
                rsid_lookup[key] = row.rsid
            elif row.rsid not in old.split(";"):
                rsid_lookup[key] = f"{old};{row.rsid}"

        print(f"Collected {len(rsid_lookup):,} rsID positions so far...")

    return rsid_lookup


def add_rsid_annotation(temp_file, output_file, rsid_lookup, chunksize):
    with gzip.open(output_file, "wt") as handle:
        wrote_header = False
        written_rows = 0
        for chunk in pd.read_csv(temp_file, sep="\t", chunksize=chunksize, low_memory=False):
            chunk["rsid"] = [
                rsid_lookup.get((str(chrom), int(pos)), "")
                for chrom, pos in zip(chunk["CHR38"], chunk["POS38"])
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
            print(f"Annotated {written_rows:,} rows so far...")


def write_coordinate_file(temp_file, coord_file, rsid_lookup, chunksize):
    wrote_header = False
    written_rows = 0

    with open(coord_file, "w") as handle:
        for chunk in pd.read_csv(temp_file, sep="\t", chunksize=chunksize, low_memory=False):
            coord = pd.DataFrame(
                {
                    "rsid": [
                        rsid_lookup.get((str(chrom), int(pos)), "")
                        for chrom, pos in zip(chunk["CHR38"], chunk["POS38"])
                    ],
                    "CHR38": chunk["CHR38"].astype(str),
                    "POS38": chunk["POS38"].astype("int64"),
                }
            )
            coord = coord[coord["rsid"] != ""].drop_duplicates()
            coord.to_csv(
                handle,
                sep="\t",
                index=False,
                header=not wrote_header,
                na_rep="",
            )
            wrote_header = True
            written_rows += len(coord)
            print(f"Wrote {written_rows:,} rsID coordinate rows so far...")


def main():
    args = parse_args()
    input_file = args.input_file
    rsid_file = args.rsid_file
    output_file = args.output_file or default_output_path(input_file)
    coord_file = args.coord_file or default_coord_path(input_file)
    output_dir = os.path.dirname(output_file) or os.getcwd()
    make_parent_dir(output_file)
    make_parent_dir(coord_file)

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".lifted.tsv",
        prefix="mi_gwas_",
        dir=output_dir,
        delete=False,
    ) as tmp:
        temp_file = tmp.name

    print(f"Input:  {input_file}")
    print(f"rsIDs:  {rsid_file}")
    if not args.coord_only:
        print(f"Output: {output_file}")
    print(f"Coords: {coord_file}")
    print(f"Temp:   {temp_file}")

    needed_positions, total_rows, lifted_rows = lift_gwas(
        input_file=input_file,
        temp_file=temp_file,
        chunksize=args.chunksize,
    )
    print(f"LiftOver kept {lifted_rows:,} of {total_rows:,} input rows.")

    rsid_lookup = build_rsid_lookup(
        rsid_file=rsid_file,
        needed_positions=needed_positions,
        chunksize=args.chunksize,
    )
    print(f"Found rsIDs for {len(rsid_lookup):,} unique GRCh38 positions.")

    write_coordinate_file(
        temp_file=temp_file,
        coord_file=coord_file,
        rsid_lookup=rsid_lookup,
        chunksize=args.chunksize,
    )

    if not args.coord_only:
        add_rsid_annotation(
            temp_file=temp_file,
            output_file=output_file,
            rsid_lookup=rsid_lookup,
            chunksize=args.chunksize,
        )

    if args.keep_temp:
        print(f"Kept temp file: {temp_file}")
    else:
        os.unlink(temp_file)

    print("Done.")


if __name__ == "__main__":
    main()
