#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = [
    "order",
    "sample",
    "genome_path",
    "mask_done",
    "model_done",
    "concat_done",
    "mask_completed_at",
    "model_completed_at",
    "concat_completed_at",
    "last_successful_step",
    "status",
    "notes",
]


GENOME_SUFFIXES = [
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
    ".fas.gz",
    ".fasta",
    ".fa",
    ".fna",
    ".fas",
]


def sample_name(path: str) -> str:
    p = Path(path)
    name = p.name

    for suffix in GENOME_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]

    return p.stem


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a workflow status TSV from an ordered genome list."
    )

    parser.add_argument(
        "--input-genomes",
        required=True,
        type=Path,
        help="Ordered input genome list.",
    )
    parser.add_argument(
        "--master-tsv",
        required=True,
        type=Path,
        help="Output workflow status TSV.",
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Back up an existing output TSV before overwriting it.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print the sample order after creating the TSV.",
    )

    return parser.parse_args()


def read_genome_list(input_genomes: Path) -> list[str]:
    if not input_genomes.exists():
        raise FileNotFoundError(f"Input genome list not found: {input_genomes}")

    genomes = []

    with input_genomes.open() as handle:
        for raw in handle:
            line = raw.strip()

            if not line or line.startswith("#"):
                continue

            genomes.append(line)

    return genomes


def build_status_rows(genomes: list[str]) -> list[dict]:
    rows = []

    for i, genome in enumerate(genomes, start=1):
        is_seed = i == 1

        rows.append(
            {
                "order": i,
                "sample": sample_name(genome),
                "genome_path": genome,
                "mask_done": "True" if is_seed else "False",
                "model_done": "False",
                "concat_done": "False",
                "mask_completed_at": "seed_no_mask" if is_seed else "",
                "model_completed_at": "",
                "concat_completed_at": "",
                "last_successful_step": "seed_no_mask" if is_seed else "0",
                "status": "awaiting_model" if is_seed else "not_started",
                "notes": (
                    "seed genome; starts with modeling without masking"
                    if is_seed
                    else ""
                ),
            }
        )

    return rows


def write_status_tsv(rows: list[dict], master_tsv: Path) -> None:
    df = pd.DataFrame(rows, columns=REQUIRED_COLUMNS)

    master_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(master_tsv, sep="\t", index=False)


def backup_existing_file(path: Path) -> None:
    if not path.exists():
        return

    backup_path = path.with_suffix(path.suffix + ".bak")
    shutil.copy2(path, backup_path)

    print(f"Backup written to: {backup_path}")


def main() -> int:
    args = parse_args()

    if args.master_tsv.exists() and args.backup:
        backup_existing_file(args.master_tsv)

    genomes = read_genome_list(args.input_genomes)
    rows = build_status_rows(genomes)

    write_status_tsv(rows, args.master_tsv)

    print(f"Created status TSV: {args.master_tsv}")
    print(f"Genomes loaded: {len(rows)}")

    if args.verbose:
        print("Sample order:")

        for row in rows:
            seed_label = " [seed; starts at model]" if row["order"] == 1 else ""
            print(f"  {row['order']}: {row['sample']}{seed_label}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())