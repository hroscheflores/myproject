#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = [
    "order",
    "sample",
    "genome_path",
    "mask_done",
    "model_done",
    "concat_done",
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

STEP_MASK = "masking"
STEP_MODEL = "model"
STEP_CONCAT = "concat"
STEP_COMPLETE = "complete"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Determine the next genome and workflow step."
    )

    parser.add_argument(
        "--input-genomes",
        type=Path,
        required=True,
        help="Path to the ordered input genome list.",
    )
    parser.add_argument(
        "--master-tsv",
        type=Path,
        required=True,
        help="Path to the workflow status TSV.",
    )
    parser.add_argument(
        "--libraries-dir",
        type=Path,
        default=None,
        help="Optional library output directory.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional workflow information.",
    )

    return parser.parse_args()


def read_input_genomes(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Input genome list not found: {path}")

    genomes = []

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()

            if not line or line.startswith("#"):
                continue

            genomes.append(line)

    return list(dict.fromkeys(genomes))


def sample_name(genome_path: str) -> str:
    name = Path(genome_path).name

    for suffix in GENOME_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]

    return Path(genome_path).stem


def read_master_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=REQUIRED_COLUMNS)

    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(
            f"Workflow status TSV is missing required columns: {', '.join(missing)}"
        )

    df["order"] = pd.to_numeric(df["order"], errors="coerce")

    for col in ["mask_done", "model_done", "concat_done"]:
        df[col] = df[col].astype(str).str.lower().isin(["true", "1", "yes"])

    return df


def determine_next_step(row: pd.Series) -> str:
    if not bool(row["mask_done"]):
        return STEP_MASK

    if not bool(row["model_done"]):
        return STEP_MODEL

    if not bool(row["concat_done"]):
        return STEP_CONCAT

    return STEP_COMPLETE


def find_last_successful_job(tsv: pd.DataFrame) -> str:
    successful_rows = tsv[
        (tsv["mask_done"] == True)
        | (tsv["model_done"] == True)
        | (tsv["concat_done"] == True)
    ].copy()

    if successful_rows.empty:
        return "none"

    successful_rows = successful_rows.sort_values("order")
    last_row = successful_rows.iloc[-1]
    last_step = determine_next_step(last_row)

    if last_step == STEP_COMPLETE:
        return f"{last_row['sample']}: {STEP_CONCAT}"

    if bool(last_row["model_done"]):
        return f"{last_row['sample']}: {STEP_MODEL}"

    if bool(last_row["mask_done"]):
        return f"{last_row['sample']}: {STEP_MASK}"

    return "none"


def main() -> int:
    args = parse_args()

    genomes_in_input = read_input_genomes(args.input_genomes)
    tsv = read_master_tsv(args.master_tsv)

    input_df = pd.DataFrame(
        {
            "genome_path": genomes_in_input,
            "sample": [sample_name(path) for path in genomes_in_input],
        }
    )

    if tsv.empty:
        if args.verbose:
            print("Completed genomes: none")
            print("Last successful job: none")

        if input_df.empty:
            print("Next genome: none")
            print("Next step: none")
            print("Status: no genomes listed in input file")
            return 0

        print(f"Next genome: {input_df.iloc[0]['sample']}")
        print(f"Next step: {STEP_MASK}")
        print("Status: no workflow status TSV found")
        return 0

    merged = input_df.merge(
        tsv,
        on=["genome_path", "sample"],
        how="left",
        suffixes=("", "_tsv"),
    )

    known_in_tsv = merged["order"].notna().sum()
    new_not_in_tsv = merged["order"].isna().sum()

    completed = tsv.loc[tsv["concat_done"] == True].copy()
    completed = completed.sort_values("order")
    completed_samples = completed["sample"].tolist()

    last_successful_job = find_last_successful_job(tsv)

    next_genome = None
    next_step = None

    for _, row in merged.iterrows():
        if pd.isna(row["order"]):
            next_genome = row["sample"]
            next_step = STEP_MASK
            break

        row_step = determine_next_step(row)

        if row_step != STEP_COMPLETE:
            next_genome = row["sample"]
            next_step = row_step
            break

    if args.verbose:
        completed_text = ", ".join(completed_samples) if completed_samples else "none"

        print(f"Completed genomes: {completed_text}")
        print(f"Last successful job: {last_successful_job}")
        print(f"Genomes in input list: {len(input_df)}")
        print(f"Genomes already tracked in TSV: {known_in_tsv}")
        print(f"New genomes not yet in TSV: {new_not_in_tsv}")

    if next_genome is None:
        print("Next genome: none")
        print("Next step: none")
        print("Status: all genomes currently listed in input appear complete")
        return 0

    print(f"Next genome: {next_genome}")
    print(f"Next step: {next_step}")

    if new_not_in_tsv > 0:
        print("Status: new genome detected in input list")
    else:
        print("Status: resume existing tracked work")

    return 0


if __name__ == "__main__":
    sys.exit(main())