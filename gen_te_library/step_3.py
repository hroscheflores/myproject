#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


STEP_MASK = "masking"
STEP_MODEL = "model"
STEP_CONCAT = "concat"
VALID_STEPS = {STEP_MASK, STEP_MODEL, STEP_CONCAT}


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


def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_status_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Status JSON not found: {path}")

    with path.open() as handle:
        return json.load(handle)


def read_master_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Workflow status TSV not found: {path}")

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


def save_master_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    out = df.copy()
    out = out.sort_values("order")
    out.to_csv(path, sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Update the workflow status TSV after a successful job."
    )

    parser.add_argument(
        "--status-json",
        type=Path,
        required=True,
        help="Path to the step status JSON.",
    )
    parser.add_argument(
        "--master-tsv",
        type=Path,
        required=True,
        help="Path to the workflow status TSV.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional update information.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    payload = read_status_json(args.status_json)

    step = payload.get("step")
    sample = payload.get("sample")
    order = payload.get("order")
    success = payload.get("success", False)

    if step not in VALID_STEPS:
        print(f"ERROR: invalid step in status JSON: {step}", file=sys.stderr)
        return 2

    if not success:
        print(
            f"ERROR: status JSON does not indicate success for {sample} {step}",
            file=sys.stderr,
        )
        return 1

    df = read_master_tsv(args.master_tsv)

    match = df[(df["order"] == int(order)) & (df["sample"] == str(sample))]

    if match.empty:
        print(
            f"ERROR: no TSV row found for order={order}, sample={sample}",
            file=sys.stderr,
        )
        return 2

    if len(match) > 1:
        print(
            f"ERROR: multiple TSV rows found for order={order}, sample={sample}",
            file=sys.stderr,
        )
        return 2

    idx = match.index[0]
    timestamp = payload.get("submitted_at") or now_iso()

    if step == STEP_MASK:
        df.at[idx, "mask_done"] = True

        if not str(df.at[idx, "mask_completed_at"]).strip():
            df.at[idx, "mask_completed_at"] = timestamp

        df.at[idx, "last_successful_step"] = STEP_MASK
        df.at[idx, "status"] = "awaiting_model"

    elif step == STEP_MODEL:
        df.at[idx, "model_done"] = True

        if not str(df.at[idx, "model_completed_at"]).strip():
            df.at[idx, "model_completed_at"] = timestamp

        df.at[idx, "last_successful_step"] = STEP_MODEL
        df.at[idx, "status"] = "awaiting_concat"

    elif step == STEP_CONCAT:
        df.at[idx, "concat_done"] = True

        if not str(df.at[idx, "concat_completed_at"]).strip():
            df.at[idx, "concat_completed_at"] = timestamp

        df.at[idx, "last_successful_step"] = STEP_CONCAT
        df.at[idx, "status"] = "complete"

    update_note = f"{now_iso()} updated from {args.status_json.name}"
    notes = str(df.at[idx, "notes"]).strip()

    if notes:
        df.at[idx, "notes"] = f"{notes} | {update_note}"
    else:
        df.at[idx, "notes"] = update_note

    save_master_tsv(df, args.master_tsv)

    print(f"Updated status TSV: {sample} -> {step}")

    if args.verbose:
        print(f"Status JSON: {args.status_json}")
        print(f"Status TSV: {args.master_tsv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())