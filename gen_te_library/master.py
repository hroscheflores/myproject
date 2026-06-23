#!/usr/bin/env python3

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


DEFAULT_PROJECT_ROOT = Path.cwd()
DEFAULT_SCRIPT_ROOT = DEFAULT_PROJECT_ROOT / "scripts"

DEFAULT_STEP_ONE = DEFAULT_SCRIPT_ROOT / "step_1.py"
DEFAULT_STEP_TWO = DEFAULT_SCRIPT_ROOT / "step_2.py"
DEFAULT_STEP_THREE = DEFAULT_SCRIPT_ROOT / "step_3.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the TE library pipeline workflow."
    )

    parser.add_argument("--step-one", type=Path, default=DEFAULT_STEP_ONE)
    parser.add_argument("--step-two", type=Path, default=DEFAULT_STEP_TWO)
    parser.add_argument("--step-three", type=Path, default=DEFAULT_STEP_THREE)

    parser.add_argument(
        "--input-genomes",
        type=Path,
        required=True,
        help="Input genome list for this run.",
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        required=True,
        help="Output directory for workflow state, runs, and libraries.",
    )

    parser.add_argument("--account", default="")
    parser.add_argument("--nodes", type=int, default=1)
    parser.add_argument("--ntasks", type=int, default=1)
    parser.add_argument("--cpus", type=int, default=12)
    parser.add_argument("--time", default="3-00:00")
    parser.add_argument("--partition", default="general")
    parser.add_argument("--mem", default="100G")

    parser.add_argument("--mask-threads", type=int, default=None)
    parser.add_argument("--model-threads", type=int, default=None)

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional workflow information.",
    )

    return parser.parse_args()


def run_cmd(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True)


def parse_step_one_output(text: str) -> dict:
    result = {"next_genome": None, "next_step": None, "status": None}

    for line in text.splitlines():
        line = line.strip()

        if line.startswith("Next genome:"):
            value = line.split(":", 1)[1].strip()
            result["next_genome"] = None if value.lower() == "none" else value

        elif line.startswith("Next step:"):
            value = line.split(":", 1)[1].strip()
            result["next_step"] = None if value.lower() == "none" else value

        elif line.startswith("Status:"):
            result["status"] = line.split(":", 1)[1].strip()

    return result


def read_master_row(sample: str, master_tsv: Path) -> dict:
    import pandas as pd

    if not master_tsv.exists():
        raise FileNotFoundError(f"Master TSV not found: {master_tsv}")

    df = pd.read_csv(master_tsv, sep="\t", dtype=str).fillna("")
    match = df[df["sample"] == sample].copy()

    if match.empty:
        raise ValueError(f"No row found in master TSV for sample: {sample}")

    if len(match) > 1:
        raise ValueError(f"Multiple rows found in master TSV for sample: {sample}")

    row = match.iloc[0].to_dict()
    row["order"] = int(row["order"])

    return row


def find_genome_path_from_input(sample: str, input_genomes: Path) -> str:
    if not input_genomes.exists():
        raise FileNotFoundError(f"Input genome list not found: {input_genomes}")

    candidates = []

    with input_genomes.open() as handle:
        for raw in handle:
            line = raw.strip()

            if not line or line.startswith("#"):
                continue

            p = Path(line)
            name = p.name

            for suffix in [
                ".fasta.gz",
                ".fa.gz",
                ".fna.gz",
                ".fas.gz",
                ".fasta",
                ".fa",
                ".fna",
                ".fas",
            ]:
                if name.endswith(suffix):
                    name = name[: -len(suffix)]
                    break
            else:
                name = p.stem

            if name == sample:
                candidates.append(line)

    if not candidates:
        raise ValueError(f"Could not find genome path in input list for sample: {sample}")

    if len(candidates) > 1:
        raise ValueError(f"Multiple genome paths matched sample {sample}: {candidates}")

    return candidates[0]


def status_json_path(runs_root: Path, order: int, sample: str, step: str) -> Path:
    return runs_root / f"run_{order:03d}_{sample}" / f"{step}.status.json"


def library_path(libraries_dir: Path, iteration: int) -> Path:
    return libraries_dir / f"ref_iter{iteration:03d}.fa.classified"


def build_step_two_command(
    sample: str,
    step: str,
    args: argparse.Namespace,
    master_tsv: Path,
    input_genomes: Path,
    libraries_dir: Path,
    runs_root: Path,
) -> list[str]:
    row = read_master_row(sample, master_tsv)
    order = row["order"]
    genome_path = row["genome_path"].strip() or find_genome_path_from_input(
        sample,
        input_genomes,
    )

    mask_threads = args.mask_threads if args.mask_threads is not None else args.cpus
    model_threads = args.model_threads if args.model_threads is not None else args.cpus

    cmd = [
        sys.executable,
        str(args.step_two),
        "--work-root",
        str(args.work_root),
        "--runs-root",
        str(runs_root),
        "--libraries-dir",
        str(libraries_dir),
        "--order",
        str(order),
        "--sample",
        sample,
        "--step",
        step,
        "--genome-path",
        genome_path,
        "--nodes",
        str(args.nodes),
        "--ntasks",
        str(args.ntasks),
        "--cpus",
        str(args.cpus),
        "--time",
        args.time,
        "--partition",
        args.partition,
        "--mem",
        args.mem,
        "--mask-threads",
        str(mask_threads),
        "--model-threads",
        str(model_threads),
    ]

    if args.account:
        cmd.extend(["--account", args.account])

    if step == "masking":
        if order > 1:
            cmd.extend(["--library-in", str(library_path(libraries_dir, order - 1))])

    elif step == "concat":
        if order > 1:
            cmd.extend(["--library-in", str(library_path(libraries_dir, order - 1))])

        cmd.extend(["--library-out", str(library_path(libraries_dir, order))])

    return cmd


def main() -> int:
    args = parse_args()

    work_root = args.work_root
    input_genomes = args.input_genomes
    state_dir = work_root / "state"
    runs_root = work_root / "runs"
    libraries_dir = work_root / "libraries"
    master_tsv = state_dir / "master_genome_status.tsv"

    for path in [args.step_one, args.step_two, args.step_three]:
        if not path.exists():
            print(f"ERROR: missing script: {path}", file=sys.stderr)
            return 2

    if not input_genomes.exists():
        print(f"ERROR: input genome list does not exist: {input_genomes}", file=sys.stderr)
        return 2

    work_root.mkdir(parents=True, exist_ok=True)
    state_dir.mkdir(parents=True, exist_ok=True)
    runs_root.mkdir(parents=True, exist_ok=True)
    libraries_dir.mkdir(parents=True, exist_ok=True)

    print("Starting workflow.")

    if args.verbose:
        print(f"Work root: {work_root}")
        print(f"Input genomes: {input_genomes}")
        print(f"Master TSV: {master_tsv}")
        print(f"Runs root: {runs_root}")
        print(f"Libraries dir: {libraries_dir}")
        print(
            "Resources: "
            f"account={args.account or 'none'} "
            f"nodes={args.nodes} "
            f"ntasks={args.ntasks} "
            f"cpus={args.cpus} "
            f"time={args.time} "
            f"partition={args.partition} "
            f"mem={args.mem}"
        )

        effective_mask_threads = (
            args.mask_threads if args.mask_threads is not None else args.cpus
        )
        effective_model_threads = (
            args.model_threads if args.model_threads is not None else args.cpus
        )

        print(
            "Tool threads: "
            f"mask={effective_mask_threads} "
            f"model={effective_model_threads}"
        )

    while True:
        step1 = run_cmd(
            [
                sys.executable,
                str(args.step_one),
                "--input-genomes",
                str(input_genomes),
                "--master-tsv",
                str(master_tsv),
                "--libraries-dir",
                str(libraries_dir),
            ]
        )

        if step1.returncode != 0:
            print("ERROR: step 1 failed.", file=sys.stderr)
            print(step1.stdout)
            print(step1.stderr, file=sys.stderr)
            return 1

        parsed = parse_step_one_output(step1.stdout)

        next_genome = parsed["next_genome"]
        next_step = parsed["next_step"]

        if not next_genome or not next_step:
            print("No further work detected.")
            return 0

        print(f"Processing sample={next_genome}, step={next_step}")

        step2_cmd = build_step_two_command(
            sample=next_genome,
            step=next_step,
            args=args,
            master_tsv=master_tsv,
            input_genomes=input_genomes,
            libraries_dir=libraries_dir,
            runs_root=runs_root,
        )

        step2 = run_cmd(step2_cmd)

        if step2.stdout.strip():
            print(step2.stdout.strip())

        if step2.returncode != 0:
            print("ERROR: step 2 failed.", file=sys.stderr)
            print(step2.stderr, file=sys.stderr)
            return 1

        row = read_master_row(next_genome, master_tsv)
        status_json = status_json_path(runs_root, row["order"], next_genome, next_step)

        if not status_json.exists():
            print(f"ERROR: expected status JSON not found: {status_json}", file=sys.stderr)
            return 1

        print("Updating workflow state.")

        step3 = run_cmd(
            [
                sys.executable,
                str(args.step_three),
                "--master-tsv",
                str(master_tsv),
                "--status-json",
                str(status_json),
            ]
        )

        if step3.stdout.strip():
            print(step3.stdout.strip())

        if step3.returncode != 0:
            print("ERROR: step 3 failed.", file=sys.stderr)
            print(step3.stderr, file=sys.stderr)
            return 1

        print("Step completed successfully.\n")


if __name__ == "__main__":
    sys.exit(main())