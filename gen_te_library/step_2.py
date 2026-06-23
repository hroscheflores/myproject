#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path


STEP_MASK = "masking"
STEP_MODEL = "model"
STEP_CONCAT = "concat"
VALID_STEPS = {STEP_MASK, STEP_MODEL, STEP_CONCAT}

SLURM_SCRIPT_ROOT: Path
RUNS_ROOT: Path
LIBRARIES_ROOT: Path


def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def clean_sample_name(sample: str) -> str:
    sample = sample.strip()
    return re.sub(r"[^A-Za-z0-9._-]+", "_", sample)


def build_sbatch_header(
    account: str | None,
    nodes: int,
    ntasks: int,
    cpus: int,
    time_limit: str,
    partition: str,
    mem: str,
    mail_user: str | None,
) -> str:
    lines = [
        "#!/usr/bin/env bash",
    ]

    if account:
        lines.append(f"#SBATCH -A {account}")

    lines.extend(
        [
            f"#SBATCH -N {nodes}",
            f"#SBATCH -n {ntasks}",
            f"#SBATCH -c {cpus}",
            f"#SBATCH -t {time_limit}",
            f"#SBATCH -p {partition}",
            f"#SBATCH --mem={mem}",
        ]
    )

    if mail_user:
        lines.extend(
            [
                "#SBATCH --mail-type=ALL",
                f"#SBATCH --mail-user={mail_user}",
            ]
        )

    return "\n".join(lines) + "\n"


def run_dir(order: int, sample: str) -> Path:
    return RUNS_ROOT / f"run_{order:03d}_{sample}"


def slurm_dir(order: int, sample: str) -> Path:
    return SLURM_SCRIPT_ROOT / f"run_{order:03d}_{sample}"


def masked_dir(order: int, sample: str) -> Path:
    return run_dir(order, sample) / "masked"


def model_dir(order: int, sample: str) -> Path:
    return run_dir(order, sample) / "model"


def logs_dir(order: int, sample: str) -> Path:
    return run_dir(order, sample) / "logs"


def state_file(order: int, sample: str, step: str) -> Path:
    return run_dir(order, sample) / f"{step}.status.json"


def script_file(order: int, sample: str, step: str) -> Path:
    return slurm_dir(order, sample) / f"{order:03d}_{sample}.{step}.slurm"


def expected_masked_genome(order: int, sample: str, genome_path: Path) -> Path:
    return masked_dir(order, sample) / f"{genome_path.name}.masked"


def ensure_dirs(order: int, sample: str) -> None:
    for directory in [
        slurm_dir(order, sample),
        run_dir(order, sample),
        masked_dir(order, sample),
        model_dir(order, sample),
        logs_dir(order, sample),
    ]:
        directory.mkdir(parents=True, exist_ok=True)


def write_status(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def build_mask_script(
    order: int,
    sample: str,
    genome_path: Path,
    library_in: Path,
    threads: int,
    sbatch_header: str,
) -> str:
    outdir = masked_dir(order, sample)
    logfile = logs_dir(order, sample) / f"{sample}.masking.log"

    return sbatch_header + f"""
set -euo pipefail

mkdir -p "{outdir}"
mkdir -p "{logs_dir(order, sample)}"

RepeatMasker -pa {threads} -lib "{library_in}" -dir "{outdir}" "{genome_path}" > "{logfile}" 2>&1
"""


def build_model_script(
    order: int,
    sample: str,
    genome_path: Path,
    threads: int,
    sbatch_header: str,
) -> str:
    if order == 1:
        model_input = genome_path
    else:
        model_input = expected_masked_genome(order, sample, genome_path)

    mdir = model_dir(order, sample)
    db_prefix = mdir / f"{sample}_db"
    logfile = logs_dir(order, sample) / f"{sample}.model.log"

    return sbatch_header + f"""
set -euo pipefail

mkdir -p "{mdir}"
mkdir -p "{logs_dir(order, sample)}"

if [[ ! -f "{model_input}" ]]; then
    echo "Missing model input genome: {model_input}" >&2
    exit 1
fi

cd "{mdir}"

BuildDatabase -name "{db_prefix}" "{model_input}" > "{logfile}" 2>&1
RepeatModeler -database "{db_prefix}" -LTRStruct -threads {threads} >> "{logfile}" 2>&1
"""


def build_concat_script(
    order: int,
    sample: str,
    library_in: Path | None,
    library_out: Path,
    sbatch_header: str,
) -> str:
    logfile = logs_dir(order, sample) / f"{sample}.concat.log"

    if order == 1:
        return sbatch_header + f"""
set -euo pipefail

mkdir -p "{library_out.parent}"
mkdir -p "{logs_dir(order, sample)}"

CONSENSI_FILE=$(find "{model_dir(order, sample)}" -type f -path "*/RM_*/consensi.fa.classified" | head -n 1)

if [[ -z "$CONSENSI_FILE" ]]; then
    echo "No consensi.fa.classified found under {model_dir(order, sample)}" >&2
    exit 1
fi

cp "$CONSENSI_FILE" "{library_out}" > "{logfile}" 2>&1
echo "Created {library_out}" >> "{logfile}"
"""

    return sbatch_header + f"""
set -euo pipefail

mkdir -p "{library_out.parent}"
mkdir -p "{logs_dir(order, sample)}"

if [[ ! -f "{library_in}" ]]; then
    echo "Missing input library: {library_in}" >&2
    exit 1
fi

CONSENSI_FILE=$(find "{model_dir(order, sample)}" -type f -path "*/RM_*/consensi.fa.classified" | head -n 1)

if [[ -z "$CONSENSI_FILE" ]]; then
    echo "No consensi.fa.classified found under {model_dir(order, sample)}" >&2
    exit 1
fi

TMP_OUT="{library_out}.tmp"

cp "{library_in}" "$TMP_OUT"
cat "$CONSENSI_FILE" >> "$TMP_OUT"
mv "$TMP_OUT" "{library_out}" > "{logfile}" 2>&1

echo "Created {library_out}" >> "{logfile}"
"""


def make_script(
    order: int,
    sample: str,
    step: str,
    genome_path: Path,
    library_in: Path | None,
    library_out: Path | None,
    mask_threads: int,
    model_threads: int,
    sbatch_header: str,
) -> str:
    if step == STEP_MASK:
        if library_in is None:
            raise ValueError("--library-in is required for masking")
        return build_mask_script(
            order,
            sample,
            genome_path,
            library_in,
            mask_threads,
            sbatch_header,
        )

    if step == STEP_MODEL:
        return build_model_script(
            order,
            sample,
            genome_path,
            model_threads,
            sbatch_header,
        )

    if step == STEP_CONCAT:
        if library_out is None:
            raise ValueError("--library-out is required for concat")
        return build_concat_script(
            order,
            sample,
            library_in,
            library_out,
            sbatch_header,
        )

    raise ValueError(f"Unsupported step: {step}")


def submit_and_wait(script_path: Path) -> tuple[int, str]:
    proc = subprocess.run(
        ["sbatch", "--wait", str(script_path)],
        capture_output=True,
        text=True,
    )

    combined = (proc.stdout or "") + (proc.stderr or "")

    return proc.returncode, combined.strip()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create and submit one SLURM job for a workflow step."
    )

    parser.add_argument("--work-root", type=Path, default=None)
    parser.add_argument("--runs-root", type=Path, default=None)
    parser.add_argument("--libraries-dir", type=Path, default=None)
    parser.add_argument("--slurm-script-root", type=Path, default=None)

    parser.add_argument("--order", type=int, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--step", required=True, choices=sorted(VALID_STEPS))
    parser.add_argument("--genome-path", required=True)

    parser.add_argument("--library-in", type=Path, default=None)
    parser.add_argument("--library-out", type=Path, default=None)

    parser.add_argument("--mask-threads", type=int, default=None)
    parser.add_argument("--model-threads", type=int, default=None)

    parser.add_argument("--account", default="")
    parser.add_argument("--nodes", type=int, default=1)
    parser.add_argument("--ntasks", type=int, default=1)
    parser.add_argument("--cpus", type=int, default=12)
    parser.add_argument("--time", dest="time_limit", default="0-12:00")
    parser.add_argument("--partition", default="general")
    parser.add_argument("--mem", default="100G")
    parser.add_argument("--mail-user", default="")

    parser.add_argument("--print-script-only", action="store_true")
    parser.add_argument("--verbose", action="store_true")

    return parser.parse_args()


def resolve_roots(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    if args.work_root is not None:
        runs_root = args.runs_root or args.work_root / "runs"
        libraries_root = args.libraries_dir or args.work_root / "libraries"
        slurm_script_root = args.slurm_script_root or args.work_root / "slurm_scripts"
        return runs_root, libraries_root, slurm_script_root

    current_dir = Path.cwd()

    runs_root = args.runs_root or current_dir / "runs"
    libraries_root = args.libraries_dir or current_dir / "libraries"
    slurm_script_root = args.slurm_script_root or current_dir / "slurm_scripts"

    return runs_root, libraries_root, slurm_script_root


def main() -> int:
    global SLURM_SCRIPT_ROOT, RUNS_ROOT, LIBRARIES_ROOT

    args = parse_args()

    RUNS_ROOT, LIBRARIES_ROOT, SLURM_SCRIPT_ROOT = resolve_roots(args)

    sample = clean_sample_name(args.sample)
    order = args.order
    step = args.step
    genome_path = Path(args.genome_path)

    library_in = args.library_in
    library_out = args.library_out

    if args.verbose:
        print(f"Runs root: {RUNS_ROOT}")
        print(f"Libraries dir: {LIBRARIES_ROOT}")
        print(f"SLURM script root: {SLURM_SCRIPT_ROOT}")

    if not genome_path.exists():
        print(f"ERROR: genome path does not exist: {genome_path}", file=sys.stderr)
        return 2

    if step == STEP_MASK and library_in is not None and not library_in.exists():
        print(f"ERROR: input library does not exist: {library_in}", file=sys.stderr)
        return 2

    if step == STEP_CONCAT and order > 1:
        if library_in is None:
            print("ERROR: --library-in is required for concat after the first genome.", file=sys.stderr)
            return 2

        if not library_in.exists():
            print(f"ERROR: input library does not exist: {library_in}", file=sys.stderr)
            return 2

    mask_threads = args.mask_threads if args.mask_threads is not None else args.cpus
    model_threads = args.model_threads if args.model_threads is not None else args.cpus

    sbatch_header = build_sbatch_header(
        account=args.account,
        nodes=args.nodes,
        ntasks=args.ntasks,
        cpus=args.cpus,
        time_limit=args.time_limit,
        partition=args.partition,
        mem=args.mem,
        mail_user=args.mail_user,
    )

    ensure_dirs(order, sample)

    script_text = make_script(
        order=order,
        sample=sample,
        step=step,
        genome_path=genome_path,
        library_in=library_in,
        library_out=library_out,
        mask_threads=mask_threads,
        model_threads=model_threads,
        sbatch_header=sbatch_header,
    )

    slurm_path = script_file(order, sample, step)
    slurm_path.write_text(script_text)
    slurm_path.chmod(0o755)

    status_path = state_file(order, sample, step)

    initial_status = {
        "order": order,
        "sample": sample,
        "step": step,
        "created_at": now_iso(),
        "submitted": False,
        "completed": False,
        "success": False,
        "script_path": str(slurm_path),
        "genome_path": str(genome_path),
        "library_in": str(library_in) if library_in else None,
        "library_out": str(library_out) if library_out else None,
        "runs_root": str(RUNS_ROOT),
        "libraries_root": str(LIBRARIES_ROOT),
        "slurm_script_root": str(SLURM_SCRIPT_ROOT),
        "account": args.account or None,
        "nodes": args.nodes,
        "ntasks": args.ntasks,
        "cpus": args.cpus,
        "time_limit": args.time_limit,
        "partition": args.partition,
        "mem": args.mem,
        "mail_user": args.mail_user or None,
        "mask_threads": mask_threads,
        "model_threads": model_threads,
    }

    write_status(status_path, initial_status)

    if args.print_script_only:
        print(f"Created script: {slurm_path}")
        print("Submission skipped.")
        return 0

    rc, submit_output = submit_and_wait(slurm_path)

    final_status = dict(initial_status)
    final_status.update(
        {
            "submitted_at": now_iso(),
            "submitted": True,
            "completed": rc == 0,
            "success": rc == 0,
            "sbatch_return_code": rc,
            "sbatch_output": submit_output,
        }
    )

    write_status(status_path, final_status)

    if rc == 0:
        print(f"SUCCESS: {sample} {step}")
        print(f"Script: {slurm_path}")
        print(f"Status file: {status_path}")
        return 0

    print(f"FAILED: {sample} {step}", file=sys.stderr)
    print(f"Script: {slurm_path}", file=sys.stderr)
    print(f"Status file: {status_path}", file=sys.stderr)

    if submit_output:
        print(submit_output, file=sys.stderr)

    return 1


if __name__ == "__main__":
    sys.exit(main())