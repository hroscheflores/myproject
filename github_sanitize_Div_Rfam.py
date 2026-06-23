#!/usr/bin/env python3

import argparse
import os
import re
import sys
from collections import defaultdict


DEFAULT_TAXONOMY_DICT = "path/to/taxonomy_dictionary.txt"


def sanitize_header_token(token: str) -> str:
    if token in ("Rclass", "Rfam"):
        return token

    return (
        token.strip()
        .rstrip(":")
        .replace("[", "")
        .replace("]", "")
        .replace(";", "_")
    )


def sanitize_to_tsv(input_path: str, output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)

    base_name = os.path.basename(input_path)
    out_name = base_name[:-4] + ".tsv" if base_name.endswith(".tab") else base_name + ".tsv"
    output_path = os.path.join(output_dir, out_name)

    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for line_number, line in enumerate(fin):
            raw = line.rstrip("\n")

            if not raw.strip() or line_number == 0:
                continue

            tokens = re.split(r"\s+", raw.strip())
            clean_tokens = (
                [sanitize_header_token(token) for token in tokens]
                if line_number == 1
                else tokens
            )

            if len(tokens) != len(clean_tokens):
                sys.exit(f"ERROR: Column count mismatch on line {line_number + 1}")

            fout.write("\t".join(clean_tokens) + "\n")

    return output_path


def load_taxonomy_dict(path: str) -> dict:
    fam_to_superfamily = {}

    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) >= 2:
                fam_to_superfamily[parts[0]] = parts[1]

    return fam_to_superfamily


def infer_superfamily(order: str, family: str, fam_to_superfamily: dict) -> str:
    if family in fam_to_superfamily:
        return fam_to_superfamily[family]

    if order == "DNA" and family == "DNA":
        return "DNA"

    if order == "LINE" and family == "LINE":
        return "LINE"

    if order == "LTR" and family in ("LTR", "Unknown"):
        return "LTR"

    if order == "Unknown" and family == "Unknown":
        return "Unknown"

    return family


def add_superfamily_column(input_path: str, taxonomy_path: str, output_dir: str) -> str:
    fam_to_superfamily = load_taxonomy_dict(taxonomy_path)
    os.makedirs(output_dir, exist_ok=True)

    base_name = os.path.basename(input_path)
    out_name = (
        base_name[:-4] + ".with_superfamily.tsv"
        if base_name.endswith(".tsv")
        else base_name + ".with_superfamily.tsv"
    )
    output_path = os.path.join(output_dir, out_name)

    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for line_number, line in enumerate(fin):
            line = line.rstrip("\n")

            if not line.strip():
                continue

            cols = line.split("\t")

            if line_number == 0:
                if len(cols) < 2:
                    sys.exit("ERROR: Header has fewer than 2 columns.")

                if cols[0] == "Rclass":
                    cols[0] = "Order"

                fout.write("\t".join([cols[0], "superfamily"] + cols[1:]) + "\n")
                continue

            if len(cols) < 2:
                continue

            order = cols[0]
            family = cols[1]
            superfamily = infer_superfamily(order, family, fam_to_superfamily)

            fout.write("\t".join([order, superfamily, family] + cols[2:]) + "\n")

    return output_path


def collapse_to_superfamily(input_path: str, output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)

    base_name = os.path.basename(input_path)
    out_name = (
        base_name[:-4] + ".collapsed_superfamily.tsv"
        if base_name.endswith(".tsv")
        else base_name + ".collapsed_superfamily.tsv"
    )
    output_path = os.path.join(output_dir, out_name)

    with open(input_path, "r") as fin:
        lines = [line.rstrip("\n") for line in fin if line.strip()]

    if not lines:
        sys.exit("ERROR: Input file is empty.")

    header = lines[0].split("\t")

    try:
        order_idx = header.index("Order")
        superfamily_idx = header.index("superfamily")
        rfam_idx = header.index("Rfam")
    except ValueError:
        sys.exit("ERROR: Required columns not found: Order, superfamily, Rfam")

    bin_start = rfam_idx + 1
    bin_headers = header[bin_start:]
    sums = defaultdict(lambda: [0] * len(bin_headers))

    for line in lines[1:]:
        cols = line.split("\t")

        if len(cols) < bin_start:
            continue

        order = cols[order_idx]
        superfamily = cols[superfamily_idx]
        values = cols[bin_start:]

        if len(values) != len(bin_headers):
            sys.stderr.write(f"WARNING: Skipping row with wrong number of values:\n{line}\n")
            continue

        acc = sums[(order, superfamily)]

        for i, value in enumerate(values):
            value = value.strip()

            if not value:
                continue

            try:
                acc[i] += int(value)
            except ValueError:
                try:
                    acc[i] += float(value)
                except ValueError:
                    sys.stderr.write(
                        f"WARNING: Non-numeric value '{value}' treated as 0.\n"
                    )

    with open(output_path, "w") as fout:
        fout.write("\t".join(["Order", "superfamily"] + bin_headers) + "\n")

        for (order, superfamily), values in sorted(sums.items()):
            fout.write(
                "\t".join([order, superfamily] + [str(value) for value in values]) + "\n"
            )

    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sanitize a RepeatMasker landscape table and collapse values by superfamily."
    )

    parser.add_argument(
        "input_file",
        help="Input .tab file from parseRM or RepeatMasker processing."
    )

    parser.add_argument(
        "-t",
        "--taxonomy-dict",
        default=DEFAULT_TAXONOMY_DICT,
        help="Tab-delimited taxonomy dictionary with family and superfamily columns."
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default="repeatmasker_landscape_processed",
        help="Directory for all output files."
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.input_file):
        sys.exit(f"ERROR: Input file not found: {args.input_file}")

    if not os.path.exists(args.taxonomy_dict):
        sys.exit(f"ERROR: Taxonomy dictionary not found: {args.taxonomy_dict}")

    sanitized_dir = os.path.join(args.output_dir, "tsv")
    superfamily_dir = os.path.join(args.output_dir, "with_superfamily")
    collapsed_dir = os.path.join(args.output_dir, "collapsed_superfamily")

    sanitized_path = sanitize_to_tsv(args.input_file, sanitized_dir)
    superfamily_path = add_superfamily_column(
        sanitized_path,
        args.taxonomy_dict,
        superfamily_dir
    )
    collapse_to_superfamily(superfamily_path, collapsed_dir)


if __name__ == "__main__":
    main()