
#!/usr/bin/env python3

"""
Assembly contiguity sensitivity analysis.

Artificially fragment genome assemblies to match the empirical contig
counts of two reference assemblies (IN and CA).

Experimental design:
    - Fragment UK and IT to match the IN contig count.
    - Fragment UK, IT, and IN to match the CA contig count.
    - Retain unchanged copies of IN and CA as reference controls.

Fragmentation is performed by repeatedly splitting the longest
available sequence at its midpoint until the target contig count
is reached.


Usage:
    python fragment_genomes.py \
        --uk path/to/uk.fasta \
        --it path/to/it.fasta \
        --in-genome path/to/in.fasta \
        --ca path/to/ca.fasta \
        --output-dir fragmented_genomes
"""

import argparse
import heapq
import shutil
from pathlib import Path


# ---------------------------------------------------------------------
# FASTA utilities
# ---------------------------------------------------------------------

def read_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file."""

    header = None
    seq_parts = []

    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip()

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)

                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

        if header is not None:
            yield header, "".join(seq_parts)


def count_contigs(path):
    """Count FASTA records without loading the genome into memory."""

    count = 0

    with open(path, "r") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1

    return count


def write_sequence(handle, sequence, width=80):
    """Write sequence using fixed-width FASTA lines."""

    for i in range(0, len(sequence), width):
        handle.write(sequence[i:i + width] + "\n")


# ---------------------------------------------------------------------
# Fragmentation
# ---------------------------------------------------------------------

def fragment_genome(input_path, output_path, target_contigs):
    """
    Fragment an assembly until it contains target_contigs records.

    The longest available fragment is repeatedly split at its
    midpoint. No nucleotide sequence is added, removed, shuffled,
    or modified.
    """

    original_contigs = count_contigs(input_path)

    print(
        f"Fragmenting {input_path.name}: "
        f"{original_contigs:,} -> {target_contigs:,} contigs"
    )

    if original_contigs > target_contigs:
        raise ValueError(
            f"{input_path.name} already contains "
            f"{original_contigs:,} contigs, exceeding "
            f"the target of {target_contigs:,}."
        )

    if original_contigs == target_contigs:
        shutil.copy2(input_path, output_path)
        return

    # Heap entries:
    # (-length, unique_id, original_header, start, end, sequence)
    #
    # Negative lengths implement a max heap.

    heap = []
    unique_id = 0

    for header, sequence in read_fasta(input_path):

        heapq.heappush(
            heap,
            (
                -len(sequence),
                unique_id,
                header,
                0,
                len(sequence),
                sequence,
            ),
        )

        unique_id += 1

    current_count = len(heap)

    while current_count < target_contigs:

        neg_length, _, header, start, end, sequence = (
            heapq.heappop(heap)
        )

        length = -neg_length

        if length < 2:
            raise RuntimeError(
                "Cannot generate additional fragments: "
                "the longest remaining sequence is shorter "
                "than 2 bp."
            )

        midpoint = length // 2

        left_seq = sequence[:midpoint]
        right_seq = sequence[midpoint:]

        genomic_midpoint = start + midpoint

        heapq.heappush(
            heap,
            (
                -len(left_seq),
                unique_id,
                header,
                start,
                genomic_midpoint,
                left_seq,
            ),
        )

        unique_id += 1

        heapq.heappush(
            heap,
            (
                -len(right_seq),
                unique_id,
                header,
                genomic_midpoint,
                end,
                right_seq,
            ),
        )

        unique_id += 1
        current_count += 1

    # Sort fragments by original sequence name and coordinate.

    fragments = sorted(
        heap,
        key=lambda x: (x[2], x[3])
    )

    with open(output_path, "w") as out:

        for _, _, header, start, end, sequence in fragments:

            # Coordinates are 1-based and inclusive.

            new_header = (
                f"{header}__fragment_{start + 1}_{end}"
            )

            out.write(f">{new_header}\n")

            write_sequence(out, sequence)


# ---------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Fragment genome assemblies to match empirical "
            "reference assembly contig counts."
        )
    )

    parser.add_argument(
        "--uk",
        required=True,
        type=Path,
        help="UK genome FASTA"
    )

    parser.add_argument(
        "--it",
        required=True,
        type=Path,
        help="IT genome FASTA"
    )

    parser.add_argument(
        "--in-genome",
        required=True,
        type=Path,
        help="IN genome FASTA"
    )

    parser.add_argument(
        "--ca",
        required=True,
        type=Path,
        help="CA genome FASTA"
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory"
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    args = parse_args()

    genomes = {
        "UK": args.uk,
        "IT": args.it,
        "IN": args.in_genome,
        "CA": args.ca,
    }

    output_dir = args.output_dir

    in_dir = output_dir / "IN_fragmentation"
    ca_dir = output_dir / "CA_fragmentation"

    # Verify input files.

    for name, path in genomes.items():

        if not path.is_file():
            raise FileNotFoundError(
                f"{name} genome not found: {path}"
            )

    # Create output directories.

    in_dir.mkdir(parents=True, exist_ok=True)
    ca_dir.mkdir(parents=True, exist_ok=True)

    # Determine empirical contig counts.

    contig_counts = {
        name: count_contigs(path)
        for name, path in genomes.items()
    }

    print("Original assembly contig counts:")

    for name in ["UK", "IT", "IN", "CA"]:
        print(f"  {name}: {contig_counts[name]:,}")

    print()

    in_target = contig_counts["IN"]
    ca_target = contig_counts["CA"]

    output_files = []

    # -------------------------------------------------------------
    # IN fragmentation experiment
    # -------------------------------------------------------------

    print("Creating IN-based fragmentation set...")

    for name in ["UK", "IT"]:

        output_path = (
            in_dir / f"{name}_fragmented_to_IN.fasta"
        )

        fragment_genome(
            genomes[name],
            output_path,
            in_target,
        )

        output_files.append(output_path)

    in_control = in_dir / "IN_unfragmented.fasta"

    shutil.copy2(
        genomes["IN"],
        in_control,
    )

    output_files.append(in_control)

    print()

    # -------------------------------------------------------------
    # CA fragmentation experiment
    # -------------------------------------------------------------

    print("Creating CA-based fragmentation set...")

    for name in ["UK", "IT", "IN"]:

        output_path = (
            ca_dir / f"{name}_fragmented_to_CA.fasta"
        )

        fragment_genome(
            genomes[name],
            output_path,
            ca_target,
        )

        output_files.append(output_path)

    ca_control = ca_dir / "CA_unfragmented.fasta"

    shutil.copy2(
        genomes["CA"],
        ca_control,
    )

    output_files.append(ca_control)

    # -------------------------------------------------------------
    # Report generated files
    # -------------------------------------------------------------

    print()
    print("=" * 70)
    print("Generated FASTA files:")
    print("=" * 70)

    for path in output_files:
        print(path.resolve())


if __name__ == "__main__":
    main()
