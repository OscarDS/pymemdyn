#!/usr/bin/env python3

import argparse


def remove_chain_id(input_file, output_file):
    """Remove chain IDs from ATOM/HETATM records in a PDB file."""

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 22:
                line = line[:21] + " " + line[22:]

            outfile.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove chain IDs from a PDB file while preserving its format."
    )
    parser.add_argument("input", help="Input PDB file")
    parser.add_argument("output", help="Output PDB file")

    args = parser.parse_args()
    remove_chain_id(args.input, args.output)

    print(f"Chain IDs removed. Output saved to {args.output}")
