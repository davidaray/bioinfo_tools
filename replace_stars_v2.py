#!/usr/bin/env python3

import argparse
import os

def parse_out_line(line):
    """Parse a RepeatMasker .out line into columns."""
    parts = line.rstrip("\n").split()
    if len(parts) < 15:
        return None
    return parts

def make_key(parts):
    """Create a stable matching key."""
    return (
        parts[0],   # SW score
        parts[1],   # div
        parts[2],   # del
        parts[3],   # ins
        parts[8],   # strand
        parts[9],   # repeat
        parts[10],  # class/family
        parts[11],  # r_begin
        parts[12],  # r_end
        parts[13],  # r_left
        parts[14],  # ID
    )

def find_original_starred(input_dir):
    """Find all lines ending with * in original .out files."""
    starred_keys = set()
    count_files = 0

    for root, _, files in os.walk(input_dir):
        for fname in files:
            if fname.endswith(".out"):
                count_files += 1
                with open(os.path.join(root, fname)) as f:
                    for line in f:
                        if line.rstrip().endswith("*"):
                            parts = parse_out_line(line)
                            if parts:
                                key = make_key(parts)
                                starred_keys.add(key)

    print(f"Identified {count_files} .out files.")
    print(f"Collected {len(starred_keys)} starred lines.")

    return starred_keys

def process_lifted_file(lifted_file, starred_keys, output_file):
    """Append '*' to lifted lines that match starred_keys."""
    with open(lifted_file) as inp, open(output_file, "w") as out:

        for line in inp:
            stripped = line.rstrip("\n")

            # Keep headers unchanged
            if stripped.startswith("SW") or stripped.startswith("score"):
                out.write(line)
                continue

            parts = parse_out_line(stripped)
            if not parts:
                out.write(line)
                continue

            key = make_key(parts)

            if key in starred_keys:
                out.write(stripped + "\t*\n")
            else:
                out.write(line)

def main():
    parser = argparse.ArgumentParser(description="Restore '*' marks to lifted RepeatMasker .out files.")
    parser.add_argument("-d", "--directory", required=True,
                        help="Directory containing original .out files")
    parser.add_argument("-i", "--input", required=True,
                        help="Lifted .out file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output .out file with '*' restored")

    args = parser.parse_args()

    starred_keys = find_original_starred(args.directory)
    process_lifted_file(args.input, starred_keys, args.output)

    print(f"Done. Output written to {args.output}")

if __name__ == "__main__":
    main()
