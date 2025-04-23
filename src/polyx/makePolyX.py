#!/usr/bin/env python
"""

"""

import os
import sys
import random
import gzip
import argparse

def main():
    args = argparse.ArgumentParser(description="Make PolyX file")
    args.add_argument(
        "-o", "--output", type=str, required=True, help="Output file"
    )
    args.add_argument(
        "-m", "--max-length", type=int, default=1000, help="Max length of the sequence"
    )
    args.add_argument(
        "-p", "--max-polytail", type=int, default=100, help="Max length of the poly tail"
    )
    args.add_argument(
        "-n", "--num-sequences", type=int, default=1000, help="Number of sequences to generate"
    )
    args = args.parse_args()

    # Output file is gzipped
    if args.output.endswith(".gz"):
        output_file = gzip.open(args.output, "wt")
    else:
        output_file = open(args.output, "w")
    
    for i in range(args.num_sequences  - 1):
        # Generate a random poly tail with a random length up to max_polytail
        polytail_length = random.randint(1, args.max_polytail)
        polytail = random.choice("ACGT") * polytail_length
        # Generate a random sequence
        seq = "".join(random.choices("ACGT", k=args.max_length-len(polytail)))
        qual = "I" * (len(seq) - 1) + "!" + "9" * len(polytail)
        # Write the sequence to the output file
        output_file.write(f"@seq_{i} poly={polytail_length} seqlen={len(seq)}\n{seq}{polytail}\n+\n{qual}\n")

    # last sequence to be a single polytail
    polytail_length = random.randint(1, args.max_polytail)
    polytail = random.choice("ACGT") * polytail_length
    qual = "9" * len(polytail)
    output_file.write(f"@seq_{args.num_sequences-1} poly={polytail_length} seqlen=0\n{polytail}\n+\n{qual}\n")
    print(f"Generated {args.num_sequences} sequences with max length {args.max_length} and max polytail {args.max_polytail} in {args.output}", file=sys.stderr)
    # Close the output file
    if args.output.endswith(".gz"):
        output_file.close()
    else:
        output_file.close()

if __name__ == "__main__":

    main()