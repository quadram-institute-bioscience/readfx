#!/usr/bin/env python
"""
A program to generate FASTQ.gz files with a number of random sequences with a poly tail.
The last sequence will be a polytail only.
The quality scores are fixed and designed to make it easy to test if the polytail is
correctly detected.
"""

import os
import sys
import random
import gzip
import argparse

def main():
    parser = argparse.ArgumentParser(description="Make PolyX file")
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output file (will be gzipped if ending in .gz)"
    )
    parser.add_argument(
        "-m", "--max-length", type=int, default=1000, help="Max length of the sequence [default: %(default)s]"
    )
    parser.add_argument(
        "-p", "--max-polytail", type=int, default=100, help="Max length of the poly tail [default: %(default)s]"
    )
    parser.add_argument(
        "-n", "--num-sequences", type=int, default=1000, help="Number of sequences to generate"
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=0, help="Random seed for reproducibility [default: %(default)s]"
    )
    args = parser.parse_args()

    # Output file is gzipped
    if args.output.endswith(".gz"):
        output_file = gzip.open(args.output, "wt")
    else:
        output_file = open(args.output, "w")
    
    for i in range(args.num_sequences - 1):
        if i % 10000 == 0:
            percent = (i / args.num_sequences) * 100
            print(f"Generating sequence {i} of {args.num_sequences} ({percent:.2f}%)", file=sys.stderr)
        
        # Generate a random poly tail with a random length up to max_polytail
        polytail_length = random.randint(1, args.max_polytail)
        polytail_base = random.choice("ACGT")
        polytail = polytail_base * polytail_length
        
        # Generate a random sequence of random length up to (max_length - polytail_length)
        seq_length = random.randint(0, args.max_length - polytail_length)
        seq = "".join(random.choices("ACGT", k=seq_length))
        
        # Combine sequence and polytail
        full_seq = seq + polytail
        
        # Generate quality string matching the sequence length
        # Using I (high quality) for main sequence and 9 (lower quality) for polytail
        # with ! as a quality marker between them
        if seq_length > 0:
            qual = "I" * (seq_length - 1) + "<" + ">" + "9" * (polytail_length - 1)
        else:
            qual = "9" * polytail_length
            
        # Write the sequence to the output file
        output_file.write(f"@seq_{i} poly={polytail_length} seqlen={seq_length}\n{full_seq}\n+\n{qual}\n")

    # Last sequence to be a single polytail
    polytail_length = random.randint(1, args.max_polytail)
    polytail_base = random.choice("ACGT")
    polytail = polytail_base * polytail_length
    qual = "9" * polytail_length
    output_file.write(f"@seq_{args.num_sequences-1} poly={polytail_length} seqlen=0 onlytail=true\n{polytail}\n+\n{qual}\n")

    print(f"Generated {args.num_sequences} sequences with max length {args.max_length} and max polytail {args.max_polytail} in {args.output}", file=sys.stderr)
    
    # Close the output file
    output_file.close()

if __name__ == "__main__":
    main()