#!/usr/bin/env python
"""
A program to generate FASTQ.gz files with a number of random sequences with a poly tail.
The last sequence will be a polytail only.
The quality scores are fixed and designed to make it easy to test if the polytail is
correctly detected.

Optimized using multithreading for faster sequence generation.
"""

import os
import sys
import random
import gzip
import argparse
import concurrent.futures
from threading import Lock
import io
from typing import List, Tuple

def generate_sequence(seq_id: int, max_length: int, max_polytail: int, is_last: bool = False) -> str:
    """Generate a single fastq sequence entry with random sequence and polytail"""
    # For the last sequence, create only a polytail
    if is_last:
        polytail_length = random.randint(1, max_polytail)
        polytail_base = random.choice("ACGT")
        polytail = polytail_base * polytail_length
        qual = "9" * polytail_length
        return f"@seq_{seq_id} poly={polytail_length} seqlen=0 onlytail=true\n{polytail}\n+\n{qual}\n"
        
    # Generate a random poly tail with a random length up to max_polytail
    polytail_length = random.randint(1, max_polytail)
    polytail_base = random.choice("ACGT")
    polytail = polytail_base * polytail_length
    
    # Generate a random sequence of random length up to (max_length - polytail_length)
    seq_length = random.randint(0, max_length - polytail_length)
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
        
    return f"@seq_{seq_id} poly={polytail_length} seqlen={seq_length}\n{full_seq}\n+\n{qual}\n"

def generate_batch(start_id: int, batch_size: int, max_length: int, max_polytail: int, 
                  last_batch: bool = False) -> List[str]:
    """Generate a batch of sequences"""
    sequences = []
    for i in range(batch_size):
        seq_id = start_id + i
        # Check if this is the very last sequence
        is_last = last_batch and i == batch_size - 1
        sequences.append(generate_sequence(seq_id, max_length, max_polytail, is_last))
        
    return sequences

def write_batch_to_buffer(batch: List[str]) -> str:
    """Join a batch of sequences into a single string"""
    return ''.join(batch)

def main():
    parser = argparse.ArgumentParser(description="Make PolyX file (multithreaded)")
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
    parser.add_argument(
        "-t", "--threads", type=int, default=os.cpu_count(), 
        help="Number of threads to use [default: number of CPU cores]"
    )
    parser.add_argument(
        "-b", "--batch-size", type=int, default=10000,
        help="Batch size for sequence generation [default: %(default)s]"
    )
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    random.seed(args.seed)
    
    # Calculate number of batches
    batch_size = min(args.batch_size, args.num_sequences)
    num_batches = (args.num_sequences + batch_size - 1) // batch_size
    
    # Open output file
    output_mode = "wt" if args.output.endswith(".gz") else "w"
    output_file = gzip.open(args.output, output_mode) if args.output.endswith(".gz") else open(args.output, output_mode)
    
    total_sequences = 0
    print(f"Generating {args.num_sequences} sequences using {args.threads} threads", file=sys.stderr)
    
    # Use ThreadPoolExecutor for parallel sequence generation
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        for batch_idx in range(num_batches):
            is_last_batch = batch_idx == num_batches - 1
            start_id = batch_idx * batch_size
            
            # Calculate actual batch size for this batch
            actual_batch_size = min(batch_size, args.num_sequences - start_id)
            
            # Submit batch generation task
            future = executor.submit(
                generate_batch, 
                start_id, 
                actual_batch_size, 
                args.max_length, 
                args.max_polytail,
                is_last_batch
            )
            
            # Get the batch results and write to file
            batch_sequences = future.result()
            batch_text = write_batch_to_buffer(batch_sequences)
            output_file.write(batch_text)
            
            total_sequences += actual_batch_size
            percent = (total_sequences / args.num_sequences) * 100
            print(f"Generated batch {batch_idx+1}/{num_batches} - {total_sequences} of {args.num_sequences} sequences ({percent:.2f}%)", file=sys.stderr)
    
    # Close the output file
    output_file.close()
    
    print(f"Successfully generated {args.num_sequences} sequences with max length {args.max_length} and max polytail {args.max_polytail} in {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()