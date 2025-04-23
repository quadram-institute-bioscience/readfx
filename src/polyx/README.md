# Polytail Trimming 

Test detection and trimming of homopolymer tails from FASTQ sequences. 

## Features

- Processes FASTQ/FASTA files with automatic compression detection
- Filters homopolymer tails of any base (A, C, G, T)
- Supports both single-threaded and multi-threaded processing
- Configurable minimum homopolymer length parameter

## Usage

```bash
./polyx input.fastq.gz > trimmed.fastq
./polyx_threads input.fastq.gz > trimmed_threads.fastq
```

## Build

```bash
make                # Build both versions
make polyx          # Build only single-threaded version
make polyx_threads  # Build only multi-threaded version
```

## Test

```bash
make test        # Generates test data and runs both implementations
```

## Performance

Multi-threaded implementation was found in the example approximately 2.2x faster than single-threaded:

| Command | Mean [s] | Relative |
|:---|---:|---:|
| `./polyx test.fq.gz` | 7.405 | 2.19 |
| `./polyx_threads test.fq.gz` | 3.385 | 1.00 |

## Test Data Generation

Generate synthetic FASTQ data with random homopolymer tails:

```bash
python3 makePolyX.py -n NUM_READS -o output.fastq.gz --max-polytail LENGTH
```