# ReadFX

ReadFX is a Nim library for parsing FASTQ/FASTA files with high performance.

[![Nim Tests](https://github.com/quadram-institute-bioscience/readfx/actions/workflows/test.yml/badge.svg)](https://github.com/quadram-institute-bioscience/readfx/actions/workflows/test.yml)

ReadFX provides efficient parsing and manipulation of FASTA and FASTQ files (FASTX). It wraps Heng Li's `kseq.h` C library for maximum throughput and also includes a native Nim implementation.

## Installation

```bash
nimble install readfx
```

## Quick Start

```nim
import readfx

# Simple string-based iteration
for record in readFQ("example.fastq.gz"):
  echo record.name, ": ", record.sequence.len

# High-performance pointer-based iteration (no string copies)
for record in readFQPtr("example.fastq.gz"):
  echo $record.name, ": ", len($record.sequence)

# Low-level buffered reader
var r: FQRecord
var f = xopen[GzFile]("example.fastq.gz")
defer: f.close()
while f.readFastx(r):
  echo r.name, ": ", r.sequence.len

# Paired-end reads
for pair in readFQPair("sample_R1.fastq.gz", "sample_R2.fastq.gz"):
  echo pair.read1.name, " / ", pair.read2.name
```

## Key Features

- FASTA and FASTQ format support (auto-detected)
- Transparent gzip decompression
- Stdin support via `"-"` as filename
- Three parsing APIs with different performance/convenience tradeoffs
- Paired-end read support with optional name validation
- Sequence utilities: reverse complement, GC content, quality trimming, subsequence extraction
- IUPAC primer matching

## Documentation

- [Methods](METHODS.md) - All available procedures and iterators
- [Data Structures](DATA_STRUCTURE.md) - Type definitions (`FQRecord`, `FQPair`, `SeqComp`, ...)
- [Parsing Methods](PARSING.md) - Comparison of `readFQ`, `readFQPtr`, `readFastx`, `readFQPair`
- [Sequence Utilities](FQRECORD_UTILS.md) - Sequence manipulation functions
- [Utility Functions Reference](UTILS.md) - Full utility function reference
- [Repository Structure](REPO_STRUCTURE.md) - Project layout
- [API Reference](API.md) - Concise API summary
