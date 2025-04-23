# ReadFX

ReadFX is a Nim library for parsing FASTQ/FASTA files with high performance.

[![Nim Tests](https://github.com/quadram-institute-bioscience/readfx/actions/workflows/test.yml/badge.svg)](https://github.com/quadram-institute-bioscience/readfx/actions/workflows/test.yml)

ReadFX is a Nim library for parsing FASTA and FASTQ files (collectively known as FASTX). It combines two approaches to sequence file parsing:

1. A Nim wrapper around Heng Li's kseq.h C library (klib)
2. A native Nim implementation (nimklib)

Both implementations provide efficient parsing with different tradeoffs in terms of performance and memory usage.

## Installation

```bash
nimble install readfx
```

## Usage

```nim
import readfx

# Using the C wrapper (klib)
for record in readFQ("example.fastq"):
  echo record.name, ": ", record.sequence.len
  
# Using pointer-based version (more efficient, reuses memory)
for record in readFQPtr("example.fastq.gz"):
  echo $record.name, ": ", $record.sequence.len
  
# Using the native Nim implementation
var r: FQRecord
var f = xopen[GzFile]("example.fastq.gz")
defer: f.close()
while f.readFastx(r):
  echo r.name, ": ", r.sequence.len
```

## Key Features

- Support for both FASTA and FASTQ formats
- Automatic handling of gzipped files
- Memory-efficient parsing options
- Support for reading from stdin using "-" as filename

## Documentation

- [Methods](METHODS.md) - Detailed documentation of all available methods
  - [Utils methods](UTILS.md) - Utility methods for manipulating sequence records
- [Data Structure](DATA_STRUCTURE.md) - Overview of internal data structures
- [Parsing algorithms](PARSING.md) - Details about the parsing implementation

- [Repository Structure](REPO_STRUCTURE.md) - Overview of the project organization
- [API Reference](API.md) - Detailed API documentation
