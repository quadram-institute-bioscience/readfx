# ReadFX - Nim FASTQ/FASTA Parser

[![Nim Tests](https://github.com/quadram-institute-bioscience/readfx/actions/workflows/test.yml/badge.svg)](https://github.com/quadram-institute-bioscience/readfx/actions/workflows/test.yml)

The Nim FASTA/FASTQ parsing library for [SeqFu](https://github.com/telatin/seqfu2).

* [API Documentation](https://quadram-institute-bioscience.github.io/readfx/readfx.html)


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

# Reading paired-end FASTQ files
for pair in readFQPair("sample_R1.fastq.gz", "sample_R2.fastq.gz"):
  echo "R1: ", pair.read1.name, " (", pair.read1.sequence.len, " bp)"
  echo "R2: ", pair.read2.name, " (", pair.read2.sequence.len, " bp)"
```

## Authors

- Original library by Heng Li (kseq.h) and Andreas Wilm ([readfq](https://github.com/andreas-wilm/nimreadfq))
- Updated and maintained by Andrea Telatin and the [Quadram Institute Bioscience](https://www.quadram.ac.uk) Core Bioinformatics team
- Co-authored by [Claude code](https://docs.anthropic.com/en/docs/agents-and-tools/claude-code/overview)