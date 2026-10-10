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
  echo $record.name, ": ", record.sequenceLen
  
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

# Pointer-based paired-end iteration (zero-copy, pointers are reused)
for pair in readFQPairPtr("sample_R1.fastq.gz", "sample_R2.fastq.gz"):
  echo "Pair length: ", pair.read1.sequenceLen + pair.read2.sequenceLen

# Pointer-based interleaved paired-end iteration
for pair in readFQInterleavedPairPtr("sample.interleaved.fastq.gz", checkNames = true):
  echo "Pair length: ", pair.read1.sequenceLen + pair.read2.sequenceLen

# String-based interleaved paired-end iteration (safe to store records)
for pair in readFQInterleavedPair("sample.interleaved.fastq.gz", checkNames = true):
  echo "R1: ", pair.read1.name, " / R2: ", pair.read2.name
```

`FQRecordPtr` now includes cached `nameLen`, `commentLen`, `sequenceLen`, and
`qualityLen` fields. These lengths exclude the trailing NUL terminator and
avoid repeated `cstring` scans in pointer-heavy code. Because `FQRecordPtr`
layout changed, downstream projects must be recompiled after upgrading to this
release.

## Paired Output

Use `pairedFastxWriter` for two files or `interleavedFastxWriter` for one stream.
Both accept `FQPair` and `FQPairPtr` through `writePair`:

```nim
import readfx

var split = pairedFastxWriter(
  fileDestination("out_R1.fastq.gz"),
  fileDestination("out_R2.fastq.gz"),
  compression = true
)
try:
  for pair in readFQPairPtr("in_R1.fastq.gz", "in_R2.fastq.gz"):
    split.writePair(pair)
finally:
  split.close()

var interleaved = interleavedFastxWriter(
  fileDestination("out.interleaved.fastq.gz"), compression = true
)
try:
  for pair in readFQPair("in_R1.fastq.gz", "in_R2.fastq.gz"):
    interleaved.writePair(pair)
finally:
  interleaved.close()
```

The default format is FASTQ; `format = fxfFasta` emits paired FASTA records.
Compression is explicit, independent of filename extensions. Buffer size,
compression level/threads, and FASTA wrapping use the same options as
`FastxWriter`. Interleaved output defaults to stdout; split output permits at
most one stdout destination and rejects identical file destinations.

Both mates are validated before either is buffered. Optional `checkNames = true`
checks matching IDs, accepting `/1` and `/2` (or ` 1` and ` 2`) mate suffixes.
Pointer data is consumed immediately; write each pair before advancing its
reader. Do not copy open writer objects.

`flush` writes application buffers; `close` finalizes gzip output and closes
both destinations. Split writes are not atomic: an I/O error may leave unequal
mate counts. After a write/flush error, subsequent writes and flushes fail;
call `close` and treat the outputs as incomplete. Validation errors leave the
writer usable. Construction opens/truncates destinations in order; if opening
R2 fails, R1 is closed but its previous contents are not restored.

## Authors

- Original library by Heng Li (kseq.h) and Andreas Wilm ([readfq](https://github.com/andreas-wilm/nimreadfq))
- Updated and maintained by Andrea Telatin and the [Quadram Institute Bioscience](https://www.quadram.ac.uk) Core Bioinformatics team
- Co-authored by [Claude code](https://docs.anthropic.com/en/docs/agents-and-tools/claude-code/overview)
