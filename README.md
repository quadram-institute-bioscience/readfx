# ReadFX - Nim FASTQ/FASTA Parser

![Nim Tests](https://github.com/telatin/readfx/workflows/Nim%20Tests/badge.svg)

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

## API Reference

### C Wrapper (klib)

#### Types

- `FQRecord*`: Main type for FASTQ/FASTA records with string fields
  - `name*`: Sequence name/ID (string)
  - `comment*`: Optional comment (string)
  - `sequence*`: The sequence (string)
  - `quality*`: Quality string for FASTQ (string)
  - `status`, `lastChar`: Internal fields

- `FQRecordPtr*`: Pointer-based version for more efficient memory usage
  - `name*`: Sequence name/ID (ptr char)
  - `comment*`: Optional comment (ptr char)
  - `sequence*`: The sequence (ptr char)
  - `quality*`: Quality string for FASTQ (ptr char)

#### Functions and Iterators

- `iterator readFQ*(path: string): FQRecord`: 
  Parses FASTQ/FASTA records from a file, converting to strings

- `iterator readFQPtr*(path: string): FQRecordPtr`: 
  Memory-efficient version using pointers (faster but requires careful handling)

- `proc `$`*(rec: FQRecord): string`: 
  Formats a record as a FASTA/FASTQ string

- `proc `$`*(rec: FQRecordPtr): string`: 
  Formats a pointer-based record as a FASTA/FASTQ string

### Native Nim Implementation (nimklib)

#### Types

- `FQRecord*`: Record type for FASTQ/FASTA data
  - `sequence`, `quality`, `name`, `comment`: String fields
  - `status`, `lastChar`: Internal fields

- `Bufio*[T]`: Buffered reader for efficient file reading
  - Generic over file type to support both regular and gzipped files

#### Functions

- `proc readFastx*[T](f: var Bufio[T], r: var FQRecord): bool`: 
  Parses a single FASTQ/FASTA record from a buffer

- `proc xopen*[T](fn: string, mode: FileMode = fmRead, sz: int = 0x10000): Bufio[T]`: 
  Opens a file with buffered reading

- `proc close*[T](f: var Bufio[T]): int`: 
  Closes a buffered file

- `proc eof*[T](f: Bufio[T]): bool`: 
  Checks if the file has reached EOF

- `proc readLine*[T](f: var Bufio[T], buf: var string): bool`: 
  Reads a line from a buffered file

### Additional Utilities

- `GzFile*`: Type for handling gzipped files
- `Interval*[S,T]`: Generic interval type with interval tree operations
- Various buffer manipulation functions

## Performance

ReadFX is designed for efficiency, with the pointer-based `readFQPtr` implementation offering the best performance for most use cases. The native Nim implementation provides more flexibility at a slight performance cost.

## Todo

[ ]  The two libraries use different definitions for FQRecord - one as an object and one as a tuple with different field orders. We
  should unify them to have a single consistent type.
[ ] 