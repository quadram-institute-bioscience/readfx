# ReadFX Methods

This document provides a detailed overview of available methods in the ReadFX library for parsing and manipulating FASTA/FASTQ files.

## Core Types

### FQRecord

Represents a FASTA/FASTQ record with string-based fields.

```nim
type FQRecord* = object
  name*: string      # Sequence identifier
  comment*: string   # Optional comment
  sequence*: string  # Nucleotide sequence
  quality*: string   # Optional quality scores (for FASTQ)
  status*, lastChar*: int # Internal state variables for parsing
```

### FQRecordPtr

Pointer-based representation for efficient parsing.

```nim
type FQRecordPtr* = object
  name*: ptr char     # Pointer to sequence identifier
  comment*: ptr char  # Pointer to optional comment
  sequence*: ptr char # Pointer to nucleotide sequence
  quality*: ptr char  # Pointer to optional quality scores
```

### Bufio[T]

Generic buffered I/O utility for low-level file reading operations.

```nim
type Bufio[T] = tuple[fp: T, buf: string, st, en, sz: int, EOF: bool]
```

## Parsing Methods

### readFQ

```nim
iterator readFQ*(path: string): FQRecord
```

High-level iterator for FASTA/FASTQ parsing.

- **Parameters**: 
  - `path: string` - File path (use "-" for stdin)
- **Returns**: FQRecord objects with string fields
- **Usage**: For clean, easy-to-use code with string manipulation

### readFQPtr

```nim
iterator readFQPtr*(path: string): FQRecordPtr
```

High-performance iterator with pointer-based records.

- **Parameters**: 
  - `path: string` - File path (use "-" for stdin)
- **Returns**: FQRecordPtr objects with C-like pointers
- **Usage**: For maximum performance and low memory usage
- **Note**: Pointers are reused on each iteration!

### readFastx

```nim
proc readFastx*[T](f: var Bufio[T], r: var FQRecord): bool
```

Low-level procedure for custom parsing workflows.

- **Parameters**: 
  - `f: var Bufio[T]` - Buffered input stream
  - `r: var FQRecord` - Record to populate
- **Returns**: Boolean indicating success
- **Usage**: For integration with custom I/O workflows

## File Handling Utilities

### xopen

```nim
proc xopen*[T](fn: string, mode: FileMode = fmRead, sz: int = 0x10000): Bufio[T]
```

Opens a file with buffered I/O.

- **Parameters**:
  - `fn: string` - File path
  - `mode: FileMode` - File mode (default: fmRead)
  - `sz: int` - Buffer size (default: 0x10000)
- **Returns**: Bufio object

### open

```nim
proc open*[T](f: var Bufio[T], fn: string, mode: FileMode = fmRead, sz: int = 0x10000): int
```

Opens a file for buffered I/O.

- **Parameters**: Same as xopen
- **Returns**: Status code (0 for success)

### close

```nim
proc close*[T](f: var Bufio[T]): int
```

Closes a buffered I/O handle.

- **Parameters**:
  - `f: var Bufio[T]` - File handle to close
- **Returns**: Status code

## Sequence Manipulation Functions

### reverseComplement

```nim
proc reverseComplement*(sequence: string): string
```

Reverse complements a DNA sequence.

- **Parameters**:
  - `sequence: string` - Input DNA sequence
- **Returns**: Reverse complemented sequence
- **Example**: `reverseComplement("ATGC")` returns "GCAT"

### reverseComplementRecord (in-place)

```nim
proc reverseComplementRecord*(record: var FQRecord)
```

Reverse complements a sequence record in place.

- **Parameters**:
  - `record: var FQRecord` - Record to modify
- **Note**: Modifies the record in place, also reversing quality if present

### reverseComplementRecord (copy)

```nim
proc reverseComplementRecord*(record: FQRecord): FQRecord
```

Creates a new record with reverse-complemented sequence.

- **Parameters**:
  - `record: FQRecord` - Input record
- **Returns**: New FQRecord with reverse-complemented sequence

### gcContent

```nim
proc gcContent*(sequence: string): float
```

Calculates GC content of a DNA sequence.

- **Parameters**:
  - `sequence: string` - DNA sequence
- **Returns**: GC content as a fraction between 0.0 and 1.0

## Quality-Based Operations

### trimQuality

```nim
proc trimQuality*(quality: string, minQual: int, offset: int = 33): string
```

Trims a quality string based on minimum quality threshold.

- **Parameters**:
  - `quality: string` - Quality string
  - `minQual: int` - Minimum quality value (0-40)
  - `offset: int` - Quality score offset (default: 33)
- **Returns**: Trimmed quality string

### qualityTrim

```nim
proc qualityTrim*(record: var FQRecord, minQual: int, offset: int = 33)
```

Trims a record based on quality scores.

- **Parameters**:
  - `record: var FQRecord` - Record to modify
  - `minQual: int` - Minimum quality value
  - `offset: int` - Quality score offset (default: 33)
- **Note**: Modifies both sequence and quality in place

### maskLowQuality

```nim
proc maskLowQuality*(record: var FQRecord, minQual: int, offset: int = 33, maskChar: char = 'N')
```

Masks sequence positions with low quality scores.

- **Parameters**:
  - `record: var FQRecord` - Record to modify
  - `minQual: int` - Minimum quality value
  - `offset: int` - Quality score offset (default: 33)
  - `maskChar: char` - Character to use for masking (default: 'N')

## Record Manipulation

### subSequence

```nim
proc subSequence*(record: FQRecord, start: int, length: int = -1): FQRecord
```

Extracts a subsequence from a record.

- **Parameters**:
  - `record: FQRecord` - Input record
  - `start: int` - Start position (0-based)
  - `length: int` - Length to extract (-1 for end of sequence)
- **Returns**: New FQRecord with extracted subsequence

### $ (string representation)

```nim
proc `$`*(rec: FQRecord): string
proc `$`*(rec: FQRecordPtr): string 
```

Converts records to string representation in FASTA/FASTQ format.

- **Parameters**:
  - `rec: FQRecord/FQRecordPtr` - Record to convert
- **Returns**: String in FASTA or FASTQ format based on record content

## Other Utilities

### Interval Operations

```nim
type Interval*[S,T] = tuple[st, en: S, data: T, max: S]
proc sort*[S,T](a: var seq[Interval[S,T]])
proc index*[S,T](a: var seq[Interval[S,T]]): int
iterator overlap*[S,T](a: seq[Interval[S,T]], st: S, en: S): Interval[S,T]
```

Interval tree implementation for genomic intervals.

- **Usage**: For interval overlap queries in genomic data

## Usage Examples

### Basic Parsing

```nim
import readfx

# Using readFQ (string-based)
for record in readFQ("sample.fastq.gz"):
  echo record.name, " has length ", record.sequence.len
  
# Using readFQPtr (pointer-based)
for record in readFQPtr("sample.fastq.gz"):
  echo $record.name, " has length ", len($record.sequence)

# Using readFastx (low-level)
var record: FQRecord
var f = xopen[GzFile]("sample.fastq.gz")
defer: f.close()
while f.readFastx(record):
  echo record.name, " has length ", record.sequence.len
```

### Record Manipulation

```nim
# GC content calculation
let gc = gcContent(record.sequence)

# Reverse complement
let rcSeq = reverseComplement(record.sequence)
let rcRecord = reverseComplementRecord(record)

# Quality trimming
qualityTrim(record, 20)  # Trim bases with quality < 20

# Masking low quality bases
maskLowQuality(record, 20)  # Mask bases with quality < 20 as 'N'

# Extract subsequence
let firstTenBases = subSequence(record, 0, 10)
```