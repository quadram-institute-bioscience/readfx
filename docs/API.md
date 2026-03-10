## API Reference

### Types

#### `FQRecordPtr`
Pointer-based record for high-performance streaming. Pointers are valid only during the current iteration.

```nim
type FQRecordPtr* = object
  name*: ptr char      # Sequence identifier
  nameLen*: int        # Cached name length (excluding trailing NUL)
  comment*: ptr char   # Optional comment
  commentLen*: int     # Cached comment length (excluding trailing NUL)
  sequence*: ptr char  # Nucleotide sequence
  sequenceLen*: int    # Cached sequence length (excluding trailing NUL)
  quality*: ptr char   # Quality scores (nil for FASTA)
  qualityLen*: int     # Cached quality length (excluding trailing NUL)
```

Nil pointers carry a cached length of `0`. Pointer fields are valid only until
the next iterator advance.

#### `FQPairPtr`
Pointer-based paired-end record containing two `FQRecordPtr` objects.

```nim
type FQPairPtr* = object
  read1*: FQRecordPtr
  read2*: FQRecordPtr
```

#### `FQRecord`
String-based record. Safe to store and manipulate after iteration.

```nim
type FQRecord* = object
  name*: string
  comment*: string
  sequence*: string
  quality*: string        # Empty for FASTA records
  status*, lastChar*: int # Internal parsing state
```

#### `FQPair`
Paired-end record containing two `FQRecord` objects.

```nim
type FQPair* = object
  read1*: FQRecord   # Forward read (R1)
  read2*: FQRecord   # Reverse read (R2)
```

#### `SeqComp`
Nucleotide composition statistics.

```nim
type SeqComp* = object
  A*, C*, G*, T*, N*, Other*: int
  GC*: float   # Fraction of G+C bases
```

#### `Bufio[T]`
Generic buffered reader. Typically used as `Bufio[GzFile]`.

#### `Interval[S,T]`
Genomic interval for use with the built-in interval tree.

---

### Parsing Iterators

#### `readFQ`
```nim
iterator readFQ*(path: string): FQRecord
```
Yields `FQRecord` objects with string fields. Use `"-"` for stdin.

#### `readFQPtr`
```nim
iterator readFQPtr*(path: string): FQRecordPtr
```
Yields pointer-based records. Faster than `readFQ` but pointers are reused on each iteration. Cached lengths are populated on every record.

#### `readFQPairPtr`
```nim
iterator readFQPairPtr*(path1, path2: string, checkNames: bool = false): FQPairPtr
```
Yields synchronized pointer-based paired-end records. Cached lengths are
available on both mates.

#### `readFQInterleavedPairPtr`
```nim
iterator readFQInterleavedPairPtr*(path: string, checkNames: bool = false): FQPairPtr
```
Yields pointer-based paired-end records from one interleaved FASTQ stream.
`read1` is scratch-backed and both mates remain valid until the next iterator
advance. FASTQ only. Raises `IOError` on incomplete trailing pairs.

#### `readFQPair`
```nim
iterator readFQPair*(path1, path2: string, checkNames: bool = false): FQPair
```
Yields synchronized paired-end records. Raises `IOError` on length mismatch; raises `ValueError` on name mismatch when `checkNames = true`.

#### `readFastx`
```nim
proc readFastx*[T](f: var Bufio[T], r: var FQRecord): bool
```
Low-level reader. Returns `false` at EOF. Used directly when managing your own file handles.

---

### File I/O

```nim
proc xopen*[T](fn: string, mode: FileMode = fmRead, sz: int = 0x10000): Bufio[T]
proc open*[T](f: var Bufio[T], fn: string, mode: FileMode = fmRead, sz: int = 0x10000): int
proc close*[T](f: var Bufio[T]): int
proc eof*[T](f: Bufio[T]): bool
proc readLine*[T](f: var Bufio[T], buf: var string): bool
```

---

### Formatting

```nim
proc `$`*(rec: FQRecord): string      # FASTA or FASTQ string
proc `$`*(rec: FQRecordPtr): string   # FASTA or FASTQ string
proc fafmt*(rec: FQRecord, width: int = 60): string  # Wrapped FASTA
```

---

### Sequence Utilities

```nim
proc revCompl*(sequence: string): string          # Reverse complement of a string
proc revCompl*(record: var FQRecord)              # In-place reverse complement
proc revCompl*(record: FQRecord): FQRecord        # Returns new reverse-complemented record
proc gcContent*(sequence: string): float          # GC fraction (0.0–1.0)
proc gcContent*(record: FQRecord): float
proc composition*(record: FQRecord): SeqComp      # Full nucleotide composition
proc subSequence*(record: FQRecord, start: int, length: int = -1): FQRecord
proc trimStart*(record: FQRecord, bases: int): FQRecord
proc trimEnd*(record: FQRecord, bases: int): FQRecord
```

### Quality Utilities

```nim
proc qualCharToInt*(c: char, offset: int = 33): int
proc qualIntToChar*(q: int, offset: int = 33): char
proc avgQuality*(record: FQRecord, offset: int = 33): float
proc avgQuality*(quality: string, offset: int = 33): float
proc trimQuality*(quality: string, minQual: int, offset: int = 33): string
proc qualityTrim*(record: var FQRecord, minQual: int, offset: int = 33)
proc maskLowQuality*(record: var FQRecord, minQual: int, offset: int = 33, maskChar: char = 'N')
```

### IUPAC Primer Matching

```nim
proc matchIUPAC*(primerBase, referenceBase: char): bool
```

---

### Interval Tree

```nim
type Interval*[S,T] = tuple[st, en: S, data: T, max: S]
proc sort*[S,T](a: var seq[Interval[S,T]])
proc index*[S,T](a: var seq[Interval[S,T]]): int
iterator overlap*[S,T](a: seq[Interval[S,T]], st, en: S): Interval[S,T]
```
