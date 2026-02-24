# ReadFX Methods Reference

## Parsing Iterators

### `readFQ`

```nim
iterator readFQ*(path: string): FQRecord
```

High-level iterator yielding `FQRecord` objects with string fields.

- `path` — file path; use `"-"` for stdin
- Supports FASTA, FASTQ, and gzipped variants

### `readFQPtr`

```nim
iterator readFQPtr*(path: string): FQRecordPtr
```

High-performance iterator yielding pointer-based records.

- `path` — file path; use `"-"` for stdin
- **Pointers are reused on each iteration** — copy data if you need to retain it

### `readFQPair`

```nim
iterator readFQPair*(path1, path2: string, checkNames: bool = false): FQPair
```

Synchronized paired-end iterator yielding `FQPair` objects.

- `path1` — R1 file; use `"-"` for stdin
- `path2` — R2 file (stdin not supported for both files)
- `checkNames` — if `true`, validates that read names match (strips `/1`/`/2` suffixes)
- Raises `IOError` if files have different lengths
- Raises `ValueError` on name mismatch when `checkNames = true`

### `readFastx`

```nim
proc readFastx*[T](f: var Bufio[T], r: var FQRecord): bool
```

Low-level reader populating `r` from a buffered stream.

- Returns `false` at EOF or on error
- `r.status` contains the sequence length (>0) or an error code on failure

---

## File I/O

### `xopen`

```nim
proc xopen*[T](fn: string, mode: FileMode = fmRead, sz: int = 0x10000): Bufio[T]
```

Opens a file and returns a `Bufio[T]`. Typically used as `xopen[GzFile](fn)`.

### `open`

```nim
proc open*[T](f: var Bufio[T], fn: string, mode: FileMode = fmRead, sz: int = 0x10000): int
```

Opens a file into an existing `Bufio[T]`. Returns 0 on success, -1 on failure.

### `close`

```nim
proc close*[T](f: var Bufio[T]): int
```

Closes a `Bufio[T]` handle. Returns the underlying close status code.

### `eof`

```nim
proc eof*[T](f: Bufio[T]): bool
```

Returns `true` if the buffer is exhausted and EOF has been reached.

### `readLine`

```nim
proc readLine*[T](f: var Bufio[T], buf: var string): bool
```

Reads one line from a `Bufio[T]`. Returns `false` at EOF.

---

## Formatting

### `$` (string conversion)

```nim
proc `$`*(rec: FQRecord): string
proc `$`*(rec: FQRecordPtr): string
```

Formats a record as FASTQ (if quality is present) or FASTA.

### `fafmt`

```nim
proc fafmt*(rec: FQRecord, width: int = 60): string
```

Formats a record as wrapped FASTA with lines of `width` characters.

---

## Sequence Operations

### `revCompl`

```nim
proc revCompl*(sequence: string): string      # Returns new reverse-complemented string
proc revCompl*(record: var FQRecord)          # In-place (also reverses quality)
proc revCompl*(record: FQRecord): FQRecord    # Returns new record
```

Reverse complements a DNA sequence. The in-place variant also reverses the quality string.

### `gcContent`

```nim
proc gcContent*(sequence: string): float
proc gcContent*(record: FQRecord): float
```

Returns the GC fraction (0.0–1.0).

### `composition`

```nim
proc composition*(record: FQRecord): SeqComp
```

Returns a `SeqComp` object with per-base counts (A, C, G, T, N, Other) and GC fraction.

### `subSequence`

```nim
proc subSequence*(record: FQRecord, start: int, length: int = -1): FQRecord
```

Extracts a subsequence starting at `start` (0-based). `length = -1` means to the end. Returns a new `FQRecord`.

### `trimStart`

```nim
proc trimStart*(record: FQRecord, bases: int): FQRecord
```

Removes `bases` characters from the 5' end. Returns a new `FQRecord`.

### `trimEnd`

```nim
proc trimEnd*(record: FQRecord, bases: int): FQRecord
```

Removes `bases` characters from the 3' end. Returns a new `FQRecord`.

---

## Quality Operations

### `qualCharToInt`

```nim
proc qualCharToInt*(c: char, offset: int = 33): int
```

Converts a quality character to its integer Phred score.

### `qualIntToChar`

```nim
proc qualIntToChar*(q: int, offset: int = 33): char
```

Converts a Phred integer score to its quality character.

### `avgQuality`

```nim
proc avgQuality*(record: FQRecord, offset: int = 33): float
proc avgQuality*(quality: string, offset: int = 33): float
```

Returns the mean Phred quality score.

### `trimQuality`

```nim
proc trimQuality*(quality: string, minQual: int, offset: int = 33): string
```

Trims trailing low-quality bases from a quality string.

### `qualityTrim`

```nim
proc qualityTrim*(record: var FQRecord, minQual: int, offset: int = 33)
```

Trims both sequence and quality in place based on a minimum quality threshold.

### `maskLowQuality`

```nim
proc maskLowQuality*(record: var FQRecord, minQual: int, offset: int = 33, maskChar: char = 'N')
```

Replaces bases with quality below `minQual` with `maskChar` (default `'N'`).

---

## IUPAC Primer Matching

### `matchIUPAC`

```nim
proc matchIUPAC*(primerBase, referenceBase: char): bool
```

Returns `true` if `primerBase` (supporting IUPAC ambiguity codes) matches `referenceBase`.

---

## Interval Tree

```nim
type Interval*[S,T] = tuple[st, en: S, data: T, max: S]

proc sort*[S,T](a: var seq[Interval[S,T]])
proc index*[S,T](a: var seq[Interval[S,T]]): int
iterator overlap*[S,T](a: seq[Interval[S,T]], st, en: S): Interval[S,T]
```

Sort and index a sequence of intervals, then query overlapping intervals efficiently.

---

## Usage Examples

### Basic parsing

```nim
import readfx

for record in readFQ("sample.fastq.gz"):
  echo record.name, " (", record.sequence.len, " bp)"
```

### Paired-end reads

```nim
for pair in readFQPair("R1.fastq.gz", "R2.fastq.gz", checkNames = true):
  echo pair.read1.name, " / ", pair.read2.name
```

### Sequence manipulation

```nim
for record in readFQ("sample.fastq.gz"):
  let gc = gcContent(record.sequence)
  let rc = revCompl(record.sequence)
  let comp = composition(record)
  echo record.name, " GC=", gc, " A=", comp.A, " N=", comp.N
```

### Quality trimming

```nim
var record: FQRecord
var f = xopen[GzFile]("sample.fastq.gz")
defer: f.close()
while f.readFastx(record):
  qualityTrim(record, 20)      # trim 3' low-quality bases
  maskLowQuality(record, 15)   # mask remaining low-quality positions
  echo $record
```
