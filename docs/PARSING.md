# FASTX Parsing in ReadFX

ReadFX provides six primary methods for parsing FASTA and FASTQ files:

1. **`readFQ`** — String-based iterator (convenient, safe)
2. **`readFQPtr`** — Pointer-based iterator (fastest, requires care)
3. **`readFastx`** — Low-level buffered reader (most flexible)
4. **`readFQPair`** — Synchronized paired-end iterator
5. **`readFQInterleavedPairPtr`** — Pointer-based interleaved paired-end iterator
6. **`readFQInterleavedPair`** — String-based interleaved paired-end iterator

## Comparison at a Glance

| Method | Memory | Speed | Ease of Use | Use case |
|---------------|----------|-----------|----------------|----------------------------|
| `readFQ` | Higher | Good | Excellent | General use |
| `readFQPtr` | Low | Excellent | Moderate | High-throughput streaming |
| `readFastx` | Custom | Excellent | Requires setup | Custom I/O workflows |
| `readFQPair` | Moderate | Good | Excellent | Separate paired-end files |
| `readFQInterleavedPairPtr` | Low | Excellent | Moderate | Interleaved paired-end streams |
| `readFQInterleavedPair` | Moderate | Good | Excellent | Interleaved paired-end streams |

---

## `readFQ`

```nim
iterator readFQ*(path: string): FQRecord
```

Yields `FQRecord` objects with Nim strings. Records are safe to store after the loop.

```nim
import readfx

for record in readFQ("sample.fastq.gz"):
  echo record.name, " (", record.sequence.len, " bp)"
```

**When to use**: General-purpose parsing where convenience matters more than raw throughput.

---

## `readFQPtr`

```nim
iterator readFQPtr*(path: string): FQRecordPtr
```

Yields pointer-based records. The underlying buffer is reused on every iteration — do not store pointers across iterations.

```nim
import readfx

for record in readFQPtr("sample.fastq.gz"):
  echo $record.name, " (", len($record.sequence), " bp)"
  # To keep data, copy to a string:
  # let name = $record.name
```

**When to use**: Processing very large files where memory allocation overhead matters.

Malformed `@` FASTQ records raise `ValueError`; valid FASTA records with `>`
headers still parse without qualities.

**Important**: Pointers in `FQRecordPtr` are invalidated on the next iteration. Cached lengths (`nameLen`, `commentLen`, `sequenceLen`, `qualityLen`) exclude the trailing NUL terminator. If you need to retain data, copy it explicitly:

```nim
var names: seq[string]
for record in readFQPtr("sample.fastq.gz"):
  names.add($record.name)
```

---

## `readFastx`

```nim
proc readFastx*[T](f: var Bufio[T], r: var FQRecord): bool
```

Low-level reader that processes one record at a time from a `Bufio` stream.

```nim
import readfx

var record: FQRecord
var f = xopen[GzFile]("sample.fastq.gz")
defer: f.close()
while f.readFastx(record):
  echo record.name, " (", record.sequence.len, " bp)"
```

**When to use**: Custom parsing workflows, interleaving reads with other I/O, or when you need fine-grained control over the parse loop.

On parse failure, `readFastx` returns `false` and stores a negative status in
the record. Status `-4` means malformed FASTQ, including an `@` record without
a `+` line or a sequence/quality length mismatch.

---

## `readFQPair`

```nim
iterator readFQPair*(path1, path2: string, checkNames: bool = false): FQPair
```

Reads two FASTQ files in lockstep, yielding an `FQPair` with `read1` and `read2` for each pair.

```nim
import readfx

for pair in readFQPair("sample_R1.fastq.gz", "sample_R2.fastq.gz"):
  echo "R1: ", pair.read1.name
  echo "R2: ", pair.read2.name
```

- If one file runs out before the other, an `IOError` is raised.
- With `checkNames = true`, the iterator strips common suffixes (`/1`, `/2`, ` 1`, ` 2`) and raises `ValueError` if names don't match.
- Stdin (`"-"`) is supported for `path1` but not for both files simultaneously.

**When to use**: Any paired-end sequencing pipeline (Illumina R1/R2 files).

---

## `readFQInterleavedPairPtr`

```nim
iterator readFQInterleavedPairPtr*(path: string, checkNames: bool = false): FQPairPtr
```

Reads one interleaved FASTQ stream and yields `FQPairPtr` values. `read1` is
copied into scratch storage before `read2` is read, so both mates remain valid
until the next iterator advance.

```nim
import readfx

for pair in readFQInterleavedPairPtr("sample.interleaved.fastq.gz", checkNames = true):
  echo pair.read1.sequenceLen + pair.read2.sequenceLen
```

- FASTQ only.
- Raises `IOError` if the stream ends with an incomplete pair.
- With `checkNames = true`, the iterator reuses the same mate-name normalization
  as `readFQPairPtr`.
- Stdin (`"-"`) is supported.

**When to use**: Single-stream paired-end workflows where allocations and repeated
`cstring` scans are measurable overhead.

---

## `readFQInterleavedPair`

```nim
iterator readFQInterleavedPair*(path: string, checkNames: bool = false): FQPair
```

String-based counterpart to `readFQInterleavedPairPtr`: reads one interleaved
FASTQ stream and yields `FQPair` records with copied data, safe to store after
the loop.

```nim
import readfx

for pair in readFQInterleavedPair("sample.interleaved.fastq.gz", checkNames = true):
  echo pair.read1.name, " / ", pair.read2.name
```

- FASTQ only; same validation and error behavior as the pointer variant.
- Stdin (`"-"`) is supported.

**When to use**: Interleaved paired-end workflows where convenience matters
more than raw throughput — the safe default for interleaved data.

---

## Implementation Notes

- `readFQ` is built on top of `readFQPtr` and converts pointers to strings on each yield (using cached field lengths, no `strlen` rescans).
- `readFQPair` and `readFQInterleavedPair` are built the same way on `readFQPairPtr` and `readFQInterleavedPairPtr`.
- `readFQPtr`, `readFQPairPtr`, and `readFQInterleavedPairPtr` use Heng Li's `kseq.h` C library directly via FFI.
- `readFastx` is a native Nim implementation in `readfx/nimklib.nim`.
- All methods support transparent gzip decompression. The interleaved readers
  are FASTQ-only; the other readers auto-detect FASTA vs FASTQ.
