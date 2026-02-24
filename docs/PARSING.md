# FASTX Parsing in ReadFX

ReadFX provides four methods for parsing FASTA and FASTQ files:

1. **`readFQ`** — String-based iterator (convenient, safe)
2. **`readFQPtr`** — Pointer-based iterator (fastest, requires care)
3. **`readFastx`** — Low-level buffered reader (most flexible)
4. **`readFQPair`** — Synchronized paired-end iterator

## Comparison at a Glance

| Method        | Memory   | Speed     | Ease of Use    | Use case                   |
|---------------|----------|-----------|----------------|----------------------------|
| `readFQ`      | Higher   | Good      | Excellent      | General use                |
| `readFQPtr`   | Low      | Excellent | Moderate       | High-throughput streaming  |
| `readFastx`   | Custom   | Excellent | Requires setup | Custom I/O workflows       |
| `readFQPair`  | Moderate | Good      | Excellent      | Paired-end reads           |

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

**Important**: Pointers in `FQRecordPtr` are invalidated on the next iteration. If you need to retain data, copy it explicitly:

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

## Implementation Notes

- `readFQ` is built on top of `readFQPtr` and converts pointers to strings on each yield.
- `readFQPtr` and `readFQPair` use Heng Li's `kseq.h` C library directly via FFI.
- `readFastx` is a native Nim implementation in `readfx/nimklib.nim`.
- All methods support both FASTA and FASTQ formats (auto-detected) and transparent gzip decompression.
