# ReadFX Data Structures

## `FQRecordPtr`

Pointer-based record for high-performance streaming via the `kseq` C library. Pointers are reused on each iteration and must not be stored without copying.

```nim
type FQRecordPtr* = object
  name*: ptr char      # Sequence identifier (null-terminated)
  comment*: ptr char   # Optional description (null-terminated)
  sequence*: ptr char  # Nucleotide sequence (null-terminated)
  quality*: ptr char   # Quality scores (null-terminated; nil for FASTA)
```

Use with `readFQPtr`. To retain data beyond the current iteration, convert to string: `$record.name`.

---

## `FQRecord`

String-based record. Safe to store and manipulate after iteration. Used by `readFQ`, `readFastx`, and as the element type inside `FQPair`.

```nim
type FQRecord* = object
  name*: string         # Sequence identifier
  comment*: string      # Optional description
  sequence*: string     # Nucleotide sequence
  quality*: string      # Quality scores (empty for FASTA)
  status*, lastChar*: int  # Internal parsing state
```

`status` codes (relevant when using `readFastx`):
- `> 0` — sequence length (normal record)
- `-1` — end of file
- `-2` — stream error
- `-3` — other parsing error
- `-4` — sequence and quality length mismatch

---

## `FQPair`

Paired-end record containing two `FQRecord` objects. Yielded by `readFQPair`.

```nim
type FQPair* = object
  read1*: FQRecord   # Forward read (R1)
  read2*: FQRecord   # Reverse read (R2)
```

```nim
for pair in readFQPair("R1.fastq.gz", "R2.fastq.gz"):
  echo pair.read1.name, " / ", pair.read2.name
```

---

## `SeqComp`

Nucleotide composition statistics, returned by `composition()`.

```nim
type SeqComp* = object
  A*, C*, G*, T*: int   # Per-base counts
  N*: int               # Ambiguous bases
  Other*: int           # Non-ACGTN characters
  GC*: float            # GC fraction (0.0–1.0)
```

```nim
let comp = composition(record)
echo "GC=", comp.GC, " Ns=", comp.N
```

---

## `Bufio[T]`

Generic buffered reader. Typically instantiated as `Bufio[GzFile]` via `xopen[GzFile](path)`.

```nim
type Bufio*[T] = tuple[fp: T, buf: string, st, en, sz: int, EOF: bool]
```

---

## `Interval[S,T]`

Genomic interval for use with the built-in interval tree.

```nim
type Interval*[S,T] = tuple[st, en: S, data: T, max: S]
```

Use `index()` after building a `seq[Interval]`, then query with `overlap()`.

---

## Choosing Between `FQRecord` and `FQRecordPtr`

| | `FQRecord` | `FQRecordPtr` |
|---|---|---|
| Memory | Allocates strings per record | Reuses a single buffer |
| Safety | Safe to store/manipulate | Valid only during current iteration |
| Iterator | `readFQ` | `readFQPtr` |
| Performance | Good | Excellent |

For most use cases `readFQ` / `FQRecord` is the right choice. Reach for `readFQPtr` when processing tens of millions of records and allocation overhead is measurable.

---

## Performance Notes

Benchmarks (`benchmark/`) show that with `--opt:speed --gc:arc`, the performance difference between object and tuple layouts becomes negligible. Recommended compile flags for production:

```bash
nim c --opt:speed --gc:arc myprogram.nim
```
