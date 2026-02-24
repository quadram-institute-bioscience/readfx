# FQRecord Utility Functions

Utility functions for manipulating FASTA/FASTQ records are exported automatically when you `import readfx`.

## DNA Sequence Operations

### `revCompl`

```nim
proc revCompl*(sequence: string): string
proc revCompl*(record: var FQRecord)        # in-place; also reverses quality
proc revCompl*(record: FQRecord): FQRecord  # returns new record
```

Reverse complements a DNA sequence.

```nim
let rc = revCompl("ATGC")  # "GCAT"

var rec: FQRecord
# ... populate rec ...
revCompl(rec)               # modify in place
let rcRec = revCompl(rec)   # get a new copy
```

### `gcContent`

```nim
proc gcContent*(sequence: string): float
proc gcContent*(record: FQRecord): float
```

Returns the GC fraction (0.0–1.0).

```nim
let gc = gcContent("ATGC")  # 0.5
```

### `composition`

```nim
proc composition*(record: FQRecord): SeqComp
```

Returns counts of A, C, G, T, N, and Other bases, plus GC fraction.

```nim
let comp = composition(record)
echo "GC=", comp.GC, " N=", comp.N
```

---

## Quality Operations

### `qualCharToInt` / `qualIntToChar`

```nim
proc qualCharToInt*(c: char, offset: int = 33): int
proc qualIntToChar*(q: int, offset: int = 33): char
```

Convert between quality characters and Phred integers. Default offset is 33 (Sanger/Illumina 1.8+).

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

Trims trailing bases below `minQual` from a quality string.

### `qualityTrim`

```nim
proc qualityTrim*(record: var FQRecord, minQual: int, offset: int = 33)
```

Trims both sequence and quality strings in place.

```nim
qualityTrim(record, 20)  # trim bases with Phred < 20
```

### `maskLowQuality`

```nim
proc maskLowQuality*(record: var FQRecord, minQual: int, offset: int = 33, maskChar: char = 'N')
```

Replaces low-quality bases with `maskChar` without trimming.

```nim
maskLowQuality(record, 20)       # replace Q<20 with 'N'
maskLowQuality(record, 20, maskChar = 'X')
```

---

## Record Manipulation

### `subSequence`

```nim
proc subSequence*(record: FQRecord, start: int, length: int = -1): FQRecord
```

Returns a new record containing a slice of the sequence (and quality).

```nim
let first50 = subSequence(record, 0, 50)
let fromPos10 = subSequence(record, 10)   # to end
```

### `trimStart` / `trimEnd`

```nim
proc trimStart*(record: FQRecord, bases: int): FQRecord
proc trimEnd*(record: FQRecord, bases: int): FQRecord
```

Remove bases from the 5' or 3' end. Return new records.

```nim
let trimmed = trimStart(record, 5)   # remove first 5 bases
```

---

## Complete Example

```nim
import readfx
import strutils

for record in readFQ("example.fastq.gz"):
  echo "Name: ", record.name

  let comp = composition(record)
  echo "GC=", (comp.GC * 100).formatFloat(ffDecimal, 1), "%  N=", comp.N

  let avg = avgQuality(record)
  echo "Mean quality: ", avg.formatFloat(ffDecimal, 1)

  # Create a modified copy
  var modified = record
  qualityTrim(modified, 20)
  maskLowQuality(modified, 15)

  let first50 = subSequence(modified, 0, 50)
  echo "First 50 bp: ", first50.sequence

  echo revCompl(record.sequence)
```

---

## Implementation

These utilities are defined in `readfx/sequtils.nim` and `readfx/oligoutils.nim`, and are re-exported by the main `readfx.nim` module.
