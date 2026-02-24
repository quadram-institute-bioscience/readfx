# ReadFX Utility Functions Reference

All functions are exported from `readfx` via `readfx/sequtils.nim`.

---

## Trimming

### `trimStart`
```nim
proc trimStart*(record: FQRecord, bases: int): FQRecord
```
Removes `bases` characters from the 5' (start) end. Returns a new record.

### `trimEnd`
```nim
proc trimEnd*(record: FQRecord, bases: int): FQRecord
```
Removes `bases` characters from the 3' (end) end. Returns a new record.

---

## Quality Conversion

### `qualCharToInt`
```nim
proc qualCharToInt*(c: char, offset: int = 33): int
```
Converts a quality character to its Phred integer value. Default offset 33 = Sanger/Illumina 1.8+.

### `qualIntToChar`
```nim
proc qualIntToChar*(q: int, offset: int = 33): char
```
Converts a Phred integer value to its quality character.

---

## Quality Statistics

### `avgQuality` (record overload)
```nim
proc avgQuality*(record: FQRecord, offset: int = 33): float
```
Returns the mean Phred quality score for a record.

### `avgQuality` (string overload)
```nim
proc avgQuality*(quality: string, offset: int = 33): float
```
Returns the mean Phred quality score for a raw quality string.

---

## Quality Trimming

### `trimQuality`
```nim
proc trimQuality*(quality: string, minQual: int, offset: int = 33): string
```
Trims trailing bases below `minQual` from a quality string.

### `qualityTrim`
```nim
proc qualityTrim*(record: var FQRecord, minQual: int, offset: int = 33)
```
Trims both sequence and quality in place.

### `maskLowQuality`
```nim
proc maskLowQuality*(record: var FQRecord, minQual: int, offset: int = 33, maskChar: char = 'N')
```
Replaces bases whose quality is below `minQual` with `maskChar`.

---

## Reverse Complement

### `revCompl` (string)
```nim
proc revCompl*(sequence: string): string
```
Returns the reverse complement of a DNA sequence string.

### `revCompl` (in-place)
```nim
proc revCompl*(record: var FQRecord)
```
Reverse complements a record in place. Also reverses the quality string.

### `revCompl` (copy)
```nim
proc revCompl*(record: FQRecord): FQRecord
```
Returns a new record with a reverse-complemented sequence.

---

## Subsequence

### `subSequence`
```nim
proc subSequence*(record: FQRecord, start: int, length: int = -1): FQRecord
```
Extracts a subsequence from position `start` (0-based) for `length` bases. `length = -1` means to the end of the sequence.

---

## Composition

### `composition`
```nim
proc composition*(record: FQRecord): SeqComp
```
Returns a `SeqComp` with per-base counts (A, C, G, T, N, Other) and GC fraction.

### `gcContent` (string)
```nim
proc gcContent*(sequence: string): float
```
Returns GC fraction (0.0–1.0) for a sequence string.

### `gcContent` (record)
```nim
proc gcContent*(record: FQRecord): float
```
Returns GC fraction for a record's sequence.
