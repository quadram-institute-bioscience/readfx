# ReadFX Utility Functions

This file documents the utility functions available in the ReadFX library.

## trimStart*

```nim
proc trimStart*(record: FQRecord, bases: int): FQRecord =
```

Remove N bases from the start (5' end) of a sequence record

Args:
  record: Input FQRecord
  bases: Number of bases to remove from the start

Returns:
  New FQRecord with trimmed sequence and quality

## trimEnd*

```nim
proc trimEnd*(record: FQRecord, bases: int): FQRecord =
```

Remove N bases from the end (3' end) of a sequence record

Args:
  record: Input FQRecord
  bases: Number of bases to remove from the end

Returns:
  New FQRecord with trimmed sequence and quality

## qualCharToInt*

```nim
proc qualCharToInt*(c: char, offset: int = 33): int =
```

Convert a quality character to its integer value

Args:
  c: Quality character
  offset: Quality score offset (default is 33 for Sanger/Illumina 1.8+)

Returns:
  Integer value of the quality character

## qualIntToChar*

```nim
proc qualIntToChar*(q: int, offset: int = 33): char =
```

Convert an integer quality value to its character representation

Args:
  q: Integer quality value
  offset: Quality score offset (default is 33 for Sanger/Illumina 1.8+)

Returns:
  Character representation of the quality value

## avgQuality*

```nim
proc avgQuality*(record: FQRecord, offset: int = 33): float =
```

Calculate the average quality of a sequence record

Args:
  record: FQRecord to analyze
  offset: Quality score offset (default is 33 for Sanger/Illumina 1.8+)

Returns:
  Average quality score as a float

## avgQuality*

```nim
proc avgQuality*(quality: string, offset: int = 33): float =
```

Calculate the average quality of a quality string

Args:
  quality: Quality string
  offset: Quality score offset (default is 33 for Sanger/Illumina 1.8+)

Returns:
  Average quality score as a float

## rc_string

```nim
proc rc_string(sequence: string): string =
```

Reverse complement a DNA sequence

Example:
  let rc = reverseComplement("ATGC")  # returns "GCAT"

## trimQuality*

```nim
proc trimQuality*(quality: string, minQual: int, offset: int = 33): string =
```

Trim a quality string based on minimum quality threshold

Args:
  quality: Quality string
  minQual: Minimum quality value (0-40)
  offset: Quality score offset (33 for Sanger/Illumina 1.8+)

Returns:
  Trimmed quality string

## qualityTrim*

```nim
proc qualityTrim*(record: var FQRecord, minQual: int, offset: int = 33) =
```

Trim a record based on quality scores

Args:
  record: FQRecord to modify
  minQual: Minimum quality value
  offset: Quality score offset (33 for Sanger/Illumina 1.8+)

## revCompl*

```nim
proc revCompl*(sequence: string): string =
```

Reverse complement a DNA sequence

Args:
  sequence: DNA sequence

Returns:
  Reverse-complemented sequence

## revCompl*

```nim
proc revCompl*(record: var FQRecord) =
```

Reverse complement a sequence record in place

Args:
  record: FQRecord to modify

## revCompl*

```nim
proc revCompl*(record: FQRecord): FQRecord =
```

Create a new record with reverse-complemented sequence

Args:
  record: Input FQRecord

Returns:
  New FQRecord with reverse-complemented sequence

## subSequence*

```nim
proc subSequence*(record: FQRecord, start: int, length: int = -1): FQRecord =
```

Extract a subsequence from a record

Args:
  record: Input FQRecord
  start: Start position (0-based)
  length: Length of subsequence to extract (-1 for end of sequence)

Returns:
  New FQRecord with extracted subsequence

## composition*

```nim
proc composition*(record: FQRecord): SeqComp =
```

Calculate composition of a DNA Sequence
Returns a SeqComp object with counts of A, C, G, T, N, and Other (int) and GC content (float)

## gcContent*

```nim
proc gcContent*(sequence: string): float =
```

Calculate GC content of a DNA sequence

Args:
  sequence: DNA sequence

Returns:
  GC content as a fraction between 0.0 and 1.0

## gcContent*

```nim
proc gcContent*(record: FQRecord): float =
```

Calculate GC content of a FQRecord (using its sequence)

## maskLowQuality*

```nim
proc maskLowQuality*(record: var FQRecord, minQual: int, offset: int = 33, maskChar: char = 'N') =
```

Mask sequence positions with low quality scores

Args:
  record: FQRecord to modify
  minQual: Minimum quality value
  offset: Quality score offset
  maskChar: Character to use for masking

