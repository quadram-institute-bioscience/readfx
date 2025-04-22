# FQRecord Utility Functions

The ReadFX library includes a utility module for manipulating FASTA/FASTQ records. This document covers how to use these utilities in your projects.

## Importing the Utilities

The utility functions are automatically exported when you import the main ReadFX module:

```nim
import readfx
# Now you have access to all fqrecordUtils functions
```

## Available Functions

### DNA Sequence Operations

- **`reverseComplement(sequence: string): string`**  
  Returns the reverse complement of a DNA sequence.

  ```nim
  let rc = reverseComplement("ATGC")  # Returns "GCAT"
  ```

- **`reverseComplementRecord(record: var FQRecord)`**  
  Reverse complements a sequence record in place.

  ```nim
  var record = readFQ("example.fastq").next
  reverseComplementRecord(record)
  ```

- **`reverseComplementRecord(record: FQRecord): FQRecord`**  
  Creates a new record with reverse-complemented sequence.

  ```nim
  let record = readFQ("example.fastq").next
  let rcRecord = reverseComplementRecord(record)
  ```

- **`gcContent(sequence: string): float`**  
  Calculates GC content as a fraction between 0.0 and 1.0.

  ```nim
  let gc = gcContent("ATGC")  # Returns 0.5
  ```

### Quality-Based Operations

- **`trimQuality(quality: string, minQual: int, offset: int = 33): string`**  
  Trims a quality string based on minimum quality threshold.

  ```nim
  let trimmedQuality = trimQuality(record.quality, 20)
  ```

- **`qualityTrim(record: var FQRecord, minQual: int, offset: int = 33)`**  
  Trims a record based on quality scores, modifying both sequence and quality.

  ```nim
  var record = readFQ("example.fastq").next
  qualityTrim(record, 20)  # Trim bases with quality < 20
  ```

- **`maskLowQuality(record: var FQRecord, minQual: int, offset: int = 33, maskChar: char = 'N')`**  
  Masks sequence positions with low quality scores with a specified character.

  ```nim
  var record = readFQ("example.fastq").next
  maskLowQuality(record, 20)  # Mask bases with quality < 20 as 'N'
  ```

### Record Manipulation

- **`subSequence(record: FQRecord, start: int, length: int = -1): FQRecord`**  
  Extracts a subsequence from a record, returning a new record.

  ```nim
  let record = readFQ("example.fastq").next
  # Extract first 10 bases
  let firstTenBases = subSequence(record, 0, 10)
  # Extract from position 10 to the end
  let fromTenToEnd = subSequence(record, 10)
  ```

## Complete Example

Here's a complete example showing how to use these utilities:

```nim
import readfx
import strutils

# Process a FASTQ file
let inputFile = "example.fastq"
for record in readFQ(inputFile):
  echo "Processing record: ", record.name
  
  # Get GC content
  let gc = gcContent(record.sequence) * 100
  echo "GC content: ", gc.formatFloat(ffDecimal, 2), "%"
  
  # Create a reverse complemented version
  let rcRecord = reverseComplementRecord(record)
  echo "Original: ", record.sequence
  echo "Reverse:  ", rcRecord.sequence
  
  # Make a modified copy we can change
  var modifiedRecord = record
  
  # Mask low quality bases (Q < 20)
  maskLowQuality(modifiedRecord, 20)
  
  # Trim low quality ends
  qualityTrim(modifiedRecord, 20)
  
  # Extract a subsequence
  let subseq = subSequence(modifiedRecord, 0, 50)
  echo "First 50 bases: ", subseq.sequence
```

For more examples, see the provided `fqrecordUtils_example.nim` file.

## Implementation Details

These utilities are defined in the `readfx/fqrecordUtils.nim` file, which is automatically imported and exported by the main `readfx.nim` module.

The utilities are compatible with both the C wrapper and native Nim implementation of ReadFX, as they operate on the unified `FQRecord` type.