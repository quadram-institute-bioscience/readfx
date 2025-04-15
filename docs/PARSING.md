# FASTX Parsing Methods in ReadFX

ReadFX provides three primary methods for parsing FASTA and FASTQ files, each with different characteristics and use cases:

1. **readFQ** - String-based high-level iterator
2. **readFQPtr** - Pointer-based high-performance iterator
3. **readFastx** - Lower-level reader for custom workflow integration

## Comparison at a Glance

| Method      | Memory Usage | Speed     | Ease of Use    | Flexibility |
|-------------|--------------|-----------|----------------|-------------|
| `readFQ`    | Higher       | Good      | Excellent      | Good        |
| `readFQPtr` | Low          | Excellent | Moderate       | Good        |
| `readFastx` | Customizable | Excellent | Requires setup | Excellent   |

## readFQ

```nim
iterator readFQ*(path: string): FQRecord
```

`readFQ` is a high-level iterator that returns `FQRecord` objects with string fields.

### How to Use

```nim
import readfx

for record in readFQ("sample.fastq.gz"):
  echo record.name, " has length ", record.sequence.len
  # Manipulate record.sequence, record.quality, etc. as strings
```

### When to Use

- When you want clean, easy-to-use code
- When you need to manipulate sequence or quality data
- When you need to store records for later use
- When working with other Nim code expecting strings

### Why Use

- Safer and more idiomatic Nim code
- No memory management concerns
- Easy string manipulation
- Records persist after iteration

## readFQPtr

```nim
iterator readFQPtr*(path: string): FQRecordPtr
```

`readFQPtr` is a high-performance iterator that returns `FQRecordPtr` objects with pointer fields for maximum efficiency.

### How to Use

```nim
import readfx

for record in readFQPtr("sample.fastq.gz"):
  echo $record.name, " has length ", $record.sequence.len
  # Be careful! These pointers are reused on each iteration
```

### When to Use

- When processing very large files (millions of records)
- When memory usage is a concern
- When maximum performance is required
- For read-only operations where you don't need to keep records

### Why Use

- Significantly lower memory usage
- Often 1.5-2x faster than `readFQ`
- No string allocation for each record
- Avoids garbage collection overhead

### Important Note

The pointers in `FQRecordPtr` are reused with each iteration! If you need to keep a record after moving to the next iteration, you must copy the data:

```nim
var savedNames: seq[string]
for record in readFQPtr("sample.fastq.gz"):
  savedNames.add($record.name)  # Make a copy
```

## readFastx

```nim
proc readFastx*[T](f: var Bufio[T], r: var FQRecord): bool
```

`readFastx` is a lower-level procedure that reads one record at a time from a buffered input stream.

### How to Use

```nim
import readfx

var record: FQRecord
var f = xopen[GzFile]("sample.fastq.gz")
defer: f.close()
while f.readFastx(record):
  echo record.name, " has length ", record.sequence.len
```

### When to Use

- When you need more control over the parsing process
- To integrate with custom I/O workflows
- When you want to manage your own file handles
- For advanced use cases requiring custom buffering

### Why Use

- Finest control over parsing behavior
- Integrates with custom file handling
- Allows for interleaving with other operations
- Can be more efficient for specialized workflows

## Implementation Details

- `readFQ` is actually built on top of `readFQPtr`, converting pointers to strings
- `readFastx` is the native Nim implementation used internally
- Both implementations support FASTA and FASTQ formats, gzipped files, and reading from stdin

## Performance Considerations

In benchmarks on large files:

- `readFQPtr` is typically the fastest but requires careful memory management
- `readFQ` is slightly slower due to string allocations but much safer
- `readFastx` performance depends on how you implement the surrounding code

Choose the method that best balances your needs for performance, safety, and code simplicity.