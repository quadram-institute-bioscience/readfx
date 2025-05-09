# ReadFX Data Structures

## FQRecord and FQRecordPtr Types

ReadFX provides two primary data structures for working with FASTA/FASTQ sequence records:

### FQRecord

`FQRecord` is a Nim object that stores sequence data with full string copies. It provides a convenient and safe way to work with sequence records in memory.

```nim
type
  FQRecord* = object
    name*: string        # Sequence identifier (required)
    comment*: string     # Optional sequence description/comment
    sequence*: string    # The actual nucleotide sequence
    quality*: string     # Quality scores (optional, present in FASTQ files)
    status*, lastChar*: int  # Internal state for parsing
```

- **Usage**: Ideal for most processing tasks where data needs to be manipulated or stored.
- **Benefits**: Safe, easy to use, fully manages its own memory.
- **Memory**: Allocates new strings for all fields.

### FQRecordPtr

`FQRecordPtr` is a lightweight structure that uses pointers to character data. It's designed for high-performance streaming applications where memory efficiency is critical.

```nim
type
  FQRecordPtr* = object
    name*: ptr char      # Pointer to sequence identifier
    comment*: ptr char   # Pointer to optional description/comment
    sequence*: ptr char  # Pointer to nucleotide sequence
    quality*: ptr char   # Pointer to quality scores (optional)
```

- **Usage**: Best for high-throughput parsing and streaming applications.
- **Benefits**: Very memory efficient, no allocation overhead during parsing.
- **Caution**: Pointers are only valid during iteration. The data is owned by the parser and will be overwritten or freed when reading the next record.

### Choosing Between Types

- Use `FQRecord` when you need to store or manipulate records.
- Use `FQRecordPtr` when parsing millions of records for streaming applications where performance is critical.

### Example

```nim
# Reading with FQRecord (safe, copies data)
for record in readFQ("sample.fastq"):
  echo record.name
  echo record.sequence.len

# Reading with FQRecordPtr (fast, uses pointers)
for record in readFQPtr("sample.fastq"):
  echo $record.name
  echo $record.sequence
```

## Performance Benchmarks and Optimization Results

> See benchmark/data_structure.nim

## Key Findings

Based on the benchmark results (see below) across different compiler optimization settings and garbage collectors:

* Creation Performance

Nearly identical creation time between tuples and objects across all configurations
Object creation was marginally faster with default GC and speed optimizations (1.02x)
With ARC, tuple creation was slightly faster (1.00x - effectively identical)

* Field Access

With default GC: Tuples were slightly faster (1.03x) without optimizations
With optimizations: Objects were faster in some tests (1.02x with default GC, 1.28x with ARC)
With ORC: Tuples were slightly faster (1.08x)

* Assignment Operations: Most significant difference: Objects are dramatically faster for assignments with default GC:

  *  4.05x faster without optimizations
  *  6.77x faster with optimizations


With ARC/ORC: Performance becomes virtually identical (difference is 1.00x-1.01x)

### Memory Usage

Nearly identical memory usage across all configurations (ratio consistently 1.000)
ARC/ORC reduced overall memory usage for both types compared to default GC

Impact of Compiler Settings
Speed Optimization (--opt:speed)

Reduced creation time by ~60% for both types
Improved access times by ~70-90%
Amplified the assignment performance gap with default GC

### Garbage Collector Options

ARC/ORC: Eliminated the assignment performance difference between tuples and objects
ARC/ORC: Reduced memory usage by ~15% compared to default GC
ARC/ORC: Dramatically improved field access times (up to 10x faster)

### Summary 
For optimal performance in bioinformatics applications:

Compile with `--opt:speed` `--gc:arc` or `--opt:speed` `--gc:orc`
With these settings, the choice between tuple or object becomes less critical

## Results 

### With 1_000_000 iterations [no compile optimizations]

```text
Running benchmarks with 1000000 iterations
========================================
Benchmarking record creation...
  Tuple creation time: 52.854600 seconds
  Object creation time: 52.623764 seconds
  Object creation is 1.00x faster
----------------------------------------
Benchmarking field access...
  Tuple access time: 0.069178 seconds
  Object access time: 0.071343 seconds
  Tuple access is 1.03x faster
----------------------------------------
Benchmarking assignment operations...
  Tuple assignment time: 1.927895 seconds
  Object assignment time: 0.476055 seconds
  Object assignment is 4.05x faster
----------------------------------------
Memory comparison through GC stats...
  Tuple memory: 43280944 bytes
  Object memory: 43264864 bytes
  Memory ratio: 1.000
========================================
Benchmarking completed.
```

### With 1_000_000 iterations and `--opt:speed`

```text
Running benchmarks with 1000000 iterations
========================================
Benchmarking record creation...
  Tuple creation time: 20.316888 seconds
  Object creation time: 19.914439 seconds
  Object creation is 1.02x faster
----------------------------------------
Benchmarking field access...
  Tuple access time: 0.019489 seconds
  Object access time: 0.019074 seconds
  Object access is 1.02x faster
----------------------------------------
Benchmarking assignment operations...
  Tuple assignment time: 0.908726 seconds
  Object assignment time: 0.134202 seconds
  Object assignment is 6.77x faster
----------------------------------------
Memory comparison through GC stats...
  Tuple memory: 43269824 bytes
  Object memory: 43253408 bytes
  Memory ratio: 1.000
========================================
Benchmarking completed.
```

### With 1_000_000 iterations and `--opt:speed` and `--gc:arc`

```text
Running benchmarks with 1000000 iterations
========================================
Benchmarking record creation...
  Tuple creation time: 21.704738 seconds
  Object creation time: 21.787406 seconds
  Tuple creation is 1.00x faster
----------------------------------------
Benchmarking field access...
  Tuple access time: 0.002041 seconds
  Object access time: 0.001592 seconds
  Object access is 1.28x faster
----------------------------------------
Benchmarking assignment operations...
  Tuple assignment time: 0.081974 seconds
  Object assignment time: 0.081971 seconds
  Object assignment is 1.00x faster
----------------------------------------
Memory comparison through GC stats...
  Tuple memory: 36924864 bytes
  Object memory: 36924864 bytes
  Memory ratio: 1.000
========================================
Benchmarking completed.
```

### With 1_000_000 iterations and `--opt:speed` and `--gc:orc`

```text
Running benchmarks with 1000000 iterations
========================================
Benchmarking record creation...
  Tuple creation time: 21.475449 seconds
  Object creation time: 21.414222 seconds
  Object creation is 1.00x faster
----------------------------------------
Benchmarking field access...
  Tuple access time: 0.002275 seconds
  Object access time: 0.002452 seconds
  Tuple access is 1.08x faster
----------------------------------------
Benchmarking assignment operations...
  Tuple assignment time: 0.081635 seconds
  Object assignment time: 0.082110 seconds
  Tuple assignment is 1.01x faster
----------------------------------------
Memory comparison through GC stats...
  Tuple memory: 36924864 bytes
  Object memory: 36924864 bytes
  Memory ratio: 1.000
========================================
Benchmarking completed.
```