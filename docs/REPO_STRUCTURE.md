# Utility structure

Suggested file structure for reorganizing the codebase:

```text
  readfx/
  ├── readfx.nim             # Main library file (exports all modules)
  ├── readfx.nimble          # Package file
  ├── src/
  │   ├── types.nim          # Type definitions
  │   ├── io.nim             # I/O operations (GzFile, Bufio)
  │   ├── parsers.nim        # Parsing functions (readFastx, readFQPtr, readFQ)
  │   ├── utils.nim          # Utility functions for FQRecord
  │   └── intervals.nim      # Interval operations
  ├── readfx/                # C bindings
  │   ├── kseq.h
  │   ├── klib/
  │   │   ├── README.md
  │   │   └── kseq.h
  │   └── bindings.nim       # C bindings for kseq.h
  └...
```

##  Import hierarchy:

  1. types.nim:
    - Imports: None (base module)
    - Contains: All type definitions (FQRecord, FQRecordPtr, etc.)
  2. io.nim:
    - Imports: types.nim
    - Contains: GzFile operations, Bufio implementation
  3. bindings.nim:
    - Imports: types.nim
    - Contains: C bindings to klib/kseq.h, kseq_init, kseq_read, etc.
  4. parsers.nim:
    - Imports: types.nim, io.nim, bindings.nim
    - Contains: readFastx, readFQPtr, readFQ implementations
  5. utils.nim:
    - Imports: types.nim
    - Contains: reverseComplement, gcContent, qualityTrim, etc.
  6. intervals.nim:
    - Imports: types.nim
    - Contains: Interval types and operations
  7. readfx.nim:
    - Imports and re-exports: types.nim, io.nim, parsers.nim, utils.nim, intervals.nim
    - Minimal glue code

  This structure creates a clear dependency direction with no circular references:

```text
  types.nim <── io.nim
      ^         ^
      |         |
      |         |
  bindings.nim  |
      ^         |
      |         |
      └── parsers.nim
           ^
           |
  utils.nim├── intervals.nim
      ^    ^
      |    |
      └────┴── readfx.nim (main module)
```