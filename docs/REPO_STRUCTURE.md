# Repository Structure

```text
readfx/
├── readfx.nim          # Main library — imports and re-exports all submodules
├── readfx.nimble       # Nimble package file
│
├── readfx/             # Submodules
│   ├── seqtypes.nim    # Type definitions: FQRecord, FQRecordPtr, FQPair, SeqComp, Strand
│   ├── sequtils.nim    # Sequence utilities: revCompl, gcContent, qualityTrim, …
│   ├── oligoutils.nim  # IUPAC primer matching: matchIUPAC
│   ├── nimklib.nim     # Native Nim FASTX parser (used by readFastx)
│   ├── kseq.h          # Heng Li's kseq C library (used by readFQPtr, readFQ, readFQPair)
│   └── klib/           # klib headers
│
├── tests/              # Unit tests (nimble test)
│   ├── tester.nim
│   ├── illumina_1.fq.gz  # Example paired-end R1 file
│   └── illumina_2.fq.gz  # Example paired-end R2 file
│
├── benchmark/          # Performance benchmarks
│   └── benchmark.nim
│
├── docs/               # GitHub Pages documentation (this directory)
│
└── src/demos/          # Demo programs
```

## Module Dependencies

```
seqtypes.nim            (no dependencies)
    ↑
sequtils.nim            (imports seqtypes)
oligoutils.nim          (imports stdlib only)
nimklib.nim             (imports seqtypes, zlib)
    ↑
readfx.nim              (imports all of the above + kseq.h via FFI)
```

## Key Files

| File | Purpose |
|------|---------|
| `readfx.nim` | Entry point — `import readfx` gives you everything |
| `readfx/seqtypes.nim` | All type definitions |
| `readfx/sequtils.nim` | Sequence manipulation utilities |
| `readfx/oligoutils.nim` | IUPAC primer/barcode matching |
| `readfx/nimklib.nim` | Native Nim buffered parser (`readFastx`) |
| `readfx/kseq.h` | C parser used by `readFQ`, `readFQPtr`, `readFQPair` |
