# nimreadfq

A Nim wrapper for [Heng Li's kseq/readfq](https://github.com/lh3/readfq/), an efficient and fast parser for FastQ and Fasta files.
nimreadfq supports reading of FastQ and Fasta files from stdin (use "-"), gzipped or flat files and is **very** fast (see benchmark below).

The main function is `readFQ()`, an iterator that yields `FQRecord(s)`. An alternative is `readFQPtr()`, which returns `FQRecordPtr(s)`. The difference is that the latter uses `ptr char` instead of strings and is thus potentially faster but memory is reused during iterations.

See `example.nim` and `tests/tester.nim` for code examples.

The initial Nim integration (and hard work) was done by [Haibao Tang](https://github.com/tanghaibao) as part of his [bio-pipeline
repo](https://github.com/tanghaibao/bio-pipeline/). Haibao generously [granted full rights to his code base](https://github.com/tanghaibao/bio-pipeline/issues/4), after which I started this separate package called [nimreadfq](https://github.com/andreas-wilm/nimreadfq) for integration into nimble.

This repository was started before Heng Li wrote his article ["Fast high-level programming languages"](https://lh3.github.io/2020/05/17/fast-high-level-programming-languages), which contains a native Nim implementation (not included here).

## Benchmark

nimreadfq is significantly faster than packages with similar functionality. Below are example timings for reading 5,682,010 sequences from `M_abscessus_HiSeq.fq` ([source](https://github.com/lh3/biofast/releases/tag/biofast-data-v1)) run on my MacBook Pro 2019:

fastq:
- readfq: 10.7s
- bioseq: 49.6s
- fastx: 47.6s

fastq.gz:
- readfq: 20.4s
- bioseq: 163.3s

How to reproduce results:

    cd ./benchmark
    nimble build
    ./benchmark
