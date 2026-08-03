import unittest
import os
import random
import strutils

import ../readfx

type RecTuple = tuple[name, comment, sequence, quality: string]

proc slurp(path: string): seq[RecTuple] =
  for rec in readFQ(path):
    result.add((rec.name, rec.comment, rec.sequence, rec.quality))

proc toTuples(records: seq[FQRecord]): seq[RecTuple] =
  for r in records:
    result.add((r.name, r.comment, r.sequence, r.quality))

proc tmpPath(name: string): string =
  getTempDir() / name

# ============================================================
# Round-trip identity: read -> write -> read must be stable
# ============================================================

test "round-trip: FASTQ plain and gzip":
  let expected = slurp("./tests/fastq_demo.fq")
  check expected.len > 0
  for compression in [false, true]:
    let outPath = tmpPath(if compression: "readfx_rt.fq.gz" else: "readfx_rt.fq")
    if fileExists(outPath):
      removeFile(outPath)
    var w = fastxWriter(
      format = fxfFastq,
      compression = compression,
      destination = fileDestination(outPath),
      bufferSize = 128
    )
    for rec in readFQ("./tests/fastq_demo.fq"):
      w.writeRecord(rec)
    w.close()
    check slurp(outPath) == expected
    removeFile(outPath)

test "round-trip: FASTA plain, gzip and forced wrapping":
  let expected = slurp("./tests/fasta_demo.fa")
  check expected.len > 0
  for (compression, width) in [(false, 60), (true, 60), (false, 5)]:
    let outPath = tmpPath(if compression: "readfx_rt.fa.gz" else: "readfx_rt.fa")
    if fileExists(outPath):
      removeFile(outPath)
    var w = fastxWriter(
      format = fxfFasta,
      compression = compression,
      destination = fileDestination(outPath),
      bufferSize = 128,
      fastaWidth = width
    )
    for rec in readFQ("./tests/fasta_demo.fa"):
      w.writeRecord(rec)
    w.close()
    let got = slurp(outPath)
    check got == expected           # wrapping must not alter parsed content
    for rec in got:
      check rec.quality == ""       # FASTA round-trip carries no qualities
    removeFile(outPath)

test "round-trip: seeded random records":
  var rng = initRand(20260802)
  var records: seq[FQRecord]
  for i in 0..<200:
    let seqLen = 1 + rng.rand(300)
    var sq = newString(seqLen)
    for c in sq.mitems:
      c = "ACGTNacgtn"[rng.rand(9)]
    var q = newString(seqLen)
    for c in q.mitems:
      c = char(33 + rng.rand(40))   # full Phred+33 printable range
    let comment = if rng.rand(1.0) < 0.5: "cmt " & $i else: ""
    records.add(FQRecord(name: "rec" & $i, comment: comment,
                         sequence: sq, quality: q))

  # FASTQ, plain and gzip: full identity including qualities
  for compression in [false, true]:
    let outPath = tmpPath(if compression: "readfx_rt_rand.fq.gz" else: "readfx_rt_rand.fq")
    if fileExists(outPath):
      removeFile(outPath)
    var w = fastxWriter(
      format = fxfFastq,
      compression = compression,
      destination = fileDestination(outPath),
      bufferSize = 256
    )
    for r in records:
      w.writeRecord(r)
    w.close()
    check slurp(outPath) == records.toTuples()
    removeFile(outPath)

  # FASTA: identity on name/comment/sequence, qualities dropped
  block:
    let outPath = tmpPath("readfx_rt_rand.fa")
    if fileExists(outPath):
      removeFile(outPath)
    var w = fastxWriter(
      format = fxfFasta,
      compression = false,
      destination = fileDestination(outPath),
      bufferSize = 256,
      fastaWidth = 47                # odd width on purpose
    )
    for r in records:
      w.writeRecord(r)
    w.close()
    let got = slurp(outPath)
    check got.len == records.len
    for i in 0..<records.len:
      check got[i].name == records[i].name
      check got[i].comment == records[i].comment
      check got[i].sequence == records[i].sequence
      check got[i].quality == ""
    removeFile(outPath)

# ============================================================
# Malformed-input corpus: pinned parser behaviors
# ============================================================

test "corpus: wrapped FASTQ sequence and quality are joined":
  check slurp("./tests/corpus/valid_wrapped.fq") == @[
    ("read1", "comment here", "ACGTACGTAC", "IIIIIIIIII"),
    ("read2", "", "TTGG", "####")
  ]

test "corpus: blank lines between records are skipped":
  check slurp("./tests/corpus/blank_lines.fq") == @[
    ("read1", "", "ACGT", "IIII"),
    ("read2", "", "TTTT", "####")
  ]

test "corpus: no final newline parses cleanly":
  check slurp("./tests/corpus/no_final_newline.fq") == @[
    ("read1", "", "ACGT", "IIII")
  ]

test "corpus: missing '+' line degrades the record to FASTA":
  let recs = slurp("./tests/corpus/missing_plus.fq")
  check recs == @[
    ("read1", "", "ACGT", ""),      # no quality: parsed as FASTA
    ("read2", "", "TTTT", "IIII")
  ]

test "corpus: malformed records stop kseq iterators silently":
  # kseq_read returns a negative code for these, which the iterators
  # treat like EOF: zero records, no exception. Documented behavior.
  for file in ["truncated_quality.fq", "truncated_after_plus.fq",
               "empty_quality.fq"]:
    var n = 0
    for rec in readFQ("./tests/corpus/" & file):
      inc n
    check n == 0

test "corpus: readFastx reports malformed records with status -4":
  # The native parser distinguishes malformed records from clean EOF.
  for file in ["truncated_quality.fq", "truncated_after_plus.fq",
               "empty_quality.fq"]:
    var r: FQRecord
    var f = xopen[GzFile]("./tests/corpus/" & file)
    check f.readFastx(r) == false
    check r.status == -4
    f.close()

# ============================================================
# Truncation fuzz: parsers must degrade, never crash
# ============================================================

test "truncation fuzz: truncated FASTQ never crashes either engine":
  let srcPath = tmpPath("readfx_fuzz_src.fq")
  if fileExists(srcPath):
    removeFile(srcPath)

  var rng = initRand(99)
  var w = fastxWriter(
    format = fxfFastq,
    compression = false,
    destination = fileDestination(srcPath)
  )
  for i in 0..<50:
    let seqLen = 1 + rng.rand(100)
    var sq = newString(seqLen)
    for c in sq.mitems:
      c = "ACGTN"[rng.rand(4)]
    w.writeRecord("r" & $i, sq, repeat('I', seqLen))
  w.close()

  let content = readFile(srcPath)
  for cut in [1, 7, 13, 50, 137, content.len div 2,
              content.len - 3, content.len - 1]:
    let truncPath = tmpPath("readfx_fuzz_trunc.fq")
    writeFile(truncPath, content[0 ..< cut])

    # kseq iterators: no exception, at most 50 records
    var n = 0
    for rec in readFQ(truncPath):
      inc n
    check n <= 50

    # readFastx: loop ends with false and a negative status, no crash
    var r: FQRecord
    var f = xopen[GzFile](truncPath)
    var m = 0
    while f.readFastx(r):
      inc m
    check m <= 50
    check r.status < 0
    f.close()
    removeFile(truncPath)

  removeFile(srcPath)

test "truncation fuzz: truncated gzip stream degrades cleanly":
  let srcPath = tmpPath("readfx_fuzz_src.fq.gz")
  if fileExists(srcPath):
    removeFile(srcPath)

  var w = fastxWriter(
    format = fxfFastq,
    compression = true,
    destination = fileDestination(srcPath)
  )
  for i in 0..<20:
    w.writeRecord("r" & $i, "ACGTACGTAC", "IIIIIIIIII")
  w.close()

  let content = readFile(srcPath)
  for cut in [1, content.len div 2, content.len - 2]:
    let truncPath = tmpPath("readfx_fuzz_trunc.fq.gz")
    writeFile(truncPath, content[0 ..< cut])

    var n = 0
    for rec in readFQ(truncPath):    # must not raise
      inc n
    check n <= 20

    var r: FQRecord
    var f = xopen[GzFile](truncPath)
    var m = 0
    while f.readFastx(r):
      inc m
    check m <= 20
    check r.status < 0               # EOF or stream error, never a crash
    f.close()
    removeFile(truncPath)

  removeFile(srcPath)
