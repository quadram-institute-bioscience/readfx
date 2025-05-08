import unittest
import md5
import os
import strutils

import ../readfx

proc checkRecord(r: FQRecord, quality = false): bool =
  var ok = true
  if r.name == "":
    ok = false
  if r.sequence == "":
    ok = false
  if quality and r.quality == "":
    ok = false
  
test "input files":
  # Check if the test files exist
  check fileExists("./tests/fastq_demo.fq")
  check fileExists("./tests/seq.txt")
  check fileExists("./tests/SRR396637_1.seqs1-2.fastq.gz")
  check fileExists("./tests/test.fasta.gz")
  check fileExists("./tests/fasta_demo.fa")

test "(1) readfq test.fasta.gz":
  var res = ""
  for rec in readfq("./tests/test.fasta.gz"):
    res = res & $rec & "\n"
  check $toMD5($res) == "21aa45c3b9110a7df328680f8b8753e8"#  gzip -dc tests/test.fasta.gz | md5sum


test "(1) readfq seq.txt":
  # tests mixed fa and fastq and messy input
  var i = 0
  for rec in readfq("./tests/seq.txt"):
    inc i
    check rec.name == $i
    if i == 1:
      check len(rec.sequence) == 15 and len(rec.quality) == 0
    elif i == 2:
      check len(rec.sequence) == 10 and len(rec.comment) > 0
    elif i == 3:
      check len(rec.quality) == len(rec.sequence)


test "(1) readfq SRR396637_1.seqs1-2.fastq.gz":
  var res = ""
  for rec in readfq("./tests/SRR396637_1.seqs1-2.fastq.gz"):
    res = res & $rec & "\n"
  check $toMD5($res) == "299882b15a2dc87f496a88173dd485ad"#  gzip -dc SRR396637_1.seqs1-2.fastq.gz | md5sum


test "(2) readFQPtr test.fasta.gz":
  var res = ""
  var recs: seq[string]
  for rec in readFQPtr("./tests/test.fasta.gz"):
    # ptr char are reused but here we convert to string on the fly
    recs.add($rec)
  res = $recs.join("\n") & "\n"
  check $toMD5($res) == "21aa45c3b9110a7df328680f8b8753e8"#  gzip -dc tests/test.fasta.gz | md5sum


test "(3) readFastx test.fasta.gz":
  var res = ""
  var r: FQRecord
  var f = xopen[GzFile]("./tests/test.fasta.gz")
  defer: f.close()
  while f.readFastx(r):
    res = res & $r & "\n"
  check $toMD5($res) == "21aa45c3b9110a7df328680f8b8753e8"#  gzip -dc tests/test.fasta.gz | md5sum


test "readFastx seq.txt":
  # tests mixed fa and fastq and messy input
  var i = 0
  var r: FQRecord
  var f = xopen[GzFile]("./tests/seq.txt")
  defer: f.close()
  while f.readFastx(r):
    inc i
    check r.name == $i
    if i == 1:
      check len(r.sequence) == 15 and len(r.quality) == 0
    elif i == 2:
      check len(r.sequence) == 10 and len(r.comment) > 0
    elif i == 3:
      check len(r.quality) == len(r.sequence)


test "readFastx SRR396637_1.seqs1-2.fastq.gz":
  var res = ""
  var r: FQRecord
  var f = xopen[GzFile]("./tests/SRR396637_1.seqs1-2.fastq.gz")
  defer: f.close()
  while f.readFastx(r):
    res = res & $r & "\n"
  check $toMD5($res) == "299882b15a2dc87f496a88173dd485ad"#  gzip -dc SRR396637_1.seqs1-2.fastq.gz | md5sum


test "FQRecord: Fasta: readfq()":
  for rec in readfq("./tests/fasta_demo.fa"):
    check len(rec.name) > 0
    check len(rec.sequence) > 0
    check len(rec.comment) > 0

test "FQRecord: Fasta: readfqptr()":
  for rec in readfqptr("./tests/fasta_demo.fa"):
    check len(rec.name) > 0
    check len(rec.sequence) > 0
    check len(rec.comment) > 0

test "FQRecord: Fasta: readfx()":
  var rec: FQRecord
  var f = xopen[GzFile]("./tests/fasta_demo.fa")
  defer: f.close()
  while f.readFastx(rec):
    check len(rec.name) > 0
    check len(rec.sequence) > 0
    check len(rec.comment) > 0
    break

#===
test "FQRecord: FASTQ: readfq()":
  for rec in readfq("./tests/fastq_demo.fq"):
    check len(rec.name) > 0
    check len(rec.sequence) > 0
    check len(rec.comment) > 0
    check len(rec.quality) == len(rec.sequence)

test "FQRecord: FASTQ: readfqptr()":
  for rec in readfqptr("./tests/fastq_demo.fq"):
    check len(rec.name) > 0
    check len(rec.sequence) > 0
    check len(rec.comment) > 0
    check len(rec.quality) == len(rec.sequence)

test "FQRecord: FASTQ: readfx()":
  var rec: FQRecord
  var f = xopen[GzFile]("./tests/fastq_demo.fq")
  defer: f.close()
  while f.readFastx(rec):
    check len(rec.name) > 0
    check len(rec.sequence) > 0
    check len(rec.comment) > 0
    check len(rec.quality) == len(rec.sequence)
    break
test "utils: revCompl()":
  var r = FQRecord()
  r.sequence = "GAAA"
  r.quality  = "IIIA"
  r.revCompl()
  check r.sequence == "TTTC"
  check r.quality == "AIII"

test "utils: maskLowQuality()":
  var r = FQRecord()
  r.sequence = "GAAATTT"
  r.quality  = "IIIA888"
  r.maskLowQuality(25)
  check r.sequence == "GAAANNN"
  check r.quality == "IIIA888"

test "utils: trimStart()":
  block:
    var r = FQRecord()
    r.sequence = "GAAATTT"
    r.quality  = "IIIA888"
    let result = trimStart(r, 3)
    check result.sequence == "ATTT"
    check result.quality == "A888"
  
  # Zero or negative bases parameter
  block:
    var r = FQRecord()
    r.sequence = "GAAATTT"
    r.quality  = "IIIA888"
    let result = trimStart(r, 0)
    check result.sequence == "GAAATTT"
    check result.quality == "IIIA888"
  
  # Trim more bases than sequence length
  block:
    var r = FQRecord()
    r.sequence = "GAAATTT"
    r.quality  = "IIIA888"
    let result = trimStart(r, 10)
    check result.sequence == ""
    check result.quality  == ""
  
  # FASTA record (no quality)
  block:
    var r = FQRecord()
    r.sequence = "GAAATTT"
    r.quality  = ""
    let result = trimStart(r, 3)
    check result.sequence == "ATTT"
    check result.quality  == ""

test "utils: trimEnd()":
  block:
    var r = FQRecord()
    r.sequence = "GAAATTT"
    r.quality  = "IIIA888"
    let result = trimEnd(r, 3)
    check result.sequence == "GAAA"
    check result.quality  == "IIIA"
  
  # Zero or negative bases parameter
  block:
    var r = FQRecord()
    r.sequence = "GAAATTT"
    r.quality  = "IIIA888"
    let result = trimEnd(r, 0)
    check result.sequence == "GAAATTT"
    check result.quality  == "IIIA888"
  
  # Trim more bases than sequence length
  block:
    var r = FQRecord()
    r.sequence = "GAAATTT"
    r.quality  = "IIIA888"
    let result = trimEnd(r, 10)
    check result.sequence == ""
    check result.quality == ""
  
  # FASTA record (no quality)
  block:
    var r = FQRecord()
    r.sequence = "GAAATTT"
    r.quality  = ""
    let result = trimEnd(r, 3)
    check result.sequence == "GAAA"
    check result.quality  == ""


test "utils: filtPolyX()":
  # Example with long poly-A tail
  block:
    var r = FQRecord()
    r.sequence = "ACGTAAAAAAAAAAAAA" # 5 bases + 12 As
    r.quality  = "IIIIIIIIIIIIIIIII"
    let result = filtPolyX(r, minLen = 10)
    # Should trim (poly-A length is 12, minLen is 10)
    check result.sequence == "ACGT"
    check result.quality == "IIII"
  
  # Example with poly-A containing mismatches
  block:
    var r = FQRecord()
    r.sequence = "ACGTAAAAACAAAAAAAA" # 5 bases + A's with C mismatch
    r.quality  = "IIIIIIIIIIIIIIIIIII"
    let result = filtPolyX(r, minLen = 10)
    # Should trim despite mismatch (12 As with 1 mismatch)
    check result.sequence == "ACGT"
    check result.quality == "IIII"
  
  # Example with lowercase bases
  block:
    var r = FQRecord()
    r.sequence = "ACGTaaaaaaaaaaa" # 5 bases + 10 lowercase a's
    r.quality  = "IIIIIIIIIIIIIII"
    let result = filtPolyX(r, minLen = 10)
    # Should trim (lowercase should be handled)
    check result.sequence == "ACGT"
    check result.quality == "IIII"
  
  # Example with poly-T instead of poly-A
  block:
    var r = FQRecord()
    r.sequence = "ACGATTTTTTTTTTTTT" # 5 bases + 12 T's
    r.quality  = "IIIIIIIIIIIIIIIII"
    let result = filtPolyX(r, minLen = 10)
    # Should trim poly-T
    check result.sequence == "ACGA"
    check result.quality == "IIII"
  
  # Example with N's in poly-A tail
  block:
    var r = FQRecord()
    r.sequence = "ACGTAAAAANAAAAAA" # 5 bases + As with N
    r.quality  = "IIIIIIIIIIIIIIIII"
    let result = filtPolyX(r, minLen = 10)
    # Should trim (N counts as A)
    check result.sequence == "ACGT"
    check result.quality == "IIII"
  
  # Example with short read
  block:
    var r = FQRecord()
    r.sequence = "AAA"
    r.quality  = "III"
    let result = filtPolyX(r, minLen = 10)
    # Shouldn't trim (read too short)
    check result.sequence == "AAA"
    check result.quality == "III"
  
  # Example with empty sequence
  block:
    var r = FQRecord()
    r.sequence = ""
    r.quality  = ""
    let result = filtPolyX(r, minLen = 10)
    # Should handle empty sequence gracefully
    check result.sequence == ""
    check result.quality == ""

test "utils: filtPolyX() minlen":
  # Basic example with simple poly-A tail
  block: 
    var r = FQRecord()
    r.sequence = "ACGTAAAAAAAAA" # 5 bases + 8 As
    r.quality  = "IIIIIIIIIIIII"
    let result = filtPolyX(r, minLen = 10)
    # Shouldn't trim (poly-A length is 8, minLen is 10)
    check result.sequence == "ACGTAAAAAAAAA"
    check result.quality == "IIIIIIIIIIIII"

test "utils: filtPolyX() only tail!":
  # Basic example with simple poly-A tail
  block: 
    var r = FQRecord()
    r.sequence = "AAAAAAAAAAAAA" # 5 bases + 8 As
    r.quality  = "IIIIIIIIIIIII"
    let result = filtPolyX(r, minLen = 10)
    # Shouldn't trim (poly-A length is 8, minLen is 10)
    check result.sequence == ""
    check result.quality == ""

test "seq content":
  var r = FQRecord()
  r.sequence = "ACGT"
  r.quality  = "IIII"
  var sc = SeqComp()
  sc = composition(r)
  check sc.A == 1
  check sc.C == 1
  check sc.G == 1
  check sc.T == 1
  check sc.GC == 0.5

  r.sequence = "aaaaatttttccggNNNnnnnnnnnnnxx"
  r.quality  = "IIIIIIIIIIIIIIIIIIIIIIIIIIIII"
  sc = composition(r)
  check sc.A == 5
  check sc.C == 2
  check sc.G == 2
  check sc.T == 5
  check sc.GC == float(4/14)
  check sc.N == 13
  check sc.Other == 2

test "GC content":
  var r = FQRecord()
  r.sequence = "GGGGGGGTA"
  r.quality  = "IIIIIII88"

  var sc = composition(r)
  var gc = gcContent(r.sequence)
  check sc.GC == gc