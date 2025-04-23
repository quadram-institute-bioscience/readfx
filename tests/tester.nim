import unittest
import md5
import strutils

import ../readfx


test "readfq test.fasta.gz":
  var res = ""
  for rec in readfq("./tests/test.fasta.gz"):
    res = res & $rec & "\n"
  check $toMD5($res) == "21aa45c3b9110a7df328680f8b8753e8"#  gzip -dc tests/test.fasta.gz | md5sum


test "readfq seq.txt":
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


test "readfq SRR396637_1.seqs1-2.fastq.gz":
  var res = ""
  for rec in readfq("./tests/SRR396637_1.seqs1-2.fastq.gz"):
    res = res & $rec & "\n"
  check $toMD5($res) == "299882b15a2dc87f496a88173dd485ad"#  gzip -dc SRR396637_1.seqs1-2.fastq.gz | md5sum


test "readFQPtr test.fasta.gz":
  var res = ""
  var recs: seq[string]
  for rec in readFQPtr("./tests/test.fasta.gz"):
    # ptr char are reused but here we convert to string on the fly
    recs.add($rec)
  res = $recs.join("\n") & "\n"
  check $toMD5($res) == "21aa45c3b9110a7df328680f8b8753e8"#  gzip -dc tests/test.fasta.gz | md5sum


test "readFastx test.fasta.gz":
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
