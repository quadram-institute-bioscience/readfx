import unittest
import md5
import os
import strutils
when defined(posix):
  import posix

import ../readfx
import gzfast

proc cstrOrEmpty(p: ptr char): string =
  if p.isNil:
    ""
  else:
    $cast[cstring](p)

proc writeGzipText(path, text: string) =
  var w = gzfast.openGzFastWriter(path)
  try:
    discard w.writeString(text)
  finally:
    w.close()
  
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

test "(2) readFQPtr missing file raises IOError":
  expect IOError:
    for rec in readFQPtr("./tests/nonexistent.fastq.gz"):
      discard rec

test "(2) readFQPtr cached lengths":
  var count = 0
  for rec in readFQPtr("./tests/fastq_demo.fq"):
    inc count
    check rec.nameLen == cstrOrEmpty(rec.name).len
    check rec.commentLen == cstrOrEmpty(rec.comment).len
    check rec.sequenceLen == cstrOrEmpty(rec.sequence).len
    check rec.qualityLen == cstrOrEmpty(rec.quality).len
  check count > 0

test "(2) readFQPtr '-' keeps stdin open":
  when defined(posix):
    var pipeFds: array[2, cint]
    check posix.pipe(pipeFds) == 0
    let savedStdin = posix.dup(0)
    check savedStdin >= 0

    let input = "@r1\nA\n+\nI\n"
    discard posix.write(pipeFds[1], cast[pointer](unsafeAddr input[0]), input.len)
    discard posix.close(pipeFds[1])
    check posix.dup2(pipeFds[0], 0) >= 0
    discard posix.close(pipeFds[0])

    var count = 0
    try:
      for rec in readFQPtr("-"):
        discard rec
        inc count
      check count == 1

      let probe = posix.dup(0)
      check probe >= 0
      if probe >= 0:
        discard posix.close(probe)
    finally:
      check posix.dup2(savedStdin, 0) >= 0
      discard posix.close(savedStdin)

test "(2) readFQPtr '-' supports gzip stdin":
  when defined(posix):
    var pipeFds: array[2, cint]
    check posix.pipe(pipeFds) == 0
    let savedStdin = posix.dup(0)
    check savedStdin >= 0

    let input = readFile("./tests/test.fasta.gz")
    discard posix.write(pipeFds[1], cast[pointer](unsafeAddr input[0]), input.len)
    discard posix.close(pipeFds[1])
    check posix.dup2(pipeFds[0], 0) >= 0
    discard posix.close(pipeFds[0])

    var res = ""
    var recs: seq[string]
    try:
      for rec in readFQPtr("-"):
        recs.add($rec)
      res = $recs.join("\n") & "\n"
      check $toMD5($res) == "21aa45c3b9110a7df328680f8b8753e8"

      let probe = posix.dup(0)
      check probe >= 0
      if probe >= 0:
        discard posix.close(probe)
    finally:
      check posix.dup2(savedStdin, 0) >= 0
      discard posix.close(savedStdin)


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
    let name = cstrOrEmpty(rec.name)
    let sequence = cstrOrEmpty(rec.sequence)
    let comment = cstrOrEmpty(rec.comment)
    check len(name) > 0
    check len(sequence) > 0
    check len(comment) > 0

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
    let name = cstrOrEmpty(rec.name)
    let sequence = cstrOrEmpty(rec.sequence)
    let comment = cstrOrEmpty(rec.comment)
    let quality = cstrOrEmpty(rec.quality)
    check len(name) > 0
    check len(sequence) > 0
    check len(comment) > 0
    check len(quality) == len(sequence)

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
test "readFastx: reuse FQRecord across two FASTA files":
  var r: FQRecord
  var f1 = xopen[GzFile]("./tests/test.fasta.gz")
  while f1.readFastx(r): discard  # exhausts f1, leaves lastChar stale without fix
  f1.close()
  var f2 = xopen[GzFile]("./tests/fasta_demo.fa")
  defer: f2.close()
  if f2.readFastx(r):
    check not r.name.startsWith(">")
    check len(r.name) > 0

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

  r.sequence = "NNNxxx"
  r.quality  = "IIIIII"
  sc = composition(r)
  check sc.A == 0
  check sc.C == 0
  check sc.G == 0
  check sc.T == 0
  check sc.N == 3
  check sc.Other == 3
  check sc.GC == 0.0

test "GC content":
  var r = FQRecord()
  r.sequence = "GGGGGGGTA"
  r.quality  = "IIIIIII88"

  var sc = composition(r)
  var gc = gcContent(r.sequence)
  check sc.GC == gc

test "fastxWriter plain FASTQ":
  let outPath = getTempDir() / "readfx_writer_plain.fastq"
  if fileExists(outPath):
    removeFile(outPath)

  var w = fastxWriter(
    format = fxfFastq,
    compression = false,
    destination = fileDestination(outPath),
    bufferSize = 32
  )
  defer:
    if w.isOpen:
      w.close()
    if fileExists(outPath):
      removeFile(outPath)

  w.writeRecord("r1", "ACGT", "IIII")
  w.writeRecord(FQRecord(name: "r2", comment: "paired", sequence: "TTAA", quality: "####"))
  w.close()

  var names: seq[string] = @[]
  var quals: seq[string] = @[]
  for rec in readFQ(outPath):
    names.add(rec.name)
    quals.add(rec.quality)

  check names == @["r1", "r2"]
  check quals == @["IIII", "####"]

test "fastxWriter gzip FASTQ":
  let outPath = getTempDir() / "readfx_writer.fastq.gz"
  if fileExists(outPath):
    removeFile(outPath)

  var w = fastxWriter(
    format = fxfFastq,
    compression = true,
    destination = fileDestination(outPath),
    bufferSize = 32,
    compressionLevel = 6,
    compressionThreads = 2
  )
  defer:
    if w.isOpen:
      w.close()
    if fileExists(outPath):
      removeFile(outPath)

  w.writeRecord("gz1", "AACCGG", "IIIIII")
  w.close()

  var n = 0
  for rec in readFQ(outPath):
    inc n
    check rec.name == "gz1"
    check rec.sequence == "AACCGG"
    check rec.quality == "IIIIII"
  check n == 1

test "fastxWriter gzip stdout":
  when defined(posix):
    var pipeFds: array[2, cint]
    check posix.pipe(pipeFds) == 0
    var readFd = pipeFds[0]
    var writeFd = pipeFds[1]
    var savedStdout = posix.dup(1)
    check savedStdout >= 0

    var compressed = ""
    let outPath = getTempDir() / "readfx_writer_stdout.fastq.gz"
    if fileExists(outPath):
      removeFile(outPath)

    try:
      flushFile(stdout)
      check posix.dup2(writeFd, 1) >= 0
      discard posix.close(writeFd)
      writeFd = -1

      var w = fastxWriter(
        format = fxfFastq,
        compression = true,
        destination = stdoutDestination(),
        bufferSize = 32,
        compressionLevel = 6
      )
      w.writeRecord("stdout_gz", "ACGT", "IIII")
      w.close()

      check posix.dup2(savedStdout, 1) >= 0
      discard posix.close(savedStdout)
      savedStdout = -1

      var chunk: array[1024, char]
      while true:
        let n = posix.read(readFd, addr chunk[0], chunk.len)
        if n <= 0:
          break
        let oldLen = compressed.len
        compressed.setLen(oldLen + int(n))
        copyMem(addr compressed[oldLen], addr chunk[0], int(n))
      discard posix.close(readFd)
      readFd = -1

      check compressed.len > 0
      writeFile(outPath, compressed)

      var n = 0
      for rec in readFQ(outPath):
        inc n
        check rec.name == "stdout_gz"
        check rec.sequence == "ACGT"
        check rec.quality == "IIII"
      check n == 1
    finally:
      if savedStdout >= 0:
        discard posix.dup2(savedStdout, 1)
        discard posix.close(savedStdout)
      if writeFd >= 0:
        discard posix.close(writeFd)
      if readFd >= 0:
        discard posix.close(readFd)
      if fileExists(outPath):
        removeFile(outPath)

test "fastxWriter FASTA ignores quality":
  let outPath = getTempDir() / "readfx_writer.fasta"
  if fileExists(outPath):
    removeFile(outPath)

  var w = fastxWriter(
    format = fxfFasta,
    compression = false,
    destination = fileDestination(outPath),
    bufferSize = 16,
    fastaWidth = 3
  )
  defer:
    if w.isOpen:
      w.close()
    if fileExists(outPath):
      removeFile(outPath)

  w.writeRecord("fa1", "ACGTAC", "!!!!!!", "comment")
  w.close()

  var n = 0
  for rec in readFQ(outPath):
    inc n
    check rec.name == "fa1"
    check rec.comment == "comment"
    check rec.sequence == "ACGTAC"
    check rec.quality == ""
  check n == 1

test "fastxWriter FASTQ validates quality":
  let outPath = getTempDir() / "readfx_writer_invalid.fastq"
  if fileExists(outPath):
    removeFile(outPath)

  var w = fastxWriter(
    format = fxfFastq,
    compression = false,
    destination = fileDestination(outPath)
  )
  defer:
    if w.isOpen:
      w.close()
    if fileExists(outPath):
      removeFile(outPath)

  expect ValueError:
    w.writeRecord("bad1", "ACGT", "")

test "fastxWriter FQRecordPtr round-trip (FASTQ)":
  let outPath = getTempDir() / "readfx_writer_ptr.fastq"
  if fileExists(outPath):
    removeFile(outPath)

  var w = fastxWriter(
    format = fxfFastq,
    compression = false,
    destination = fileDestination(outPath),
    bufferSize = 64
  )
  defer:
    if w.isOpen:
      w.close()
    if fileExists(outPath):
      removeFile(outPath)

  # zero-copy pipeline: pointer records written directly
  for rec in readFQPtr("./tests/fastq_demo.fq"):
    w.writeRecord(rec)
  w.close()

  # written output must match a direct string-based read of the source
  var expected: seq[tuple[name, sequence, quality: string]]
  for rec in readFQ("./tests/fastq_demo.fq"):
    expected.add((rec.name, rec.sequence, rec.quality))

  var got: seq[tuple[name, sequence, quality: string]]
  for rec in readFQ(outPath):
    got.add((rec.name, rec.sequence, rec.quality))

  check got == expected
  check got.len > 0

test "fastxWriter FQRecordPtr round-trip (FASTA)":
  let outPath = getTempDir() / "readfx_writer_ptr.fasta"
  if fileExists(outPath):
    removeFile(outPath)

  var w = fastxWriter(
    format = fxfFasta,
    compression = false,
    destination = fileDestination(outPath),
    bufferSize = 64,
    fastaWidth = 60
  )
  defer:
    if w.isOpen:
      w.close()
    if fileExists(outPath):
      removeFile(outPath)

  for rec in readFQPtr("./tests/fasta_demo.fa"):
    w.writeRecord(rec)
  w.close()

  var expected: seq[tuple[name, sequence: string]]
  for rec in readFQ("./tests/fasta_demo.fa"):
    expected.add((rec.name, rec.sequence))

  var got: seq[tuple[name, sequence: string]]
  for rec in readFQ(outPath):
    check rec.quality == ""        # FASTA output carries no qualities
    got.add((rec.name, rec.sequence))

  check got == expected
  check got.len > 0

test "fastxWriter FQRecordPtr FASTQ validates quality":
  let outPath = getTempDir() / "readfx_writer_ptr_invalid.fastq"
  if fileExists(outPath):
    removeFile(outPath)

  var w = fastxWriter(
    format = fxfFastq,
    compression = false,
    destination = fileDestination(outPath)
  )
  defer:
    if w.isOpen:
      w.close()
    if fileExists(outPath):
      removeFile(outPath)

  # a FASTA record (no qualities) must be rejected in FASTQ mode
  expect ValueError:
    for rec in readFQPtr("./tests/fasta_demo.fa"):
      w.writeRecord(rec)

#===
# Bufio low-level API
#===

test "Bufio: readUntil field/line modes, readLine, readByte, eof":
  let path = getTempDir() / "readfx_bufio.txt"
  writeFile(path, "field1 field2\tfield3\nline two\r\nlast line")
  defer:
    if fileExists(path):
      removeFile(path)

  var f = xopen[GzFile](path)
  defer: f.close()
  var buf = ""
  var dret: char

  # field mode (-2) stops at space
  check f.readUntil(buf, dret, -2) == 6
  check buf == "field1"
  check dret == ' '
  # field mode (-2) stops at tab
  check f.readUntil(buf, dret, -2) == 6
  check buf == "field2"
  check dret == '\t'
  # line mode (-1) reads to newline
  check f.readUntil(buf, dret, -1) == 6
  check buf == "field3"
  check dret == '\n'
  # line mode (-1) strips a trailing CR
  check f.readUntil(buf, dret, -1) == 8
  check buf == "line two"
  # readLine at EOF without trailing newline
  check f.readLine(buf)
  check buf == "last line"
  # EOF signaling
  check f.readByte() == -1
  check f.eof()

test "Bufio: readUntil custom delimiter":
  let path = getTempDir() / "readfx_bufio2.txt"
  writeFile(path, "abc:def:ghi")
  defer:
    if fileExists(path):
      removeFile(path)

  var f = xopen[GzFile](path)
  defer: f.close()
  var buf = ""
  var dret: char

  check f.readUntil(buf, dret, int(':')) == 3
  check buf == "abc"
  check dret == ':'
  check f.readUntil(buf, dret, int(':')) == 3
  check buf == "def"
  check f.readUntil(buf, dret, int(':')) == 3  # last field, ends at EOF
  check buf == "ghi"

test "GzFile: direct reads support non-FASTX gzip content":
  let path = getTempDir() / "readfx_arbitrary.txt.gz"
  let text = "alpha,beta\nplain text\nnot fastx\n"
  writeGzipText(path, text)
  defer:
    if fileExists(path):
      removeFile(path)

  var f: GzFile
  f.open(path)
  defer: f.close()

  var got = ""
  var chunk = ""
  while true:
    let n = f.read(chunk, 5)
    check n >= 0
    if n == 0:
      break
    got.add(chunk)
  check got == text

test "readGzChunks reads arbitrary gzip and plain content":
  let gzPath = getTempDir() / "readfx_chunks.txt.gz"
  let plainPath = getTempDir() / "readfx_chunks.txt"
  let text = "one|two|three|four"
  writeGzipText(gzPath, text)
  writeFile(plainPath, text)
  defer:
    if fileExists(gzPath):
      removeFile(gzPath)
    if fileExists(plainPath):
      removeFile(plainPath)

  var gzJoined = ""
  for chunk in readGzChunks(gzPath, chunkSize = 4):
    check chunk.len <= 4
    gzJoined.add(chunk)
  check gzJoined == text

  var plainJoined = ""
  for chunk in readGzChunks(plainPath, chunkSize = 6):
    check chunk.len <= 6
    plainJoined.add(chunk)
  check plainJoined == text

test "readGzLines reads arbitrary gzip lines":
  let path = getTempDir() / "readfx_lines.txt.gz"
  writeGzipText(path, "header\tvalue\r\nx\t1\nlast")
  defer:
    if fileExists(path):
      removeFile(path)

  var lines: seq[string]
  for line in readGzLines(path):
    lines.add(line)
  check lines == @["header\tvalue", "x\t1", "last"]

test "readGzChunks rejects invalid chunk size":
  expect ValueError:
    for chunk in readGzChunks("./tests/seq.txt", chunkSize = 0):
      discard chunk

#===
# readFastx edge cases
#===

test "readFastx: malformed FASTQ (quality shorter than sequence) sets status -4":
  let path = getTempDir() / "readfx_badqual.fq"
  writeFile(path, "@r1\nACGT\n+\nII\n")
  defer:
    if fileExists(path):
      removeFile(path)

  var r: FQRecord
  var f = xopen[GzFile](path)
  defer: f.close()
  check f.readFastx(r) == false
  check r.status == -4

test "readFastx: wrapped multi-line FASTA and EOF status":
  let path = getTempDir() / "readfx_wrapped.fa"
  writeFile(path, ">seq1 comment here\nACGT\nACGT\nAC\n>seq2\nTT\nGG\n")
  defer:
    if fileExists(path):
      removeFile(path)

  var r: FQRecord
  var f = xopen[GzFile](path)
  defer: f.close()
  check f.readFastx(r)
  check r.name == "seq1"
  check r.comment == "comment here"
  check r.sequence == "ACGTACGTAC"
  check f.readFastx(r)
  check r.name == "seq2"
  check r.sequence == "TTGG"
  check f.readFastx(r) == false
  check r.status == -1

test "readFastx: CRLF line endings":
  let path = getTempDir() / "readfx_crlf.fq"
  # r1 has a comment; r2 has no comment (name directly followed by CR)
  writeFile(path, "@r1 cmt\r\nACGT\r\n+\r\nIIII\r\n@r2\r\nAC\r\n+\r\nII\r\n")
  defer:
    if fileExists(path):
      removeFile(path)

  var r: FQRecord
  var f = xopen[GzFile](path)
  defer: f.close()
  check f.readFastx(r)
  check r.name == "r1"
  check r.comment == "cmt"
  check r.sequence == "ACGT"
  check r.quality == "IIII"
  check f.readFastx(r)
  check r.name == "r2"          # no trailing CR in the name
  check r.comment == ""
  check r.sequence == "AC"
  check r.quality == "II"
