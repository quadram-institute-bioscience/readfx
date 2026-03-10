import unittest
import os
when defined(posix):
  import posix

import ../readfx

proc cstrOrEmpty(p: ptr char): string =
  if p.isNil:
    ""
  else:
    $cast[cstring](p)

proc checkPtrLengths(rec: FQRecordPtr) =
  check rec.nameLen == cstrOrEmpty(rec.name).len
  check rec.commentLen == cstrOrEmpty(rec.comment).len
  check rec.sequenceLen == cstrOrEmpty(rec.sequence).len
  check rec.qualityLen == cstrOrEmpty(rec.quality).len

proc writeTempInterleavedGzip(srcPath: string): string =
  result = getTempDir() / "readfx_interleaved.fastq.gz"
  if fileExists(result):
    removeFile(result)

  var w = fastxWriter(
    format = fxfFastq,
    compression = true,
    destination = fileDestination(result)
  )
  defer:
    if w.isOpen:
      w.close()

  for rec in readFQ(srcPath):
    w.writeRecord(rec)
  w.close()

suite "Paired-end reading tests":

  test "readFQPairPtr basic functionality (uncompressed synthetic data)":
    var count = 0
    for pair in readFQPairPtr("./tests/test_R1.fq", "./tests/test_R2.fq"):
      count += 1
      let seq1 = cstrOrEmpty(pair.read1.sequence)
      let seq2 = cstrOrEmpty(pair.read2.sequence)
      let qual1 = cstrOrEmpty(pair.read1.quality)
      let qual2 = cstrOrEmpty(pair.read2.quality)
      check seq1.len + seq2.len == 10
      check cstrOrEmpty(pair.read1.name).len > 0
      check cstrOrEmpty(pair.read2.name).len > 0
      check seq1.len > 0
      check seq2.len > 0
      check qual1.len > 0
      check qual2.len > 0
      check seq1.len == qual1.len
      check seq2.len == qual2.len
      checkPtrLengths(pair.read1)
      checkPtrLengths(pair.read2)
    check count > 0

  test "readFQPairPtr basic functionality":
    var count = 0
    var totalLen1, totalLen2 = 0

    for pair in readFQPairPtr("./tests/illumina_1.fq.gz", "./tests/illumina_2.fq.gz"):
      count += 1
      let seq1 = cstrOrEmpty(pair.read1.sequence)
      let seq2 = cstrOrEmpty(pair.read2.sequence)
      let qual1 = cstrOrEmpty(pair.read1.quality)
      let qual2 = cstrOrEmpty(pair.read2.quality)
      totalLen1 += seq1.len
      totalLen2 += seq2.len

      check cstrOrEmpty(pair.read1.name).len > 0
      check cstrOrEmpty(pair.read2.name).len > 0
      check seq1.len > 0
      check seq2.len > 0
      check qual1.len > 0
      check qual2.len > 0
      check seq1.len == qual1.len
      check seq2.len == qual2.len
      checkPtrLengths(pair.read1)
      checkPtrLengths(pair.read2)

    check count == 7
    check totalLen1 > 0
    check totalLen2 > 0

  test "readFQPairPtr name checking":
    var count = 0
    for pair in readFQPairPtr("./tests/illumina_1.fq.gz", "./tests/illumina_2.fq.gz", checkNames = true):
      discard pair
      count += 1
    check count == 7

  test "readFQPairPtr error handling for non-existent files":
    expect IOError:
      for pair in readFQPairPtr("./tests/nonexistent1.fq", "./tests/nonexistent2.fq"):
        discard pair

  test "readFQPairPtr error handling for single non-existent file":
    expect IOError:
      for pair in readFQPairPtr("./tests/illumina_1.fq.gz", "./tests/nonexistent.fq"):
        discard pair

  test "readFQInterleavedPairPtr basic functionality":
    var count = 0
    for pair in readFQInterleavedPairPtr("./tests/test_interleaved.fq"):
      inc count
      check cstrOrEmpty(pair.read1.name) == "read" & $count
      check cstrOrEmpty(pair.read2.name) == "read" & $count
      check pair.read1.sequenceLen + pair.read2.sequenceLen == 10
      check pair.read1.sequenceLen == pair.read1.qualityLen
      check pair.read2.sequenceLen == pair.read2.qualityLen
      checkPtrLengths(pair.read1)
      checkPtrLengths(pair.read2)
    check count == 3

  test "readFQInterleavedPairPtr name checking":
    var count = 0
    for pair in readFQInterleavedPairPtr("./tests/test_interleaved.fq", checkNames = true):
      discard pair
      inc count
    check count == 3

  test "readFQInterleavedPairPtr gzip input":
    let path = writeTempInterleavedGzip("./tests/test_interleaved.fq")
    defer:
      if fileExists(path):
        removeFile(path)

    var count = 0
    for pair in readFQInterleavedPairPtr(path, checkNames = true):
      discard pair
      inc count
    check count == 3

  test "readFQInterleavedPairPtr odd record count raises IOError":
    expect IOError:
      for pair in readFQInterleavedPairPtr("./tests/test_interleaved_odd.fq"):
        discard pair

  test "readFQInterleavedPairPtr mismatched names raise ValueError":
    let path = getTempDir() / "readfx_interleaved_bad_names.fastq"
    writeFile(path, "@r1/1\nAAAA\n+\nIIII\n@other/2\nTTTT\n+\n####\n")
    defer:
      if fileExists(path):
        removeFile(path)

    expect ValueError:
      for pair in readFQInterleavedPairPtr(path, checkNames = true):
        discard pair

  test "readFQInterleavedPairPtr supports stdin":
    when defined(posix):
      var pipeFds: array[2, cint]
      check posix.pipe(pipeFds) == 0
      let savedStdin = posix.dup(0)
      check savedStdin >= 0

      let input = readFile("./tests/test_interleaved.fq")
      discard posix.write(pipeFds[1], cast[pointer](input.cstring), input.len)
      discard posix.close(pipeFds[1])
      check posix.dup2(pipeFds[0], 0) >= 0
      discard posix.close(pipeFds[0])

      var count = 0
      try:
        for pair in readFQInterleavedPairPtr("-", checkNames = true):
          discard pair
          inc count
        check count == 3

        let probe = posix.dup(0)
        check probe >= 0
        if probe >= 0:
          discard posix.close(probe)
      finally:
        check posix.dup2(savedStdin, 0) >= 0
        discard posix.close(savedStdin)

  test "readFQInterleavedPairPtr pointers stay valid until next yield":
    iterator nextInterleavedPair(): FQPairPtr {.closure.} =
      for pair in readFQInterleavedPairPtr("./tests/test_interleaved.fq", checkNames = true):
        yield pair

    var nextPair = nextInterleavedPair
    let first = nextPair()
    check cstrOrEmpty(first.read1.name) == "read1"
    check cstrOrEmpty(first.read1.sequence) == "AAAA"
    check cstrOrEmpty(first.read2.sequence) == "AAAAAA"

    let second = nextPair()
    check cstrOrEmpty(second.read1.name) == "read2"
    check cstrOrEmpty(second.read1.sequence) == "AA"
    check cstrOrEmpty(second.read2.sequence) == "AAAAAAAA"

  test "readFQPair basic functionality (uncompressed synthetic data)":
    var count = 0
    
    for pair in readFQPair("./tests/test_R1.fq", "./tests/test_R2.fq"):
      count += 1
      # this dataset is made so that the total length of both reads is 10
      check pair.read1.sequence.len  + pair.read2.sequence.len == 10
    
      
      check pair.read1.name.len > 0
      check pair.read2.name.len > 0
      check pair.read1.sequence.len > 0
      check pair.read2.sequence.len > 0
      check pair.read1.quality.len > 0
      check pair.read2.quality.len > 0
      
  test "readFQPair basic functionality":
    var count = 0
    var totalLen1, totalLen2 = 0
    
    for pair in readFQPair("./tests/illumina_1.fq.gz", "./tests/illumina_2.fq.gz"):
      count += 1
      totalLen1 += pair.read1.sequence.len
      totalLen2 += pair.read2.sequence.len
      
      check pair.read1.name.len > 0
      check pair.read2.name.len > 0
      check pair.read1.sequence.len > 0
      check pair.read2.sequence.len > 0
      check pair.read1.quality.len > 0
      check pair.read2.quality.len > 0
      
      # Check that quality and sequence lengths match
      check pair.read1.sequence.len == pair.read1.quality.len
      check pair.read2.sequence.len == pair.read2.quality.len
    
    check count == 7  # Should have 7 read pairs
    check totalLen1 > 0
    check totalLen2 > 0
  
  test "readFQPair name checking":
    # This should work without exceptions since the names match (base names)
    var count = 0
    for pair in readFQPair("./tests/illumina_1.fq.gz", "./tests/illumina_2.fq.gz", checkNames = true):
      count += 1
    check count == 7
  
  test "FQPair type properties":
    for pair in readFQPair("./tests/illumina_1.fq.gz", "./tests/illumina_2.fq.gz"):
      # Test that we can access all fields
      discard pair.read1.name
      discard pair.read1.comment
      discard pair.read1.sequence
      discard pair.read1.quality
      discard pair.read2.name
      discard pair.read2.comment
      discard pair.read2.sequence
      discard pair.read2.quality
      break  # Just test first pair
  
  test "Error handling for non-existent files":
    expect IOError:
      for pair in readFQPair("./tests/nonexistent1.fq", "./tests/nonexistent2.fq"):
        discard pair
  
  test "Error handling for single non-existent file":
    expect IOError:
      for pair in readFQPair("./tests/illumina_1.fq.gz", "./tests/nonexistent.fq"):
        discard pair
