## ReadFX: A Nim library for bioinformatics sequence file parsing
## ===============================================================
##
## This module provides efficient parsing and manipulation of FASTA/FASTQ format
## sequence files commonly used in bioinformatics.
##
## Features:
## * Fast FASTA/FASTQ sequence parsing (supports gzipped files)
## * Buffered I/O for efficient file reading
## * Interval tree implementation for genomic interval operations
##
## Example:
##
## ```nim
## import readfx
##
## # Read a FASTQ file
## for record in readFQ("sample.fastq.gz"):
##   echo "Sequence name: ", record.name
##   echo "Sequence: ", record.sequence
##   echo "Quality: ", record.quality
##
## ```
##
import strutils
when defined(posix):
  import posix
import readfx/seqtypes
export seqtypes

import readfx/sequtils
export sequtils

import readfx/oligoutils
export oligoutils

import readfx/nimklib
export nimklib

import readfx/writer
export writer
# https://forum.nim-lang.org/t/2668
from os import splitPath
const kseqh = currentSourcePath().splitPath.head & "/readfx/kseq.h"

# https://github.com/nim-lang/nimble/issues/157
{.passL: "-lz".}


type
  kstring_t {.importc, header: kseqh.} = object
    # C definition: `size_t l, m; char *s;`
    l: csize_t
    m: csize_t
    s: ptr char
  kstream_t {.importc, header: kseqh.} = object
    # Only ever used as an opaque pointer (`ptr kstream_t`), so exact
    # field layout is irrelevant; C fields are `int` -> cint.
    begin: cint
    endd: cint
    is_eof: cint
  kseq_t {.importc, header: kseqh.} = object
    name: kstring_t
    comment: kstring_t
    seq: kstring_t
    qual: kstring_t
    last_char: cint
    f: ptr kstream_t
  gzFile = pointer

# ------------------------------------------------------------------
# Private zlib bindings for the kseq C wrapper.
# The native Nim parser (`readFastx`, `Bufio`, `Interval`) and its own
# GzFile handling live in readfx/nimklib.nim and are re-exported above.
# ------------------------------------------------------------------

when defined(windows):
  const libz = "zlib1.dll"
elif defined(macosx):
  const libz = "libz.dylib"
else:
  const libz = "libz.so.1"

proc gzopen(path: cstring, mode: cstring): gzFile{.cdecl, dynlib: libz,
    importc: "gzopen".}
proc gzdopen(fd: int32, mode: cstring): gzFile{.cdecl, dynlib: libz,
    importc: "gzdopen".}
proc gzclose(thefile: gzFile): int32{.cdecl, dynlib: libz, importc: "gzclose".}

## Initialize a kseq parser handle from an open gzFile stream.
##
## Args:
##   fp: Open gzip/plain-text stream handle
##
## Returns:
##   Pointer to an initialized parser state
proc kseq_init*(fp: gzFile): ptr kseq_t {.header: kseqh, importc: "kseq_init".}


## Reset parser state to the beginning of the input stream.
##
## Args:
##   seq: Parser state previously created with `kseq_init`
proc kseq_rewind*(seq: ptr kseq_t) {.header: kseqh, importc: "kseq_rewind".}


## Read the next FASTA/FASTQ record from the parser stream.
##
## Args:
##   seq: Parser state previously created with `kseq_init`
##
## Returns:
##   Record sequence length on success, or a negative status code on EOF/error
proc kseq_read*(seq: ptr kseq_t): cint {.header: kseqh, importc: "kseq_read".}

proc openGzForRead(path: string): GzFile =
  if path == "-":
    when defined(posix):
      let stdinDup = posix.dup(0)
      if stdinDup < 0:
        raise newException(IOError, "Cannot duplicate stdin for reading")
      result = gzdopen(int32(stdinDup), "r")
      if result == nil:
        discard posix.close(stdinDup)
        raise newException(IOError, "Cannot open stdin for reading")
    else:
      result = gzdopen(0, "r")
      if result == nil:
        raise newException(IOError, "Cannot open stdin for reading")
  else:
    result = gzopen(path, "r")
    if result == nil:
      raise newException(IOError, "Cannot open file: " & path)

proc cstrOrEmpty(p: ptr char): string {.inline.} =
  if p.isNil:
    ""
  else:
    $cast[cstring](p)

proc normalizePairName(name: string, mate: int): string {.inline.} =
  result = name
  if mate == 1 and (result.endsWith("/1") or result.endsWith(" 1")):
    result = result[0..^3]
  elif mate == 2 and (result.endsWith("/2") or result.endsWith(" 2")):
    result = result[0..^3]

proc ptrFieldLen(field: kstring_t): int {.inline.} =
  if field.s.isNil:
    0
  else:
    # field.l is a C size_t; record fields never approach int range limits
    int(field.l)

proc setRecordPtrFields(
    record: var FQRecordPtr,
    name: ptr char, nameLen: int,
    comment: ptr char, commentLen: int,
    sequence: ptr char, sequenceLen: int,
    quality: ptr char, qualityLen: int
  ) {.inline.} =
  record.name = name
  record.nameLen = nameLen
  record.comment = comment
  record.commentLen = commentLen
  record.sequence = sequence
  record.sequenceLen = sequenceLen
  record.quality = quality
  record.qualityLen = qualityLen

proc setRecordPtrFields(record: var FQRecordPtr, rec: ptr kseq_t) {.inline.} =
  setRecordPtrFields(
    record,
    rec.name.s, ptrFieldLen(rec.name),
    rec.comment.s, ptrFieldLen(rec.comment),
    rec.seq.s, ptrFieldLen(rec.seq),
    rec.qual.s, ptrFieldLen(rec.qual)
  )

proc copyScratchField(dst: var string, src: ptr char, len: int): ptr char =
  if src.isNil:
    dst.setLen(0)
    return nil

  dst.setLen(len + 1)
  if len > 0:
    copyMem(addr dst[0], src, len)
  dst[len] = '\0'
  addr dst[0]

proc requireFastqRecord(rec: ptr kseq_t, path: string, pairNumber: int, mate: int) =
  if rec.qual.s.isNil:
    raise newException(
      ValueError,
      "readFQInterleavedPairPtr requires FASTQ input; missing qualities at pair " &
      $pairNumber & ", mate " & $mate & " in " & path
    )

## Iterator for reading FASTQ files, returning pointers to record data
##
## Note: The pointers are reused between iterations, so don't store them.
## For stdin input, use "-" as the path parameter.
## Cached lengths are available on each `FQRecordPtr` field and exclude the
## trailing NUL terminator.
##
## Args:
##   path: Path to the FASTQ file (supports gzipped files)
##
## Returns:
##   An iterator yielding FQRecordPtr objects
##
## Raises:
##   IOError: If the input stream cannot be opened
##
## Example:
##
## ```nim
## for rec in readFQPtr("sample.fastq.gz"):
##   echo $rec.name
##   echo $rec.sequence
## ```
iterator readFQPtr*(path: string): FQRecordPtr =
  # - ptr char will be reused on next iteration
  # - for stdin use "-" as path
  # - gz[d]open default even for flat file format
  var result: FQRecordPtr# 'result' not implicit in iterators
  let fp = openGzForRead(path)
  let rec = kseq_init(fp)
  while true:
    if kseq_read(rec) < 0:
      break
    result.setRecordPtrFields(rec)
    yield result
  discard gzclose(fp)

## Iterator for reading FASTQ files, returning copies of record data
##
## This iterator creates copies of the strings, unlike readFQPtr which
## returns pointers to the underlying data.
##
## Args:
##   path: Path to the FASTQ file (supports gzipped files)
##
## Returns:
##   An iterator yielding FQRecord objects with copied data
##
## Raises:
##   IOError: Propagated from `readFQPtr` if the input stream cannot be opened
##
## Example:
##
## ```nim
## for rec in readFQ("sample.fastq.gz"):
##   echo rec.name
##   echo rec.sequence
## ```
iterator readFQ*(path: string): FQRecord =
  var result: FQRecord  # 'result' not implicit in iterators
  for rec in readFQPtr(path):
    result.name = cstrOrEmpty(rec.name)
    result.comment = cstrOrEmpty(rec.comment)
    result.sequence = cstrOrEmpty(rec.sequence)
    result.quality = cstrOrEmpty(rec.quality)
    yield result

## Iterator for reading paired-end FASTQ files synchronously with pointers
##
## Reads two FASTQ files in parallel, yielding pairs of corresponding records.
## Pointer fields are reused between iterations; convert to strings if data
## must be retained after the next `yield`.
## Cached lengths are available on both `read1` and `read2`.
##
## Args:
##   path1: Path to the first FASTQ file (R1, forward reads)
##   path2: Path to the second FASTQ file (R2, reverse reads)
##   checkNames: Whether to verify that read names match (default: false)
##
## Returns:
##   An iterator yielding FQPairPtr objects with synchronized reads
##
## Raises:
##   IOError: If files cannot be opened or have mismatched lengths
##   ValueError: If checkNames is true and read names don't match
iterator readFQPairPtr*(path1: string, path2: string, checkNames: bool = false): FQPairPtr =
  var fp1, fp2: GzFile

  fp1 = openGzForRead(path1)

  if path2 == "-":
    discard gzclose(fp1)
    raise newException(IOError, "Cannot use stdin for both paired files")
  else:
    fp2 = gzopen(path2, "r")
  if fp2 == nil:
    discard gzclose(fp1)
    raise newException(IOError, "Cannot open file: " & path2)

  let rec1 = kseq_init(fp1)
  let rec2 = kseq_init(fp2)
  var pair: FQPairPtr
  var count = 0

  try:
    while true:
      let ret1 = kseq_read(rec1)
      let ret2 = kseq_read(rec2)

      if ret1 < 0 and ret2 < 0:
        break

      if ret1 < 0:
        raise newException(IOError, "File " & path1 & " ended prematurely after " & $count & " sequences")
      if ret2 < 0:
        raise newException(IOError, "File " & path2 & " ended prematurely after " & $count & " sequences")

      count += 1

      pair.read1.setRecordPtrFields(rec1)
      pair.read2.setRecordPtrFields(rec2)

      if checkNames:
        let rawName1 = cstrOrEmpty(pair.read1.name)
        let rawName2 = cstrOrEmpty(pair.read2.name)
        let name1 = normalizePairName(rawName1, 1)
        let name2 = normalizePairName(rawName2, 2)
        if name1 != name2:
          raise newException(ValueError, "Sequence name mismatch at record " & $count & ": '" &
                           rawName1 & "' != '" & rawName2 & "'")

      yield pair

  finally:
    discard gzclose(fp1)
    discard gzclose(fp2)

## Iterator for reading interleaved paired-end FASTQ files with pointers.
##
## Reads one interleaved FASTQ stream and yields each adjacent R1/R2 record pair
## as `FQPairPtr`. `read1` is copied into scratch storage so that both records
## remain valid until the next `yield`. Scratch-backed pointers are always
## NUL-terminated.
##
## Args:
##   path: Path to the interleaved FASTQ file (supports gzipped files; `"-"`
##     reads from stdin)
##   checkNames: Whether to verify that read names match after mate suffix
##     normalization (default: false)
##
## Returns:
##   An iterator yielding `FQPairPtr` objects with cached field lengths
##
## Raises:
##   IOError: If the input stream cannot be opened or ends with an incomplete pair
##   ValueError: If input is not FASTQ or `checkNames` detects a mismatch
iterator readFQInterleavedPairPtr*(path: string, checkNames: bool = false): FQPairPtr =
  let fp = openGzForRead(path)
  let rec = kseq_init(fp)
  var pair: FQPairPtr
  var pairCount = 0

  var read1NameScratch = ""
  var read1CommentScratch = ""
  var read1SequenceScratch = ""
  var read1QualityScratch = ""

  try:
    while true:
      let ret1 = kseq_read(rec)
      if ret1 < 0:
        break

      let pairNumber = pairCount + 1
      requireFastqRecord(rec, path, pairNumber, 1)

      let read1Name = copyScratchField(read1NameScratch, rec.name.s, ptrFieldLen(rec.name))
      let read1Comment = copyScratchField(read1CommentScratch, rec.comment.s, ptrFieldLen(rec.comment))
      let read1Sequence = copyScratchField(read1SequenceScratch, rec.seq.s, ptrFieldLen(rec.seq))
      let read1Quality = copyScratchField(read1QualityScratch, rec.qual.s, ptrFieldLen(rec.qual))
      pair.read1.setRecordPtrFields(
        read1Name, ptrFieldLen(rec.name),
        read1Comment, ptrFieldLen(rec.comment),
        read1Sequence, ptrFieldLen(rec.seq),
        read1Quality, ptrFieldLen(rec.qual)
      )

      let ret2 = kseq_read(rec)
      if ret2 < 0:
        raise newException(IOError, "Interleaved file " & path &
          " ended prematurely after " & $pairCount & " complete pairs")

      requireFastqRecord(rec, path, pairNumber, 2)
      pair.read2.setRecordPtrFields(rec)

      if checkNames:
        let rawName1 = cstrOrEmpty(pair.read1.name)
        let rawName2 = cstrOrEmpty(pair.read2.name)
        let name1 = normalizePairName(rawName1, 1)
        let name2 = normalizePairName(rawName2, 2)
        if name1 != name2:
          raise newException(ValueError, "Sequence name mismatch at record " & $pairNumber & ": '" &
                           rawName1 & "' != '" & rawName2 & "'")

      pairCount = pairNumber
      yield pair

  finally:
    discard gzclose(fp)

## Iterator for reading paired-end FASTQ files synchronously
##
## Reads two FASTQ files in parallel, yielding pairs of corresponding records.
## The files must have the same number of sequences in the same order.
##
## Args:
##   path1: Path to the first FASTQ file (R1, forward reads)
##   path2: Path to the second FASTQ file (R2, reverse reads)
##   checkNames: Whether to verify that read names match (default: false)
##
## Returns:
##   An iterator yielding FQPair objects with synchronized reads
##
## Raises:
##   IOError: If files cannot be opened or have mismatched lengths
##   ValueError: If checkNames is true and read names don't match
##
## Example:
##
## ```nim
## for pair in readFQPair("sample_R1.fastq.gz", "sample_R2.fastq.gz"):
##   echo "Forward: ", pair.read1.name
##   echo "Reverse: ", pair.read2.name
##   processReadPair(pair.read1.sequence, pair.read2.sequence)
## ```
iterator readFQPair*(path1: string, path2: string, checkNames: bool = false): FQPair =
  var pair: FQPair
  for rec in readFQPairPtr(path1, path2, checkNames = checkNames):
    pair.read1.name = cstrOrEmpty(rec.read1.name)
    pair.read1.comment = cstrOrEmpty(rec.read1.comment)
    pair.read1.sequence = cstrOrEmpty(rec.read1.sequence)
    pair.read1.quality = cstrOrEmpty(rec.read1.quality)

    pair.read2.name = cstrOrEmpty(rec.read2.name)
    pair.read2.comment = cstrOrEmpty(rec.read2.comment)
    pair.read2.sequence = cstrOrEmpty(rec.read2.sequence)
    pair.read2.quality = cstrOrEmpty(rec.read2.quality)
    yield pair

## Formats a sequence record as a FASTA or FASTQ string
##
## Args:
##   name: Sequence name/identifier
##   comment: Sequence comment (optional)
##   sequence: The sequence string
##   quality: Quality scores (empty for FASTA format)
##
## Returns:
##   Formatted FASTA/FASTQ string
proc fqfmt(name: string, comment: string, sequence: string, quality: string): string =
  var fastq = false
  var header = ">"
  if len(sequence) == 0:
    return ""
  if len(quality) > 0:
    fastq = true
    header = "@"
  result = header & name
  if comment != "":
    result = result & " " & comment
  result = result & "\n" & sequence
  if fastq:
    result = result & "\n+\n" & quality

## Print FASTA record splitting the sequence into lines of N characters
## Returns:
##  Formatted FASTA string
proc fafmt(name: string, comment: string, sequence: string, width: int = 60): string =
  result = ">" & name
  if comment != "":
    result = result & " " & comment
  result = result & "\n"
  for i in 0 ..< sequence.len:
    if i > 0 and i mod width == 0:
      result = result & "\n"
    result = result & sequence[i]
  return result

## Print FASTA record splitting the sequence into lines of N characters
## Returns:
## Formatted FASTA string
proc fafmt*(rec: FQRecord, width: int = 60): string =
  return fafmt(rec.name, rec.comment, rec.sequence, width)

# Convert a FQRecord to a string (FASTA or FASTQ format)
##
## Returns:
##   Formatted FASTA/FASTQ string
proc `$`*(rec: FQRecord): string =
  return fqfmt(rec.name, rec.comment, rec.sequence, rec.quality)

## Convert a FQRecordPtr to a string (FASTA or FASTQ format)
##
## Returns:
##   Formatted FASTA/FASTQ string
proc `$`*(rec: FQRecordPtr): string =
  let nameStr = cstrOrEmpty(rec.name)
  let commentStr = cstrOrEmpty(rec.comment)
  let sequenceStr = cstrOrEmpty(rec.sequence)
  let qualityStr = cstrOrEmpty(rec.quality)

  return fqfmt(nameStr, commentStr, sequenceStr, qualityStr)
