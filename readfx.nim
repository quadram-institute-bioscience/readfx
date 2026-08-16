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
    curr_char: cint
    f: ptr kstream_t
  KseqGzFile = pointer
  KseqGzHandleObj = object
    file: GzFile
    path: string
    lastError: string
    closed: bool
  KseqGzHandle = ref KseqGzHandleObj

# ------------------------------------------------------------------
# Private gzfast-backed input shim for the kseq C wrapper.
# The native Nim parser (`readFastx`, `Bufio`, `Interval`) and its own
# GzFile handling live in readfx/nimklib.nim and are re-exported above.
# ------------------------------------------------------------------

proc asKseqFile(handle: KseqGzHandle): KseqGzFile {.inline.} =
  cast[KseqGzFile](handle)

proc readfx_gzfast_read(fp: KseqGzFile; buf: pointer; len: cint): cint
    {.cdecl, exportc.} =
  let handle = cast[ptr KseqGzHandleObj](fp)
  if handle.isNil or handle.closed:
    return -1
  if len <= 0:
    return 0
  try:
    let n = handle.file.readRaw(buf, int(len))
    if n < 0:
      handle.lastError = "error reading " & handle.path
      return -1
    result = cint(n)
  except CatchableError as error:
    handle.lastError = error.msg
    result = -1

## Initialize a kseq parser handle from an open gzfast-backed input stream.
##
## Args:
##   fp: Open gzip/plain-text stream handle adapted for kseq
##
## Returns:
##   Pointer to an initialized parser state
proc kseq_init*(fp: KseqGzFile): ptr kseq_t {.header: kseqh, importc: "kseq_init".}


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


## Destroy a kseq parser handle and release its internal buffers.
##
## Args:
##   seq: Parser state previously created with `kseq_init`
proc kseq_destroy*(seq: ptr kseq_t) {.header: kseqh, importc: "kseq_destroy".}

proc openGzForRead(path: string): KseqGzHandle =
  result = KseqGzHandle(path: path)
  try:
    result.file.open(path)
  except CatchableError as error:
    raise newException(IOError, "Cannot open file: " & path & ": " & error.msg)

proc closeGzForRead(handle: KseqGzHandle): int {.discardable.} =
  if handle.isNil or handle.closed:
    return 0
  result = handle.file.close()
  handle.closed = true

proc cstrOrEmpty(p: ptr char): string {.inline.} =
  if p.isNil:
    ""
  else:
    $cast[cstring](p)

proc strFromPtr(p: ptr char, len: int): string {.inline.} =
  ## Build a Nim string from a kseq field pointer and its cached length,
  ## avoiding the strlen rescan done by cstrOrEmpty. Nil pointers map to "".
  if p.isNil or len <= 0:
    return ""
  result = newString(len)
  copyMem(addr result[0], p, len)

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

proc raiseKseqReadError(path: string, status: cint, recordNumber: int;
    handle: KseqGzHandle = nil) {.noReturn.} =
  case status
  of -2:
    raise newException(
      ValueError,
      "Malformed FASTQ in " & path & " at record " & $recordNumber &
      ": truncated or mismatched quality string"
    )
  of -3:
    if not handle.isNil and handle.lastError.len > 0:
      raise newException(
        IOError,
        "Error reading " & path & " at record " & $recordNumber &
        ": " & handle.lastError
      )
    raise newException(
      IOError,
      "Error reading " & path & " at record " & $recordNumber
    )
  else:
    raise newException(
      IOError,
      "kseq_read failed for " & path & " at record " & $recordNumber &
      " with status " & $status
    )

proc requireQualityForFastqHeader(rec: ptr kseq_t, path: string, recordNumber: int) =
  if rec.curr_char == cint('@') and rec.qual.s.isNil:
    raise newException(
      ValueError,
      "Malformed FASTQ in " & path & " at record " & $recordNumber &
      ": '@' record is missing '+' line and qualities"
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
##   ValueError: If an `@` FASTQ record is malformed
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
  # - the gzfast-backed shim accepts both gzip and plain FASTX streams
  var result: FQRecordPtr# 'result' not implicit in iterators
  let fp = openGzForRead(path)
  let rec = kseq_init(fp.asKseqFile())
  var recordNumber = 0
  try:
    while true:
      let ret = kseq_read(rec)
      if ret == -1:
        break
      if ret < -1:
        raiseKseqReadError(path, ret, recordNumber + 1, fp)

      inc recordNumber
      requireQualityForFastqHeader(rec, path, recordNumber)
      result.setRecordPtrFields(rec)
      yield result
  finally:
    kseq_destroy(rec)
    discard closeGzForRead(fp)

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
##   ValueError: Propagated from `readFQPtr` if an `@` FASTQ record is malformed
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
    result.name = strFromPtr(rec.name, rec.nameLen)
    result.comment = strFromPtr(rec.comment, rec.commentLen)
    result.sequence = strFromPtr(rec.sequence, rec.sequenceLen)
    result.quality = strFromPtr(rec.quality, rec.qualityLen)
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
  var fp1, fp2: KseqGzHandle

  fp1 = openGzForRead(path1)

  if path2 == "-":
    discard closeGzForRead(fp1)
    raise newException(IOError, "Cannot use stdin for both paired files")
  else:
    try:
      fp2 = openGzForRead(path2)
    except CatchableError:
      discard closeGzForRead(fp1)
      raise

  let rec1 = kseq_init(fp1.asKseqFile())
  let rec2 = kseq_init(fp2.asKseqFile())
  var pair: FQPairPtr
  var count = 0

  try:
    while true:
      let ret1 = kseq_read(rec1)
      let ret2 = kseq_read(rec2)

      if ret1 < -1:
        raiseKseqReadError(path1, ret1, count + 1, fp1)
      if ret2 < -1:
        raiseKseqReadError(path2, ret2, count + 1, fp2)

      if ret1 == -1 and ret2 == -1:
        break

      if ret1 == -1:
        raise newException(IOError, "File " & path1 & " ended prematurely after " & $count & " sequences")
      if ret2 == -1:
        raise newException(IOError, "File " & path2 & " ended prematurely after " & $count & " sequences")

      count += 1

      requireQualityForFastqHeader(rec1, path1, count)
      requireQualityForFastqHeader(rec2, path2, count)
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
    kseq_destroy(rec1)
    kseq_destroy(rec2)
    discard closeGzForRead(fp1)
    discard closeGzForRead(fp2)

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
  let rec = kseq_init(fp.asKseqFile())
  var pair: FQPairPtr
  var pairCount = 0

  var read1NameScratch = ""
  var read1CommentScratch = ""
  var read1SequenceScratch = ""
  var read1QualityScratch = ""

  try:
    while true:
      let ret1 = kseq_read(rec)
      if ret1 == -1:
        break
      if ret1 < -1:
        raiseKseqReadError(path, ret1, pairCount * 2 + 1, fp)

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
      if ret2 == -1:
        raise newException(IOError, "Interleaved file " & path &
          " ended prematurely after " & $pairCount & " complete pairs")
      if ret2 < -1:
        raiseKseqReadError(path, ret2, pairNumber * 2, fp)

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
    kseq_destroy(rec)
    discard closeGzForRead(fp)

## Iterator for reading interleaved paired-end FASTQ files
##
## Reads one interleaved FASTQ stream and yields each adjacent R1/R2 record
## pair as `FQPair` with copied string data (safe to store after the loop).
## Built on `readFQInterleavedPairPtr`; use that instead when allocation
## overhead matters.
##
## Args:
##   path: Path to the interleaved FASTQ file (supports gzipped files; `"-"`
##     reads from stdin)
##   checkNames: Whether to verify that read names match after mate suffix
##     normalization (default: false)
##
## Returns:
##   An iterator yielding `FQPair` objects with synchronized reads
##
## Raises:
##   IOError: If the input stream cannot be opened or ends with an incomplete pair
##   ValueError: If input is not FASTQ or `checkNames` detects a mismatch
##
## Example:
##
## ```nim
## for pair in readFQInterleavedPair("sample.interleaved.fastq.gz"):
##   echo "R1: ", pair.read1.name
##   echo "R2: ", pair.read2.name
## ```
iterator readFQInterleavedPair*(path: string, checkNames: bool = false): FQPair =
  var pair: FQPair
  for rec in readFQInterleavedPairPtr(path, checkNames = checkNames):
    pair.read1.name = strFromPtr(rec.read1.name, rec.read1.nameLen)
    pair.read1.comment = strFromPtr(rec.read1.comment, rec.read1.commentLen)
    pair.read1.sequence = strFromPtr(rec.read1.sequence, rec.read1.sequenceLen)
    pair.read1.quality = strFromPtr(rec.read1.quality, rec.read1.qualityLen)

    pair.read2.name = strFromPtr(rec.read2.name, rec.read2.nameLen)
    pair.read2.comment = strFromPtr(rec.read2.comment, rec.read2.commentLen)
    pair.read2.sequence = strFromPtr(rec.read2.sequence, rec.read2.sequenceLen)
    pair.read2.quality = strFromPtr(rec.read2.quality, rec.read2.qualityLen)
    yield pair

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
    pair.read1.name = strFromPtr(rec.read1.name, rec.read1.nameLen)
    pair.read1.comment = strFromPtr(rec.read1.comment, rec.read1.commentLen)
    pair.read1.sequence = strFromPtr(rec.read1.sequence, rec.read1.sequenceLen)
    pair.read1.quality = strFromPtr(rec.read1.quality, rec.read1.qualityLen)

    pair.read2.name = strFromPtr(rec.read2.name, rec.read2.nameLen)
    pair.read2.comment = strFromPtr(rec.read2.comment, rec.read2.commentLen)
    pair.read2.sequence = strFromPtr(rec.read2.sequence, rec.read2.sequenceLen)
    pair.read2.quality = strFromPtr(rec.read2.quality, rec.read2.qualityLen)
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
