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
import zip/zlib
import algorithm
import strutils
import readfx/seqtypes
export seqtypes

import readfx/sequtils
export sequtils

import readfx/oligoutils
export oligoutils
# https://forum.nim-lang.org/t/2668
from os import splitPath
const kseqh = currentSourcePath().splitPath.head & "/readfx/kseq.h"

# https://github.com/nim-lang/nimble/issues/157
{.passL: "-lz".}


type
  kstring_t {.importc, header: kseqh.} = object
    l: int  # Use the C definition as is
    m: int
    s: ptr char
  kstream_t {.importc, header: kseqh.} = object
    begin: int
    endd: int
    is_eof: int
  kseq_t {.importc, header: kseqh.} = object
    name: kstring_t
    comment: kstring_t
    seq: kstring_t
    qual: kstring_t
    last_char: int
    f: ptr kstream_t
  gzFile = pointer

proc kseq_init*(fp: gzFile): ptr kseq_t {.header: kseqh, importc: "kseq_init".}


proc kseq_rewind*(seq: ptr kseq_t) {.header: kseqh, importc: "kseq_rewind".}


proc kseq_read*(seq: ptr kseq_t): int {.header: kseqh, importc: "kseq_read".}

## Iterator for reading FASTQ files, returning pointers to record data
##
## Note: The pointers are reused between iterations, so don't store them.
## For stdin input, use "-" as the path parameter.
##
## Args:
##   path: Path to the FASTQ file (supports gzipped files)
##
## Returns:
##   An iterator yielding FQRecordPtr objects
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
  var fp: GzFile
  if path == "-":
    fp = gzdopen(0, "r")
  else:
    fp = gzopen(path, "r")

  doAssert fp != nil
  let rec = kseq_init(fp)
  while true:
    if kseq_read(rec) < 0:
      break
    result.name = rec.name.s
    result.comment = rec.comment.s
    result.sequence = rec.seq.s
    result.quality = rec.qual.s
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
    # Explicitly convert ptr char to string with proper nil checks
    result.name = if rec.name.isNil: "" else: $cast[cstring](rec.name)
    result.comment = if rec.comment.isNil: "" else: $cast[cstring](rec.comment)
    result.sequence = if rec.sequence.isNil: "" else: $cast[cstring](rec.sequence)
    result.quality = if rec.quality.isNil: "" else: $cast[cstring](rec.quality)
    yield result

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
  var fp1, fp2: GzFile
  
  # Open first file
  if path1 == "-":
    fp1 = gzdopen(0, "r")
  else:
    fp1 = gzopen(path1, "r")
  if fp1 == nil:
    raise newException(IOError, "Cannot open file: " & path1)
  
  # Open second file  
  if path2 == "-":
    raise newException(IOError, "Cannot use stdin for both paired files")
  else:
    fp2 = gzopen(path2, "r")
  if fp2 == nil:
    discard gzclose(fp1)
    raise newException(IOError, "Cannot open file: " & path2)
  
  let rec1 = kseq_init(fp1)
  let rec2 = kseq_init(fp2)
  var pair: FQPair
  var count = 0
  
  try:
    while true:
      let ret1 = kseq_read(rec1)
      let ret2 = kseq_read(rec2)
      
      # Check if both files ended simultaneously
      if ret1 < 0 and ret2 < 0:
        break
      
      # Check for premature end of either file
      if ret1 < 0:
        raise newException(IOError, "File " & path1 & " ended prematurely after " & $count & " sequences")
      if ret2 < 0:
        raise newException(IOError, "File " & path2 & " ended prematurely after " & $count & " sequences")
      
      count += 1
      
      # Convert records to FQRecord format
      pair.read1.name = if rec1.name.s.isNil: "" else: $cast[cstring](rec1.name.s)
      pair.read1.comment = if rec1.comment.s.isNil: "" else: $cast[cstring](rec1.comment.s)
      pair.read1.sequence = if rec1.seq.s.isNil: "" else: $cast[cstring](rec1.seq.s)
      pair.read1.quality = if rec1.qual.s.isNil: "" else: $cast[cstring](rec1.qual.s)
      
      pair.read2.name = if rec2.name.s.isNil: "" else: $cast[cstring](rec2.name.s)
      pair.read2.comment = if rec2.comment.s.isNil: "" else: $cast[cstring](rec2.comment.s)
      pair.read2.sequence = if rec2.seq.s.isNil: "" else: $cast[cstring](rec2.seq.s)
      pair.read2.quality = if rec2.qual.s.isNil: "" else: $cast[cstring](rec2.qual.s)
      
      # Optional name checking
      if checkNames:
        var name1 = pair.read1.name
        var name2 = pair.read2.name
        
        # Remove common suffixes like /1, /2, or trailing spaces
        if name1.endsWith("/1"):
          name1 = name1[0..^3]
        elif name1.endsWith(" 1"):
          name1 = name1[0..^3]
        
        if name2.endsWith("/2"):
          name2 = name2[0..^3]
        elif name2.endsWith(" 2"):
          name2 = name2[0..^3]
        
        if name1 != name2:
          raise newException(ValueError, "Sequence name mismatch at record " & $count & ": '" & 
                           pair.read1.name & "' != '" & pair.read2.name & "'")
      
      yield pair
      
  finally:
    discard gzclose(fp1)
    discard gzclose(fp2)

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
  var result = ">" & name
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
  # Explicitly convert ptr char to string first
  let nameStr = if rec.name.isNil: "" else: $cast[cstring](rec.name)
  let commentStr = if rec.comment.isNil: "" else: $cast[cstring](rec.comment)
  let sequenceStr = if rec.sequence.isNil: "" else: $cast[cstring](rec.sequence)
  let qualityStr = if rec.quality.isNil: "" else: $cast[cstring](rec.quality)
  
  return fqfmt(nameStr, commentStr, sequenceStr, qualityStr)



# -----------------------
# gzip file I/O #
# -----------------------

when defined(windows):
  const libz = "zlib1.dll"
elif defined(macosx):
  const libz = "libz.dylib"
  const libc = "libc.dylib"
else:
  const libz = "libz.so.1"
  const libc = "libc.so.6"

type
  GzFile* = gzFile
proc gzopen(path: cstring, mode: cstring): gzFile{.cdecl, dynlib: libz,
    importc: "gzopen".}
proc gzdopen(fd: int32, mode: cstring): gzFile{.cdecl, dynlib: libz,
    importc: "gzdopen".}
proc gzread(thefile: gzFile, buf: pointer, length: int): int32{.cdecl,
    dynlib: libz, importc: "gzread".}
proc gzclose(thefile: gzFile): int32{.cdecl, dynlib: libz, importc: "gzclose".}



proc open(f: var GzFile, fn: string,
    mode: FileMode = fmRead): int {.discardable.} =
  assert(mode == fmRead or mode == fmWrite)
  result = 0
  if fn == "-" or fn == "":
    if mode == fmRead: f = gzdopen(0, cstring("r"))
    elif mode == fmWrite: f = gzdopen(1, cstring("w"))
  else:
    if mode == fmRead: f = gzopen(cstring(fn), cstring("r"))
    elif mode == fmWrite: f = gzopen(cstring(fn), cstring("w"))
  if f == nil:
    result = -1
    raise newException(IOError, "error opening " & fn)

proc close(f: var GzFile): int {.discardable.} =
  if f != nil:
    result = int(gzclose(f))
    f = nil
  else: result = 0

proc read(f: var GzFile, buf: var string, sz: int, offset: int = 0):
    int {.discardable.} =
  if buf.len < offset + sz: buf.setLen(offset + sz)
  result = gzread(f, buf[offset].addr, buf.len)
  buf.setLen(result)

# -----------------------
# Buffered reader #
# -----------------------

type
  Bufio*[T] = tuple[fp: T, buf: string, st, en, sz: int, EOF: bool]

proc open*[T](f: var Bufio[T], fn: string, mode: FileMode = fmRead,
    sz: int = 0x10000): int {.discardable.} =
  assert(mode == fmRead) # only fmRead is supported for now
  result = f.fp.open(fn, mode)
  (f.st, f.en, f.sz, f.EOF) = (0, 0, sz, false)
  f.buf.setLen(sz)

proc xopen*[T](fn: string, mode: FileMode = fmRead,
    sz: int = 0x10000): Bufio[T] =
  var f: Bufio[T]
  f.open(fn, mode, sz)
  return f

proc close*[T](f: var Bufio[T]): int {.discardable.} =
  return f.fp.close()

proc eof*[T](f: Bufio[T]): bool {.noSideEffect.} =
  result = (f.EOF and f.st >= f.en)

proc readByte*[T](f: var Bufio[T]): int =
  if f.EOF and f.st >= f.en: return -1
  if f.st >= f.en:
    (f.st, f.en) = (0, f.fp.read(f.buf, f.sz))
    if f.en == 0: f.EOF = true; return -1
    if f.en < 0: f.EOF = true; return -2
  result = int(f.buf[f.st])
  f.st += 1

proc read*[T](f: var Bufio[T], buf: var string, sz: int,
    offset: int = 0): int {.discardable.} =
  if f.EOF and f.st >= f.en: return 0
  buf.setLen(offset)
  var off = offset
  var rest = sz
  while rest > f.en - f.st:
    if f.en > f.st:
      let l = f.en - f.st
      if buf.len < off + l: buf.setLen(off + l)
      copyMem(buf[off].addr, f.buf[f.st].addr, l)
      rest -= l
      off += l
    (f.st, f.en) = (0, f.fp.read(f.buf, f.sz))
    if f.en < f.sz: f.EOF = true
    if f.en == 0: return off - offset
  if buf.len < off + rest: buf.setLen(off + rest)
  copyMem(buf[off].addr, f.buf[f.st].addr, rest)
  f.st += rest
  return off + rest - offset

proc memchr(buf: pointer, c: cint, sz: csize_t): pointer {.cdecl, dynlib: libc,
    importc: "memchr".}

proc readUntil*[T](f: var Bufio[T], buf: var string, dret: var char,
    delim: int = -1, offset: int = 0): int {.discardable.} =
  if f.EOF and f.st >= f.en: return -1
  buf.setLen(offset)
  var off = offset
  var gotany = false
  while true:
    if f.en < 0: return -3
    if f.st >= f.en: # buffer is empty
      if not f.EOF:
        (f.st, f.en) = (0, f.fp.read(f.buf, f.sz))
        if f.en < f.sz: f.EOF = true
        if f.en == 0: break
        if f.en < 0:
          f.EOF = true
          return -2
      else: break
    var x: int = f.en
    if delim == -1: # read a line
      #for i in f.st..<f.en:
      #  if f.buf[i] == '\n': x = i; break
      var p = memchr(f.buf[f.st].addr, cint(0xa), csize_t(f.en - f.st))
      if p != nil: x = cast[int](p) - cast[int](f.buf[0].addr)
    elif delim == -2: # read a field
      for i in f.st..<f.en:
        if f.buf[i] == '\t' or f.buf[i] == ' ' or f.buf[i] == '\n':
          x = i; break
    else: # read to other delimitors
      #for i in f.st..<f.en:
      #  if f.buf[i] == char(delim): x = i; break
      var p = memchr(f.buf[f.st].addr, cint(delim), csize_t(f.en - f.st))
      if p != nil: x = cast[int](p) - cast[int](f.buf[0].addr)
    gotany = true
    if x > f.st: # something to write to buf[]
      let l = x - f.st
      if buf.len < off + l: buf.setLen(off + l)
      copyMem(buf[off].addr, f.buf[f.st].addr, l)
      off += l
    f.st = x + 1
    if x < f.en: dret = f.buf[x]; break
  if not gotany and f.eof(): return -1
  if delim == -1 and off > 0 and buf[off - 1] == '\r':
    off -= 1
    buf.setLen(off)
  return off - offset

proc readLine*[T](f: var Bufio[T], buf: var string): bool {.discardable.} =
  var dret: char
  var ret = readUntil(f, buf, dret)
  return if ret >= 0: true else: false

# -----------------------
# Fastx Reader #
# -----------------------


proc readFastx*[T](f: var Bufio[T], r: var FQRecord): bool {.discardable.} =
  var x: int
  var c: char
  if r.lastChar == 0: # the header character hasn't been read yet
    while true:       # look for the header character '>' or '@'
      x = f.readByte()
      if x < 0 or x == int('>') or x == int('@'): break
    if x < 0: r.status = x; return false # end-of-file or stream error
    r.lastChar = x
  r.sequence.setLen(0); r.quality.setLen(0); r.comment.setLen(0)
  x = f.readUntil(r.name, c, -2)
  if x < 0: r.status = x; return false       # EOF or stream error
  if c != '\n': f.readUntil(r.comment, c) # read FASTA/Q comment
  while true:         # read sequence
    x = f.readByte()  # read the first char on a line
    if x < 0 or x == int('>') or x == int('+') or x == int('@'): break
    if x == int('\n'): continue
    r.sequence.add(char(x))
    f.readUntil(r.sequence, c, -1, r.sequence.len)  # read the rest of the seq line
  r.status = r.sequence.len   # for normal records, this keeps the sequence length
  if x == int('>') or x == int('@'): r.lastChar = x
  if x != int('+'): return true
  while true:         # skip the rest of the "+" line
    x = f.readByte()
    if x < 0 or x == int('\n'): break
  if x < 0: r.status = x; return false  # error: no quality
  while true:         # read quality
    x = f.readUntil(r.quality, c, -1, r.quality.len)
    if x < 0 or r.quality.len >= r.sequence.len: break
  if x == -3: r.status = -3; return false # other stream error
  r.lastChar = 0
  if r.sequence.len != r.quality.len: r.status = -4; return false
  return true

# -----------------------
# Intervals #
# -----------------------

type
  Interval*[S,T] = tuple[st, en: S, data: T, max: S]

proc sort*[S,T](a: var seq[Interval[S,T]]) =
  a.sort do (x, y: Interval[S,T]) -> int:
    if x.st < y.st: -1
    elif x.st > y.st: 1
    else: 0

proc index*[S,T](a: var seq[Interval[S,T]]): int {.discardable.} =
  if a.len == 0: return 0
  var is_srt = true
  for i in 1..<a.len:
    if a[i-1].st > a[i].st:
      is_srt = false; break
  if not is_srt: a.sort()
  var last_i: int
  var last: S
  for i in countup(0, a.len-1, 2): # leaves (i.e. at level 0)
    (last_i, last, a[i].max) = (i, a[i].en, a[i].en)
  var k = 1
  while 1 shl k <= a.len: # process internal nodes in the bottom-up order
    let x = 1 shl (k - 1)
    let i0 = (x shl 1) - 1 # the first node at level k
    let step = x shl 2
    for i in countup(i0, a.len - 1, step): # traverse nodes at level k
      let el = a[i - x].max  # max value of the left child
      let er = if i + x < a.len: a[i + x].max else: last # of the right child
      var e = a[i].en
      if e < el: e = el
      if e < er: e = er
      a[i].max = e
    # point last_i to the parent of the original last_i
    last_i = if ((last_i shr k) and 1) != 0: last_i - x else: last_i + x
    if last_i < a.len and a[last_i].max > last: # update last accordingly
      last = a[last_i].max
    k += 1
  return k - 1

iterator overlap*[S,T](a: seq[Interval[S,T]], st: S, en: S): Interval[S,T] {.noSideEffect.} =
  var h: int = 0
  while 1 shl h <= a.len: h += 1
  h -= 1 # h is the height of the tree
  var stack: array[64, tuple[k, x, w:int]] # 64 is the max possible tree height
  var t: int = 0
  stack[t] = (h, (1 shl h) - 1, 0); t += 1 # push the root
  while t > 0: # the following guarantees sorted "yield"
    t -= 1
    let (k, x, w) = stack[t] # pop from the stack
    if k <= 3: # in a small subtree, traverse everything
      let i0 = (x shr k) shl k
      var i1 = i0 + (1 shl (k + 1)) - 1
      if i1 >= a.len: i1 = a.len
      for i in countup(i0, i1 - 1):
        if a[i].st >= en: break  # out of range; no need to proceed
        if st < a[i].en: yield a[i] # overlap! yield
    elif w == 0: # the left child not processed
      let y = x - (1 shl (k - 1)) # the left child of z.x; y may >=a.len
      stack[t] = (k, x, 1); t += 1
      if y >= a.len or a[y].max > st:
        stack[t] = (k-1, y, 0); t += 1 # add left child
    elif x < a.len and a[x].st < en: # need to push the right child
      if st < a[x].en: yield a[x] # test if x overlaps the query
      stack[t] = (k - 1, x + (1 shl (k - 1)), 0); t += 1