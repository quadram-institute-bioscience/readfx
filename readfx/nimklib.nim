## Native Nim FASTA/FASTQ parser and buffered I/O
## ==========================================================
## Based on the "BioFast" repository by Heng Li (https://github.com/lh3/biofast)
##
## This module is the canonical home of:
## * `GzFile`: gzip/plain input handle used for transparent FASTX I/O
## * `Bufio[T]`: generic buffered reader
## * `readFastx`: native Nim FASTA/FASTQ record parser
## * `Interval[S,T]`: implicit interval tree for genomic overlap queries
##
## It is imported and re-exported by `readfx.nim`; prefer `import readfx`
## over importing this module directly.

import os, algorithm
import std/[streams, strutils]
import gzfast except open
import seqtypes

const klibVer* = "0.2"

################
# FASTX file I/O #
################

when defined(windows):
  const libc = "msvcrt.dll"  # memchr lives in the MS C runtime
elif defined(macosx):
  const libc = "libc.dylib"
else:
  const libc = "libc.so.6"

const GzipMagic = "\x1f\x8b"

const
  # Compile-time knobs for benchmarking the gzfast-backed native reader.
  readfxGzfastThreads* {.intdefine.}: int = 0
  readfxGzfastDecodedChunkSize* {.intdefine.}: int = 4 * 1024 * 1024
  readfxGzfastInputPageSize* {.intdefine.}: int = 1024 * 1024

type
  GzFileKind = enum
    gfClosed,
    gfPlain,
    gfGzip

  PrefixFileStream = ref object of StreamObj
    file: File
    prefix: string
    prefixPos: int
    closed: bool
    ownsFile: bool

  GzFile* = object
    kind: GzFileKind
    path: string
    plainFile: File
    plainPrefix: string
    plainPrefixPos: int
    ownsPlainFile: bool
    gzipStream: GzFastStream
    gzipSpan: DecodedSpan
    gzipSpanPos: int
    gzipFinished: bool

proc readPrefix(file: File, maxLen: int): string =
  if maxLen <= 0:
    return ""
  result = newString(maxLen)
  let n = file.readBuffer(addr result[0], maxLen)
  if n < 0:
    raise newException(IOError, "error reading input prefix")
  result.setLen(n)

proc hasGzipMagic(prefix: string): bool {.inline.} =
  prefix.len >= 2 and prefix[0] == GzipMagic[0] and prefix[1] == GzipMagic[1]

proc fileHasGzipMagic(path: string): bool =
  var file: File
  if not open(file, path, fmRead):
    raise newException(IOError, "error opening " & path)
  try:
    file.readPrefix(2).hasGzipMagic()
  finally:
    file.close()

proc gzfastConfig(bufferSize: int): GzFastConfig =
  result = defaultGzFastConfig()
  result.threads = readfxGzfastThreads
  result.decodedChunkSize = max(readfxGzfastDecodedChunkSize, max(bufferSize, 4096))
  result.inputPageSize = max(readfxGzfastInputPageSize, 4096)
  if result.maxSpeculativeOutput < result.decodedChunkSize:
    result.maxSpeculativeOutput = result.decodedChunkSize
  result.validate()

proc pfsReadData(s: Stream; buffer: pointer; bufLen: int): int
    {.nimcall, gcsafe, tags: [ReadIOEffect], raises: [IOError].} =
  let ps = PrefixFileStream(s)
  if ps.closed:
    raise newException(IOError, "input stream is closed")
  if bufLen <= 0:
    return 0

  var written = 0
  if ps.prefixPos < ps.prefix.len:
    let n = min(bufLen, ps.prefix.len - ps.prefixPos)
    copyMem(buffer, addr ps.prefix[ps.prefixPos], n)
    ps.prefixPos += n
    written += n

  if written < bufLen:
    let dst = cast[pointer](cast[uint](buffer) + written.uint)
    let n = ps.file.readBuffer(dst, bufLen - written)
    if n < 0:
      raise newException(IOError, "error reading input stream")
    written += n

  result = written

proc pfsReadDataStr(s: Stream; buffer: var string; slice: Slice[int]): int
    {.nimcall, gcsafe, tags: [ReadIOEffect], raises: [IOError].} =
  if slice.b < slice.a:
    return 0
  if buffer.len <= slice.b:
    buffer.setLen(slice.b + 1)
  result = pfsReadData(s, addr buffer[slice.a], slice.b - slice.a + 1)

proc pfsAtEnd(s: Stream): bool
    {.nimcall, gcsafe, tags: [ReadIOEffect], raises: [].} =
  let ps = PrefixFileStream(s)
  if ps.closed:
    return true
  if ps.prefixPos < ps.prefix.len:
    return false
  result = ps.file.endOfFile()

proc pfsClose(s: Stream)
    {.nimcall, gcsafe, tags: [WriteIOEffect], raises: [].} =
  let ps = PrefixFileStream(s)
  if not ps.closed:
    ps.closed = true
    if ps.ownsFile:
      ps.file.close()

proc pfsUnsupported(s: Stream; pos: int)
    {.nimcall, gcsafe, tags: [], raises: [IOError].} =
  raise newException(IOError, "input stream is forward-only")

proc pfsUnsupportedPos(s: Stream): int
    {.nimcall, gcsafe, tags: [], raises: [IOError].} =
  raise newException(IOError, "input stream is forward-only")

proc pfsUnsupportedPeek(s: Stream; buffer: pointer; bufLen: int): int
    {.nimcall, gcsafe, tags: [ReadIOEffect], raises: [IOError].} =
  raise newException(IOError, "input stream does not support peeking")

proc pfsUnsupportedWrite(s: Stream; buffer: pointer; bufLen: int)
    {.nimcall, gcsafe, tags: [WriteIOEffect], raises: [IOError].} =
  raise newException(IOError, "input stream is read-only")

proc pfsFlush(s: Stream)
    {.nimcall, gcsafe, tags: [WriteIOEffect], raises: [].} =
  discard

proc newPrefixFileStream(file: File; prefix: string;
    ownsFile: bool): PrefixFileStream =
  result = PrefixFileStream(file: file, prefix: prefix, ownsFile: ownsFile)
  result.readDataImpl = pfsReadData
  result.readDataStrImpl = pfsReadDataStr
  result.atEndImpl = cast[typeof(result.atEndImpl)](pfsAtEnd)
  result.closeImpl = pfsClose
  result.setPositionImpl = pfsUnsupported
  result.getPositionImpl = pfsUnsupportedPos
  result.peekDataImpl = pfsUnsupportedPeek
  result.writeDataImpl = pfsUnsupportedWrite
  result.flushImpl = pfsFlush

proc openPlain(f: var GzFile; file: File; prefix: string; ownsFile: bool;
    path: string) =
  f.kind = gfPlain
  f.path = path
  f.plainFile = file
  f.plainPrefix = prefix
  f.plainPrefixPos = 0
  f.ownsPlainFile = ownsFile
  f.gzipStream = nil
  f.gzipSpan = DecodedSpan(data: nil, len: 0)
  f.gzipSpanPos = 0
  f.gzipFinished = false

proc openGzip(f: var GzFile; path: string; bufferSize: int) =
  f.kind = gfGzip
  f.path = path
  f.plainFile = nil
  f.plainPrefix.setLen(0)
  f.plainPrefixPos = 0
  f.ownsPlainFile = false
  f.gzipStream = gzfast.open(initGzFastDecoder(gzfastConfig(bufferSize)), path)
  f.gzipSpan = DecodedSpan(data: nil, len: 0)
  f.gzipSpanPos = 0
  f.gzipFinished = false

proc openGzip(f: var GzFile; file: File; prefix: string; ownsFile: bool;
    path: string; bufferSize: int) =
  f.kind = gfGzip
  f.path = path
  f.plainFile = nil
  f.plainPrefix.setLen(0)
  f.plainPrefixPos = 0
  f.ownsPlainFile = false
  f.gzipStream = openGzFastSequential(
    newPrefixFileStream(file, prefix, ownsFile),
    gzfastConfig(bufferSize)
  )
  f.gzipSpan = DecodedSpan(data: nil, len: 0)
  f.gzipSpanPos = 0
  f.gzipFinished = false

proc open(f: var GzFile, fn: string,
    mode: FileMode = fmRead, bufferSize: int = 0x10000): int {.discardable.} =
  assert(mode == fmRead)
  result = 0
  if fn == "-" or fn == "":
    let prefix = stdin.readPrefix(2)
    if prefix.hasGzipMagic():
      f.openGzip(stdin, prefix, ownsFile = false, path = "stdin", bufferSize)
    else:
      f.openPlain(stdin, prefix, ownsFile = false, path = "stdin")
  else:
    let gzipInput = fn.endsWith(".gz") or fn.fileHasGzipMagic()
    if gzipInput:
      f.openGzip(fn, bufferSize)
    else:
      var file: File
      if not open(file, fn, fmRead):
        result = -1
        raise newException(IOError, "error opening " & fn)
      f.openPlain(file, "", ownsFile = true, path = fn)

proc close(f: var GzFile): int {.discardable.} =
  result = 0
  case f.kind
  of gfPlain:
    if f.ownsPlainFile and f.plainFile != nil:
      f.plainFile.close()
  of gfGzip:
    if not f.gzipStream.isNil:
      f.gzipStream.close()
  of gfClosed:
    discard
  f.kind = gfClosed
  f.path.setLen(0)
  f.plainFile = nil
  f.plainPrefix.setLen(0)
  f.plainPrefixPos = 0
  f.ownsPlainFile = false
  f.gzipStream = nil
  f.gzipSpan = DecodedSpan(data: nil, len: 0)
  f.gzipSpanPos = 0
  f.gzipFinished = false

proc read(f: var GzFile, buf: var string, sz: int, offset: int = 0):
    int {.discardable, inline.} =
  if sz <= 0:
    if buf.len < offset:
      buf.setLen(offset)
    return 0
  if buf.len < offset + sz:
    buf.setLen(offset + sz)

  case f.kind
  of gfPlain:
    var written = 0
    if f.plainPrefixPos < f.plainPrefix.len:
      let n = min(sz, f.plainPrefix.len - f.plainPrefixPos)
      copyMem(addr buf[offset], addr f.plainPrefix[f.plainPrefixPos], n)
      f.plainPrefixPos += n
      written += n
    if written < sz:
      let n = f.plainFile.readBuffer(addr buf[offset + written], sz - written)
      if n < 0:
        buf.setLen(offset + written)
        return -1
      written += n
    result = written
    buf.setLen(offset + result)
  of gfGzip:
    result = f.gzipStream.readData(addr buf[offset], sz)
    buf.setLen(offset + result)
    if result == 0 and not f.gzipFinished:
      discard f.gzipStream.finish()
      f.gzipFinished = true
  of gfClosed:
    raise newException(IOError, "input stream is closed")

proc finishGzipAtEnd(f: var GzFile) {.inline.} =
  if not f.gzipFinished:
    discard f.gzipStream.finish()
    f.gzipFinished = true

proc currentGzipSpan(f: var GzFile): DecodedSpan {.inline.} =
  let available = f.gzipSpan.len - f.gzipSpanPos
  if available <= 0:
    return DecodedSpan(data: nil, len: 0)
  DecodedSpan(
    data: cast[ptr UncheckedArray[byte]](
      cast[uint](f.gzipSpan.data) + uint(f.gzipSpanPos)
    ),
    len: available
  )

proc clearGzipSpan(f: var GzFile) {.inline.} =
  f.gzipSpan = DecodedSpan(data: nil, len: 0)
  f.gzipSpanPos = 0

proc peekGzipDecoded(f: var GzFile): DecodedSpan {.inline.} =
  if f.kind != gfGzip or f.gzipStream.isNil:
    raise newException(IOError, "gzip stream is not open")
  result = f.currentGzipSpan()
  if result.len > 0:
    return
  if f.gzipFinished:
    return DecodedSpan(data: nil, len: 0)
  f.gzipSpan = f.gzipStream.peekDecoded()
  f.gzipSpanPos = 0
  result = f.currentGzipSpan()
  if result.len == 0:
    f.clearGzipSpan()
    f.finishGzipAtEnd()

proc consumeGzipDecoded(f: var GzFile; n: int) {.inline.} =
  let available = f.gzipSpan.len - f.gzipSpanPos
  if n < 0 or n > available:
    raise newException(ValueError,
      "consumeDecoded count exceeds available decoded bytes")
  f.gzipSpanPos += n
  if f.gzipSpanPos == f.gzipSpan.len and f.gzipSpan.len > 0:
    f.gzipStream.consumeDecoded(f.gzipSpan.len)
    f.clearGzipSpan()

###################
# Buffered reader #
###################

type
  ## Generic buffered reader state used by `readFastx` and helper readers.
  ##
  ## Fields:
  ##   fp: Underlying file/stream handle
  ##   buf: Internal byte buffer
  ##   st: Current read index inside buffer
  ##   en: End index of valid data in buffer
  ##   sz: Buffer size in bytes
  ##   EOF: EOF flag for the underlying stream
  Bufio*[T] = tuple[fp: T, buf: string, st, en, sz: int, EOF: bool]

## Open a buffered reader over a gzip/plain file handle.
##
## Args:
##   f: Buffer object to initialize
##   fn: Input filename (`"-"` means stdin)
##   mode: File mode (only `fmRead` is supported)
##   sz: Internal buffer size in bytes
##
## Returns:
##   `0` on success; may raise on open failure
proc open*[T](f: var Bufio[T], fn: string, mode: FileMode = fmRead,
    sz: int = 0x10000): int {.discardable.} =
  assert(mode == fmRead) # only fmRead is supported for now
  when T is GzFile:
    result = f.fp.open(fn, mode, sz)
  else:
    result = f.fp.open(fn, mode)
  (f.st, f.en, f.sz, f.EOF) = (0, 0, sz, false)
  f.buf.setLen(sz)

## Create and open a buffered reader in one call.
##
## Args:
##   fn: Input filename (`"-"` means stdin)
##   mode: File mode (only `fmRead` is supported)
##   sz: Internal buffer size in bytes
##
## Returns:
##   Initialized `Bufio[T]` instance
proc xopen*[T](fn: string, mode: FileMode = fmRead,
    sz: int = 0x10000): Bufio[T] =
  var f: Bufio[T]
  f.open(fn, mode, sz)
  return f

## Close a buffered reader and its underlying file handle.
##
## Args:
##   f: Buffered reader to close
##
## Returns:
##   Close status from the underlying file handle
proc close*[T](f: var Bufio[T]): int {.discardable.} =
  return f.fp.close()

proc refill[T](f: var Bufio[T]): int {.inline.} =
  f.st = 0
  f.en = f.fp.read(f.buf, f.sz)
  if f.en <= 0:
    f.EOF = true
  result = f.en

## Report whether buffered reader has reached end-of-file.
##
## Args:
##   f: Buffered reader
##
## Returns:
##   `true` when no more bytes can be read
proc eof*[T](f: Bufio[T]): bool {.noSideEffect, inline.} =
  result = (f.EOF and f.st >= f.en)

## Read a single byte from the buffered reader.
##
## Args:
##   f: Buffered reader
##
## Returns:
##   Byte value (0..255) on success, or:
##   - `-1`: EOF
##   - `-2`: stream read error
proc readByte*[T](f: var Bufio[T]): int {.inline.} =
  when T is GzFile:
    if f.fp.kind == gfGzip:
      if f.EOF: return -1
      let span = f.fp.peekGzipDecoded()
      if span.len == 0:
        f.EOF = true
        return -1
      result = int(span.data[0])
      f.fp.consumeGzipDecoded(1)
      return

  if f.EOF and f.st >= f.en: return -1
  if f.st >= f.en:
    let n = f.refill()
    if n == 0: return -1
    if n < 0: return -2
  result = int(f.buf[f.st])
  f.st += 1

## Read up to `sz` bytes from the buffered reader into `buf`.
##
## Args:
##   f: Buffered reader
##   buf: Destination buffer
##   sz: Number of bytes to read
##   offset: Write position inside `buf`
##
## Returns:
##   Number of bytes written into `buf`
proc read*[T](f: var Bufio[T], buf: var string, sz: int,
    offset: int = 0): int {.discardable.} =
  when T is GzFile:
    if f.fp.kind == gfGzip:
      if f.EOF:
        return 0
      buf.setLen(offset)
      var off = offset
      var rest = sz
      while rest > 0:
        let span = f.fp.peekGzipDecoded()
        if span.len == 0:
          f.EOF = true
          break
        let n = min(rest, span.len)
        if buf.len < off + n: buf.setLen(off + n)
        copyMem(addr buf[off], span.data, n)
        f.fp.consumeGzipDecoded(n)
        off += n
        rest -= n
      return off - offset

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
    let n = f.refill()
    if n == 0: return off - offset
    if n < 0:
      if off > offset:
        return off - offset
      return -2
  if buf.len < off + rest: buf.setLen(off + rest)
  copyMem(buf[off].addr, f.buf[f.st].addr, rest)
  f.st += rest
  return off + rest - offset

proc memchr(buf: pointer, c: cint, sz: csize_t): pointer {.cdecl, dynlib: libc,
    importc: "memchr".}

proc findDelimiter(span: DecodedSpan; delim: int): int {.inline.} =
  result = span.len
  if span.len <= 0:
    return
  if delim == -1:
    let p = memchr(span.data, cint(0x0a), csize_t(span.len))
    if p != nil:
      result = int(cast[uint](p) - cast[uint](span.data))
  elif delim == -2:
    for i in 0..<span.len:
      let b = span.data[i]
      if b == 0x09'u8 or b == 0x20'u8 or b == 0x0a'u8 or b == 0x0d'u8:
        return i
  else:
    let p = memchr(span.data, cint(delim), csize_t(span.len))
    if p != nil:
      result = int(cast[uint](p) - cast[uint](span.data))

proc appendSpan(buf: var string; off: var int; span: DecodedSpan; n: int)
    {.inline.} =
  if n <= 0:
    return
  if buf.len < off + n:
    buf.setLen(off + n)
  copyMem(addr buf[off], span.data, n)
  off += n

proc readUntilGzip(f: var Bufio[GzFile], buf: var string, dret: var char;
    delim: int = -1, offset: int = 0): int {.discardable, inline.} =
  if f.EOF:
    return -1
  buf.setLen(offset)
  var off = offset
  var gotany = false
  while true:
    let span = f.fp.peekGzipDecoded()
    if span.len == 0:
      f.EOF = true
      break
    let x = findDelimiter(span, delim)
    gotany = true
    appendSpan(buf, off, span, x)
    if x < span.len:
      dret = char(span.data[x])
      f.fp.consumeGzipDecoded(x + 1)
      break
    f.fp.consumeGzipDecoded(span.len)
  if not gotany and f.EOF:
    return -1
  if delim == -1 and off > 0 and buf[off - 1] == '\r':
    off -= 1
    buf.setLen(off)
  return off - offset

## Read from buffered stream until a delimiter or EOF.
##
## Args:
##   f: Buffered reader
##   buf: Destination buffer
##   dret: Delimiter character found
##   delim: Delimiter mode:
##     - `-1`: read line
##     - `-2`: read field (space/tab/newline)
##     - other: read until the given byte value
##   offset: Write position inside `buf`
##
## Returns:
##   Number of bytes appended, or negative status code:
##   - `-1`: EOF
##   - `-2`: stream read error
##   - `-3`: internal buffered-state error
proc readUntil*[T](f: var Bufio[T], buf: var string, dret: var char,
    delim: int = -1, offset: int = 0): int {.discardable, inline.} =
  when T is GzFile:
    if f.fp.kind == gfGzip:
      return readUntilGzip(f, buf, dret, delim, offset)

  if f.EOF and f.st >= f.en: return -1
  buf.setLen(offset)
  var off = offset
  var gotany = false
  while true:
    if f.en < 0: return -3
    if f.st >= f.en: # buffer is empty
      if not f.EOF:
        let n = f.refill()
        if n == 0: break
        if n < 0: return -2
      else: break
    var x: int = f.en
    if delim == -1: # read a line
      var p = memchr(f.buf[f.st].addr, cint(0xa), csize_t(f.en - f.st))
      if p != nil: x = cast[int](p) - cast[int](f.buf[0].addr)
    elif delim == -2: # read a field
      # CR is a delimiter too: matches kseq's isspace() behavior and
      # keeps CRLF names clean ("@r1\r\n" -> "r1", not "r1\r")
      for i in f.st..<f.en:
        if f.buf[i] == '\t' or f.buf[i] == ' ' or f.buf[i] == '\n' or
            f.buf[i] == '\r':
          x = i; break
    else: # read to other delimitors
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

## Read one line from buffered stream.
##
## Args:
##   f: Buffered reader
##   buf: Line destination buffer
##
## Returns:
##   `true` if a line was read, `false` on EOF/error
proc readLine*[T](f: var Bufio[T], buf: var string): bool {.discardable, inline.} =
  var dret: char
  var ret = readUntil(f, buf, dret)
  return if ret >= 0: true else: false

################
# Fastx Reader #
################

## Read one FASTA/FASTQ record from a buffered stream.
##
## Args:
##   f: Buffered reader
##   r: Record object reused for output
##
## Returns:
##   `true` if one record was parsed, `false` on EOF/error
##
## Notes:
##   On failure, `r.status` is set to a negative code:
##   - `-1`: EOF
##   - `-2`: stream read error
##   - `-3`: parser/stream state error
##   - `-4`: malformed FASTQ (`@` record without `+`, or quality length mismatch)
proc readFastx*[T](f: var Bufio[T], r: var FQRecord): bool {.discardable.} =
  var x: int
  var c: char
  if r.lastChar == 0: # the header character hasn't been read yet
    while true:       # look for the header character '>' or '@'
      x = f.readByte()
      if x < 0 or x == int('>') or x == int('@'): break
    if x < 0: r.status = x; return false # end-of-file or stream error
    r.lastChar = x
  let headerChar = r.lastChar
  r.sequence.setLen(0); r.quality.setLen(0); r.comment.setLen(0)
  x = f.readUntil(r.name, c, -2)
  if x < 0: r.status = x; r.lastChar = 0; return false  # EOF or stream error
  if c != '\n': f.readUntil(r.comment, c) # read FASTA/Q comment
  while true:         # read sequence
    x = f.readByte()  # read the first char on a line
    if x < 0 or x == int('>') or x == int('+') or x == int('@'): break
    if x == int('\n'): continue
    r.sequence.add(char(x))
    f.readUntil(r.sequence, c, -1, r.sequence.len)  # read the rest of the seq line
  r.status = r.sequence.len   # for normal records, this keeps the sequence length
  if x == int('>') or x == int('@'): r.lastChar = x
  if x != int('+'):
    if headerChar == int('@'):
      r.status = -4
      return false
    return true
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

#############
# Intervals #
#############

type
  ## Interval record used by the implicit interval-tree utilities.
  ##
  ## Fields:
  ##   st: Interval start coordinate
  ##   en: Interval end coordinate
  ##   data: Payload associated with interval
  ##   max: Cached subtree maximum end coordinate (set by `index`)
  Interval*[S,T] = tuple[st, en: S, data: T, max: S]

## Sort intervals in-place by their start coordinate.
##
## Args:
##   a: Interval sequence to sort
proc sort*[S,T](a: var seq[Interval[S,T]]) =
  a.sort do (x, y: Interval[S,T]) -> int:
    if x.st < y.st: -1
    elif x.st > y.st: 1
    else: 0

## Build interval-tree auxiliary maxima for overlap queries.
##
## Args:
##   a: Interval sequence (sorted automatically if needed)
##
## Returns:
##   Height of the implicit interval tree
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

## Iterate intervals that overlap the query range `[st, en)`.
##
## Args:
##   a: Indexed interval sequence
##   st: Query start coordinate (inclusive)
##   en: Query end coordinate (exclusive)
##
## Yields:
##   Intervals that overlap the query range
iterator overlap*[S,T](a: seq[Interval[S,T]], st: S, en: S): Interval[S,T] {.noSideEffect.} =
  if a.len > 0: # guard against negative tree height on empty input
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
