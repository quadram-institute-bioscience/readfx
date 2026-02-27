## Fastx writer utilities for buffered FASTA/FASTQ output.
## =======================================================
##
## This module provides an efficient `FastxWriter` abstraction that supports:
## - FASTA or FASTQ formatting
## - Plain-text or gzip-compressed output
## - File or stdout destination
## - Buffered chunked writes for throughput
##
## Typical usage:
##
## ```nim
## import readfx/writer
## import readfx/seqtypes
##
## var w = fastxWriter(
##   format = fxfFastq,
##   compression = true,
##   destination = fileDestination("out.fastq.gz"),
##   bufferSize = 4 * 1024 * 1024,
##   compressionLevel = 6
## )
## defer: w.close()
##
## w.writeRecord(FQRecord(name: "r1", sequence: "ACGT", quality: "IIII"))
## ```
import zip/zlib
when defined(posix):
  import posix

import ../readfx/seqtypes

const
  DefaultBufferSize* = 4 * 1024 * 1024
  DefaultCompressionLevel* = 6
  DefaultFastaWidth* = 60

type
  ## Output format emitted by `FastxWriter`.
  FastxFormat* = enum
    fxfFasta, ## Emit FASTA records
    fxfFastq  ## Emit FASTQ records

  ## Destination kind used by `FastxDestination`.
  FastxDestinationKind* = enum
    fxdStdout, ## Write to stdout
    fxdFile    ## Write to file path

  ## Destination configuration for `FastxWriter`.
  FastxDestination* = object
    case kind*: FastxDestinationKind
    of fxdStdout:
      discard
    of fxdFile:
      path*: string

  FastxBackend = enum
    fxbPlain,
    fxbGzip

  ## Buffered writer for FASTA/FASTQ output.
  ##
  ## Use `fastxWriter` to initialize this object, then call:
  ## - `writeRecord` repeatedly
  ## - `flush` optionally
  ## - `close` once at the end
  FastxWriter* = object
    format*: FastxFormat
    compression*: bool
    destination*: FastxDestination
    bufferSize*: int
    compressionLevel*: int
    fastaWidth*: int
    isOpen*: bool
    backend: FastxBackend
    plainFile: File
    gzFile: GzFile
    ownsPlainHandle: bool
    outputBuffer: string

proc stdoutDestination*(): FastxDestination =
  ## Create a stdout destination descriptor.
  ##
  ## Returns:
  ##   `FastxDestination(kind: fxdStdout)`
  result = FastxDestination(kind: fxdStdout)

proc fileDestination*(path: string): FastxDestination =
  ## Create a file destination descriptor.
  ##
  ## Args:
  ##   path: Output file path
  ##
  ## Returns:
  ##   `FastxDestination(kind: fxdFile, path: path)`
  result = FastxDestination(kind: fxdFile, path: path)

proc appendWrappedSequence(dst: var string, sequence: string, width: int) =
  if width <= 0:
    dst.add(sequence)
    dst.add('\n')
    return

  var i = 0
  while i < sequence.len:
    let j = min(i + width, sequence.len)
    dst.add(sequence[i ..< j])
    dst.add('\n')
    i = j

proc ensureOpen(w: FastxWriter) =
  if not w.isOpen:
    raise newException(IOError, "FastxWriter is closed")

proc writeChunkPlain(w: var FastxWriter, data: string) =
  if data.len == 0:
    return
  let written = writeBuffer(w.plainFile, cast[pointer](data.cstring), data.len)
  if written != data.len:
    raise newException(IOError, "Short write on plain output stream")

proc zlibErrorMessage(w: FastxWriter): string =
  var errNo: cint
  let msg = gzerror(w.gzFile, errNo)
  if msg != nil:
    result = $cast[cstring](msg)
  else:
    result = "zlib error code " & $errNo

proc writeChunkGzip(w: var FastxWriter, data: string) =
  if data.len == 0:
    return
  let written = gzwrite(w.gzFile, cast[pointer](data.cstring), cuint(data.len))
  if written < 0 or written != cint(data.len):
    raise newException(IOError, "gzwrite failed: " & zlibErrorMessage(w))

proc flush*(w: var FastxWriter) =
  ## Flush pending buffered bytes to the destination stream.
  ##
  ## Notes:
  ##   This writes the current application buffer. For gzip streams, final trailer
  ##   bytes are emitted during `close`.
  ensureOpen(w)
  if w.outputBuffer.len == 0:
    return

  case w.backend
  of fxbPlain:
    writeChunkPlain(w, w.outputBuffer)
  of fxbGzip:
    writeChunkGzip(w, w.outputBuffer)

  w.outputBuffer.setLen(0)

proc maybeFlushBuffer(w: var FastxWriter) =
  if w.outputBuffer.len >= w.bufferSize:
    w.flush()

proc appendFastqRecord(dst: var string, record: FQRecord) =
  dst.add('@')
  dst.add(record.name)
  if record.comment.len > 0:
    dst.add(' ')
    dst.add(record.comment)
  dst.add('\n')
  dst.add(record.sequence)
  dst.add("\n+\n")
  dst.add(record.quality)
  dst.add('\n')

proc appendFastaRecord(dst: var string, record: FQRecord, width: int) =
  dst.add('>')
  dst.add(record.name)
  if record.comment.len > 0:
    dst.add(' ')
    dst.add(record.comment)
  dst.add('\n')
  appendWrappedSequence(dst, record.sequence, width)

proc validateRecord(w: FastxWriter, record: FQRecord) =
  if record.name.len == 0:
    raise newException(ValueError, "Record name cannot be empty")

  case w.format
  of fxfFastq:
    if record.quality.len == 0:
      raise newException(ValueError, "FASTQ output requires a non-empty quality string")
    if record.quality.len != record.sequence.len:
      raise newException(ValueError, "FASTQ quality length must match sequence length")
  of fxfFasta:
    discard

proc writeRecord*(w: var FastxWriter, record: FQRecord) =
  ## Append one record to the writer.
  ##
  ## Args:
  ##   w: Writer instance
  ##   record: Record to write
  ##
  ## Raises:
  ##   IOError: If the writer is closed or I/O fails
  ##   ValueError: If record validation fails for selected format
  ensureOpen(w)
  validateRecord(w, record)

  case w.format
  of fxfFastq:
    appendFastqRecord(w.outputBuffer, record)
  of fxfFasta:
    appendFastaRecord(w.outputBuffer, record, w.fastaWidth)

  maybeFlushBuffer(w)

proc writeRecord*(
    w: var FastxWriter,
    name: string,
    sequence: string,
    quality: string = "",
    comment: string = ""
  ) =
  ## Convenience overload to write from scalar fields.
  ##
  ## Args:
  ##   w: Writer instance
  ##   name: Record identifier
  ##   sequence: Sequence string
  ##   quality: FASTQ quality string (required for FASTQ mode)
  ##   comment: Optional header comment
  writeRecord(w, FQRecord(name: name, sequence: sequence, quality: quality, comment: comment))

proc close*(w: var FastxWriter) =
  ## Flush and close the writer stream.
  ##
  ## It is safe to call `close` multiple times.
  if not w.isOpen:
    return

  w.flush()

  case w.backend
  of fxbPlain:
    if w.ownsPlainHandle and w.plainFile != nil:
      close(w.plainFile)
    elif w.plainFile != nil:
      flushFile(w.plainFile)
  of fxbGzip:
    if w.gzFile != nil:
      let rc = gzclose(w.gzFile)
      if rc != Z_OK:
        raise newException(IOError, "gzclose failed with code " & $rc)

  w.plainFile = nil
  w.gzFile = nil
  w.isOpen = false
  w.outputBuffer.setLen(0)

proc compressionMode(level: int): cstring =
  case level
  of 0: "wb0"
  of 1: "wb1"
  of 2: "wb2"
  of 3: "wb3"
  of 4: "wb4"
  of 5: "wb5"
  of 6: "wb6"
  of 7: "wb7"
  of 8: "wb8"
  of 9: "wb9"
  else:
    "wb6"

proc fastxWriter*(
    format: FastxFormat,
    compression: bool = false,
    destination: FastxDestination = stdoutDestination(),
    bufferSize: int = DefaultBufferSize,
    compressionLevel: int = DefaultCompressionLevel,
    fastaWidth: int = DefaultFastaWidth
  ): FastxWriter =
  ## Create and open a configured FASTA/FASTQ writer.
  ##
  ## Args:
  ##   format: Output format (`fxfFasta` or `fxfFastq`)
  ##   compression: Enable gzip compression for output stream
  ##   destination: Output destination (stdout or file)
  ##   bufferSize: Application-level output buffer size in bytes
  ##   compressionLevel: zlib compression level (0..9)
  ##   fastaWidth: FASTA sequence line width (used only in FASTA mode)
  ##
  ## Returns:
  ##   Initialized `FastxWriter`
  ##
  ## Raises:
  ##   IOError: If destination stream cannot be opened
  ##   ValueError: If configuration values are invalid
  if bufferSize <= 0:
    raise newException(ValueError, "bufferSize must be > 0")
  if compressionLevel < 0 or compressionLevel > 9:
    raise newException(ValueError, "compressionLevel must be between 0 and 9")
  if format == fxfFasta and fastaWidth <= 0:
    raise newException(ValueError, "fastaWidth must be > 0 in FASTA mode")

  result.format = format
  result.compression = compression
  result.destination = destination
  result.bufferSize = bufferSize
  result.compressionLevel = compressionLevel
  result.fastaWidth = fastaWidth
  result.outputBuffer = newStringOfCap(min(bufferSize, 128 * 1024))

  if compression:
    result.backend = fxbGzip
    let modeC = compressionMode(compressionLevel)
    case destination.kind
    of fxdFile:
      if destination.path.len == 0:
        raise newException(ValueError, "Destination path cannot be empty")
      var pathCopy = destination.path & '\0'
      let pathC = cast[cstring](addr pathCopy[0])
      result.gzFile = gzopen(pathC, modeC)
      if result.gzFile == nil:
        raise newException(IOError, "Cannot open gzip output: " & destination.path)
    of fxdStdout:
      when defined(posix):
        let outDup = posix.dup(1)
        if outDup < 0:
          raise newException(IOError, "Cannot duplicate stdout for gzip output")
        result.gzFile = gzdopen(cint(outDup), modeC)
        if result.gzFile == nil:
          discard posix.close(outDup)
          raise newException(IOError, "Cannot open gzip stdout stream")
      else:
        raise newException(IOError, "Gzip stdout is only supported on POSIX platforms")
  else:
    result.backend = fxbPlain
    case destination.kind
    of fxdFile:
      if destination.path.len == 0:
        raise newException(ValueError, "Destination path cannot be empty")
      if not open(result.plainFile, destination.path, fmWrite):
        raise newException(IOError, "Cannot open output file: " & destination.path)
      result.ownsPlainHandle = true
    of fxdStdout:
      result.plainFile = stdout
      result.ownsPlainHandle = false

  result.isOpen = true
