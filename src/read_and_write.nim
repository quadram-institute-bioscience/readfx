## read_and_write: rename FASTA/FASTQ records while copying input to output
## =======================================================================
##
## Usage:
##   read_and_write -i INPUT_FILE [-o OUTPUT_FILE] [-t THREADS]
##
## Behavior:
## - Record names are rewritten to `seq.<n>` (1-based progressive counter)
## - If `-o/--output` is omitted, output is written to stdout
## - If output filename ends with `.gz`, output is gzip-compressed
## - `-t/--threads` controls compression threads (pigz when available)
## - Output format (FASTA/FASTQ) is inferred from the first input record
import std/[os, strutils]
import ../readfx

type
  CliOptions = object
    input: string
    output: string
    hasOutput: bool
    threads: int

proc usage() =
  stderr.writeLine("Usage: read_and_write -i INPUT_FILE [-o OUTPUT_FILE] [-t THREADS]")
  stderr.writeLine("  -i, --input   Input FASTA/FASTQ file (plain or .gz)")
  stderr.writeLine("  -o, --output  Output file path (optional; defaults to stdout)")
  stderr.writeLine("                If it ends with .gz, output is compressed")
  stderr.writeLine("  -t, --threads Compression threads (default: 1)")

proc parseCli(): CliOptions =
  result.threads = 1
  let args = commandLineParams()
  var i = 0
  while i < args.len:
    let a = args[i]
    case a
    of "-h", "--help":
      usage()
      quit(0)
    of "-i", "--input":
      inc i
      if i >= args.len:
        raise newException(ValueError, "Missing value for --input")
      result.input = args[i]
    of "-o", "--output":
      inc i
      if i >= args.len:
        raise newException(ValueError, "Missing value for --output")
      result.output = args[i]
      result.hasOutput = true
    of "-t", "--threads":
      inc i
      if i >= args.len:
        raise newException(ValueError, "Missing value for --threads")
      result.threads = parseInt(args[i])
    else:
      if a.startsWith("--input="):
        result.input = a.split("=", maxsplit = 1)[1]
      elif a.startsWith("--output="):
        result.output = a.split("=", maxsplit = 1)[1]
        result.hasOutput = true
      elif a.startsWith("--threads="):
        result.threads = parseInt(a.split("=", maxsplit = 1)[1])
      else:
        raise newException(ValueError, "Unknown argument: " & a)
    inc i

  if result.input.len == 0:
    raise newException(ValueError, "Missing required --input")
  if result.threads <= 0:
    raise newException(ValueError, "--threads must be > 0")

proc openWriter(opts: CliOptions, format: FastxFormat): FastxWriter =
  let compress = opts.hasOutput and opts.output.endsWith(".gz")
  let dest = if opts.hasOutput: fileDestination(opts.output) else: stdoutDestination()
  result = fastxWriter(
    format = format,
    compression = compress,
    destination = dest,
    bufferSize = DefaultBufferSize,
    compressionLevel = DefaultCompressionLevel,
    compressionThreads = opts.threads
  )

when isMainModule:
  try:
    let opts = parseCli()
    var inputStream = xopen[GzFile](opts.input)
    defer:
      discard inputStream.close()

    var record: FQRecord
    if not inputStream.readFastx(record):
      # Empty input: create empty output stream only if output file was requested.
      if opts.hasOutput:
        var w = openWriter(opts, fxfFasta)
        w.close()
      quit(0)

    let format = if record.quality.len > 0: fxfFastq else: fxfFasta
    var writer = openWriter(opts, format)
    defer:
      if writer.isOpen:
        writer.close()

    var n = 1
    record.name = "seq." & $n
    writer.writeRecord(record)

    while inputStream.readFastx(record):
      inc n
      record.name = "seq." & $n
      writer.writeRecord(record)

    writer.close()
  except ValueError as e:
    stderr.writeLine("Error: ", e.msg)
    usage()
    quit(2)
  except IOError as e:
    stderr.writeLine("I/O error: ", e.msg)
    quit(1)
