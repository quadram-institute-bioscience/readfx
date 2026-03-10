import std/[algorithm, parseopt, strformat, strutils, times]
import ../readfx

type
  BenchResult = tuple[units: int, bases: int]
  BenchProc = proc(primaryPath: string, secondaryPath: string, loops: int): BenchResult {.nimcall.}

proc median(values: seq[float]): float =
  if values.len == 0:
    return 0.0
  var sorted = values
  sorted.sort()
  let mid = sorted.len div 2
  if (sorted.len mod 2) == 1:
    sorted[mid]
  else:
    (sorted[mid - 1] + sorted[mid]) / 2.0

proc countReadFQPtr(primaryPath: string, secondaryPath: string, loops: int): BenchResult =
  discard secondaryPath
  for _ in 0..<loops:
    for rec in readFQPtr(primaryPath):
      inc result.units
      result.bases += rec.sequenceLen

proc countReadFQ(primaryPath: string, secondaryPath: string, loops: int): BenchResult =
  discard secondaryPath
  for _ in 0..<loops:
    for rec in readFQ(primaryPath):
      inc result.units
      result.bases += rec.sequence.len

proc countReadFastx(primaryPath: string, secondaryPath: string, loops: int): BenchResult =
  discard secondaryPath
  for _ in 0..<loops:
    var record: FQRecord
    var f = xopen[GzFile](primaryPath)
    defer: f.close()
    while f.readFastx(record):
      inc result.units
      result.bases += record.sequence.len

proc countReadFQPairPtr(primaryPath: string, secondaryPath: string, loops: int): BenchResult =
  for _ in 0..<loops:
    for pair in readFQPairPtr(primaryPath, secondaryPath):
      inc result.units
      result.bases += pair.read1.sequenceLen + pair.read2.sequenceLen

proc countReadFQPair(primaryPath: string, secondaryPath: string, loops: int): BenchResult =
  for _ in 0..<loops:
    for pair in readFQPair(primaryPath, secondaryPath):
      inc result.units
      result.bases += pair.read1.sequence.len + pair.read2.sequence.len

proc countReadFQInterleavedPairPtr(primaryPath: string, secondaryPath: string, loops: int): BenchResult =
  discard secondaryPath
  for _ in 0..<loops:
    for pair in readFQInterleavedPairPtr(primaryPath):
      inc result.units
      result.bases += pair.read1.sequenceLen + pair.read2.sequenceLen

proc countReadFQStringRegroup(primaryPath: string, secondaryPath: string, loops: int): BenchResult =
  discard secondaryPath
  for _ in 0..<loops:
    var pending: FQRecord
    var havePending = false
    for rec in readFQ(primaryPath):
      if not havePending:
        pending = rec
        havePending = true
      else:
        inc result.units
        result.bases += pending.sequence.len + rec.sequence.len
        havePending = false

proc runBench(
    label: string,
    dataset: string,
    unitLabel: string,
    bench: BenchProc,
    primaryPath: string,
    secondaryPath: string,
    loops: int,
    repeats: int
  ) =
  var seconds: seq[float] = @[]
  var firstRun: BenchResult

  for rep in 0..<repeats:
    let start = epochTime()
    let benchResult = bench(primaryPath, secondaryPath, loops)
    let elapsed = epochTime() - start
    seconds.add(elapsed)
    if rep == 0:
      firstRun = benchResult

  let medianSec = median(seconds)
  let unitsPerSec = if medianSec > 0: firstRun.units.float / medianSec else: 0.0
  echo label & "\t" & dataset & "\t" & unitLabel & "\t" & primaryPath & "\t" &
       secondaryPath & "\t" & $loops & "\t" & $repeats & "\t" &
       $firstRun.units & "\t" & $firstRun.bases & "\t" &
       fmt"{medianSec:.6f}" & "\t" & fmt"{unitsPerSec:.2f}"

proc main() =
  var
    singlePath = "tests/illumina_1.fq.gz"
    pairPath1 = "tests/illumina_1.fq.gz"
    pairPath2 = "tests/illumina_2.fq.gz"
    interleavedPath = "tests/test_interleaved.fq"
    loops = 2000
    repeats = 5

  var p = initOptParser()
  while true:
    p.next()
    case p.kind
    of cmdEnd:
      break
    of cmdLongOption, cmdShortOption:
      if p.key.len == 0:
        continue
      case p.key
      of "single-path":
        singlePath = p.val
      of "pair-path1":
        pairPath1 = p.val
      of "pair-path2":
        pairPath2 = p.val
      of "interleaved-path":
        interleavedPath = p.val
      of "loops":
        loops = parseInt(p.val)
      of "repeats":
        repeats = parseInt(p.val)
      else:
        raise newException(ValueError, "Unknown option: " & p.key)
    of cmdArgument:
      discard

  if loops <= 0:
    raise newException(ValueError, "--loops must be > 0")
  if repeats <= 0:
    raise newException(ValueError, "--repeats must be > 0")

  echo "name\tdataset\tunit\tprimary_path\tsecondary_path\tloops\trepeats\ttotal_units\ttotal_bases\tmedian_seconds\tunits_per_sec"
  runBench("readFQPtr", "single", "records", countReadFQPtr, singlePath, "", loops, repeats)
  runBench("readFQ", "single", "records", countReadFQ, singlePath, "", loops, repeats)
  runBench("readFastx", "single", "records", countReadFastx, singlePath, "", loops, repeats)
  runBench("readFQPairPtr", "paired", "pairs", countReadFQPairPtr, pairPath1, pairPath2, loops, repeats)
  runBench("readFQPair", "paired", "pairs", countReadFQPair, pairPath1, pairPath2, loops, repeats)
  runBench("readFQInterleavedPairPtr", "interleaved", "pairs", countReadFQInterleavedPairPtr, interleavedPath, "", loops, repeats)
  runBench("readFQ-string-regroup", "interleaved", "pairs", countReadFQStringRegroup, interleavedPath, "", loops, repeats)

when isMainModule:
  main()
