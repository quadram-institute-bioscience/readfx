import std/[parseopt, strformat, strutils, times]
import ../readfx

type
  BenchResult = tuple[pairs: int, bases: int]
  BenchProc = proc(path: string, loops: int): BenchResult {.nimcall.}

proc countFlatPtr(path: string, loops: int): BenchResult =
  var records = 0
  for _ in 0..<loops:
    for rec in readFQPtr(path):
      inc records
      result.bases += rec.sequenceLen
  result.pairs = records div 2

proc countInterleavedPairPtr(path: string, loops: int): BenchResult =
  for _ in 0..<loops:
    for pair in readFQInterleavedPairPtr(path):
      inc result.pairs
      result.bases += pair.read1.sequenceLen + pair.read2.sequenceLen

proc countStringRegroup(path: string, loops: int): BenchResult =
  for _ in 0..<loops:
    var pending: FQRecord
    var havePending = false
    for rec in readFQ(path):
      if not havePending:
        pending = rec
        havePending = true
      else:
        inc result.pairs
        result.bases += pending.sequence.len + rec.sequence.len
        havePending = false

proc runBench(label: string, bench: BenchProc, path: string, loops: int) =
  let memBefore = getOccupiedMem()
  let start = epochTime()
  let resultCounts = bench(path, loops)
  let elapsed = epochTime() - start
  let memAfter = getOccupiedMem()
  let pairsPerSec = if elapsed > 0: resultCounts.pairs.float / elapsed else: 0.0

  echo label & "\t" & path & "\t" & $loops & "\t" &
       $resultCounts.pairs & "\t" & $resultCounts.bases & "\t" &
       fmt"{elapsed:.6f}" & "\t" & fmt"{pairsPerSec:.2f}" & "\t" &
       $(memAfter - memBefore)

proc main() =
  var
    path = "tests/test_interleaved.fq"
    loops = 50000

  var p = initOptParser()
  while true:
    p.next()
    case p.kind
    of cmdEnd:
      break
    of cmdLongOption, cmdShortOption:
      case p.key
      of "path":
        path = p.val
      of "loops":
        loops = parseInt(p.val)
      else:
        raise newException(ValueError, "Unknown option: " & p.key)
    of cmdArgument:
      discard

  if loops <= 0:
    raise newException(ValueError, "--loops must be > 0")

  echo "name\tpath\tloops\ttotal_pairs\ttotal_bases\tseconds\tpairs_per_sec\talloc_delta"
  runBench("readFQPtr-flat", countFlatPtr, path, loops)
  runBench("readFQInterleavedPairPtr", countInterleavedPairPtr, path, loops)
  runBench("readFQ-string-regroup", countStringRegroup, path, loops)

when isMainModule:
  main()
