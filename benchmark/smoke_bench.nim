import std/[algorithm, parseopt, strformat, strutils, times]
import ../readfx

type
  CounterProc = proc(path: string): int {.nimcall.}

proc countReadFQPtr(path: string): int =
  for rec in readFQPtr(path):
    discard rec
    inc result

proc countReadFQ(path: string): int =
  for rec in readFQ(path):
    discard rec
    inc result

proc countReadFastx(path: string): int =
  var r: FQRecord
  var f = xopen[GzFile](path)
  defer: f.close()
  while f.readFastx(r):
    inc result

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

proc runCounter(label: string, counter: CounterProc, path: string, loops: int, repeats: int) =
  var seconds: seq[float] = @[]
  var recordsPerRun = 0

  for rep in 0..<repeats:
    let start = epochTime()
    var total = 0
    for i in 0..<loops:
      total += counter(path)
    let elapsed = epochTime() - start
    seconds.add(elapsed)
    if rep == 0:
      recordsPerRun = total

  let medianSec = median(seconds)
  let recordsPerSec = if medianSec > 0: recordsPerRun.float / medianSec else: 0.0
  echo label & "\t" & path & "\t" & $loops & "\t" & $repeats & "\t" &
       $recordsPerRun & "\t" & fmt"{medianSec:.6f}" & "\t" & fmt"{recordsPerSec:.2f}"

proc main() =
  var
    path = "tests/illumina_1.fq.gz"
    loops = 2000
    repeats = 5

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

  echo "name\tpath\tloops\trepeats\ttotal_records\tmedian_seconds\trecords_per_sec"
  runCounter("readFQPtr", countReadFQPtr, path, loops, repeats)
  runCounter("readFQ", countReadFQ, path, loops, repeats)
  runCounter("readFastx", countReadFastx, path, loops, repeats)

when isMainModule:
  main()
