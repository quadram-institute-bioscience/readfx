import times
import strutils
import zip/zlib
import zip/gzipfiles
import os
import ../readfx
#import nimbioseq
import fastx_reader



# https://stackoverflow.com/questions/36577570/how-to-benchmark-few-lines-of-code-in-nim
template benchmark(benchmarkName: string, code: untyped) =
  block:
    let t0 = epochTime()
    code
    let elapsed = (epochTime() - t0 ) * 1000
    let elapsedStr = elapsed.formatFloat(format = ffDecimal, precision = 3) 
    echo  benchmarkName, "\t", elapsedStr, "ms"


proc readfqPtr_count(path: string): int =
  for rec in readfx.readFQPtr(path):
    inc result

proc readfq_count(path: string): int =
  for rec in readfx.readFQ(path):
    inc result

proc readFastq_count(path: string): int =
  var i = 0
  #for rec in readfastq(path):
  #  inc result

proc klib_count(path: string): int =
  var r: FQRecord
  var n = 0
  var f = xopen[readfx.GzFile](path)
  defer: f.close()
  while f.readFastx(r):
    inc result

proc fastq_reader_count(path: string): int =
  for name, sequence, quality in fastq_reader(open(path)):
    inc result


when isMainModule:
  echo "Start"
  # see https://github.com/lh3/biofast/releases/tag/biofast-data-v1)
  var fq = "./hiseq.fq"

  benchmark "readfq count":
    echo "# n=" & $readfq_count(fq)

  benchmark "readfqPtr count":
    echo "# n=" & $readfqptr_count(fq)

  benchmark "klib count":
    echo "# n=" & $klib_count(fq)

#  benchmark "bioseq count":
#    echo "# n=" & $readFastq_count(fq)

  benchmark "fastx count":
    echo "# n=" & $fastq_reader_count(fq)

  # see https://github.com/lh3/biofast/releases/tag/biofast-data-v1)
  var fqgz = "./hiseq.fq.gz"

  benchmark "readfq gz count":
    echo "# n=" & $readfq_count(fqgz)

  benchmark "klib gz count":
    echo "# n=" & $klib_count(fqgz)

#  benchmark "bioseq gz count":
#    echo "# n=" & $readFastq_count(fqgz)

  var fa = "./genome.fa"

  benchmark "readfq count (fa)":
    echo "# n=" & $readfq_count(fa)

  benchmark "readfqPtr count (fa)":
    echo "# n=" & $readfqptr_count(fa)

  benchmark "klib count (fa)":
    echo "# n=" & $klib_count(fa)

  var fagz = "./genome.fa.gz"

  benchmark "readfq gz count (fa.gz)":
    echo "# n=" & $readfq_count(fagz)

  benchmark "klib gz count (fa.gz)":
    echo "# n=" & $klib_count(fagz)