import strutils
import unittest

const ProductionSources = [
  "readfx.nim",
  "readfx/kseq.h",
  "readfx/nimklib.nim",
  "readfx/writer.nim"
]

const ForbiddenFragments = [
  "import zip/zlib",
  "from zip/zlib",
  "passL: \"-lz\"",
  "passL:\"-lz\"",
  "importc: \"gzopen\"",
  "importc: \"gzdopen\"",
  "importc: \"gzread\"",
  "importc: \"gzclose\"",
  "KSEQ_INIT(gzFile, gzread)"
]

suite "system libz independence":
  test "production source does not bind system zlib":
    for path in ProductionSources:
      let source = readFile(path)
      for fragment in ForbiddenFragments:
        check source.find(fragment) == -1

  test "kseq is wired to the gzfast-backed read callback":
    let source = readFile("readfx/kseq.h")
    check source.contains("readfx_gzfast_read")
    check source.contains("KSEQ_INIT(readfx_gzfast_file, readfx_gzfast_read)")
