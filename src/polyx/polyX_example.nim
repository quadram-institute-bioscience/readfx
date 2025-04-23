import os

import ../../readfx
import strformat, times

when isMainModule:
  let args = commandLineParams()
  if len(args) == 0:
    stderr.writeLine "Missing input parameter [FILENAME]"
    quit(1)
  var
    inputFile = args[0]
  stderr.writeLine("Input file: ", inputFile)

    
  for rec in readfq(inputFile):
    echo filtPolyX(rec, minLen = 10)