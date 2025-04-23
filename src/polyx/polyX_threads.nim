import os, strformat, times
import std/cpuinfo
import std/threadpool
import ../../readfx

proc processChunk(chunk: seq[FQRecord], minLen: int = 10): seq[FQRecord] =
  result = newSeq[FQRecord](chunk.len)
  for i in 0 ..< chunk.len:
    result[i] = filtPolyX(chunk[i], minLen = minLen)

when isMainModule:
  let args = commandLineParams()
  if len(args) == 0:
    stderr.writeLine "Missing input parameter [FILENAME]"
    quit(1)
  
  let inputFile = args[0]
  stderr.writeLine("Input file: ", inputFile)
  
  # Load all records first
  var records: seq[FQRecord] = @[]
  for rec in readfq(inputFile):
    records.add(rec)
  
  stderr.writeLine(fmt"Processing {records.len} records")
  
  # Set thread pool size based on CPU count
  let numThreads = countProcessors()
  stderr.writeLine(fmt"Using {numThreads} threads")
  setMaxPoolSize(numThreads)
  
  # Determine chunk size and number of chunks
  let chunkSize = max(1, records.len div numThreads)
  let numChunks = (records.len + chunkSize - 1) div chunkSize
  
  # Create chunks of records
  var chunks: seq[seq[FQRecord]] = newSeq[seq[FQRecord]](numChunks)
  for i in 0 ..< numChunks:
    let start = i * chunkSize
    let endIdx = min(start + chunkSize, records.len)
    chunks[i] = records[start ..< endIdx]
  
  # Process chunks in parallel
  var futures = newSeq[FlowVar[seq[FQRecord]]](numChunks)
  for i in 0 ..< numChunks:
    futures[i] = spawn processChunk(chunks[i], minLen = 10)
  
  # Collect and output results
  for i in 0 ..< numChunks:
    let chunkResults = ^futures[i]
    for result in chunkResults:
      echo result
  
  # Make sure to wait for all threads to complete
  sync()