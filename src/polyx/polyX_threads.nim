import os, strformat, times
import std/cpuinfo
import std/threadpool
import ../../readfx

# PreAllocate result to avoid reallocation inside processChunk
const ChunkSize = 50_000

# Use a lock-free queue for collecting results
type ResultQueue = object
  data: seq[FQRecord]
  count: int

proc createResultQueue(): ptr ResultQueue =
  result = cast[ptr ResultQueue](allocShared0(sizeof(ResultQueue)))
  # Allocate more space than ChunkSize to handle larger results
  result.data = newSeq[FQRecord](ChunkSize * 4)
  result.count = 0

proc destroyResultQueue(queue: ptr ResultQueue) =
  queue.data = @[]
  deallocShared(queue)

proc processChunk(chunk: seq[FQRecord], minLen: int = 10, results: ptr ResultQueue): int =
  # Process each record in the chunk and return count of valid results
  var resultCount = 0
  for i in 0 ..< chunk.len:
    let processed = filtPolyX(chunk[i], minLen = minLen)
    if len(processed.sequence) > 0:
      # Store directly in the result queue with bounds check
      if resultCount < results.data.len:
        results.data[resultCount] = processed
        resultCount += 1
      else:
        stderr.writeLine("Warning: Result queue overflow, record skipped")
  
  # Atomically update the count in the queue
  atomicInc(results.count, resultCount)
  return resultCount

when isMainModule:
  let args = commandLineParams()
  if len(args) == 0:
    stderr.writeLine "Missing input parameter [FILENAME]"
    quit(1)
  
  let inputFile = args[0]
  stderr.writeLine("Input file: ", inputFile)
  
  let startTime = cpuTime()
  
  # Set thread pool size based on CPU count
  let numThreads = countProcessors()
  stderr.writeLine(fmt"Using {numThreads} threads")
  setMaxPoolSize(numThreads)
  
  # Initialize shared result queues
  var resultQueues = newSeq[ptr ResultQueue](numThreads)
  for i in 0 ..< numThreads:
    resultQueues[i] = createResultQueue()
  
  # Use a buffer to collect records
  var buffer = newSeq[FQRecord]()
  var threadIndex = 0
  var futures = newSeq[FlowVar[int]]()
  
  # Process FASTQ file in chunks without loading everything into memory
  for rec in readfq(inputFile):
    buffer.add(rec)
    
    # If buffer reaches chunk size, process it in a thread
    if buffer.len >= ChunkSize:
      let currentBuffer = buffer
      buffer = newSeq[FQRecord]()
      futures.add(spawn processChunk(currentBuffer, minLen = 10, resultQueues[threadIndex]))
      threadIndex = (threadIndex + 1) mod numThreads
  
  # Process any remaining records
  if buffer.len > 0:
    futures.add(spawn processChunk(buffer, minLen = 10, resultQueues[threadIndex]))
  
  # Wait for all tasks to complete and output results
  var totalResults = 0
  for future in futures:
    totalResults += ^future
  
  # Output all results
  for q in resultQueues:
    for i in 0 ..< q.count:
      echo q.data[i]
  
  # Clean up
  for q in resultQueues:
    destroyResultQueue(q)
  
  let endTime = cpuTime()
  stderr.writeLine(fmt"Processed {totalResults} records in {endTime - startTime:.2f} seconds")