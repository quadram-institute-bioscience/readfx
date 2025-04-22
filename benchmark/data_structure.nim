import times
import strutils
import random

# Define both types for comparison
type
  FQRecordTuple* = tuple[sequence, quality, name, comment: string, status, lastChar: int]

  FQRecordObject* = object
    name*: string
    comment*: string  # optional
    sequence*: string
    quality*: string  # optional
    status, lastChar: int

# Generate random DNA sequence of specified length
proc randomDNA(length: int): string =
  const bases = ['A', 'C', 'G', 'T']
  result = newString(length)
  for i in 0..<length:
    result[i] = bases[rand(3)]

# Generate random quality string of specified length
proc randomQuality(length: int): string =
  result = newString(length)
  for i in 0..<length:
    result[i] = char(33 + rand(40))  # ASCII 33-73 for quality scores

# Benchmark creation operations
proc benchmarkCreation(iterations: int) =
  echo "Benchmarking record creation..."
  
  # Time tuple creation
  let tupleStartTime = cpuTime()
  var tuples = newSeq[FQRecordTuple](iterations)
  for i in 0..<iterations:
    let seqLen = 100 + rand(900)  # Random sequence length between 100-1000
    tuples[i] = (
      sequence: randomDNA(seqLen),
      quality: randomQuality(seqLen),
      name: "read_" & $i,
      comment: "comment_" & $i,
      status: i mod 5,
      lastChar: i mod 256
    )
  let tupleDuration = cpuTime() - tupleStartTime
  echo "  Tuple creation time: ", formatFloat(tupleDuration, ffDecimal, 6), " seconds"
  
  # Time object creation
  let objectStartTime = cpuTime()
  var objects = newSeq[FQRecordObject](iterations)
  for i in 0..<iterations:
    let seqLen = 100 + rand(900)  # Same sequence length as tuple test
    objects[i] = FQRecordObject(
      sequence: randomDNA(seqLen),
      quality: randomQuality(seqLen),
      name: "read_" & $i,
      comment: "comment_" & $i,
      status: i mod 5,
      lastChar: i mod 256
    )
  let objectDuration = cpuTime() - objectStartTime
  echo "  Object creation time: ", formatFloat(objectDuration, ffDecimal, 6), " seconds"
  
  # Report speed comparison
  let ratio = tupleDuration / objectDuration
  if ratio > 1.0:
    echo "  Object creation is ", formatFloat(ratio, ffDecimal, 2), "x faster"
  else:
    echo "  Tuple creation is ", formatFloat(1.0/ratio, ffDecimal, 2), "x faster"

# Benchmark field access operations
proc benchmarkAccess(iterations: int) =
  echo "Benchmarking field access..."
  
  # Create test data first
  var tuples = newSeq[FQRecordTuple](iterations)
  var objects = newSeq[FQRecordObject](iterations)
  
  for i in 0..<iterations:
    let seqLen = 100 + rand(900)
    let seq = randomDNA(seqLen)
    let qual = randomQuality(seqLen)
    
    tuples[i] = (
      sequence: seq,
      quality: qual,
      name: "read_" & $i,
      comment: "comment_" & $i,
      status: i mod 5,
      lastChar: i mod 256
    )
    
    objects[i] = FQRecordObject(
      sequence: seq,
      quality: qual,
      name: "read_" & $i,
      comment: "comment_" & $i,
      status: i mod 5,
      lastChar: i mod 256
    )
  
  # Benchmark tuple access
  var tupleSum = 0
  var tupleSeqLen = 0
  let tupleStartTime = cpuTime()
  for i in 0..<iterations:
    # Access all fields to ensure fair comparison
    tupleSum += tuples[i].status + tuples[i].lastChar
    tupleSeqLen += tuples[i].sequence.len + tuples[i].quality.len
    discard tuples[i].name
    discard tuples[i].comment
  let tupleDuration = cpuTime() - tupleStartTime
  
  # Benchmark object access
  var objectSum = 0
  var objectSeqLen = 0
  let objectStartTime = cpuTime()
  for i in 0..<iterations:
    # Access all fields to ensure fair comparison
    objectSum += objects[i].status + objects[i].lastChar
    objectSeqLen += objects[i].sequence.len + objects[i].quality.len
    discard objects[i].name
    discard objects[i].comment
  let objectDuration = cpuTime() - objectStartTime
  
  # Report results
  echo "  Tuple access time: ", formatFloat(tupleDuration, ffDecimal, 6), " seconds"
  echo "  Object access time: ", formatFloat(objectDuration, ffDecimal, 6), " seconds"
  
  # Report speed comparison
  let ratio = tupleDuration / objectDuration
  if ratio > 1.0:
    echo "  Object access is ", formatFloat(ratio, ffDecimal, 2), "x faster"
  else:
    echo "  Tuple access is ", formatFloat(1.0/ratio, ffDecimal, 2), "x faster"

# Benchmark assignment operations
proc benchmarkAssignment(iterations: int) =
  echo "Benchmarking assignment operations..."
  
  # Create test data first
  let seqLen = 500  # Fixed size for consistent comparison
  let seq = randomDNA(seqLen)
  let qual = randomQuality(seqLen)
  
  let tupleTemplate = (
    sequence: seq,
    quality: qual,
    name: "read_template",
    comment: "comment_template",
    status: 0,
    lastChar: 0
  )
  
  let objectTemplate = FQRecordObject(
    sequence: seq,
    quality: qual,
    name: "read_template",
    comment: "comment_template",
    status: 0,
    lastChar: 0
  )
  
  # Benchmark tuple assignment (value semantics - creates copies)
  let tupleStartTime = cpuTime()
  for i in 0..<iterations:
    var tupleCopy = tupleTemplate  # Creates a full copy
    tupleCopy.status = i
    tupleCopy.lastChar = i mod 256
    tupleCopy.sequence = "MODIFIED" & tupleCopy.sequence[8..^1]
  let tupleDuration = cpuTime() - tupleStartTime
  
  # Benchmark object assignment with new instance
  let objectStartTime = cpuTime()
  for i in 0..<iterations:
    var objectCopy = objectTemplate  # Also creates a full copy due to string content
    objectCopy.status = i
    objectCopy.lastChar = i mod 256
    objectCopy.sequence = "MODIFIED" & objectCopy.sequence[8..^1]
  let objectDuration = cpuTime() - objectStartTime
  
  # Report results
  echo "  Tuple assignment time: ", formatFloat(tupleDuration, ffDecimal, 6), " seconds"
  echo "  Object assignment time: ", formatFloat(objectDuration, ffDecimal, 6), " seconds"
  
  # Report speed comparison
  let ratio = tupleDuration / objectDuration
  if ratio > 1.0:
    echo "  Object assignment is ", formatFloat(ratio, ffDecimal, 2), "x faster"
  else:
    echo "  Tuple assignment is ", formatFloat(1.0/ratio, ffDecimal, 2), "x faster"

# Benchmark memory usage through allocation counts
proc benchmarkMemory(iterations: int) =
  echo "Memory comparison through GC stats..."
  
  # Benchmark tuples
  GC_fullCollect()
  let tupleAllocsBefore = getOccupiedMem()
  var tuples = newSeq[FQRecordTuple](iterations)
  for i in 0..<iterations:
    let seqLen = 100  # Fixed size for consistent comparison
    tuples[i] = (
      sequence: randomDNA(seqLen),
      quality: randomQuality(seqLen),
      name: "read_" & $i,
      comment: "comment_" & $i,
      status: i mod 5,
      lastChar: i mod 256
    )
  let tupleAllocsAfter = getOccupiedMem()
  let tupleMemUsage = tupleAllocsAfter - tupleAllocsBefore
  
  # Clear tuples to free memory
  tuples = @[]
  GC_fullCollect()
  
  # Benchmark objects
  let objectAllocsBefore = getOccupiedMem()
  var objects = newSeq[FQRecordObject](iterations)
  for i in 0..<iterations:
    let seqLen = 100  # Same fixed size
    objects[i] = FQRecordObject(
      sequence: randomDNA(seqLen),
      quality: randomQuality(seqLen),
      name: "read_" & $i,
      comment: "comment_" & $i,
      status: i mod 5,
      lastChar: i mod 256
    )
  let objectAllocsAfter = getOccupiedMem()
  let objectMemUsage = objectAllocsAfter - objectAllocsBefore
  
  # Report results
  echo "  Tuple memory: ", tupleMemUsage, " bytes"
  echo "  Object memory: ", objectMemUsage, " bytes"
  echo "  Memory ratio: ", formatFloat(float(objectMemUsage) / float(tupleMemUsage), ffDecimal, 3)

# Run the benchmarks
when isMainModule:
  randomize()  # Initialize random number generator
  
  let iterations = 1_000_000
  echo "Running benchmarks with ", iterations, " iterations"
  echo "========================================"
  
  benchmarkCreation(iterations)
  echo "----------------------------------------"
  
  benchmarkAccess(iterations)
  echo "----------------------------------------"
  
  benchmarkAssignment(iterations div 10)  # Use fewer iterations for assignment test
  echo "----------------------------------------"
  
  benchmarkMemory(iterations div 10)
  echo "========================================"
  
  echo "Benchmarking completed."