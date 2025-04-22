import os
import strformat
import strutils
import ../readfx  # This automatically imports fqrecordUtils

when isMainModule:
  let args = commandLineParams()
  
  if len(args) == 0:
    stderr.writeLine "Usage: fqrecordUtils_example [FASTQ file]"
    quit(1)
  
  let inputFile = args[0]
  echo "Processing file: ", inputFile
  
  # Example 1: Reverse complement a sequence
  echo "\n== Example 1: Reverse complement a sequence =="
  let dnaSeq = "ATGCTAGCTAGCTA"
  let rc = reverseComplement(dnaSeq)
  echo fmt"Original sequence: {dnaSeq}"
  echo fmt"Reverse complement: {rc}"
  
  # Example 2: Process a FASTQ file and apply various operations
  echo "\n== Example 2: Process FASTQ and apply operations =="
  
  # Read FASTQ and apply operations to each record
  var recordCount = 0
  for record in readFQ(inputFile):
    inc recordCount
    if recordCount > 3:  # Limit to first 3 records for demo
      break
    
    echo "\nRecord #", recordCount, ":"
    echo "Name: ", record.name
    echo "Original sequence: ", record.sequence[0..min(29, record.sequence.high)], "..."
    echo "Length: ", record.sequence.len
    echo "GC content: ", (gcContent(record.sequence) * 100).formatFloat(ffDecimal, 2), "%"
    
    # Create a reverse complemented version
    let rcRecord = reverseComplementRecord(record)
    echo "Reverse complement: ", rcRecord.sequence[0..min(29, rcRecord.sequence.high)], "..."
    
    # Extract a subsequence (first 10 bases)
    let subseq = subSequence(record, 0, 10)
    echo "First 10 bases: ", subseq.sequence
    
    # Make a modified copy that we can change
    var modifiedRecord = record
    
    # Mask low quality bases if the record has quality scores
    if modifiedRecord.quality.len > 0:
      maskLowQuality(modifiedRecord, 20)  # Mask bases with quality < 20
      echo "After masking low quality: ", modifiedRecord.sequence[0..min(29, modifiedRecord.sequence.high)], "..."
    
    # Trim low quality ends if the record has quality scores
    if modifiedRecord.quality.len > 0:
      let originalLen = modifiedRecord.sequence.len
      qualityTrim(modifiedRecord, 20)
      echo "After quality trimming: ", modifiedRecord.sequence[0..min(29, modifiedRecord.sequence.high)], "..."
      echo "New length: ", modifiedRecord.sequence.len, " (trimmed ", originalLen - modifiedRecord.sequence.len, " bases)"