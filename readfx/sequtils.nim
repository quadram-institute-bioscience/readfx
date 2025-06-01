## SeqFu common methods
## ==========================================================
## This module provides common methods for manipulating
## and analyzing sequence data in FASTA/FASTQ format.
## It includes functions for trimming, quality scoring,
## reverse complementing, and calculating GC content.
##
import seqtypes
import std/strutils

proc trimStart*(record: FQRecord, bases: int): FQRecord =
  ## Remove N bases from the start (5' end) of a sequence record
  ## 
  ## Args:
  ##   record: Input FQRecord
  ##   bases: Number of bases to remove from the start
  ## 
  ## Returns:
  ##   New FQRecord with trimmed sequence and quality
  result = record
  let seqLen = record.sequence.len
  
  if bases <= 0 or seqLen == 0:
    return result
    
  let trimBases = min(bases, seqLen)
  result.sequence = record.sequence[trimBases..^1]
  
  if record.quality.len > 0:
    result.quality = record.quality[trimBases..^1]
  
  return result

proc trimEnd*(record: FQRecord, bases: int): FQRecord =
  ## Remove N bases from the end (3' end) of a sequence record
  ## 
  ## Args:
  ##   record: Input FQRecord
  ##   bases: Number of bases to remove from the end
  ## 
  ## Returns:
  ##   New FQRecord with trimmed sequence and quality
  result = record
  let seqLen = record.sequence.len
  
  if bases <= 0 or seqLen == 0:
    return result
    
  let trimBases = min(bases, seqLen)
  let newLen = seqLen - trimBases
  
  result.sequence = record.sequence[0..<newLen]
  
  if record.quality.len > 0:
    result.quality = record.quality[0..<newLen]
  
  return result

proc qualCharToInt*(c: char, offset: int = 33): int =
  ## Convert a quality character to its integer value
  ## 
  ## Args:
  ##   c: Quality character
  ##   offset: Quality score offset (default is 33 for Sanger/Illumina 1.8+)
  ## 
  ## Returns:
  ##   Integer value of the quality character
  return ord(c) - offset

proc qualIntToChar*(q: int, offset: int = 33): char =
    ## Convert an integer quality value to its character representation
    ## 
    ## Args:
    ##   q: Integer quality value
    ##   offset: Quality score offset (default is 33 for Sanger/Illumina 1.8+)
    ## 
    ## Returns:
    ##   Character representation of the quality value
    return chr(q + offset)

proc avgQuality*(record: FQRecord, offset: int = 33): float =
  ## Calculate the average quality of a sequence record
  ## 
  ## Args:
  ##   record: FQRecord to analyze
  ##   offset: Quality score offset (default is 33 for Sanger/Illumina 1.8+)
  ## 
  ## Returns:
  ##   Average quality score as a float
  if record.quality.len == 0 or record.quality.len != record.sequence.len:
    return -1.0 # Invalid record
  
  var totalQuality = 0
  for c in record.quality:
    totalQuality += qualCharToInt(c, offset)
  
  return totalQuality / record.sequence.len

proc avgQuality*(quality: string, offset: int = 33): float =
    ## Calculate the average quality of a quality string
    ## 
    ## Args:
    ##   quality: Quality string
    ##   offset: Quality score offset (default is 33 for Sanger/Illumina 1.8+)
    ## 
    ## Returns:
    ##   Average quality score as a float
    if quality.len == 0:
        return -1.0 # Invalid quality string
    
    var totalQuality = 0
    for c in quality:
        totalQuality += qualCharToInt(c, offset)
    
    return totalQuality / quality.len


proc filtPolyX*(read: FQRecord, minLen: int = 10, 
               trackStats: bool = false): FQRecord =
  ## Removes homopolymer runs (poly-X) from the 3' end of a read.
  ##
  ## Parameters:
  ##   read: The input sequence read
  ##   minLen: Minimum length of homopolymer run to trigger trimming (default: 10)
  ##   trackStats: Whether to track statistics about trimming (not implemented in this version)
  ##
  ## Returns:
  ##   A new FQRecord with poly-X tails trimmed if needed
  
  # Hard-coded parameters (as per original algorithm)
  const 
    allowOneMismatchForEach = 8
    maxMismatch = 5
  
  # Initialize result as a copy of the input
  result = read
  
  let readLength = read.sequence.len
  if readLength == 0:
    return result
  
  # Initialize counters for A, T, C, G
  var baseCount = [0, 0, 0, 0]
  
  # Scan from 3' end
  var position = 0
  while position < readLength:
    # Get base at current position from the 3' end
    let base = read.sequence[readLength - position - 1].toUpperAscii
    
    # Increment the counter for this base
    case base:
      of 'A': baseCount[0] += 1
      of 'T': baseCount[1] += 1
      of 'C': baseCount[2] += 1
      of 'G': baseCount[3] += 1
      of 'N': 
        # N could be any base, increment all
        baseCount[0] += 1
        baseCount[1] += 1
        baseCount[2] += 1
        baseCount[3] += 1
      else:
        discard  # Handle unexpected characters
    
    # Calculate allowed mismatches
    let current = position + 1
    let allowedMismatch = min(maxMismatch, current div allowOneMismatchForEach)
    
    # Check if any base meets homopolymer criteria
    var needToBreak = true
    for b in 0..3:
      if (current - baseCount[b]) <= allowedMismatch:
        needToBreak = false
        break
    
    # Stop scanning if no base meets criteria or we've reached our stopping conditions
    if needToBreak and position >= minLen:
      break
      
    position += 1
  
  # Now determine if we need to trim based on the longest homopolymer run
  var 
    maxBase = 0
    maxCount = baseCount[0]
  
  for b in 1..3:
    if baseCount[b] > maxCount:
      maxCount = baseCount[b]
      maxBase = b
  
  # Handle the case where the entire sequence is a homopolymer
  if maxCount >= minLen and position >= readLength:
    # If the entire sequence matches the criteria, return an empty record
    result.sequence = ""
    if read.quality.len > 0:
      result.quality = ""
    return result
  
  # Only trim if the dominant base count meets or exceeds minLen and we need to trim
  if maxCount >= minLen and position > 0:
    # Get the dominant base
    let polyBase = "ATCG"[maxBase]
    
    # Find the exact position to trim at (avoid indexing outside the string)
    var trimPos = position
    var foundValidPos = false
    var finalTrimPos = 0
    
    # Search for the last occurrence of polyBase within the window
    while trimPos > 0:
      let idx = readLength - trimPos
      if idx >= 0 and idx < readLength and read.sequence[idx].toUpperAscii == polyBase:
        foundValidPos = true
        finalTrimPos = idx
        break
      trimPos -= 1
    
    # Only trim if we found a valid position
    if foundValidPos:
      if finalTrimPos > 0:
        result.sequence = read.sequence[0 ..< finalTrimPos]
        if read.quality.len > 0:
          result.quality = read.quality[0 ..< finalTrimPos]
  
  return result

proc rc_string(sequence: string): string =
  ## Reverse complement a DNA sequence
  ## 
  ## Example:
  ##   let rc = reverseComplement("ATGC")  # returns "GCAT"
  result = newString(sequence.len)
  for i in 0 ..< sequence.len:
    let c = sequence[sequence.len - 1 - i]
    result[i] = case c
      of 'A', 'a': 'T'
      of 'U', 'u', 'T', 't': 'A'
      of 'G', 'g': 'C'
      of 'C', 'c': 'G'
      of 'N', 'n': 'N'
      of 'R', 'r': 'Y'
      of 'Y', 'y': 'R'
      of 'S', 's': 'S'
      of 'W', 'w': 'W'
      of 'K', 'k': 'M'
      of 'M', 'm': 'K'
      of 'B', 'b': 'V'
      of 'V', 'v': 'B'
      of 'D', 'd': 'H'
      of 'H', 'h': 'D'
      else: c

proc trimQuality*(quality: string, minQual: int, offset: int = 33): string =
  ## Trim a quality string based on minimum quality threshold
  ## 
  ## Args:
  ##   quality: Quality string
  ##   minQual: Minimum quality value (0-40)
  ##   offset: Quality score offset (33 for Sanger/Illumina 1.8+)
  ## 
  ## Returns:
  ##   Trimmed quality string
  var endPos = quality.high
  while endPos >= 0 and (ord(quality[endPos]) - offset) < minQual:
    dec endPos
  
  if endPos < 0:
    return ""
  
  return quality[0..endPos]

proc qualityTrim*(record: var FQRecord, minQual: int, offset: int = 33) =
  ## Trim a record based on quality scores
  ## 
  ## Args:
  ##   record: FQRecord to modify
  ##   minQual: Minimum quality value
  ##   offset: Quality score offset (33 for Sanger/Illumina 1.8+)
  if record.quality.len == 0:
    return # FASTA has no quality, nothing to do
  
  let newQuality = trimQuality(record.quality, minQual, offset)
  let newLen = newQuality.len
  
  record.quality = newQuality
  if newLen < record.sequence.len:
    record.sequence = record.sequence[0 ..< newLen]

proc revCompl*(sequence: string): string =
  ## Reverse complement a DNA sequence
  ## 
  ## Args:
  ##   sequence: DNA sequence
  ## 
  ## Returns:
  ##   Reverse-complemented sequence
  result = rc_string(sequence)

proc revCompl*(record: var FQRecord) =
  ## Reverse complement a sequence record in place
  ## 
  ## Args:
  ##   record: FQRecord to modify
  record.sequence = revCompl(record.sequence)
  if record.quality.len > 0:
    # For FASTQ records, also reverse the quality string
    var reversed = ""
    for i in countdown(record.quality.high, 0):
      reversed.add(record.quality[i])
    record.quality = reversed

proc revCompl*(record: FQRecord): FQRecord =
  ## Create a new record with reverse-complemented sequence
  ## 
  ## Args:
  ##   record: Input FQRecord
  ## 
  ## Returns:
  ##   New FQRecord with reverse-complemented sequence
  result = record
  revCompl(result)

proc subSequence*(record: FQRecord, start: int, length: int = -1): FQRecord =
  ## Extract a subsequence from a record
  ## 
  ## Args:
  ##   record: Input FQRecord
  ##   start: Start position (0-based)
  ##   length: Length of subsequence to extract (-1 for end of sequence)
  ## 
  ## Returns:
  ##   New FQRecord with extracted subsequence
  result = record
  let actualLength = if length < 0: record.sequence.len - start else: length
  let endPos = min(start + actualLength, record.sequence.len)
  
  if start >= record.sequence.len or actualLength <= 0:
    result.sequence = ""
    result.quality = ""
    return result
    
  result.sequence = record.sequence[start ..< endPos]
  if record.quality.len > 0:
    result.quality = record.quality[start ..< endPos]

proc composition*(record: FQRecord): SeqComp =
  ## Calculate composition of a DNA Sequence
  ## Returns a SeqComp object with counts of A, C, G, T, N, and Other (int) and GC content (float)
  
  for c in record.sequence:
    case c:
      of 'A', 'a': inc result.A
      of 'C', 'c': inc result.C
      of 'G', 'g': inc result.G
      of 'T', 't': inc result.T
      of 'N', 'n': inc result.N
      else: inc result.Other
  
  result.GC = float(result.G + result.C) / float(record.sequence.len - result.N - result.Other)



proc gcContent*(sequence: string): float =
  ## Calculate GC content of a DNA sequence
  ## 
  ## Args:
  ##   sequence: DNA sequence
  ## 
  ## Returns:
  ##   GC content as a fraction between 0.0 and 1.0
  if sequence.len == 0:
    return 0.0
    
  var gcCount = 0
  for c in sequence:
    if c in {'G', 'g', 'C', 'c'}:
      inc gcCount
      
  return gcCount / sequence.len

proc gcContent*(record: FQRecord): float =
  ## Calculate GC content of a FQRecord (using its sequence)
  return gcContent(record.sequence)

proc maskLowQuality*(record: var FQRecord, minQual: int, offset: int = 33, maskChar: char = 'N') =
  ## Mask sequence positions with low quality scores
  ## 
  ## Args:
  ##   record: FQRecord to modify
  ##   minQual: Minimum quality value
  ##   offset: Quality score offset
  ##   maskChar: Character to use for masking
  if record.quality.len == 0 or record.quality.len != record.sequence.len:
    return # FASTA has no quality, or invalid record
    
  for i in 0 ..< record.quality.len:
    if ord(record.quality[i]) - offset < minQual:
      record.sequence[i] = maskChar

#====