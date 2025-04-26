## This module defines core data structures for bioinformatics sequence handling:
## * Record types for FASTA/FASTQ format sequence data
## * Nucleotide composition analysis
## * Strand orientation enumeration
##
## These types are designed to efficiently represent and manipulate biological
## sequence data in memory with both pointer-based and string-based implementations.


## FQRecordPtr
## ==========================================================
## 
## **Pointer-based representation of FASTA/FASTQ records**
##
## This structure uses raw pointers for memory efficiency when processing
## large sequence files. The pointers are typically managed by the kseq
## library and should not be manipulated directly.
##
## Note: These pointers are typically volatile and may be invalidated
## when reading the next record. Use FQRecord when you need to retain data.
##
## Example:
##
## ```nim
## for rec in readFQPtr("sample.fastq"):
##   # Access data through pointers
##   echo $rec.name       # Convert ptr char to string
##   echo $rec.sequence
## ```

type
  FQRecordPtr* = object
    name*: ptr char      ## Sequence name/identifier (null-terminated)
    comment*: ptr char   ## Optional sequence description/comment (null-terminated)
    sequence*: ptr char  ## Nucleotide sequence (null-terminated)
    quality*: ptr char   ## Optional quality scores (null-terminated, empty for FASTA)


## FQRecord
## ==========================================================
## **String-based representation of FASTA/FASTQ records**
## This structure uses Nim strings to represent sequence data, making it
## more convenient but potentially less memory efficient than FQRecordPtr.
## Use this type when you need to retain sequence data across iterations.
##
## Fields:
##   name: Sequence identifier (without '>' or '@' prefix)
##   comment: Optional description/comment after sequence name
##   sequence: Nucleotide sequence
##   quality: Quality scores (empty for FASTA format)
##   status: Status code from the last read operation:
##     - >0: Length of sequence for normal records
##     - -1: End of file
##     - -2: Stream error
##     - -3: Other parsing error
##     - -4: Sequence and quality length mismatch
##   lastChar: The last character read from the file (used internally)
##
## Example:
##
## ```nim
## for rec in readFQ("sample.fastq"):
##   echo rec.name
##   echo rec.sequence
##   if rec.quality.len > 0:
##     echo "Average quality: ", rec.quality.mapIt(ord(it) - 33).sum / rec.quality.len
## ```
type
  FQRecord* = object
    name*: string         ## Sequence name/identifier 
    comment*: string      ## Optional sequence description/comment
    sequence*: string     ## Nucleotide sequence
    quality*: string      ## Optional quality scores (empty for FASTA)
    status*, lastChar*: int  ## Status code and internal parsing state

## SeqComp
## =========================================================
## **Nucleotide composition statistics**
##
## Counts of different nucleotides in a sequence and calculated GC content.
## Useful for sequence analysis and quality control.
##
## Example:
##
## ```nim
## proc calculateComposition(seq: string): SeqComp =
##   result = SeqComp()
##   for base in seq:
##     case base.toUpperAscii:
##     of 'A': result.A += 1
##     of 'C': result.C += 1
##     of 'G': result.G += 1
##     of 'T': result.T += 1
##     of 'N': result.N += 1
##     else: result.Other += 1
##   
##   let total = result.A + result.C + result.G + result.T + result.N + result.Other
##   if total > 0:
##     result.GC = (result.G + result.C) / total.float
##
## let composition = calculateComposition("ACGTACGTNNN")
## echo "GC content: ", composition.GC
## ```
type
  SeqComp* = object
    A*: int        ## Count of adenine (A) nucleotides
    C*: int        ## Count of cytosine (C) nucleotides
    G*: int        ## Count of guanine (G) nucleotides
    T*: int        ## Count of thymine (T) nucleotides
    GC*: float     ## GC content (ratio of G+C to total nucleotides)
    N*: int        ## Count of ambiguous (N) nucleotides
    Other*: int    ## Count of other characters (non-ACGTN)

  ## DNA strand orientation
  ##
  ## Represents the directionality of a DNA sequence.
  ## In molecular biology, DNA is double-stranded with complementary bases.
  ## The two strands run in opposite directions (5' to 3' and 3' to 5').
  ##
  ## Example:
  ##