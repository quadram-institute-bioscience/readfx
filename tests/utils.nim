import unittest
import os
import strutils

import ../readfx


  
test "input files":
  # Check if the test files exist
  check fileExists("./tests/fastq_demo.fq")
  check fileExists("./tests/seq.txt")
  check fileExists("./tests/SRR396637_1.seqs1-2.fastq.gz")
  check fileExists("./tests/test.fasta.gz")
  check fileExists("./tests/fasta_demo.fa")


test "matchIUPAC: match":

  # check matches
  check matchIUPAC('A', 'A')
  check matchIUPAC('C', 'C')
  check matchIUPAC('g', 'G')
  check matchIUPAC('t', 't')

  # N in primer
  check matchIUPAC('N', 'A')
  check matchIUPAC('N', 'g')


test "matchIUPAC: mismatch":
  # check matches
  check matchIUPAC('A', 'C') == false
  check matchIUPAC('C', 'T') == false
  check matchIUPAC('g', 't') == false
  check matchIUPAC('t', 'N') == false

  # N in reference
  check matchIUPAC('A', 'N') == false
  check matchIUPAC('a', 'n') == false


test "findOligoMatches: exact match non degenerate":

  let target =    "GCGTACGATCGTACGTACAGCTGATCGTACTGCTAGCTGTC"
  let primer1 =   "GCGTACGAT"
  let primer2 = "GAGCGTACGAT"

  let matches1 = findOligoMatches(target, primer1, 0.8, max_mismatches=1, min_matches=6)
  let matches2 = findOligoMatches(target, primer2, 0.8, max_mismatches=1, min_matches=6)

  check matches1.len == 1
  check matches1[0] == 0

  check matches2.len == 1
  check matches2[0] == -2

test "findOligoMatches: exact match degenerate":

  let target =     "GCGTACGATCGTACGTACAGCTGATCGTACTGCTAGCTGTC"
  let primer1 =    "NCGWACSAT"
  let primer2 = "GWANCGWACSAT"

  let matches1 = findOligoMatches(target, primer1, 0.9, max_mismatches=0, min_matches=6)
  let matches2 = findOligoMatches(target, primer2, 0.9, max_mismatches=0, min_matches=6)

  check matches1.len == 1
  check matches1[0] == 0

  check matches2.len == 1
  check matches2[0] == -3