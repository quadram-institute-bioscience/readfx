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

test "findOligoMatches: empty primer returns no matches":
  let matches = findOligoMatches("ATCG", "", 0.8, max_mismatches=0, min_matches=0)
  check matches.len == 0

test "findPrimerMatches: searches both strands":
  let target = "GCGTACGATCGTACGTACAGCTGATCGTACTGCTAGCTGTC"
  let matches = findPrimerMatches(target, "GCGTACGAT", 0.8, max_mismatches=1, min_matches=6)
  check matches.len == 2      # forward and reverse strand result sets
  check matches[0].len >= 1   # forward hit

test "utils: quality char/int conversion":
  check qualCharToInt('I') == 40
  check qualCharToInt('!') == 0
  check qualIntToChar(40) == 'I'
  check qualIntToChar(0) == '!'
  # custom offset (Illumina 1.3+)
  check qualCharToInt('h', 64) == 40
  check qualIntToChar(40, 64) == 'h'

test "utils: avgQuality()":
  var r = FQRecord()
  r.sequence = "ACGT"
  r.quality = "IIII"          # 'I' = Q40
  check avgQuality(r) == 40.0
  # string overload: '!' = Q0, 'I' = Q40
  check avgQuality("!I") == 20.0
  # invalid inputs return -1
  check avgQuality("") == -1.0
  r.quality = "II"            # length mismatch
  check avgQuality(r) == -1.0

test "utils: revCompl() string and IUPAC codes":
  check revCompl("ATGC") == "GCAT"
  check revCompl("atgc") == "gcat"   # case is preserved
  check revCompl("AtGcN") == "NgCaT" # mixed case preserved per base
  check revCompl("AU") == "AT"       # U complements to A
  check revCompl("ARY") == "RYT"     # degenerate bases are complemented
  check revCompl("ryswkmbvdh") == "dhbvkmwsry" # lowercase IUPAC codes
  check revCompl("ACGT-") == "-ACGT" # unknown characters pass through
  # double application is the identity
  check revCompl(revCompl("AcGtNrY-")) == "AcGtNrY-"

test "utils: revComplInPlace()":
  # FASTQ: sequence complemented, quality reversed (never complemented)
  block:
    var s = "GAAA"
    var q = "IIIA"
    revComplInPlace(s, q)
    check s == "TTTC"
    check q == "AIII"
  # odd length: middle base complemented in place
  block:
    var s = "ACGTN"
    var q = "12345"
    revComplInPlace(s, q)
    check s == "NACGT"
    check q == "54321"
  # case preserved per base
  block:
    var s = "AtGcN"
    var q = "abcde"
    revComplInPlace(s, q)
    check s == "NgCaT"
    check q == "edcba"
  # FASTA: empty quality is allowed
  block:
    var s = "AAGT"
    var q = ""
    revComplInPlace(s, q)
    check s == "ACTT"
    check q == ""
  # single base and empty sequence
  block:
    var s = "c"
    var q = "!"
    revComplInPlace(s, q)
    check s == "g"
    check q == "!"
  block:
    var s = ""
    var q = ""
    revComplInPlace(s, q)
    check s == ""
  # double application is the identity
  block:
    var s = "AcGtNrY"
    var q = "1234567"
    let origS = s
    let origQ = q
    revComplInPlace(s, q)
    revComplInPlace(s, q)
    check s == origS
    check q == origQ

test "utils: gcContent() over valid bases":
  check gcContent("ATGC") == 0.5
  check gcContent("GGGG") == 1.0
  check gcContent("AAAA") == 0.0
  check gcContent("acgtn") == 0.5   # case-insensitive
  # N and other symbols are excluded from the denominator
  check gcContent("GGN") == 1.0
  check gcContent("ACGTN") == 0.5
  check gcContent("ACGTXX") == 0.5
  # no valid bases -> 0.0
  check gcContent("") == 0.0
  check gcContent("NNN") == 0.0
  # consistent with composition() (canonical definition)
  var r = FQRecord()
  r.sequence = "ACGTNXX"
  check gcContent(r.sequence) == composition(r).GC
  check gcContent(r) == composition(r).GC

test "utils: subSequence()":
  var r = FQRecord()
  r.name = "t"
  r.sequence = "AACCGGTT"
  r.quality = "12345678"
  # middle slice
  var s = subSequence(r, 2, 4)
  check s.sequence == "CCGG"
  check s.quality == "3456"
  # from position to end
  s = subSequence(r, 4)
  check s.sequence == "GGTT"
  check s.quality == "5678"
  # negative start counts from the end
  s = subSequence(r, -2)
  check s.sequence == "TT"
  check s.quality == "78"
  # out of range yields empty record
  s = subSequence(r, 100)
  check s.sequence == ""
  check s.quality == ""

test "utils: trimQuality() and qualityTrim()":
  check trimQuality("III!!!", 20) == "III"
  check trimQuality("!!!", 20) == ""
  var r = FQRecord()
  r.sequence = "AACCGG"
  r.quality = "IIII!!"
  qualityTrim(r, 20)
  check r.sequence == "AACC"
  check r.quality == "IIII"
  # FASTA record (no quality) is left untouched
  var fa = FQRecord()
  fa.sequence = "AACCGG"
  qualityTrim(fa, 20)
  check fa.sequence == "AACCGG"

test "interval tree: index and overlap":
  var intervals: seq[Interval[int, string]] = @[
    (st: 10, en: 20, data: "a", max: 0),
    (st: 15, en: 25, data: "b", max: 0),
    (st: 30, en: 40, data: "c", max: 0),
    (st: 50, en: 60, data: "d", max: 0)
  ]
  index(intervals)

  var hits: seq[string]
  for iv in overlap(intervals, 12, 18):
    hits.add(iv.data)
  check hits == @["a", "b"]   # overlapping hits, yielded in sorted order

  hits.setLen(0)
  for iv in overlap(intervals, 35, 45):
    hits.add(iv.data)
  check hits == @["c"]

  hits.setLen(0)
  for iv in overlap(intervals, 26, 29):
    hits.add(iv.data)
  check hits.len == 0         # query in a gap

  hits.setLen(0)
  for iv in overlap(intervals, 0, 100):
    hits.add(iv.data)
  check hits == @["a", "b", "c", "d"]

test "interval tree: unsorted input is sorted by index":
  var intervals: seq[Interval[int, string]] = @[
    (st: 50, en: 60, data: "d", max: 0),
    (st: 10, en: 20, data: "a", max: 0),
    (st: 30, en: 40, data: "c", max: 0)
  ]
  index(intervals)
  check intervals[0].data == "a"
  var hits: seq[string]
  for iv in overlap(intervals, 55, 70):
    hits.add(iv.data)
  check hits == @["d"]

test "interval tree: edge cases":
  # empty input
  var empty: seq[Interval[int, string]] = @[]
  check index(empty) == 0
  var n = 0
  for iv in overlap(empty, 0, 100):
    inc n
  check n == 0

  # single interval
  var single: seq[Interval[int, string]] = @[(st: 5, en: 10, data: "x", max: 0)]
  index(single)
  n = 0
  for iv in overlap(single, 0, 100):
    check iv.data == "x"
    inc n
  check n == 1

  # half-open boundaries: touching at st/en does not overlap
  var ivs: seq[Interval[int, string]] = @[(st: 10, en: 20, data: "a", max: 0)]
  index(ivs)
  n = 0
  for iv in overlap(ivs, 20, 30): inc n   # query starts exactly at en
  check n == 0
  n = 0
  for iv in overlap(ivs, 0, 10): inc n    # query ends exactly at st
  check n == 0
  n = 0
  for iv in overlap(ivs, 19, 20): inc n   # overlaps the last base
  check n == 1

  # nested intervals are all reported
  var nested: seq[Interval[int, string]] = @[
    (st: 0, en: 100, data: "outer", max: 0),
    (st: 10, en: 20, data: "inner", max: 0)
  ]
  index(nested)
  var hits: seq[string]
  for iv in overlap(nested, 15, 16):
    hits.add(iv.data)
  check hits.len == 2

test "utils: fafmt() wrapping":
  var r = FQRecord()
  r.name = "seq1"
  r.comment = "a comment"
  r.sequence = "ACGTACGTAC"   # 10 bp
  check fafmt(r, 4) == ">seq1 a comment\nACGT\nACGT\nAC"
  check fafmt(r, 60) == ">seq1 a comment\nACGTACGTAC"
  r.comment = ""
  check fafmt(r, 5) == ">seq1\nACGTA\nCGTAC"

test "utils: revCompl() copy overload keeps original":
  # Note: the copy overload only resolves for immutable bindings —
  # with a mutable var, `revCompl(r)` picks the in-place var overload.
  let r = FQRecord(name: "rec", sequence: "GAAA", quality: "IIIA")
  let rc = revCompl(r)
  check rc.sequence == "TTTC"
  check rc.quality == "AIII"
  check rc.name == "rec"
  # original record is untouched
  check r.sequence == "GAAA"
  check r.quality == "IIIA"
