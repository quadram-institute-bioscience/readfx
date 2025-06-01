
import strutils
import sequtils

proc matchIUPAC*(primerBase, referenceBase: char): bool =
  ## Test if a primer base matches a reference base using IUPAC nucleotide codes.
  ## 
  ## Performs case-insensitive matching between primer and reference bases,
  ## with full support for IUPAC ambiguous nucleotide codes in the primer.
  ## 
  ## Parameters:
  ## - `primerBase`: Base from primer sequence (supports IUPAC ambiguous codes)
  ## - `referenceBase`: Base from reference/target sequence (typically A,T,G,C,N)
  ## 
  ## Returns:
  ##   `true` if the bases are compatible according to IUPAC rules, `false` otherwise
  ## 
  ## Matching Rules:
  ## - Exact matches (A-A, T-T, etc.) always return `true`
  ## - `N` in primer matches any reference base
  ## - `N` in reference never matches (considered unknown/low quality)
  ## - IUPAC ambiguous codes in primer match their corresponding base sets:
  ##   - `Y` (pyrimidine): matches C, T
  ##   - `R` (purine): matches A, G  
  ##   - `S` (strong): matches G, C
  ##   - `W` (weak): matches A, T
  ##   - `K` (keto): matches T, G
  ##   - `M` (amino): matches A, C
  ##   - `B` (not A): matches C, G, T
  ##   - `D` (not C): matches A, G, T
  ##   - `H` (not G): matches A, C, T
  ##   - `V` (not T): matches A, C, G
  ## 
  ## Note:
  ##   Function is asymmetric - IUPAC codes are only interpreted in the primer,
  ##   not the reference. This reflects typical primer design workflows.
  ## 
  ## Example:
  ## ```nim
  ## assert matchIUPAC('A', 'A') == true    # Exact match
  ## assert matchIUPAC('Y', 'C') == true    # Y matches C,T
  ## assert matchIUPAC('Y', 'A') == false   # Y doesn't match A
  ## assert matchIUPAC('N', 'G') == true    # N matches anything
  ## assert matchIUPAC('A', 'N') == false   # N in reference = no match
  ## ```
  let a = primerBase.toupperascii()
  let b = referenceBase.toupperascii()
  # a=primer; b=read
  if b == 'N': return false
  if a == b or a == 'N': return true
  
  # Direct case matching for IUPAC codes
  case a:
    of 'Y': b in {'C', 'T'}
    of 'R': b in {'A', 'G'}
    of 'S': b in {'G', 'C'}
    of 'W': b in {'A', 'T'}
    of 'K': b in {'T', 'G'}
    of 'M': b in {'A', 'C'}
    of 'B': b != 'A'  # Everything except A
    of 'D': b != 'C'  # Everything except C
    of 'H': b != 'G'  # Everything except G
    of 'V': b != 'T'  # Everything except T
    else: false


proc findOligoMatches*(sequence, primer: string, threshold: float, max_mismatches = 0, min_matches = 6): seq[int] =
  ## Find all binding sites for an oligonucleotide primer in a DNA sequence.
  ## 
  ## Uses IUPAC-aware matching to identify potential primer binding sites with
  ## configurable mismatch tolerance and scoring thresholds. Handles primer
  ## overhangs by padding the target sequence.
  ## 
  ## Parameters:
  ## - `sequence`: Target DNA sequence to search
  ## - `primer`: Oligonucleotide primer sequence (supports IUPAC ambiguous bases)
  ## - `threshold`: Minimum match score (matches/effective_primer_length, 0.0-1.0)
  ## - `max_mismatches`: Maximum allowed mismatches (default: 0)
  ## - `min_matches`: Minimum required matching bases (default: 6)
  ## 
  ## Returns:
  ##   Sequence of starting positions where primer binds. Positions are relative
  ##   to the original sequence:
  ##   - Negative values: primer overhangs before sequence start
  ##   - Values ≥ sequence.len: primer overhangs past sequence end
  ##   - 0 to sequence.len-1: primer starts within sequence bounds
  ## 
  ## Scoring:
  ## - Match score = (matching_bases) / (primer_length_excluding_gap_padding)
  ## - Only non-gap positions count toward primer length
  ## - Early termination if mismatches exceed `max_mismatches`
  ## 
  ## IUPAC Support:
  ##   Handles ambiguous nucleotide codes (Y, R, S, W, K, M, B, D, H, V, N)
  ## 
  ## Example:
  ## ```nim
  ## let positions = findOligoMatches("ATCGATCG", "ATCG", 0.75, 1, 3)
  ## # Returns positions where ATCG matches with ≥75% identity
  ## ```
  let
    dna = ('-'.repeat(len(primer) - 1) & sequence & '-'.repeat(len(primer) - 1)).toUpper()
    primer = primer.toUpper()

  for pos in 0..len(dna)-len(primer):
    let query = dna[pos..<pos+len(primer)]
    var
      matches = 0
      mismatches = 0
      primerRealLen = 0

    for c in 0..<len(query):
      if matchIUPAC(primer[c], query[c]):
        matches += 1
        primerRealLen += 1
      elif query[c] != '-':
        mismatches += 1
        primerRealLen += 1

      if mismatches > max_mismatches:
        break
 
    let
      score = float(matches) / float(primerRealLen)
    if score >= threshold and mismatches <= max_mismatches and matches >= min_matches:
      result.add(pos-len(primer)+1)

proc findPrimerMatches*(sequence, primer: string, threshold: float, max_mismatches = 0, min_matches = 6): seq[seq[int]] =
  ## Find primer binding sites on both strands of a DNA sequence.
  ## 
  ## Searches for matches of the primer and its reverse complement against the target
  ## sequence, allowing for IUPAC ambiguous bases and configurable matching criteria.
  ## 
  ## Parameters:
  ## - `sequence`: Target DNA sequence to search
  ## - `primer`: Primer sequence (supports IUPAC codes)
  ## - `threshold`: Minimum match score (matches/primer_length, 0.0-1.0)
  ## - `max_mismatches`: Maximum allowed mismatches (default: 0)
  ## - `min_matches`: Minimum required matching bases (default: 6)
  ## 
  ## Returns:
  ##   A sequence containing two elements:
  ##   - `result[0]`: Forward strand match positions
  ##   - `result[1]`: Reverse strand match positions
  ## 
  ## Position values are relative to the input sequence:
  ## - Negative positions indicate primer overhangs before sequence start
  ## - Positions ≥ sequence length indicate overhangs past sequence end
  ## 
  ## Example:
  ## ```nim
  ## let matches = findPrimerMatches("ATCGATCG", "ATCG", 0.75)
  ## # matches[0] = forward hits, matches[1] = reverse complement hits
  ## ```
  let
    forMatches    = findOligoMatches(sequence, primer, threshold, max_mismatches, min_matches)
    primerReverse = revcompl(primer)
    revMatches    = findOligoMatches(sequence, primerReverse, threshold, max_mismatches, min_matches)

  result = @[forMatches, revMatches]
