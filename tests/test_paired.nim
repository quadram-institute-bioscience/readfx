import unittest, ../readfx

suite "Paired-end reading tests":
  
  test "readFQPair basic functionality (uncompressed synthetic data)":
    var count = 0
    var totalLen1, totalLen2 = 0
    
    for pair in readFQPair("./tests/test_R1.fq", "./tests/test_R2.fq"):
      count += 1
      # this dataset is made so that the total length of both reads is 10
      check pair.read1.sequence.len  + pair.read2.sequence.len == 10
    
      
      check pair.read1.name.len > 0
      check pair.read2.name.len > 0
      check pair.read1.sequence.len > 0
      check pair.read2.sequence.len > 0
      check pair.read1.quality.len > 0
      check pair.read2.quality.len > 0
      
  test "readFQPair basic functionality":
    var count = 0
    var totalLen1, totalLen2 = 0
    
    for pair in readFQPair("./tests/illumina_1.fq.gz", "./tests/illumina_2.fq.gz"):
      count += 1
      totalLen1 += pair.read1.sequence.len
      totalLen2 += pair.read2.sequence.len
      
      check pair.read1.name.len > 0
      check pair.read2.name.len > 0
      check pair.read1.sequence.len > 0
      check pair.read2.sequence.len > 0
      check pair.read1.quality.len > 0
      check pair.read2.quality.len > 0
      
      # Check that quality and sequence lengths match
      check pair.read1.sequence.len == pair.read1.quality.len
      check pair.read2.sequence.len == pair.read2.quality.len
    
    check count == 7  # Should have 7 read pairs
    check totalLen1 > 0
    check totalLen2 > 0
  
  test "readFQPair name checking":
    # This should work without exceptions since the names match (base names)
    var count = 0
    for pair in readFQPair("./tests/illumina_1.fq.gz", "./tests/illumina_2.fq.gz", checkNames = true):
      count += 1
    check count == 7
  
  test "FQPair type properties":
    for pair in readFQPair("./tests/illumina_1.fq.gz", "./tests/illumina_2.fq.gz"):
      # Test that we can access all fields
      discard pair.read1.name
      discard pair.read1.comment
      discard pair.read1.sequence
      discard pair.read1.quality
      discard pair.read2.name
      discard pair.read2.comment
      discard pair.read2.sequence
      discard pair.read2.quality
      break  # Just test first pair
  
  test "Error handling for non-existent files":
    expect IOError:
      for pair in readFQPair("./tests/nonexistent1.fq", "./tests/nonexistent2.fq"):
        discard pair
  
  test "Error handling for single non-existent file":
    expect IOError:
      for pair in readFQPair("./tests/illumina_1.fq.gz", "./tests/nonexistent.fq"):
        discard pair