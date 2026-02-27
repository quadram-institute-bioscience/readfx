import os
import re
import docopt
import readfx

# Version information for the program
proc version(): string = "1.0.0"

# Print a FASTQ record to file or stdout
proc print_seq(record: FQRecord, outFile: File = stdout) =
  ## Outputs a FASTQ record in standard format
  if outFile == stdout:
    echo "@", record.name
    echo record.sequence  
    echo "+"
    echo record.quality
  else:
    outFile.writeLine("@", record.name)
    outFile.writeLine(record.sequence)
    outFile.writeLine("+") 
    outFile.writeLine(record.quality)

# Main program logic
proc main() =
  var argv = commandLineParams()
  
  let args = docopt("""
ilv: interleave FASTQ files

  Usage: ilv [options] -1 <forward-pair> [-2 <reverse-pair>]

  -f --for-tag <tag-1>       string identifying forward files [default: auto]
  -r --rev-tag <tag-2>       string identifying reverse files [default: auto]
  -o --output <outputfile>   save file to <out-file> instead of STDOUT
  -c --check                 enable careful mode (check sequence names and numbers)
  -v --verbose               print verbose output

  -s --strip-comments        skip comments
  -p --prefix "string"       rename sequences (append a progressive number)

guessing second file:
  by default <forward-pair> is scanned for _R1. and substitute with _R2.
  if this fails, the patterns _1. and _2. are tested.

example:

    ilv -1 file_R1.fq > interleaved.fq
  
  """, version=version(), argv=argv)

  var
    file_R1 = $args["<forward-pair>"]
    file_R2: string
    pattern_R1  = $args["--for-tag"]
    pattern_R2  = $args["--rev-tag"]
    output_file = $args["--output"]
    prefix      = $args["--prefix"]

  # Parse boolean flags from docopt
  let check = args["--check"]
  let stripComments = args["--strip-comments"]
  let verbose = args["--verbose"]
    
  # Get R2 file if provided
  if args["<reverse-pair>"]:
    file_R2 = $args["<reverse-pair>"]

  # Auto-detect R2 file if not provided
  if file_R2 == "":
    if pattern_R1 == "auto" and pattern_R2 == "auto":
        # Automatic pattern detection
        if match(file_R1, re".+_R1\..+"):           
            file_R2 = file_R1.replace(re"_R1\.", "_R2.")
        elif match(file_R1, re".+_1\..+"):            
            file_R2 = file_R1.replace(re"_1\.", "_2.")
        else:
            echo "Unable to detect --for-tag (_R1. or _1.) in <", file_R1, ">"
            quit(1)
    else:
        # User-defined patterns
        if match(file_R1, re(".+" & pattern_R1 & ".+") ):
            file_R2 = file_R1.replace(re(pattern_R1), pattern_R2)
        else:
            echo "Unable to find pattern <", pattern_R1, "> in file <", file_R1, ">"
            quit(1)

  # Setup output file
  var outFile: File
  if output_file != "nil":
      outFile = open(output_file, fmWrite)
      
  # Print verbose information
  if verbose:
    stderr.writeLine("- file1:\t", file_R1)
    stderr.writeLine("- file2:\t", file_R2)
    stderr.writeLine("- patterns:\t", pattern_R1,';',pattern_R2)
    if output_file != "nil":
        stderr.writeLine("- output:\t", output_file)

  # Validate input files
  if file_R1 == file_R2:
    echo "FATAL ERROR: First file and second file are equal."
    quit(1)

  if not fileExists(file_R1):
    echo "FATAL ERROR: First pair (", file_R1 , ") not found."
    quit(1)

  if not fileExists(file_R2):
    echo "FATAL ERROR: Second pair (", file_R2 , ") not found."
    quit(1)

  # Open FASTQ files using readfx
  var fq1 = xopen[GzFile](file_R1)
  defer: fq1.close()
  var fq2 = xopen[GzFile](file_R2)
  defer: fq2.close()

  # Process records
  var R1: FQRecord
  var R2: FQRecord
  var c = 0

  # Read and interleave sequences
  while fq1.readFastx(R1):
    c += 1
    if not fq2.readFastx(R2):
      stderr.writeLine("File R2 ended prematurely after ", c, " sequences.")
      if check:
        quit(1)
      else: 
        quit(0)

    # Check sequence name consistency if requested
    if check and R1.name != R2.name:
        echo "Sequence error [seq ", c, "], name mismatch"
        echo R1.name, " != ", R2.name
        quit(3)

    # Apply prefix if specified
    if prefix != "nil":
        R1.name = prefix & $c
        R2.name = prefix & $c
        
    # Output interleaved sequences
    echo R1
    echo R2

  # Check if R2 has remaining sequences
  if fq2.readFastx(R2):
    stderr.writeLine("File R1 ended prematurely after ", c, " sequences.")
    if check:
      quit(1)
    else: 
      quit(0)

  # Print summary if verbose
  if verbose:
    stderr.writeLine("Printed ", c, " sequence pairs from ", file_R1, " and ", file_R2)

  # Close output file if needed
  if output_file != "nil":
    outFile.close()

# Run the main program
when isMainModule:
  main()
