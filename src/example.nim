import os

import ../readfx
import strformat, times

var messageCounter = 0

proc debugMessage*(message: string, indent: bool = false) =
  let prefix = if indent: "##" 
              else: "#"
  # Increment the counter only if not indented
  if not indent:
    inc messageCounter
  
  # Get current timestamp
  let timestamp = now().format("yyyy-MM-dd HH:mm:ss")
  
  # ANSI color codes
  const 
    reset = "\e[0m"
    bold = "\e[1m"
    green = "\e[32m"
    yellow = "\e[33m"
    cyan = "\e[36m"
  
  # Indentation formatting
  let 
    indentStr = if indent: "    " else: ""
    counterStr = if indent: "  " else: fmt"[ {messageCounter} ]"
  
  # Format and print the message with counter, timestamp and colors
  if indent == false:
    echo fmt"{reset}{cyan}---------------{cyan}{timestamp}{reset}--------------- {reset}"
  echo fmt"{bold}{yellow}{prefix}  {counterStr}  {bold}{green}{message}{reset}{bold}{yellow} {reset}"


when isMainModule:
  let args = commandLineParams()
  # Check if the number of arguments is less than 1
  if args.len < 1:
    stderr.writeLine "Missing input parameter [FILENAME]"
    quit(1)

  var newSeq = FQRecord()

  debugMessage "Testing Record"
  newSeq.name     = "NAME_1"
  newSeq.comment  = "this is a comment"
  newSeq.sequence = "CAGATATATATATATATATATATATATATAT"
  
  echo "Outputting a FASTA record:"
  echo $newSeq
  echo "--- split by 10 ---"
  echo fafmt(newSeq, 10)

  debugMessage "Testing Record: with quality"
  newSeq.quality  = "IIIIBIHHHEEEEEEHHHHHHHHHHHA9987"
  echo $newSeq
  if len(args) == 0:
    stderr.writeLine "Missing input parameter [FILENAME]"
    quit(1)


  debugMessage "Testing `readfq` for each input file"
  for path in args:
    debugMessage("Running readfq: " & path, true)
    var 
      c = 0
      t = 0
    for rec in readfq(path):
      c += 1
      t += len(rec.sequence)
    echo path, "\t", c, "\t", t
    c = 0
    t = 0
    debugMessage("Running readFQPtr: " & path, true)
    for rec in readFQPtr(path):
      c += 1
      t += len(rec.sequence)
    echo path, "\t", c, "\t", t


    debugMessage "Testing `readFastx` for each input file"
    for path in args:
      debugMessage("Running readFastx: " & path, true)  
      c = 0
      t = 0
      var r: FQRecord
      var f = xopen[GzFile](path)
      defer: f.close()
      while f.readFastx(r):
        c += 1
        t += len(r.sequence)
      echo path, "\t", c, "\t", t