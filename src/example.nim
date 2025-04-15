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


  var newSeq = FQRecord()

  debugMessage "Testing Record"
  newSeq.name     = "name"
  newSeq.comment  = "comment"
  newSeq.sequence = "SEQUENCE"
  
  echo $newSeq

  debugMessage "Testing Record: with quality"
  newSeq.quality  = "Quality!"
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