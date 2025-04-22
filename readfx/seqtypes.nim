type
  # convenience type for FastQ or Fasta records
  FQRecordPtr* = object
    name*: ptr char
    comment*: ptr char# optional
    sequence*: ptr char
    quality*: ptr char# optional
  FQRecord* = object
    name*: string
    comment*: string# optional
    sequence*: string
    quality*: string# optional
    status*, lastChar*: int
  StrandDirection* = enum
    Forward, Reverse
