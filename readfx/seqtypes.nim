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
  SeqComp* = object
    A*: int
    C*: int
    G*: int
    T*: int
    GC*: float
    N*: int
    Other*: int
  StrandDirection* = enum
    Forward, Reverse
