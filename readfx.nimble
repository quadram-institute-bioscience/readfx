# Package

version       = "0.2.1"
author        = "Andreas Wilm and SeqFu team"
description   = "Parse FASTQ and FASTA files, using Heng Li's Klib"
license       = "MIT"

requires "nim >= 1.0", "zip >= 0.2.1"

skipDirs = @["tests"]

# Build the library?

task test, "run the tests":
  exec "nim c -r tests/tester"

# DOCUMENTATION
# ============================================================================
# Get Nim version by executing a nim command
proc getNimVersionStr(): string =
  # Run nim -v and get the first line
  let (output, _) = gorgeEx("nim -v")
  let firstLine = output.splitLines()[0]
  # Extract version number
  let parts = firstLine.split(" ")
  for part in parts:
    if part[0].isDigit:
      return part
  return ""

# Parse version string into a sequence of integers
proc parseVersionStr(version: string): seq[int] =
  result = @[]
  for part in version.split("."):
    try:
      result.add(parseInt(part))
    except ValueError:
      # Handle cases like "1.6.0-rc1" by ignoring non-numeric parts
      let numPart = part.split("-")[0]
      if numPart.len > 0:
        result.add(parseInt(numPart))
      else:
        result.add(0)
  # Ensure we have at least 3 components (major, minor, patch)
  while result.len < 3:
    result.add(0)

# Compare version numbers
proc versionAtLeast(version: seq[int], major, minor, patch: int): bool =
  if version.len < 3:
    return false
  if version[0] > major:
    return true
  if version[0] == major and version[1] > minor:
    return true
  if version[0] == major and version[1] == minor and version[2] >= patch:
    return true
  return false

# Get current Nim version
let nimVersionStr = getNimVersionStr()
let nimVersion = parseVersionStr(nimVersionStr)
let isNim2OrLater = versionAtLeast(nimVersion, 2, 0, 0)
task docs, "Generate documentation":
 
  if isNim2OrLater:
    exec "nim doc --git.url:https://github.com/quadram-institute-bioscience/readfx/ --git.commit:main --index:on --project --out:docs ./readfx.nim"
