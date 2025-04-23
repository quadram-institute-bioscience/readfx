
## API Reference

> This page is incomplete and a placeholder for the API reference. It will be updated with more details and examples.

### C Wrapper (klib)

#### Types

- `FQRecord*`: Main type for FASTQ/FASTA records with string fields
  - `name*`: Sequence name/ID (string)
  - `comment*`: Optional comment (string)
  - `sequence*`: The sequence (string)
  - `quality*`: Quality string for FASTQ (string)
  - `status`, `lastChar`: Internal fields

- `FQRecordPtr*`: Pointer-based version for more efficient memory usage
  - `name*`: Sequence name/ID (ptr char)
  - `comment*`: Optional comment (ptr char)
  - `sequence*`: The sequence (ptr char)
  - `quality*`: Quality string for FASTQ (ptr char)

#### Functions and Iterators

- `iterator readFQ*(path: string): FQRecord`: 
  Parses FASTQ/FASTA records from a file, converting to strings

- `iterator readFQPtr*(path: string): FQRecordPtr`: 
  Memory-efficient version using pointers (faster but requires careful handling)

- `proc `$`*(rec: FQRecord): string`: 
  Formats a record as a FASTA/FASTQ string

- `proc `$`*(rec: FQRecordPtr): string`: 
  Formats a pointer-based record as a FASTA/FASTQ string

### Native Nim Implementation (nimklib)

#### Types

- `FQRecord*`: Record type for FASTQ/FASTA data
  - `sequence`, `quality`, `name`, `comment`: String fields
  - `status`, `lastChar`: Internal fields

- `Bufio*[T]`: Buffered reader for efficient file reading
  - Generic over file type to support both regular and gzipped files

#### Functions

- `proc readFastx*[T](f: var Bufio[T], r: var FQRecord): bool`: 
  Parses a single FASTQ/FASTA record from a buffer

- `proc xopen*[T](fn: string, mode: FileMode = fmRead, sz: int = 0x10000): Bufio[T]`: 
  Opens a file with buffered reading

- `proc close*[T](f: var Bufio[T]): int`: 
  Closes a buffered file

- `proc eof*[T](f: Bufio[T]): bool`: 
  Checks if the file has reached EOF

- `proc readLine*[T](f: var Bufio[T], buf: var string): bool`: 
  Reads a line from a buffered file

### Additional Utilities

- `GzFile*`: Type for handling gzipped files
- `Interval*[S,T]`: Generic interval type with interval tree operations
- Various buffer manipulation functions
