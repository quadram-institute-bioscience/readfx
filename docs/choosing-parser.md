# Choosing a Parser in ReadFX

This guide helps you choose between `readFQ`, `readFQPtr`, `readFastx`, and `readFQPair`.

## Short answer

- Use `readFQ` by default.
- Use `readFQPtr` when profiling shows string allocation is a bottleneck.
- Use `readFastx` when you need full control over stream handling and parser status.
- Use `readFQPair` for paired-end FASTQ files.

## Decision matrix

| Need | Best choice | Why |
|---|---|---|
| Safe, simple single-file parsing | `readFQ` | Returns `FQRecord` with Nim strings; easy to keep records after iteration. |
| Maximum single-file throughput | `readFQPtr` | Pointer-based records avoid per-record string allocations. |
| Custom I/O loop and explicit parse status | `readFastx` | You control `xopen/close`, loop behavior, and can inspect `record.status`. |
| Synchronized R1/R2 parsing | `readFQPair` | Reads two files in lockstep and yields `FQPair`. |
| Validate paired names while reading | `readFQPair(..., checkNames = true)` | Normalizes `/1` `/2` and ` 1` ` 2` suffixes before comparison. |

## What each method optimizes for

### `readFQ(path)`

Best default for most applications.

- Tradeoff: easiest API, but extra allocations/copies per record.
- Record lifetime: safe to store and pass around.
- Implementation note: wraps `readFQPtr` and converts pointer fields to strings.

Use when:

- You are building tools where maintainability matters more than maximum raw speed.
- You need to keep records in memory after the loop.

### `readFQPtr(path)`

Highest-throughput streaming API.

- Tradeoff: fastest and allocation-light, but pointer lifetime is short.
- Record lifetime: pointers are reused on each iteration.
- Rule: convert fields immediately (`$record.name`, `$record.sequence`) if you must keep them.

Use when:

- You process very large files in a tight loop.
- You can consume each record immediately and avoid storing raw pointers.

Avoid when:

- You need to retain record data without explicit copying.

### `readFastx(f, record)` with `xopen`

Low-level API for custom workflows.

- Tradeoff: more boilerplate, but full control of stream lifecycle and loop behavior.
- Record lifetime: you reuse one `FQRecord` object and overwrite it each read.
- Error signaling: on failure/EOF, `record.status` carries codes (`-1` EOF, `-2` stream error, `-3` parser state error, `-4` FASTQ length mismatch).

Use when:

- You need custom buffering/control around reads.
- You want explicit status handling and tighter control of failure behavior.

### `readFQPair(path1, path2, checkNames = false)`

Purpose-built paired-end reader.

- Tradeoff: convenient and safe for paired workflows; not pointer-based.
- Behavior: throws `IOError` if one file ends before the other.
- Validation: with `checkNames = true`, throws `ValueError` on mismatched read names.
- Stdin note: `"-"` is supported for `path1`, but `path2` cannot be `"-"`.

Use when:

- You need lockstep R1/R2 processing.
- You want optional read-name validation in the same loop.

## Practical recommendations by workload

- Command-line utility or prototype: start with `readFQ`.
- Production high-throughput counting/filtering: benchmark `readFQPtr` first.
- Stream orchestration (custom open/close, status-aware parsing): use `readFastx`.
- Any paired-end pipeline: use `readFQPair`.

## Performance snapshot (repository benchmark)

From `benchmark/2026-02-27*.txt` on `benchmark/M_abscessus_HiSeq.fq`:

- `readFQPtr`: about 8.8M-9.4M records/sec
- `readFastx`: about 4.35M-4.55M records/sec
- `readFQ`: about 3.0M-3.3M records/sec

Approximate ratios in these runs:

- `readFQPtr` is about 2.8x faster than `readFQ`
- `readFQPtr` is about 2.0x faster than `readFastx`
- `readFastx` is about 1.4x faster than `readFQ`

Treat these as directional; your workload and I/O environment can shift absolute numbers.

## Suggested default policy

If you are unsure:

1. Start with `readFQ`.
2. Profile.
3. Switch only hot paths to `readFQPtr`.
4. Use `readFastx` only where control/status handling is required.
5. Use `readFQPair` whenever two FASTQ files must stay synchronized.
