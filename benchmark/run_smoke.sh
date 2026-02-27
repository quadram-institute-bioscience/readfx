#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PATH_ARG="${1:-tests/illumina_1.fq.gz}"
LOOPS_ARG="${2:-2000}"
REPEATS_ARG="${3:-5}"

echo "Building smoke benchmark (release, ORC)..." >&2
nim c --nimcache:/tmp/nimcache_readfx_smoke -d:release --mm:orc --threads:on -o:benchmark/smoke_bench benchmark/smoke_bench.nim >/dev/null

echo "Running smoke benchmark..." >&2
echo "  path=$PATH_ARG loops=$LOOPS_ARG repeats=$REPEATS_ARG" >&2
benchmark/smoke_bench --path:"$PATH_ARG" --loops:"$LOOPS_ARG" --repeats:"$REPEATS_ARG"
