#!/bin/bash
# One job of e28_reference_tables.py. Argument: a tab-separated line "sample target endpoint raw_csv".
# Run from the repository root; E28_REF is the output folder (default paper/rev07/e28_ref).
set -euo pipefail
IFS=$'\t' read -r s t e f <<< "$1"
ref=${E28_REF:-paper/rev07/e28_ref}
out=$ref/$s/${t}_${e}
mkdir -p "$out"
PYTHONPATH=src ${PYTHON:-python3} paper/rev07/e28_reference_tables.py "$f" "$t" "$e" "$out" > "$out/log.txt" 2>&1 || echo "FAILED $s $t $e" >> "$ref/failures.txt"
