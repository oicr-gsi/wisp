#!/bin/bash
# Summarise a provisioned output directory into a form a baseline can be compared against.
#
# A bare `ls` only proves files appeared. What matters is that the same sites were called and
# the same purity was fitted, so the variant counts, the PURPLE fit and the WISP estimate are
# recorded too. Nothing sample-specific that is not deterministic is included -- no paths, no
# timestamps, no tool version lines -- so the baseline does not drift between runs.

set -o nounset
set -o pipefail

cd "$1" || exit 1

# File inventory, without the sizes that move with compression level.
find . -type f | sed 's|^\./||' | sort

# Variant counts per FILTER value: what changes when a caller behaves differently.
for vcf in *.vcf.gz; do
    [ -e "$vcf" ] || continue
    echo "## ${vcf}"
    zcat "$vcf" | grep -v '^#' | awk -F'\t' '{ print $7 }' | sort | uniq -c | sort -k2,2
done

# The numbers the assay reports.
for tsv in *.purple.purity.tsv *.wisp.summary.tsv; do
    [ -e "$tsv" ] || continue
    echo "## ${tsv}"
    cat "$tsv"
done

# What the preflight decided: mode, platforms, sample ids, primary provenance.
for log in *.validation.log; do
    [ -e "$log" ] || continue
    echo "## ${log}"
    cat "$log"
done

# Archive members only; their contents are covered above where they matter.
for tar in *.tar.gz; do
    [ -e "$tar" ] || continue
    echo "## ${tar}"
    tar -tzf "$tar" | sed 's|^\./||' | sort
done
