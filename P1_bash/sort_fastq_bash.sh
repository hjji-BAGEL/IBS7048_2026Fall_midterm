#!/bin/bash
# Problem 1: copy raw/*.fastq into sorted_bash/<condition>/<read>/
# Run inside IBS7048_2026:  bash sort_fastq_bash.sh

for condition in Control Treated; do          # 2 conditions
  for read in R1 R2; do                       # 2 read directions
    mkdir -p sorted_bash/$condition/$read     # make folder (no error if it exists)
    cp raw/*_${condition}_*_${read}.fastq sorted_bash/$condition/$read/   # copy by filename pattern
  done
done

# Report the number of files in each folder
echo "== File counts =="
echo "raw: $(ls raw | wc -l)"
for d in sorted_bash/*/*; do
  echo "$d: $(ls $d | wc -l)"
done
echo "sorted_bash total: $(find sorted_bash -name '*.fastq' | wc -l)"
