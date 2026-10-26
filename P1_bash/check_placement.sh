# Problem 1 extra check: every file is in the folder that matches its own name
# Run inside IBS7048_2026:  bash check_placement.sh
n=0
for f in sorted_bash/*/*/*.fastq; do
  name=$(basename "$f" .fastq)
  cond=$(echo "$name" | cut -d_ -f2)
  rd=$(echo "$name" | cut -d_ -f4)
  [ "$f" = "sorted_bash/$cond/$rd/$name.fastq" ] || { echo "Misplaced: $f"; n=$((n+1)); }
done
echo "Misplaced files: $n"
