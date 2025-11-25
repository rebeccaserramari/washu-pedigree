#!/bin/bash

#PBS -N identify_rDNA_hg06807_A_array
#PBS -l select=1:ncpus=32:mem=64gb:scratch_local=200gb
#PBS -l walltime=12:00:00

#PBS -j oe
#PBS -o /storage/praha5-elixir/projects/bioinf-fi/polakova/logs/identify_rDNA_hg06807_A_array.log

#PBS -m abe
#PBS -M 550200@mail.muni.cz

set -euo pipefail
trap 'clean_scratch' TERM EXIT

source "/storage/praha5-elixir/projects/bioinf-fi/polakova/venv/bioinf-fi/bin/activate"

THREADS=32

REF="$1"
READS_DIR="$2"
OUTPUT_FILE_PATH="$3"  # final output txt file path

SCRATCH_OUTPUT="$SCRATCHDIR/rdna_ids.txt"
> "$SCRATCH_OUTPUT"


cd "$SCRATCHDIR" || exit 1

# Copy reference
cp "$REF" ref.fa
# Copy FASTQ files to scratch
mkdir reads
find "$READS_DIR" -maxdepth 1 -name "*.fastq" -exec cp {} reads/ \;

# Index in scratch
minimap2 -d ref.fa.mmi ref.fa

# Process FASTQ files
for file in reads/*.fastq; do
  fname=$(basename "$file")
  basename="${fname%%.*}"
  echo "Processing $basename..."

  minimap2 -t $THREADS -a -x map-ont -y --MD -Y ref.fa.mmi "$file" \
    | awk -v basename="$basename" '
        BEGIN {FS="\t"}
        !/^@/ {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^AS:i:/) {
                    split($i, as, ":")
                    if (as[3] >= 3000) print basename ":" $1
                }
            }
        }
      ' >> "$SCRATCH_OUTPUT"
done

# Copy result back
mkdir -p "$(dirname "$OUTPUT_FILE_PATH")"
cp "$SCRATCH_OUTPUT" "$OUTPUT_FILE_PATH"

echo "Done."
