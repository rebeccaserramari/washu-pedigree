#!/bin/bash
#PBS -N prepare_pod5_nonrdna_WGS
#PBS -l select=1:ncpus=1:mem=64gb:scratch_local=500gb
#PBS -l walltime=4:00:00

#PBS -j oe
#PBS -o /storage/praha5-elixir/projects/bioinf-fi/polakova/logs/prepare_pod5_nonrdna_WGS.log

#PBS -m abe
#PBS -M 550200@mail.muni.cz

venv="/storage/praha5-elixir/projects/bioinf-fi/polakova/venv/bioinf-fi/bin/activate"  # path to venv with pod5 CLI

trap 'clean_scratch' TERM EXIT

IN_DIR="$1" # folder with pod5 files in $IN_DIR/pod5
RDNA_IDs="$2"  # file with ids of reads containing rdna units (e.g. output of identify_rdna_reads.array.sh); line format -> <filename>:<read_id>
OUTPUT_FILE_PATH="$3"
LOG_FILE="$4"

# ------------------------------------------------------------
# Step 1: Copy only relevant POD5 files (files with rDNA reads) into scratch
# ------------------------------------------------------------

cd "$SCRATCHDIR" || exit 1
mkdir -p pod5_selected
cd pod5_selected || exit 1

while IFS= read -r line; do
  base="${line%%:*}"
  basename="$base.pod5"

  if [ -f "$basename" ]; then  # skip if already copied
    continue
  fi

  if [ -f "$IN_DIR/$basename" ]; then
    cp "$IN_DIR/$basename" .
    echo "$basename" >> "$LOG_FILE/raw_ont_rdna_copied.txt"
  else
    echo "missing $IN_DIR/$basename"
  fi
done < "$RDNA_IDs"

cd ..

# ------------------------------------------------------------
# Step 2: Extract only rDNA reads from copied POD5 files
# ------------------------------------------------------------

OUT_DIR="raw_ont"
TMP="$OUT_DIR/__tmp"
POD5_DIR="pod5_selected"          # all pod5 files are in scratch/pod5_selected

source "$venv"
mkdir -p "$OUT_DIR" "$TMP"

while read -r filename; do
  echo "[$filename] processing..."

  # Extract read IDs for this POD5 file
  grep "^$filename:" "$RDNA_IDs" | cut -d':' -f2 > "${TMP}/${filename}_reads.txt"

  if [[ ! -f "${POD5_DIR}/${filename}.pod5" ]]; then
    echo "Warning: Pod5 file ${filename}.pod5 not found, skipping."
    continue
  fi

  if [[ ! -s "${TMP}/${filename}_reads.txt" ]]; then
    echo "Warning: No reads listed for $filename, skipping."
    continue
  fi

  # Prepare subset CSV (mapping output filename to read IDs)
  awk -v fname="${filename}_extracted.pod5" '{ print fname "," $0 }' "${TMP}/${filename}_reads.txt" > "${TMP}/${filename}_subset.csv"

  # Extract reads
  if ! pod5 subset --missing-ok --csv "${TMP}/${filename}_subset.csv" -o "${OUT_DIR}" "${POD5_DIR}/${filename}.pod5"; then
    echo "Error: pod5 subset failed for $filename, skipping."
    continue
  fi

done < <(cut -d':' -f1 "$RDNA_IDs" | sort -u)

# ------------------------------------------------------------
# Step 3: Merge all extracted POD5 files
# ------------------------------------------------------------

pod5 merge $(find "$OUT_DIR" -name "*_extracted.pod5") -o "${OUT_DIR}/__all_reads.pod5" --force-overwrite

# ------------------------------------------------------------
# Step 4: Copy merged file back to permanent storage
# ------------------------------------------------------------

cp "${OUT_DIR}/__all_reads.pod5" "$OUTPUT_FILE_PATH"
