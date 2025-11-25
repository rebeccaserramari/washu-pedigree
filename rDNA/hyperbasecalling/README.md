# Scripts for rDNA Read Identification and Extraction

This folder contains two PBS-ready helper scripts used for processing Oxford Nanopore data in order to:
1. `identify_rdna_reads.array.sh`
   - **Identify reads containing rDNA units** within basecalled FASTQ files.
2. `prepare_rdna_pod5.sh`
    - **Extract only rDNA reads from POD5 files** and merge them into a single POD5 file for downstream analysis.

These scripts are optimized for the **MetaCentrum cluster environment** (use **scratch storage** for performance and clean up automatically), and rely on **minimap2** and the **POD5 CLI**.

---

1. `identify_rdna_reads.array.sh`
   ```bash
   identify_rdna_reads.array.sh <reference.fasta> <fastq_directory> <output_ids.txt>
   ```
    **Variables**
    - `<reference.fasta>` - **reference rDNA sequence** (array) used for alignment in fasta format
    - `<fastq_directory>` - directory containing **basecalled FASTQ files**, from which the **rDNA-containing reads will be identified** and extracted.
    - `<output_ids.txt>` - **path** where the list of detected rDNA read IDs will be written
   
    **Dependencies**
    - `minimap2` (available in the provided virtual environment)
   
    **Output**
    - **text file** listing **read IDs** that **contain rDNA sequences**
    - line format: `<fastq_filename>:<read_id>`

---

2. `prepare_rdna_pod5.sh`
   ```bash
   prepare_rdna_pod5.sh <pod5_input_directory> <rdna_ids.txt> <output_merged.pod5> <log_directory>
   ```
   
   **Variables**  
   - `<pod5_input_directory>` - directory containing the **original POD5 files** (needs to have subdirectory <pod5_input_directory>/pod5/)
   - `<rdna_ids.txt>` - text file listing rDNA reads (line format: `<fastq_filename>:<read_id>`)
   - `<output_merged.pod5>` - path where the final merged POD5 file containing only rDNA reads will be saved
   - `<log_directory>` - directory where logs will be written

   **Dependencies**
   - `POD5 CLI` tools (pod5 subset, pod5 merge)

   **Output**
   - **merged POD5 file** containing only **reads that include the rDNA reference array**.

---

## Notes
- **POD5 and FASTQ files must share the same base filenames.**  
  The pipeline assumes that a read listed as coming from `sample123.fastq` will be found in `sample123.pod5`. If filenames do not match, reads cannot be extracted correctly.

## Example workflow
### **Step 1**
- Identify rDNA reads with 
``` bash
identify_rdna_reads.array.sh reference.fasta fastq_dir rdna_read_ids.txt
```
- output: `rdna_read_ids.txt`

### **Step 2**
- Extract rDNA reads from POD5 files with 
``` bash
prepare_rdna_pod5.sh pod5_input rdna_read_ids.txt extracted_rdna_reads.pod5 logs/
```
- output: `extracted_rdna_reads.pod5` **containing only reads with rDNA**
