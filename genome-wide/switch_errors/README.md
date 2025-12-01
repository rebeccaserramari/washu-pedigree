# T2T Pedigree Switch-Error Detection Pipeline

This repository contains a multi‑step pipeline for detecting **switch errors** in T2T pedigree assemblies. The workflow processes *transmitted genomic blocks* across three generations (grandparent → mother → daughter) and identifies variants whose allelic behavior across liftover suggests switch‑error‑prone regions.

The pipeline consists of the following stages:

1. Generating block query metadata
2. Extracting block sequences
3. Variant calling inside transmitted blocks
4. Building BED files
5. Liftover to opposite haplotypes
6. Extracting lifted‑over sequence context
7. Summarizing and clustering switch‑error‑prone variants

Each stage is controlled by a dedicated script described below.

---

## 0. Generate Block Query Metadata

**Script:** `0_generate_blocks_query.py`

This script creates the initial *results.json* file used by the entire pipeline. It loads transmitted block definitions from TSV files in a directory.

### Input format

Each TSV row contains three tab‑separated genomic intervals:

```
PAN010.chr1.haplotype1:169073210-170484256	PAN027.chr1.maternal:172620692-174031726	PAN028.chr1.haplotype2:172160018-173571062
```

These represent:

* Grandparent block
* Mother block
* Daughter block

### Arguments

* `--input` (default: `./blocks`): Directory with transmitted‑block TSV files.
* `--data` (default: `./blocks/results.json`): Output path for generated metadata.

This script produces **results.json**, which manages all block metadata for subsequent steps.

---

## 1. Extract Block Sequences

**Script:** `1_extract_blocks.py`

Extracts FASTA sequences corresponding to all transmitted blocks.

### Expected genome layout

Assemblies must be stored as:

```
./genomes/{assembly}/{haplotype}/{assembly}.{chrom}.{haplotype}.fasta
```

### Arguments

* `--input` (default: `blocks/results.json`)
* `--out_dir` (default: `./extracted_blocks`)

Uses `samtools faidx` to extract the sequences for:

* grandparent blocks
* mother blocks
* daughter blocks

Outputs the extracted FASTA files for all transmitted blocks.

---

## 2. Variant Calling

**Script:** `2_call_variants.py`

Performs variant calling within each transmitted block.

### Arguments

* `--data` (default: `blocks/results.json`)
* `--out-dir` (default: `./variants`)
* `--workers` (default: `25`): multiprocessing workers
* `--ignore-variants-in-regions` (optional BED): regions to exclude (grandparent‑origin)
* `--ignore-variants-in-regions2` (optional BED): regions to exclude (mother‑origin)

### Description

Variant calling is done in *mother‑based coordinate system*:

* **grandparent → mother**
* **daughter → mother**

We then subtract granddaughter variants from grandparent variants, keeping only variants that:

* appear in the grandparent
* are present in the mother
* and are inherited by the daughter

Then variants originating in problematic regions are filtered using the provided BED masks.

Outputs VCF files for each transmitted block.

---

## 3. Build BED Files for Liftover

**Script:** `3_build_beds.py`

Converts VCF files into BED‑formatted variant records.

### Arguments

* `--data` (default: `blocks/results.json`)
* `--out-dir` (default: `./beds_pre_liftover`)
* `--indels` (flag): include INDELs (default: off)

Outputs pre‑liftover BED files.

---

## 4. Liftover to Opposite Haplotypes

**Script:** `4_liftover.py`

Generates chain files and performs liftover of variant BEDs to the opposite haplotype.

### Arguments

* `--data` (default: `blocks/results.json`)
* `--threads-chains` (default: `4`)
* `--threads-liftover` (default: `8`)
* `--only-chains` (flag): only build chains
* `--lift-indels` (flag)
* `--genomes-dir` (default: `./genomes`)
* `--out-chains` (default: `./chains`)
* `--out-beds` (default: `./beds_post_liftover`)

### Description

For each assembly × chromosome × haplotype, chain files are built to enable:

```
{haplotypeA} → {haplotypeB}
```

Then all variant BEDs from step 3 are lifted over to the opposite haplotype.

Outputs lifted BED files.

---

## 5. Extract Lifted‑Over Sequence Context

**Script:** `5_extract_lifted.py`

Extracts sequence context for lifted‑over variant regions.

### Arguments

* `--data` (default: `blocks/results.json`)
* `--genomes-dir` (default: `./genomes`)
* `--out-beds` (default: `./beds_post_extraction`)
* `--threads` (default: `20`)

This step mirrors step 1 but for lifted‑over variant intervals.

Outputs FASTA or region files depending on implementation.

---

## 6. Summarize and Identify Switch‑Error‑Prone Variants

**Script:** `6_summarize.py`

This is the main analytical step: it determines which variants are consistent with switch errors and clusters them into switch‑error‑prone regions.

### Arguments

* `--data` (default: `blocks/results.json`)
* `--threads` (default: `20`)
* `--indels` (flag): include trustworthy INDELs
* `--plots-dir` (default: `./plots`)
* `--genomes-dir` (default: `./genomes`)
* `--out_vcf` (default: `plots/sw_prone_variants.vcf`)
* `--skip-clustering` (flag)
* `--out_bed` (default: `plots/sw_prone_regions.bed`)
* `--cluster-max-dist` (default: `20000`)
* `--cluster-min-members` (default: `8`)
* `--all-variants-vcf` (default: `src/variants.vcf`)

### Switch‑error logic

For each variant, the script compares:

* reference allele in original block
* query allele in original block
* reference allele after liftover
* query allele after liftover

If the liftover‑transformed alleles explain the mutation in a way that is characteristic of haplotype switch errors, the variant is labeled **switch‑error‑prone**.

### Clustering

Switch‑error‑prone variants are clustered per contig using:

* maximum cluster distance (`--cluster-max-dist`)
* minimum cluster size (`--cluster-min-members`)

Outputs:

* **VCF** of all switch‑error‑prone variants
* **BED** of clustered switch‑error regions
* Plots showing variant density and switch‑error density in sliding windows

---

## Summary of Output Directories

| Step | Output                           |
| ---- | -------------------------------- |
| 0    | `blocks/results.json`            |
| 1    | `extracted_blocks/`              |
| 2    | `variants/`                      |
| 3    | `beds_pre_liftover/`             |
| 4    | `chains/`, `beds_post_liftover/` |
| 5    | `beds_post_extraction/`          |
| 6    | `plots/`, VCF + BED outputs      |

