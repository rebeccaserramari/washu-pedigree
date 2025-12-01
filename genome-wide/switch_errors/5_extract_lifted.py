import json
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse
from pathlib import Path
import datetime
from tqdm import tqdm
from enum import Enum

def run_cmd(cmd: str) -> str | None:
    """
    Safely executes a shell command

    Parameters
    ----------
    cmd - to be executed in shell

    Returns
    -------
    command's STDOUT if successful, None otherwise
    """

    try:
        res = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)
        return res.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {cmd}\n{e.stderr}")
        return None

def extract_seq(fasta_path: Path, chrom: str, start: int, end: int) -> str:
    """
    Extract sequence from fasta using samtools faidx

    Parameters
    ----------
    fasta_path - path to fasta file
    chrom      - contig from which to extract sequence
    start      - start position of sequence
    end        - end position of sequence

    Returns
    -------
    extracted sequence or None on failure
    """

    region: str = f"{chrom}:{start}-{end}"
    cmd: str = f"samtools faidx {fasta_path} {region} | grep -v '^>'"
    return run_cmd(cmd)

def process_block(block: dict) -> tuple[dict, int, int]:
    """
    Processes one transmitted block by extracting sequence from lifted-over beds (for either mother and grandparent
    reference beds)

    Parameters
    ----------
    block - transmitted block definition from results.json

    Returns
    -------
    updated transmitted block definition, as well as number of failed blocks for either missing data, or failed samtools
    faidx
    """

    error_missing_data = error_failed_extraction = 0

    for block_path, ident in [('lifted_bed_mother', 'mother'), ('lifted_bed_grandparent', 'grandparent')]:

        lifted_bed = block.get(block_path)
        if not lifted_bed or not os.path.exists(lifted_bed) or os.stat(lifted_bed).st_size == 0:
            error_missing_data += 1
            continue  # not fatal error, one successful liftover can still show switch error

        out_bed: Path = BEDS_DIR / f"{block['id']}_{ident}.bed"

        with open(lifted_bed) as fin, open(out_bed, "w") as fout:
            for line in fin:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue

                if len(parts) > 4:  # indels
                    ...

                lifted_chr, start, end, info, *_ = parts
                start, end = int(start), int(end)
                llen = end - start

                asm, chrom, hap = lifted_chr.split(".")  # PAN010.chr1.haplotype1
                lifted_fa: Path = GENOME_DIR / asm / hap / f"{asm}.{chrom}.{hap}.fasta"

                lifted_allele = extract_seq(lifted_fa, lifted_chr, start + 1, start + 1 + llen - 1)
                if not lifted_allele:
                    print(f"[WARN] Failed extraction for {block['id']} {line.strip()}")
                    error_failed_extraction += 1
                    continue

                new_info: str = info + f",{lifted_chr}:{start + 1}-{start + 1 + llen - 1}:{lifted_allele}"
                fout.write(f"{lifted_chr}\t{start}\t{end}\t{new_info}\n")

                block[f"extracted_bed_{ident}"] = str(out_bed)
    return block, error_missing_data, error_failed_extraction

def main(input_json: Path, workers: int = 4) -> None:
    """
    Main worker for extraction of lifted-over sequences from previous step

    Parameters
    ----------
    input_json - transmitted blocks results json file (results.json)
    workers    - number of parallel workers to use

    Returns
    -------
    None (updates results.json file)
    """

    missing_data = failed_extraction = 0

    with open(input_json) as f:
        blocks = json.load(f)["data"]

    updated: list[dict] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(process_block, block) for block in blocks]
        for fut in tqdm(as_completed(futures), total=len(blocks)):
            res, missing, failed = fut.result()
            missing_data += missing
            failed_extraction += failed
            updated.append(res)

    with open(input_json, "w") as f:
        json.dump({
            "from_step": "5_extract_lifted",
            "executed_at": str(datetime.datetime.now()),
            "data": updated
        }, f, indent=4)
    print(f"Updated JSON written to {input_json}")
    print("Finished!")
    print(f"Blocks missing lifted bed file (from either mother or grandparent): {missing_data}")
    print(f"Sequences failed to be extracted: {failed_extraction}")

parser = argparse.ArgumentParser()
parser.add_argument("--data", default="blocks/results.json", help="Path to pipeline's data file")
parser.add_argument("--genomes-dir", default="./genomes", help="Directory to genomic data")
parser.add_argument("--out-beds", default="./beds_post_extraction", help="Output directory for extracted beds")
parser.add_argument("--threads", type=int, default=20, help="Number of parallel workers")
args = parser.parse_args()

BLOCKS_FILE = Path(args.data).resolve()
GENOME_DIR = Path(args.genomes_dir).resolve()
BEDS_DIR = Path(args.out_beds).resolve()

if not os.path.exists(BEDS_DIR):
    os.makedirs(BEDS_DIR, exist_ok=True)

if __name__ == "__main__":
    print("\n" * 3)
    print(" Switch error detection pipeline ".center(50, "="))
    print("=" + "[STEP 5] Extract lifted".center(48, " ") + "=")
    print("".center(50, "="))
    print()

    with open(BLOCKS_FILE, 'r') as f:
        data = json.load(f)
        if data["from_step"] not in ("4_liftover", "5_extract_lifted"):
            print(f"The data file ({BLOCKS_FILE}) does not immediately precede previous step! Aborting.")
            assert False

    main(BLOCKS_FILE, args.threads)
