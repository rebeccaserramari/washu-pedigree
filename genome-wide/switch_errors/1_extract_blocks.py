import json
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
from pathlib import Path
import datetime
from subprocess import CompletedProcess

from tqdm import tqdm
from rich import print

def build_fasta_path(assembly_str: str) -> str:
    """
    Transforms transmitted contig (e.g., 'PAN028.chr1.haplotype2') into path for chromosomal assembly at
    ./genomes/<assembly>/<haplotype>/<assembly>.<chromosome>.<haplotype>.fasta

    Parameters
    ----------
    assembly_str - transmitted block contig

    Returns
    -------
    path to expected chromosomal assembly
    """

    assembly, chrom, haplo = assembly_str.split(".")
    return f"./genomes/{assembly}/{haplo}/{assembly}.{chrom}.{haplo}.fasta"

def extract_region(block: dict[str, str | int], role: str) -> str | None:
    """
    Extracts sequence from chromosomal assembly given transmitted role (grandparent, parent, or daughter) and region
    defined by transmitted block

    Parameters
    ----------
    block - results.json dictionary, defining transmitted block from previous script
    role  - role for which to extract sequence

    Returns
    -------
    Path to extracted fasta file
    """
    global failures, successes

    asm: str = block[f"{role}_assembly"]
    start: str = block[f"{role}_from"]
    end: str = block[f"{role}_to"]
    fasta: str = build_fasta_path(asm)

    out_path: str = os.path.join(WORKING_DIR, f"{block['id']}_{role}.fasta")

    os.makedirs(WORKING_DIR, exist_ok=True)

    region: str = f"{asm}:{start}-{end}"
    cmd: str = f"samtools faidx {fasta} {region} > {out_path}"

    res: CompletedProcess | None = None
    try:
        res = subprocess.run(cmd, shell=True, check=True, capture_output=True, text = True)
        successes += 1
    except subprocess.CalledProcessError:
        print(f"samtools faidx {fasta} {region} failed: {res.stderr if res is not None else ""}")
        failures += 1
        return None
    return out_path

def process_block(block: dict[str, str | int | dict]) -> dict[str, str | int | dict]:
    """
    Extracts sequences for all three roles (e.g., grandparent, parent, and daughter) for transmitted block

    Parameters
    ----------
    block - transmitted block definition from results.json

    Returns
    -------
    Updated transmitted block definition, with paths to extracted fasta files
    """

    result_paths: dict[str, str] = {}
    for role in ["grandparent", "parent", "daughter"]:
        result_paths[role] = extract_region(block, role)
    block["extracted_files"] = result_paths
    return block

def main(input_json: Path) -> None:
    """
    Extracts sequences for transmitted blocks, defined by previous script

    Parameters
    ----------
    input_json - path to results.json file generated from previous script

    Returns
    -------
    None
    """
    with open(str(input_json)) as f:
        blocks: dict[str, str | int | dict] = json.load(f)["data"]

    updated: list[dict[str, str | int]] = []
    for block in tqdm(blocks, desc="Extracting sequences for transmitted blocks..."):
        updated.append(process_block(block))

    with open(str(input_json), "w") as f:
        json.dump({
            "from_step": "1_extract_blocks",
            "executed_at": str(datetime.datetime.now()),
            "data": updated
        }, f, indent=4)

    print(f"Updated JSON written to {str(input_json)}")

parser = argparse.ArgumentParser()
parser.add_argument("--data", default="blocks/results.json", help="Path to pipeline's data file")
parser.add_argument("--out_dir", default="./extracted_blocks", help="Path to output directory with extracted sequences from transmitted blocks")
args = parser.parse_args()

BLOCKS_FILE: Path = Path(args.data).resolve()
WORKING_DIR: Path = Path(args.out_dir).resolve()

if not os.path.exists(WORKING_DIR):
    os.makedirs(WORKING_DIR, exist_ok=True)

successes: int = 0
failures: int = 0

if __name__ == "__main__":

    print("\n" * 3)
    print('[cyan bold]' + " Switch error detection pipeline ".center(50, "=") + '[/cyan bold]')
    print('[cyan bold]' + "=" + "[STEP 1] Extraction of transmitted sequences".center(48, " ") + "=" + '[/cyan bold]')
    print('[cyan bold]' + "".center(50, "=") + '[/cyan bold]')
    print()

    with open(BLOCKS_FILE) as f:
        data = json.load(f)
        if data["from_step"] not in ("0_generate_blocks_query", "1_extract_blocks"):
            print(f"[red bold]The data file ({BLOCKS_FILE}) does not immediately precede previous step! Aborting.[/red bold]")
            assert False

    main(BLOCKS_FILE)

    print("[green bold]Finished![/green bold]")
    print(f"Extracted {successes} sequences, failed on {failures} ones...")
