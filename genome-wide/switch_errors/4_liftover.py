import json
import os
import subprocess
from collections import defaultdict
from pathlib import Path
import argparse
from tqdm import tqdm
import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from enum import Enum
from rich import print

class LIFTOVER_STATUS(Enum):
    ERROR_MISSING_DATA = 0
    ERROR_MISSING_CHAIN_FILE = 1
    ERROR_LIFTOVER_FAILED = 2
    ERROR_COUNTING = 3
    OK = 4

def run_cmd(cmd: str) -> bool:
    """
    Safely executes shell command

    Parameters
    ----------
    cmd - shell command to be executed

    Returns
    -------
    True if command was executed successfully (exit code 0), False otherwise
    """

    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        print(f"[ERROR] Command failed: {cmd}")
        return False

def inverse_haplotype(hap: str) -> str:
    """
    Returns the inverse haplotype, as per PAN010, PAN011, and PAN027 assemblies

    Parameters
    ----------
    hap - haplotype to be inversed

    Returns
    -------
    inverse (the other) haplotype
    """

    if hap == "haplotype1":
        return "haplotype2"
    if hap == "haplotype2":
        return "haplotype1"
    if hap == "maternal":
        return "paternal"
    if hap == "paternal":
        return "maternal"
    assert False, "unknown haplotype"

def count_bed_entries(bed_path: Path) -> int | None:
    """
    Count number of entries in a BED file

    Parameters
    ----------
    bed_path - path to BED file

    Returns
    -------
    number of entries in BED file, None on failure
    """

    try:
        count: int = 0
        with open(bed_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                count += 1
        return count
    except Exception as e:
        print(f"[ERROR] Could not read {bed_path}: {e}")
        return None

def liftover_bed(block: dict, chain_dict: dict, lift_indels: bool) -> tuple[dict, LIFTOVER_STATUS]:
    """
    Performs liftover of coordinates (variants in BED files) using generated chain files and UCSC LiftOver tool

    Parameters
    ----------
    block       - transmitted block definition from results.json file, containing paths to bed files
    chain_dict  - dictionary of available chain files
    lift_indels - if True, LiftOver is set to process indels, i.e., allow partial liftovers

    Returns
    -------
    updated transmitted block definition, and liftover status
    """

    lifted_variants_mother = unlifted_variants_mother = lifted_variants_grandparent = unlifted_variants_grandparent = 0

    for bed_block, ident in [('mother_liftover', 'mother'), ('grandparent_liftover', 'grandparent')]:
        if block.get(bed_block) is None:
            return block, LIFTOVER_STATUS.ERROR_MISSING_DATA

        bed_in: Path = Path(block.get(bed_block))
        if not bed_in.exists() or bed_in.stat().st_size == 0:
            print(f"[SKIP] Missing or empty BED file for block {block['id']}")
            return block, LIFTOVER_STATUS.ERROR_MISSING_DATA

        with open(bed_in) as f: # Read first line to get origin haplotype
            first_line = f.readline().strip().split("\t")
        origin_hap: str = first_line[0]  # e.g., PAN010.chr1.haplotype1
        src_asm: str = origin_hap.split(".")[0]
        chrom: str = origin_hap.split(".")[1]
        hap_from: str = origin_hap.split(".")[2]
        hap_to: str = inverse_haplotype(hap_from)

        chain_key: str = f"{src_asm}.{chrom}.{hap_from}->{hap_to}"
        chain_file: str = chain_dict.get(chain_key)
        if not chain_file or not Path(chain_file).exists():
            print(f"[ERROR] Missing chain file for {chain_key}")
            return block, LIFTOVER_STATUS.ERROR_MISSING_CHAIN_FILE

        lifted_out: Path = LIFTED_BED_DIR / f"{block['id']}_{ident}.lifted.bed"
        unmapped_out: Path = LIFTED_BED_DIR / f"{block['id']}_{ident}.unmapped.bed"

        cmd = [
            "liftOver",
            str(bed_in),
            str(chain_file),
            str(lifted_out),
            str(unmapped_out),
        ]

        if lift_indels:
            cmd.append("-multiple")
            cmd.append("-minMatch=0.5")

        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print(f"[ERROR] liftOver failed for block {block['id']}")
            return block, LIFTOVER_STATUS.ERROR_LIFTOVER_FAILED

        lifted_bed_cnt: int = count_bed_entries(lifted_out)
        unmapped_bed_cnt: int = count_bed_entries(unmapped_out)
        if lifted_bed_cnt is None or unmapped_bed_cnt is None:
            return block, LIFTOVER_STATUS.ERROR_COUNTING

        if ident == "mother":
            lifted_variants_mother += lifted_bed_cnt
            unlifted_variants_mother += unmapped_bed_cnt
        elif ident == "grandparent":
            lifted_variants_grandparent += lifted_bed_cnt
            unlifted_variants_grandparent += unmapped_bed_cnt
        else:
            assert False, "unknown ident"

        block[f"lifted_bed_{ident}"] = str(lifted_out)
        block[f"unmapped_bed_{ident}"] = str(unmapped_out)

    block["lifted_mother_count"] = lifted_variants_mother
    block["unlifted_mother_variants"] = unlifted_variants_mother
    block["lifted_grandparent_count"] = lifted_variants_grandparent
    block["unlifted_grandparent_variants"] = unlifted_variants_grandparent
    return block, LIFTOVER_STATUS.OK

def main(input_json: Path, chain_dict: dict, workers: int = 4, lift_indels: bool = False) -> dict:
    """
    Main script for liftover of beds onto the other haplotype of the contig

    Parameters
    ----------
    input_json  - transmitted blocks file (results.json)
    chain_dict  - dictionary of available chain files
    workers     - number of threads to run in parallel
    lift_indels - if True, liftover will be optimized for indels

    Returns
    -------
    dictionary of basic statistics
    """

    debug: defaultdict[str, int] = defaultdict(int)
    lifted_total_m = unlifted_total_m = lifted_total_gp = unlifted_total_gp = 0

    with open(input_json) as f:
        blocks = json.load(f)["data"]

    updated = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(liftover_bed, block, chain_dict, lift_indels): block for block in blocks}
        for future in tqdm(as_completed(futures), total=len(blocks), desc="Lifting over..."):
            res, status = future.result()

            if status != LIFTOVER_STATUS.OK:
                debug[status.name] += 1
            else:
                lifted_total_m += res["lifted_mother_count"]
                unlifted_total_m += res["unlifted_mother_variants"]
                lifted_total_gp += res["lifted_grandparent_count"]
                unlifted_total_gp += res["unlifted_grandparent_variants"]
            updated.append(res)

    with open(str(input_json), "w") as f:
        json.dump({
            "from_step": "4_liftover",
            "executed_at": str(datetime.datetime.now()),
            "data": updated
        }, f, indent=4)

    debug["Lifted variants mother"] = lifted_total_m
    debug["Unlifted variants mother"] = unlifted_total_m
    debug["Lifted variants grandparent"] = lifted_total_gp
    debug["Unlifted variants grandparent"] = unlifted_total_gp
    print(f"Updated JSON written to {input_json}")
    return debug

def build_chain(src_fa: str, tgt_fa: str, chain_path: Path, threads: int = 8) -> Path:
    """
    Builds chain file from two fasta files

    Parameters
    ----------
    src_fa      - source fasta file
    tgt_fa      - target fasta file
    chain_path  - output path for chain file
    threads     - number of threads for minimap2

    Returns
    -------
    Path to created chain file
    """

    chain_path: Path = Path(chain_path)
    if chain_path.exists() and chain_path.stat().st_size > 0:
        print(f"[INFO] Chain already exists: {chain_path}")
        return chain_path
    print(f"[INFO] Building chain: {chain_path}")
    raw_paf: Path = chain_path.with_suffix(".raw.paf")
    sorted_paf: Path = chain_path.with_suffix(".sorted.paf")

    print(f"[INFO] Creating raw PAF: {raw_paf}")
    subprocess.run(f"minimap2 -cx asm5 --cs -t {threads} {src_fa} {tgt_fa} > {raw_paf}", shell=True, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30000)

    print(f"[INFO] Sorting PAF -> {sorted_paf}")
    subprocess.run(f"sort -k6,6 -k8,8n {raw_paf} > {sorted_paf}", shell=True, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30000)

    print(f"[INFO] Producing chain: {chain_path}")
    subprocess.run(f"paf2chain -i {sorted_paf} > {chain_path}", shell=True, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30000)

    return chain_path

def build_chain_task(args: tuple) -> tuple[str, str | None]:
    """
    Creates chain file for key (defining chain file)

    Parameters
    ----------
    args - arguments for defining chain file (source and target fasta files, output path, parallelism, and chain key)

    Returns
    -------
    chain key and path to chain file if successful, otherwise chain key and None
    """

    src_fa, tgt_fa, chain_path, threads, key = args
    try:
        chain_file = build_chain(src_fa, tgt_fa, chain_path, threads)
        return key, str(chain_file)
    except Exception as e:
        print(f"[ERROR] Failed building chain for {key}: {e}")
        return key, None

def build_chains(workers: int) -> tuple[dict, int, int]:
    """
    Builds all necessary chain files for liftover over PAN010, PAN011, and PAN027 assemblies (over both haplotypes)

    Parameters
    ----------
    workers - number of threads for minimap2

    Returns
    -------
    dictionary of chain keys and chain files, number of successfully built chain files,
    and number of unsuccessfully built chain files
    """

    chain_dict: dict[str, str] = {}
    succeeded = failed = 0

    CHAIN_TASKS = {
        "PAN010": ["haplotype1", "haplotype2"],
        "PAN011": ["haplotype1", "haplotype2"],
        "PAN027": ["maternal", "paternal"],
        "PAN028": ["haplotype1", "haplotype2"],
    }

    tasks = []
    for asm in CHAIN_TASKS.keys():
        asm_haplotypes: list[str] = CHAIN_TASKS[asm]
        chromosomes_h1: list[str] = [chrom_fa.stem.split(".")[1] for chrom_fa in sorted((GENOME_DIR / asm / asm_haplotypes[0]).glob(f"{asm}.chr*.fasta"))]
        chromosomes_h2: list[str] = [chrom_fa.stem.split(".")[1] for chrom_fa in sorted((GENOME_DIR / asm / asm_haplotypes[1]).glob(f"{asm}.chr*.fasta"))]
        chromosomes: set[str] = set(chromosomes_h1).union(set(chromosomes_h2))

        for chrom in chromosomes:
            for hap_from, hap_to in [(asm_haplotypes[0], asm_haplotypes[1]), (asm_haplotypes[1], asm_haplotypes[0])]:

                if asm == "PAN011":  # fix for only male assembly
                    if hap_from == "haplotype1" and chrom == "chrX":
                        src_fa: Path = GENOME_DIR / asm / hap_from / f"{asm}.{chrom}.{hap_from}.fasta"
                        tgt_fa: Path = GENOME_DIR / asm / hap_to / f"{asm}.chrY.{hap_to}.fasta"
                        chain_path: Path = CHAIN_DIR / f"{asm}.{chrom}.{hap_from}_to_{hap_to}.chain"
                        key: str = f"{asm}.{chrom}.{hap_from}->{hap_to}"
                        tasks.append((src_fa, tgt_fa, chain_path, workers, key))
                        continue
                    if hap_from == "haplotype2" and chrom == "chrY":
                        src_fa: Path = GENOME_DIR / asm / hap_from / f"{asm}.{chrom}.{hap_from}.fasta"
                        tgt_fa: Path = GENOME_DIR / asm / hap_to / f"{asm}.chrX.{hap_to}.fasta"
                        chain_path: Path = CHAIN_DIR / f"{asm}.{chrom}.{hap_from}_to_{hap_to}.chain"
                        key: str = f"{asm}.{chrom}.{hap_from}->{hap_to}"
                        tasks.append((src_fa, tgt_fa, chain_path, workers, key))
                        continue

                src_fa: Path = GENOME_DIR / asm / hap_from / f"{asm}.{chrom}.{hap_from}.fasta"
                tgt_fa: Path = GENOME_DIR / asm / hap_to / f"{asm}.{chrom}.{hap_to}.fasta"
                chain_path: Path = CHAIN_DIR / f"{asm}.{chrom}.{hap_from}_to_{hap_to}.chain"
                key: str = f"{asm}.{chrom}.{hap_from}->{hap_to}"
                tasks.append((src_fa, tgt_fa, chain_path, workers, key))

    print(f"[italic]Created {len(tasks)} tasks for chain files...[/italic]")

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(build_chain_task, t) for t in tasks]
        for fut in tqdm(as_completed(futures), total=len(tasks), desc="Building chain files..."):
            key, chain_file = fut.result()
            if chain_file is not None:
                chain_dict[key] = chain_file
                succeeded += 1
            else:
                failed += 1

    return chain_dict, succeeded, failed

parser = argparse.ArgumentParser()
parser.add_argument("--data", default="blocks/results.json", help="Path to pipeline's data file")
parser.add_argument("--threads-chains", type=int, default=4, help="Number of workers for chain building")
parser.add_argument("--threads-liftover", type=int, default=8, help="Number of workers for liftover")
parser.add_argument("--only-chains", action='store_true', help="Only generate chain files")
parser.add_argument("--lift-indels", action='store_true', help="Also lift indels")
parser.add_argument("--genomes-dir", default="./genomes", help="Directory to genomic data")
parser.add_argument("--out-chains", default="./chains", help="Output directory for chain files")
parser.add_argument("--out-beds", default="./beds_post_liftover", help="Output directory for lifted beds")
args = parser.parse_args()

BLOCKS_FILE = Path(args.data).resolve()
GENOME_DIR = Path(args.genomes_dir).resolve()
CHAIN_DIR = Path(args.out_chains).resolve()
LIFTED_BED_DIR = Path(args.out_beds).resolve()

if not os.path.exists(CHAIN_DIR):
    os.makedirs(CHAIN_DIR, exist_ok=True)

if not os.path.exists(LIFTED_BED_DIR):
    os.makedirs(LIFTED_BED_DIR, exist_ok=True)

if __name__ == "__main__":
    print("\n" * 3)
    print('[cyan bold]' + " Switch error detection pipeline ".center(50, "=") + '[/cyan bold]')
    print('[cyan bold]' + "=" + "[STEP 4] Liftover".center(48, " ") + "=" + '[/cyan bold]')
    print('[cyan bold]' + "".center(50, "=") + '[/cyan bold]')
    print()

    with open(BLOCKS_FILE, 'r') as f:
        data = json.load(f)
        if data["from_step"] not in ("3_build_beds", "4_liftover"):
            print(f"[red bold]The data file ({BLOCKS_FILE}) does not immediately precede previous step! Aborting.[/red bold]")
            assert False

    chain_dict, chain_success, chain_fail = build_chains(args.threads_chains)
    print(f"[green bold]Finished building chain files! Success: {chain_success} [/green bold][red bold]Failed: {chain_fail}[/red bold]")
    if args.only_chains:
        print("[orange bold][WARN] Only generating chain files[/orange bold]")
        exit(0)

    debug = main(BLOCKS_FILE, chain_dict, workers=args.threads_liftover, lift_indels=args.lift_indels)
    print("[green bold]Finished![/green bold]")
    for k, v in debug.items():
        print(f"\t{k}: {v}")
