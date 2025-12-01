#!/usr/bin/env python3
import json
import os
from collections import defaultdict
from pathlib import Path
import argparse
from tqdm import tqdm
import datetime
from rich import print

stats_variants: dict[str, int] = {
    "maternal_snv": 0,
    "maternal_indel": 0,
    "paternal_snv": 0,
    "paternal_indel": 0,
}

stats_variants_per_chr: defaultdict[bool, defaultdict[str, defaultdict[bool, int]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

stats_blocks: dict[str, int] = {
    "maternal_len": 0,
    "paternal_len": 0,
}

stats_lengths_per_chr: defaultdict[bool, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))

def parse_variant_line(line: str) -> dict:
    """
    Parse a VCF line into reference + query variant coordinates in BED coordinate system.
    Note: variants called in previous step are transmitted-block-bound (i.e., their position is defined WITHIN the
           transmitted block), this script transforms these coordinates into assembly-bound ones.

    Parameters
    ----------
    line - variant line from .vcf file (from previous step)

    Returns
    -------
    dictionary defining the variant by reference and query coordinates, as well as other metadata
    """

    global stats_variants

    parts: list[str] = line.strip().split("\t")
    ref_field: str = parts[0]           # e.g. PAN027.chr1.maternal:144216130-144939231
    ref_pos_rel: int = int(parts[1])    # relative to block start
    # id = parts[2]
    ref_allele: str = parts[3]
    alt_allele: str = parts[4]
    # quality = parts[5]
    # filter = parts[6]
    info_field: str = parts[7]

    ref_contig, coords = ref_field.split(":")
    ref_block_start, ref_block_end = map(int, coords.split("-"))

    ref_abs_start = ref_block_start + ref_pos_rel - 1
    ref_abs_end = ref_block_start + ref_pos_rel + len(ref_allele) - 2

    # Parse query information from INFO field
    qname_field: str = [f for f in info_field.split(";") if f.startswith("QNAME=")][0]
    qname: str = qname_field.split("=")[1]  # e.g. PAN010.chr1.haplotype2:144470190-145193296
    query_contig, qcoords = qname.split(":")
    q_block_start, q_block_end = map(int, qcoords.split("-"))

    qstart_field = [f for f in info_field.split(";") if f.startswith("QSTART=")][0]
    qstart_rel = int(qstart_field.split("=")[1])

    # absolute query position (1-based -> 0-based BED)
    query_abs_start = q_block_start + qstart_rel - 2
    query_abs_end = q_block_start + qstart_rel + len(alt_allele) - 3

    is_indel: bool = len(ref_allele) > 1 or len(alt_allele) > 1
    is_maternal: bool = "maternal" in ref_contig
    ref_chromosome: str = ref_contig.split('.')[1]
    stats_variants_per_chr[is_maternal][ref_chromosome][is_indel] += 1

    if not is_indel and     is_maternal: stats_variants["maternal_snv"] += 1
    if     is_indel and     is_maternal: stats_variants["maternal_indel"] += 1
    if not is_indel and not is_maternal: stats_variants["paternal_snv"] += 1
    if     is_indel and not is_maternal: stats_variants["paternal_indel"] += 1

    return {
        "ref_contig": ref_contig,
        "ref_start": ref_abs_start,
        "ref_end": ref_abs_end,
        "ref_allele": ref_allele,
        "query_contig": query_contig,
        "query_start": query_abs_start,
        "query_end": query_abs_end,
        "query_allele": alt_allele,
    }


def write_beds(block_id: int, variants: list[dict], BED_DIR: Path, indels: bool) -> tuple[str, str]:
    """
    Creates two bed files for each transmitted block (defining variants in parameters), one where mother (PAN027)
    is the reference, another one where grandparents (PAN010/PAN011) are the reference. These bed files are then
    directly used in liftover.

    Parameters
    ----------
    block_id - transmitted block ID (as in results.json file)
    variants - variants extracted from .vcf file from previous step
    BED_DIR  - directory where to save bed files
    indels   - if True, indels are included. Otherwise, only SNPs are considered.

    Returns
    -------
    path to PAN027-reference and grandparent-reference BED files
    """

    mother_path: Path = BED_DIR / f"{block_id}.mother_liftover.bed"
    gp_path: Path = BED_DIR / f"{block_id}.grandparent_liftover.bed"

    with mother_path.open("w") as mother_out, gp_path.open("w") as gp_out:
        for v in variants:

            ref_size = v['ref_end'] - v['ref_start']
            query_size = v['query_end'] - v['query_start']

            if not indels and (ref_size > 1 or len(v['ref_allele']) > 1 or query_size > 1 or len(v['query_allele']) > 1):
                continue

            # INFO string
            info = (
                f"{v['ref_contig']}:{v['ref_start']}-{v['ref_end']}:{v['ref_allele']},"
                f"{v['query_contig']}:{v['query_start']}-{v['query_end']}:{v['query_allele']}"
            )

            # PAN027 BED
            mother_out.write(
                f"{v['ref_contig']}\t{v['ref_start'] - 1}\t{v['ref_start'] - 1 + len(v['ref_allele'])}\t{info}\n"
            )

            # grandparent BED
            gp_out.write(
                f"{v['query_contig']}\t{v['query_start'] - 1}\t{v['query_start'] - 1 + len(v['query_allele'])}\t{info}\n"
            )

    return str(mother_path), str(gp_path)


def main(input_json: Path, output_dir: Path, indels: bool) -> dict:
    """
    Main script for creating bed files (mother-reference and grandparent-reference) from variant calls generated
    in previous step. Bed files are called for each transmitted block separately.

    Parameters
    ----------
    input_json - path to the results.json file.
    output_dir - path to directory where to save bed files
    indels     - if True, indels are included, otherwise only SNPs are considered.

    Returns
    -------
    dictionary of encountered errors
    """

    global stats_blocks

    debug: dict[str, int] = {
        "Error missing data": 0,
        "Warn no variants in block": 0,
    }

    with open(str(input_json)) as f:
        blocks = json.load(f)["data"]

    for block in tqdm(blocks, desc="Building BEDs for liftover..."):
        block_id: int = block["id"]

        contig, contig_from, contig_to = block["parent_assembly"], block["parent_from"], block["parent_to"]
        contig_sz = contig_to - contig_from
        contig_chr: str = contig.split('.')[1]
        contig_maternal: bool = "maternal" in contig
        if contig_maternal: stats_blocks["maternal_len"] += contig_sz
        else: stats_blocks["paternal_len"] += contig_sz

        stats_lengths_per_chr[contig_maternal][contig_chr] += contig_sz

        vcf_path = block.get("variants")
        if not vcf_path or not os.path.exists(vcf_path):
            print(f"[SKIP] Variant file missing for block {block_id}")
            debug["Error missing data"] += 1
            continue

        variants = []  # ref is always PAN027
        with open(vcf_path) as vcf_in:
            for line in vcf_in:
                if line.startswith("#") or not line.strip():
                    continue
                try:
                    v = parse_variant_line(line)
                    variants.append(v)
                except Exception as e:
                    print(f"[WARN] Block {block_id}: failed to parse line: {line.strip()} ({e})")

        if not variants:
            debug["Warn no variants in block"] += 1
            continue

        ref_bed, gp_bed = write_beds(block_id, variants, output_dir, indels)
        block["mother_liftover"] = ref_bed
        block["grandparent_liftover"] = gp_bed

    with open(str(input_json), "w") as f:
        json.dump({
            "from_step": "3_build_beds",
            "executed_at": str(datetime.datetime.now()),
            "data": blocks
        }, f, indent=4)

    return debug


parser = argparse.ArgumentParser()
parser.add_argument("--data", default="blocks/results.json", help="Path to pipeline's data file")
parser.add_argument("--out-dir", default="./beds_pre_liftover", help="Output directory for variant calling")
parser.add_argument("--indels", action='store_true', help="Also include indels")

args = parser.parse_args()

BEDS_DIR = Path(args.out_dir).resolve()
if not os.path.exists(str(BEDS_DIR)):
    os.makedirs(str(BEDS_DIR), exist_ok=True)

BLOCKS_FILE = Path(args.data).resolve()

if __name__ == "__main__":

    print("\n" * 3)
    print('[cyan bold]' + " Switch error detection pipeline ".center(50, "=") + '[/cyan bold]')
    print('[cyan bold]' + "=" + "[STEP 3] Build pre-liftover .bed's".center(48, " ") + "=" + '[/cyan bold]')
    print('[cyan bold]' + "".center(50, "=") + '[/cyan bold]')
    print()


    with open(BLOCKS_FILE, 'r') as f:
        data = json.load(f)
        if data["from_step"] not in ("2_call_variants", "3_build_beds"):
            print(f"[red bold]The data file ({BLOCKS_FILE}) does not immediately precede previous step! Aborting.[/red bold]")
            assert False

    if args.indels:
        print("[WARN] Indels are BEING INCLUDED!")
    else:
        print("[WARN] Indels are *NOT* BEING INCLUDED!")
    print()

    debug = main(BLOCKS_FILE, BEDS_DIR, args.indels)
    print("[green bold]Finished![/green bold]")
    for k, v in debug.items():
        print(f"\t{k}: {v}")
