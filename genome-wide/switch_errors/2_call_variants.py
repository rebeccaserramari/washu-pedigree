import json
import os
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import tempfile
from os.path import join
from subprocess import CompletedProcess
import gzip
from typing import Dict, List, Tuple
from tqdm import tqdm
import datetime
from pathlib import Path
from enum import Enum
from rich import print

class VARIANTCALL_RESULT(Enum):
    ERR_MISSING_DATA = 0
    ERR_VARIANT_CALL = 1
    ERR_VARIANT_COUNT = 2
    ERR_ISEC = 3
    OK = 4
    ERR_VARIANT_COUNT_SNPINDEL = 5
    ERR_ZIP = 6

def run_cmd(cmd: str) -> bool:
    """
    Safely execute a shell command

    Parameters
    ----------
    cmd - comand to be executed in shell

    Returns
    -------
    True if command executed successfully (exit code 0), False otherwise
    """

    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {cmd}")
        return False

def call_variants(ref_path: str, query_path: str, result_path: str) -> str | None:
    """
    Calls variants using minimap2 -cx asm5 -cs and paftools.js

    Parameters
    ----------
    ref_path    - reference fasta
    query_path  - query fasta
    result_path - resulting .vcf.gz path

    Returns
    -------
    result path if successful, None otherwise
    """

    with tempfile.TemporaryDirectory() as tmpdirname:
        #print(f'Calling variants ref::{ref_path} ~> query::{query_path}')

        aln_path = join(tmpdirname, "aln.paf")
        aln_sorted_path = join(tmpdirname, "aln.sorted.paf")
        vcf_path = join(tmpdirname, "out.vcf")
        vcf_gz_path = join(tmpdirname, "out.vcf.gz")

        cmd = (
            f"minimap2 -cx asm5 --cs {ref_path} {query_path} > {aln_path} && "
            f"sort -k6,6 -k8,8n {aln_path} > {aln_sorted_path} && "
            f"paftools.js call -f {ref_path} {aln_sorted_path} > {vcf_path} && "
            f"bgzip -c {vcf_path} > {vcf_gz_path} && "
            f"cp {vcf_gz_path} {result_path} && "
            f"tabix -p vcf {result_path}"
        )

        success = run_cmd(cmd)

    return result_path if success else None

def call_isec(keep_path: str, remove_path: str, result_path: str, complement: bool = True) -> str | None:
    """
    Calls bcftools isec given:

    Parameters
    ----------
    keep_path    - original .vcf.gz file
    remove_path  - .vcf.gz file to remove
    result_path  - output .vcf.gz file
    complement   - if True, variants in keep_path but not in remove_path are outputted.
                 - if False, variants from keep_path that have been removed outputted.

    Returns
    -------
    result path if successful, None otherwise
    """

    if complement:
        cmd = (
            f"bcftools isec -C -w1 {keep_path} {remove_path} -o {result_path}"
        )
    else:
        cmd = (
            f"bcftools isec -w1 -n=2 {keep_path} {remove_path} -o {result_path}"
        )
    success: bool = run_cmd(cmd)
    return result_path if success else None

def call_unzip(vcf_gz_path: str, vcf_path: str) -> str | None:
    """
    Unzips .vcf.gz file to .vcf file

    Parameters
    ----------
    vcf_gz_path - path to input .vcf.gz file
    vcf_path    - path to output .vcf file

    Returns
    -------
    vcf_path if successful, None otherwise
    """

    cmd = (
        f"gunzip -c {vcf_gz_path} > {vcf_path}"
    )
    success: bool = run_cmd(cmd)
    return vcf_path if success else None

def call_bgzip(vcf_path: str, vcf_gz_path: str) -> str | None:
    """
    Zips .vcf file to .vcf.gz file

    Parameters
    ----------
    vcf_path    - path to input .vcf file
    vcf_gz_path - path to output .vcf.gz file

    Returns
    -------
    vcf_gz_path if successful, None otherwise
    """

    cmd = {
        f"bgzip -c {vcf_path} > {vcf_gz_path} && tabix -p vcf {vcf_gz_path}"
    }
    success: bool = run_cmd(cmd)
    return vcf_gz_path if success else None

def call_index(vcf_gz_path: str) -> str | None:
    """
    Indexes .vcf.gz file using tabix

    Parameters
    ----------
    vcf_gz_path - vcf.gz file to index

    Returns
    -------
    vcf_gz_path if successful, None otherwise
    """

    cmd = {
        f"tabix -p vcf {vcf_gz_path}"
    }
    success: bool = run_cmd(cmd)
    return vcf_gz_path if success else None

def count_snp_indel(vcf_path: str) -> tuple[int | None, int | None]:
    """
    Counts SNPs and indels in a VCF or VCF.GZ file using bcftools

    Parameters
    ----------
    vcf_path - path to .vcf.gz file

    Returns
    -------
    Number of SNPs and indels if successful, (None, None) otherwise
    """

    try:
        # SNPs
        cmd_snps: str = f"bcftools view -v snps {vcf_path} | grep -vc '^#'"
        snp_result: CompletedProcess = subprocess.run(cmd_snps, shell=True, check=False,
                                                      capture_output=True, text=True)
        num_snps: int = int(snp_result.stdout.strip()) if snp_result.stdout.strip() else 0

        # indels
        cmd_indels: str = f"bcftools view -v indels {vcf_path} | grep -vc '^#'"
        indel_result: CompletedProcess = subprocess.run(cmd_indels, shell=True, check=False,
                                                        capture_output=True, text=True)
        num_indels: int = int(indel_result.stdout.strip()) if indel_result.stdout.strip() else 0

        return num_snps, num_indels

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] bcftools failed on {vcf_path}: {e}. {e.output}")
        return None, None
    except ValueError:
        print(f"[WARN] Could not parse counts from {vcf_path}")
        return None, None

def count_variants(vcf_path: str) -> int | None:
    """
    Count number of variants in a .vcf.gz file using bcftools

    Parameters
    ----------
    vcf_path - path to vcf.gz file

    Returns
    -------
    Number of variants (SNPs and indels) if successful, None otherwise
    """

    cmd: str = f"bcftools view -H {vcf_path} | wc -l"
    try:
        result: CompletedProcess = subprocess.run(cmd, shell=True, check=True,
                                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return int(result.stdout.strip())
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {cmd}")
        print(e.stderr)
        return None
    except ValueError:
        print(f"[ERROR] Could not parse output from {cmd}")
        return None

def load_bed(bed_path: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Load BED file into a dict: contig -> list of (start, end) intervals.

    Parameters
    ----------
    bed_path - to the bed file

    Returns
    -------
    dictionary of contig -> sorted list of (start, end) intervals
    """

    intervals: dict = {}
    open_fn = gzip.open if bed_path.endswith(".gz") else open

    with open_fn(bed_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)
            intervals.setdefault(chrom, []).append((start, end))

    for chrom in intervals:
        intervals[chrom].sort()

    return intervals

def pos_in_intervals(intervals: List[Tuple[int, int]], pos: int) -> bool:
    """
    Check whether 1-based position lies in any (start,end) BED interval.

    Parameters
    ----------
    intervals - bed intervals from load_bed(...)
    pos - position to be checked to be inside any interval

    Returns
    -------
    True if position is inside any interval, False otherwise
    """

    p: int = pos - 1
    for s, e in intervals:
        if p < s:
            return False
        if s <= p < e:
            return True
    return False

def filter_ignore(vcf_path: str, bed_path: str, result_path: str, is_ref: bool = False) -> None:
    """
    Filters variants in .vcf.gz file using regions defined by bed file

    Parameters
    ----------
    vcf_path     - variant file to be filtered
    bed_path     - bed file defining intervals in which to filter out variants
    result_path  - path to output .vcf.gz file
    is_ref       - if True, variant reference is taken into account, otherwise variant query is used

    Returns
    -------
    None (creates output .vcf.gz file)
    """

    bed: Dict[str, List[Tuple[int, int]]] = load_bed(bed_path)

    open_vcf = gzip.open if vcf_path.endswith(".gz") else open
    open_out = gzip.open if result_path.endswith(".gz") else open

    kept: int = 0
    total: int = 0
    skipped_missing_info: int = 0
    skipped_no_match: int = 0

    with open_vcf(vcf_path, "rt") as inp, open_out(result_path, "wt") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue

            total += 1
            fields: list[str] = line.strip().split("\t")
            info = fields[7]

            # parse QNAME and QSTART
            qname: str | None = None
            qstart: int | None = None

            if is_ref:
                ref_name = fields[0]
                ref_start, ref_end = fields[1], fields[2]
                qname: str = f"{ref_name}:{ref_start}-{ref_end}"
                qstart = int(ref_start)
            else:
                for entry in info.split(";"):
                    if entry.startswith("QNAME="):
                        qname = entry.split("=", 1)[1]
                    elif entry.startswith("QSTART="):
                        qstart = int(entry.split("=", 1)[1])

            if qname is None or qstart is None:
                skipped_missing_info += 1
                continue

            contig: str = qname.split(":", 1)[0]

            if contig not in bed:
                skipped_no_match += 1
                continue

            if pos_in_intervals(bed[contig], qstart):
                out.write(line)
                kept += 1

def process_block(block: dict, output_folder: str,
                  ignore_variants_in_regions: str, ignore_variants_in_regions2: str) -> tuple[dict, VARIANTCALL_RESULT]:
    """
    Performs variant calling and filtering for one transmitted block by:
     - calling variants grandparent -> mother
     - calling variants granddaughter -> mother
     - filtering out variants in ignore_variants_in_regions bed (query regions)
     - filtering out variants in ignore_variants_in_regions2 bed (reference regions)
     - keeping variants from grandparent -> mother which are not present (i.e., are preserved) in granddaughter -> mother

    Parameters
    ----------
    block                       - transmitted block definition from previous script
    output_folder               - folder where to output .vcf.gz files
    ignore_variants_in_regions  - bed file to define regions to ignore in variant query (grandparents)
    ignore_variants_in_regions2 - bed file to define regions to ignore in variant reference (mother)

    Returns
    -------
    updated transmitted block definition with variant calling results (if successful, otherwise None), and result of this
    step
    """

    extracted: dict = block.get("extracted_files", {})
    if not all(extracted.get(k) for k in (["grandparent", "parent", "daughter"])):
        print(f"[ERROR] Block {block['id']} does not contain all necessary extracted files")
        return block, VARIANTCALL_RESULT.ERR_MISSING_DATA

    block_id: str = block["id"]
    ref_file: str = extracted["parent"]
    gp_file: str = extracted["grandparent"]
    d_file: str = extracted["daughter"]

    variants: str = join(output_folder, f"{block_id}.vcf")

    # primary variants: grandparent -> parent
    variants_gp_to_p: str = join(output_folder, f"{block_id}_gp_to_p.vcf.gz")
    if call_variants(ref_file, gp_file, variants_gp_to_p) is None:
        return block, VARIANTCALL_RESULT.ERR_VARIANT_CALL

    variant_cnt_gp_to_p: int = count_variants(variants_gp_to_p)
    if variant_cnt_gp_to_p is None:
        return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT

    snvs_gp_to_p, indels_gp_to_p = count_snp_indel(variants_gp_to_p)
    if snvs_gp_to_p is None or indels_gp_to_p is None:
        return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT_SNPINDEL

    # ignoring variants in query regions (if set)
    variants_ignored_cnt = variants_ignored_snp = variants_ignored_indel = None
    if ignore_variants_in_regions is not None:
        try:
            variants_ignored: str = join(output_folder, f"{block_id}.ignored.vcf")
            filter_ignore(variants_gp_to_p, ignore_variants_in_regions, variants_ignored)
            variants_ignored: str = call_bgzip(variants_ignored, variants_ignored + ".gz")
        except:
            return block, VARIANTCALL_RESULT.ERR_ISEC

        variants_ignored_cnt: int | None = count_variants(variants_ignored)
        if variants_ignored_cnt is None:
            return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT

        variants_ignored_snp, variants_ignored_indel = count_snp_indel(variants_ignored)
        if variants_ignored_snp is None or variants_ignored_indel is None:
            return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT_SNPINDEL

        variants_new: str = join(output_folder, f"{block_id}.pure.vcf.gz")
        variants_new: str = call_isec(variants_gp_to_p, variants_ignored, variants_new)
        if variants_new is None:
            return block, VARIANTCALL_RESULT.ERR_ISEC
        variants_new: str = call_index(variants_new)
        if variants_new is None:
            return block, VARIANTCALL_RESULT.ERR_ZIP

        variants_gp_to_p = variants_new

    # ignoring variants in reference regions (if set)
    variants_ignored2_cnt = variants_ignored2_snp = variants_ignored2_indel = None
    if ignore_variants_in_regions2 is not None:
        try:
            variants_ignored: str = join(output_folder, f"{block_id}.ignored2.vcf")
            filter_ignore(variants_gp_to_p, ignore_variants_in_regions2, variants_ignored, is_ref=True)
            variants_ignored: str = call_bgzip(variants_ignored, variants_ignored + ".gz")
        except:
            return block, VARIANTCALL_RESULT.ERR_ISEC

        variants_ignored2_cnt: int | None = count_variants(variants_ignored)
        if variants_ignored_cnt is None:
            return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT

        variants_ignored2_snp, variants_ignored2_indel = count_snp_indel(variants_ignored)
        if variants_ignored_snp is None or variants_ignored_indel is None:
            return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT_SNPINDEL

        variants_new: str = join(output_folder, f"{block_id}.pure2.vcf.gz")
        variants_new: str = call_isec(variants_gp_to_p, variants_ignored, variants_new)
        if variants_new is None:
            return block, VARIANTCALL_RESULT.ERR_ISEC
        variants_new: str = call_index(variants_new)
        if variants_new is None:
            return block, VARIANTCALL_RESULT.ERR_ZIP
        variants_gp_to_p = variants_new

    # secondary variants: granddaughter -> parent (if pedigree)
    variants_d_to_p: str = join(output_folder, f"{block_id}_d_to_p.vcf.gz")
    if call_variants(ref_file, d_file, variants_d_to_p) is None:
        return block, VARIANTCALL_RESULT.ERR_VARIANT_CALL

    variant_cnt_d_to_p: int = count_variants(variants_d_to_p)
    if variant_cnt_d_to_p is None:
        return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT

    snvs_d_to_p, indels_d_to_p = count_snp_indel(variants_d_to_p)
    if (snvs_d_to_p is None or indels_d_to_p is None):
        return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT_SNPINDEL

    # analysis of variands in gp -> p but preserved in daughter
    variants_isec: str = call_isec(variants_gp_to_p, variants_d_to_p, variants)
    if variants_isec is None:
        return block, VARIANTCALL_RESULT.ERR_ISEC

    variant_cnt_filtered: int = count_variants(variants_isec)
    if variant_cnt_filtered is None:
        return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT

    filtered_snp, filtered_indels = count_snp_indel(variants_isec)
    if (filtered_snp is None or filtered_indels is None):
        return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT_SNPINDEL

    variants_removed: str = join(output_folder, f"{block_id}.removed.vcf")
    variants_isec_removed: str = call_isec(variants_gp_to_p, variants_d_to_p, variants_removed, complement=False)
    if variants_isec_removed is None:
        return block, VARIANTCALL_RESULT.ERR_ISEC

    variant_cnt_removed: int = count_variants(variants_isec_removed)
    if variant_cnt_removed is None:
        return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT

    removed_snp, removed_indels = count_snp_indel(variants_isec_removed)
    if (removed_snp is None or removed_indels is None):
        return block, VARIANTCALL_RESULT.ERR_VARIANT_COUNT_SNPINDEL

    variants_path = variants_isec

    block["variants"] = variants_path
    block["variants_gp_to_p"] = variant_cnt_gp_to_p
    block["variants_gp_to_p_snv"] = snvs_gp_to_p
    block["variants_gp_to_p_indel"] = indels_gp_to_p

    block["variants_d_to_p"] = variant_cnt_d_to_p
    block["variants_d_to_p_snv"] = snvs_d_to_p
    block["variants_d_to_p_indel"] = indels_d_to_p

    block["variants_gp_to_p_filtered"] = variant_cnt_filtered
    block["variants_gp_to_p_filtered_snvs"] = filtered_snp
    block["variants_gp_to_p_filtered_indels"] = filtered_indels

    block["variants_removed"] = variant_cnt_removed
    block["variants_removed_snvs"] = removed_snp
    block["variants_removed_indels"] = removed_indels

    block["variants_ignored"] = variants_ignored_cnt
    block["variants_ignored_snvs"] = variants_ignored_snp
    block["variants_ignored_indels"] = variants_ignored_indel

    block["variants_ignored_027"] = variants_ignored2_cnt
    block["variants_ignored_snvs_027"] = variants_ignored2_snp
    block["variants_ignored_indels_027"] = variants_ignored2_indel
    return block, VARIANTCALL_RESULT.OK

def main(output_folder: Path, input_json: Path, workers: int,
         ignore_variants_in_regions: str, ignore_variants_in_regions2: str) -> tuple[dict, dict]:
    """
    Main worker for variant calling, which processes each transmitted block (see process_block) and aggregates statistics

    Parameters
    ----------
    output_folder               - path to the folder for saving variant files
    input_json                  - path to the input results.json file
    workers                     - number of workers for parallel processing
    ignore_variants_in_regions  - bed file to define regions to ignore in variant query (grandparents)
    ignore_variants_in_regions2 - bed file to define regions to ignore in variant reference (mother)

    Returns
    -------
    debug dictionary of errors encountered, and statistics of SNPs/indels in various stage of the variant calling process
    """

    debug: dict[str, int] = {
        "Error missing input files": 0,
        "Error during variant calling": 0,
        "Error during variant filtering": 0,
        "Error during variant counting": 0,
        "Error during variant counting (snp/indel)": 0,
        "Error during zipping": 0,
    }

    fields: defaultdict[str, int] = defaultdict(int)

    with open(input_json) as f:
        blocks = json.load(f)["data"]

    updated: list[dict] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(process_block, block, output_folder, ignore_variants_in_regions, ignore_variants_in_regions2): block for block in blocks}
        for future in tqdm(as_completed(futures), total=len(blocks), desc="Calling variants..."):
            new_block, status = future.result()

            if status == VARIANTCALL_RESULT.ERR_MISSING_DATA:
                debug["Error missing input files"] += 1
            elif status == VARIANTCALL_RESULT.ERR_VARIANT_CALL:
                debug["Error during variant calling"] += 1
            elif status == VARIANTCALL_RESULT.ERR_ISEC:
                debug["Error during variant filtering"] += 1
            elif status == VARIANTCALL_RESULT.ERR_VARIANT_COUNT:
                debug["Error during variant counting"] += 1
            elif status == VARIANTCALL_RESULT.ERR_VARIANT_COUNT_SNPINDEL:
                debug["Error during variant counting (snp/indel)"] += 1
            elif status == VARIANTCALL_RESULT.ERR_ZIP:
                debug["Error during zipping"] += 1
            else:
                assert status == VARIANTCALL_RESULT.OK, f"unexpected variant calling status '{status}'"
                field_values = ["variants_gp_to_p", "variants_gp_to_p_snv", "variants_gp_to_p_indel",
                                "variants_d_to_p", "variants_d_to_p_snv", "variants_d_to_p_indel",
                                "variants_gp_to_p_filtered", "variants_gp_to_p_filtered_snvs", "variants_gp_to_p_filtered_indels",
                                "variants_removed", "variants_removed_snvs", "variants_removed_indels", "variants_ignored",
                                "variants_ignored_snvs", "variants_ignored_indels",
                                "variants_ignored_027", "variants_ignored_snvs_027", "variants_ignored_indels_027"]

                for field in field_values:
                    if field in new_block and new_block[field] is not None:
                        fields[field] += new_block[field]

            updated.append(new_block)

    with open(input_json, "w") as f:
        json.dump({
            "from_step": "2_call_variants",
            "executed_at": str(datetime.datetime.now()),
            "data": updated
        }, f, indent=4)

    print(f"Updated JSON written to {str(input_json)}")
    return debug, fields


parser = argparse.ArgumentParser()
parser.add_argument("--data", default="blocks/results.json", help="Path to pipeline's data file")
parser.add_argument("--out-dir", default="./variants", help="Output directory for variant calling")
parser.add_argument("--workers", type=int, default=25, help="Number of worker *processes* for variant calling")
parser.add_argument("--ignore-variants-in-regions", default=None, help="Variants to remove")
parser.add_argument("--ignore-variants-in-regions2", default=None, help="Variants to remove")

args = parser.parse_args()

VARIANT_DIR = Path(args.out_dir).resolve()
if not os.path.exists(str(VARIANT_DIR)):
    os.makedirs(str(VARIANT_DIR), exist_ok=True)

BLOCKS_FILE = Path(args.data).resolve()


if __name__ == "__main__":

    print("\n" * 3)
    print('[cyan bold]' + " Switch error detection pipeline ".center(50, "=") + '[/cyan bold]')
    print('[cyan bold]' + "=" + "[STEP 2] Variant calling and filtering".center(48, " ") + "=" + '[/cyan bold]')
    print('[cyan bold]' + "".center(50, "=") + '[/cyan bold]')
    print()

    with open(BLOCKS_FILE, 'r') as f:
        data = json.load(f)
        if data["from_step"] not in ("1_extract_blocks", "2_call_variants"):
            print(f"[red bold]The data file ({BLOCKS_FILE}) does not immediately precede previous step! Aborting.[/red bold]")
            assert False

    debug, fields = main(VARIANT_DIR, BLOCKS_FILE, args.workers, args.ignore_variants_in_regions, args.ignore_variants_in_regions2)
    print("[green bold]Finished![/green bold]")
    print("Errors:")
    for k, v in debug.items():
        print(f"\t{k}: {v}")

    print()
    print("Statistics:")
    for k, v in fields.items():
        print(f"\t{k}: {v}")
