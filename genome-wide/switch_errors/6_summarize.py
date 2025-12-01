#!/usr/bin/env python3
import os
import subprocess
import argparse
from pathlib import Path
from tqdm import tqdm
import json
from collections import defaultdict
from bisect import bisect_left

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

def fetch_variants_from_vcf(vcf_path: str) -> list[tuple]:
    """
    Reads variants from .vcf file

    Parameters
    ----------
    vcf_path - path to .vcf file

    Returns
    -------
    list of variants (defined by contig, region start, region end, and identifier)
    """

    result: list[tuple] = []
    with open(vcf_path, 'r') as fd:
        for line in fd:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue

            contig, pos_from, pos_to, ident, *_ = parts
            result.append((contig, pos_from, pos_to, ident))
    return result

def collect_variants(input_json: Path) -> tuple[list[tuple], int, int]:
    """
    Collect all lifted-over variants from mother-reference bed files (grandparent-reference will be added in the next step)

    Parameters
    ----------
    input_json - transmitted blocks definition json (results.json)

    Returns
    -------
    list of all variants, and number of blocks with and without variants
    """

    with open(input_json, "r") as in_f:
        data = json.load(in_f)["data"]

    bed_files: list[str] = []
    for entry in data:
        if "mother_liftover" in entry and entry["mother_liftover"]:
            bed_files.append(entry["mother_liftover"])
    blocks_with_variants, blocks_without_variants = len(bed_files), len(data) - len(bed_files)

    variants = []
    for file in tqdm(bed_files, desc="Collecting variants..."):
        variants.extend(fetch_variants_from_vcf(file))
    return variants, blocks_with_variants, blocks_without_variants

def build_variant_dict(variants: list[tuple]) -> tuple[dict[str, str | int], int]:
    """
    Transforms variant list (from fetch_variants_from_vcf) into variant dictionary, removing duplicates

    Parameters
    ----------
    variants - variant list from fetch_variants_from_vcf

    Returns
    -------
    dictionary of variants, as well as number of duplicates removed
    """

    variants_dict: dict = {
        info: {
            "contig": contig,
            "pos_from": pos_from,
            "pos_to": pos_to,
            "mother_contig": contig,
            "mother_region": info.split(',')[0].split(':')[1],
            "mother_allele": info.split(',')[0].split(':')[2],
            "grandparent_contig": info.split(',')[1].split(':')[0],
            "grandparent_region": info.split(',')[1].split(':')[1],
            "grandparent_allele": info.split(',')[1].split(':')[2],
        } for (contig, pos_from, pos_to, info) in variants
    }

    exists: set[tuple] = set()
    duplicates: int = 0
    for (contig, pos_from, pos_to, info) in variants:
        if info in exists:
            duplicates += 1
        else:
            exists.add(info)

    return variants_dict, duplicates

def extract_lifted_alleles(data: dict, variants_dict: dict) -> tuple[int, int]:
    """
    Associates grandparent-reference lifted-over variants where available

    Parameters
    ----------
    data            - transmitted blocks definition json (results.json)
    variants_dict   - variant list from build_variant_dict

    Returns
    -------
    number of failed associations (liftover bed is missing) and number of failed associations (grandparents contain
    variant not defined in mother-reference liftover)
    """

    error_missing_data = failed_missing_key = 0
    for entry in tqdm(data, desc="Extracting lifted alleles..."):
        for liftover in ('extracted_bed_mother', 'extracted_bed_grandparent'):
            lifted_bed = entry.get(liftover)
            if not lifted_bed or not os.path.exists(lifted_bed) or os.stat(lifted_bed).st_size == 0:
                error_missing_data += 1
                continue  # not fatal error, one successful liftover can still show switch error

            with open(lifted_bed, "r") as in_f:
                for line in in_f:
                    if line.startswith("#"): continue
                    parts: list[str] = line.strip().split('\t')
                    if len(parts) < 4: continue

                    lifted_contig, pos_from, pos_to, info, *_ = parts
                    parts: list[str] = info.strip().split(',')
                    if len(parts) < 3: continue
                    info_mother, info_grandparent, info_lifted = info.split(',')
                    variant_ident = info_mother + ',' + info_grandparent

                    if variant_ident not in variants_dict:
                        failed_missing_key += 1
                        continue

                    is_mother_other = "PAN027" in lifted_contig
                    extract_field = "mother_other" if is_mother_other else "grandparent_other"
                    variants_dict[variant_ident][f"{extract_field}_contig"] = info_lifted.split(':')[0]
                    variants_dict[variant_ident][f"{extract_field}_region"] = info_lifted.split(':')[1]
                    variants_dict[variant_ident][f"{extract_field}_allele"] = info_lifted.split(':')[2]
    return error_missing_data, failed_missing_key

def associate_lifted_alleles(input_json: Path, variants: list[tuple]) -> tuple[dict, list[tuple], int, int, int, int, int]:
    """
    Associates grandparent-reference lifted alleles onto already-loaded mother-reference lifted alleles

    Parameters
    ----------
    input_json - transmitted blocks definition (results.json)
    variants   - existing variants extracted from mother-reference bed (variant liftover) files

    Returns
    -------
    dictionary of associated variants, list of variants in input, number of fully/partially associated variants,
    number of variants where key is missing, number of variants where data is missing, and number of duplicates
    """

    with open(input_json, "r") as in_f:
        data = json.load(in_f)["data"]

    variants_dict, duplicates = build_variant_dict(variants)
    error_missing_data = failed_missing_key = extract_lifted_alleles(data, variants_dict)

    fully_assoc = partial_assoc = 0
    for ident, data in variants_dict.items():
        if "mother_other_contig" in data and "grandparent_other_contig" in data:
            fully_assoc += 1
        elif "mother_other_contig" in data or "grandparent_other_contig" in data:
            partial_assoc += 1

    return variants_dict, variants, fully_assoc, partial_assoc, failed_missing_key, error_missing_data, duplicates

def select_switch_error_variants(variants: dict, process_indels: bool):
    """
    Processes variants from association, detecting possible switch errors if variant can be explained by using
    mother's other haplotype (switch error in mother) or grandparents' other haplotype (switch error in grandparent)

    Parameters
    ----------
    variants        - associated variants dictionary
    process_indels  - if True, process indels

    Returns
    -------
    updated variant dictionary with switch-error-prone variants, and statistics on:
     - number of variants prone to be switch error, explained by mother's other haplotype
     - number of variants prone to be switch error, explained by grandparents' other haplotype
     - number of variants prone to be switch error, explained by mother's and grandparents' other haplotype (cross)
     - number of variants not prone to be switch error
     ... breakdown of SNVs and indels for each above
    """

    result: dict = {}
    mother_switch = grandparent_switch = both_switch = none_switch = 0
    m_snv = m_ind = gp_snv = gp_ind = b_snv = b_ind = n_snv = n_ind = 0

    for ident, data in variants.items():
        m_allele, gp_allele, m_other_allele, gp_other_allele = (data["mother_allele"],
                                                                data["grandparent_allele"],
                                                                data["mother_other_allele"] if "mother_other_allele" in data else None,
                                                                data["grandparent_other_allele"] if "grandparent_other_allele" in data else None)

        is_indel: bool = (len(m_allele) > 1
                         or len(gp_allele) > 1
                         or (m_other_allele is not None and len(m_other_allele) > 1)
                         or (gp_other_allele is not None and len(gp_other_allele) > 1))

        if is_indel and not process_indels:
            continue

        is_switch_error, switch_error_kind = False, None
        if m_other_allele is not None and m_other_allele == gp_allele:
            is_switch_error, switch_error_kind = True, "SWITCH ERROR IN MOTHER"
            mother_switch += 1
            if is_indel:
                m_ind += 1
            else:
                m_snv += 1
        elif gp_other_allele is not None and gp_other_allele == m_allele:
            is_switch_error, switch_error_kind = True, "SWITCH ERROR IN GRANDPARENT"
            grandparent_switch += 1
            if is_indel:
                gp_ind += 1
            else:
                gp_snv += 1
        elif gp_other_allele is not None and m_other_allele is not None and gp_other_allele == m_other_allele:
            is_switch_error, switch_error_kind = True, "SWITCH ERROR IN BOTH"
            both_switch += 1
            if is_indel:
                b_ind += 1
            else:
                b_snv += 1
        else:
            none_switch += 1
            if is_indel:
                n_ind += 1
            else:
                n_snv += 1

        if not is_switch_error:
            continue

        data["switch_error_kind"] = switch_error_kind
        data["indel"] = is_indel
        result[ident] = data.copy()
    return result, mother_switch, grandparent_switch, both_switch, none_switch, m_snv, m_ind, gp_snv, gp_ind, b_snv, b_ind, n_snv, n_ind

def save_sw_prone_variants(sw_variants_dict: dict, out_path: Path) -> None:
    """
    Saves switch-error-prone variants into vcf file

    Parameters
    ----------
    sw_variants_dict - switch-error-prone variants dictionary (from select_switch_error_variants)
    out_path - path to the output vcf file

    Returns
    -------
    None
    """

    with open(out_path, "w") as vcf:
        # VCF header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write('##INFO=<ID=MOTHER,Number=0,Type=Flag,Description="Variant at mother">\n')
        vcf.write('##INFO=<ID=INHERITANCE,Number=0,Type=Flag,Description="Variant inherited from">\n')
        vcf.write('##INFO=<ID=ERROR_TYPE,Number=0,Type=Flag,Description="Type of possible switch error">\n')
        vcf.write('##INFO=<ID=REASON,Number=0,Type=Flag,Description="Reason for switch error">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for entry in sw_variants_dict.values():
            contig = entry["contig"]
            pos = entry["pos_from"]
            ref = entry["mother_allele"]
            alt = entry["grandparent_allele"]

            mother_line = entry["mother_contig"] + ":" + entry["mother_region"] + ":" + entry["mother_allele"]
            inherited_line = entry["grandparent_contig"] + ":" + entry["grandparent_region"] + ":" + entry[
                "grandparent_allele"]
            error_type = entry["switch_error_kind"]

            if error_type == "SWITCH ERROR IN MOTHER":
                reason = entry["mother_other_contig"] + ":" + entry["mother_other_region"] + ":" + entry[
                    "mother_other_allele"] + "_explains_mutation"
            elif error_type == "SWITCH ERROR IN GRANDPARENT":
                reason = entry["grandparent_other_contig"] + ":" + entry["grandparent_other_region"] + ":" + entry[
                    "grandparent_other_allele"] + "_explains_mutation"
            elif error_type == "SWITCH ERROR IN BOTH":
                reason = "liftover_of_both_haplotypes_explains_mutation"
            else:
                assert False, "Unknown error type"

            vcf.write(
                f"{contig}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tMOTHER={mother_line};INHERITANCE={inherited_line};ERROR_TYPE={error_type.replace(' ', '_')};REASON={reason}\n")

def cluster_positions(positions: list[int], max_dist: int, min_members: int) -> list[tuple[int, int, int]]:
    """
    Cluster positions that are within max_dist of each other

    Parameters
    ----------
    positions   - list of positions
    max_dist    -
    min_members

    Returns
    -------
    list of (start, end, count)
    """

    if not positions:
        return []

    # sort positions by genomic coordinate
    positions: list[int] = sorted(positions)
    clusters: list[tuple[int, int, int]] = []
    cluster_start: int = positions[0]
    cluster_end: int = positions[0]
    count: int = 1

    for pos in positions[1:]:
        if pos - cluster_end <= max_dist:
            cluster_end = pos
            count += 1
        else:
            if count >= min_members: # finalize current cluster if large enough
                clusters.append((cluster_start, cluster_end, count))
            cluster_start = pos
            cluster_end = pos
            count = 1

    if count >= min_members: # finalize last cluster
        clusters.append((cluster_start, cluster_end, count))

    return clusters

def cluster_switch_variants(sw_variants_dict: dict, bed_out_path: str, max_dist: int = 20000, min_members: int = 10):
    """
    Cluster switch-error-prone variants by contig and haplotype.

    Parameters
    ----------
    sw_variants_dict    - Dict of entries with keys:
                            - contig: "PAN027.chr1.maternal"
                            - pos_from, pos_to: int
    bed_out_path        - Path to save resulting BED file
    max_dist            - Maximum distance (bp) between consecutive variants to be clustered
    min_members         - Minimum number of variants per cluster to keep

    Returns
    -------
    clusters
    """

    grouped: defaultdict[str, list[int]] = defaultdict(list)
    for entry in sw_variants_dict.values():
        contig: str = entry["contig"]
        start: int = int(entry["pos_from"])
        end: int = int(entry["pos_to"])
        pos: int = (start + end) // 2
        grouped[contig].append(pos)

    clusters_all: list[tuple[str, int, int, int]] = []

    for contig, positions in grouped.items():
        asm, chrom, hap = contig.split(".")
        clusters = cluster_positions(positions, max_dist, min_members)
        for start, end, count in clusters:
            clusters_all.append((contig, start, end, count))

    with open(bed_out_path, "w") as bed:
        bed.write(f"#Clusters of switch error-prone SNPs "
                  f"(clustering of SNPs within {max_dist // 1000}kb, "
                  f"min {min_members} variants per cluster)\n")
        for contig, start, end, count in clusters_all:
            bed.write(f"{contig}\t{start - 1}\t{start - 1 + (end - start)}\t\"{count} switch error-prone variants\"\n")

    print(f"[INFO] Wrote {len(clusters_all)} clusters to {bed_out_path}")
    return clusters_all

def summarize_variants_in_clusters(variants_filtered: list, clusters: list[tuple[str, int, int, int]]):
    """
    Count how many variants fall within or outside of clusters, separated by switch_error_kind

    Parameters
    ----------
    variants_filtered - variants to analyze with regards to clusters
    clusters - output of cluster_switch_variants

    Returns
    -------
    statistics of variants within/outside clusters for mother, grandparent, and cross switch-error-events
    """

    cluster_dict: dict = {}
    for contig, start, end, count in clusters:
        cluster_dict.setdefault(contig, []).append((start, end))
    for contig in cluster_dict:
        cluster_dict[contig].sort()

    # prepare results structure
    kinds: list[str] = ["SWITCH ERROR IN MOTHER", "SWITCH ERROR IN GRANDPARENT", "SWITCH ERROR IN BOTH"]
    results: dict = {k: {"inside": 0, "outside": 0} for k in kinds}
    results2: dict = {k: {True: {"inside": 0, "outside": 0}, False: {"inside": 0, "outside": 0}} for k in kinds}

    for entry in variants_filtered:
        contig = entry["contig"]
        pos = int(entry["pos_from"])
        kind = entry["switch_error_kind"]
        indel = entry["indel"]
        if kind not in results:
            continue

        in_cluster = False
        if contig in cluster_dict:
            clusters_list = cluster_dict[contig]
            i = bisect_left(clusters_list, (pos, pos))
            for j in (i - 1, i):
                if 0 <= j < len(clusters_list):
                    start, end = clusters_list[j]
                    if start <= pos <= end:
                        in_cluster = True
                        break

        if in_cluster:
            results[kind]["inside"] += 1
            results2[kind][indel]["inside"] += 1
        else:
            results[kind]["outside"] += 1
            results2[kind][indel]["outside"] += 1

    return results, results2


parser = argparse.ArgumentParser()
parser.add_argument("--data", default="blocks/results.json", help="Path to pipeline's data file")
parser.add_argument("--threads", type=int, default=20, help="Number of parallel workers")
parser.add_argument("--indels", action='store_true', help="Include trustworthy indels")
parser.add_argument("--plots-dir", default="./plots", help="Path to output plots")
parser.add_argument("--genomes-dir", default="./genomes", help="Directory to genomic data")
parser.add_argument("--out_vcf", default="plots/sw_prone_variants.vcf", help="Path to save filtered variants")

parser.add_argument("--skip-clustering", action='store_true', help="Skip clustering step")
parser.add_argument("--out_bed", default="plots/sw_prone_regions.bed", help="Path to save clustered regions")
parser.add_argument("--cluster-max-dist", type=int, default=20000, help="Clustering distance for sw-prone variants")
parser.add_argument("--cluster-min-members", type=int, default=8, help="Minimum number of variants per cluster")

args = parser.parse_args()

BLOCKS_FILE = Path(args.data).resolve()
GENOME_DIR = Path(args.genomes_dir).resolve()
PLOTS = Path(args.plots_dir).resolve()
if not os.path.exists(PLOTS):
    os.makedirs(PLOTS, exist_ok=True)
VARIANTS_OUT_VCF = Path(args.out_vcf).resolve()
REGIONS_OUT_BED = Path(args.out_bed).resolve()

if __name__ == "__main__":

    print("\n" * 3)
    print(" Switch error detection pipeline ".center(50, "="))
    print("=" + "[STEP 6] Summarize".center(48, " ") + "=")
    print("".center(50, "="))
    print()

    with open(BLOCKS_FILE, 'r') as f:
        data = json.load(f)
        if data["from_step"] != "5_extract_lifted":
            print(f"The data file ({BLOCKS_FILE}) does not immediately precede previous step! Aborting.")
            assert False

    print("Fetching all variants...")
    variants, blocks_with_variants, blocks_without_variants = collect_variants(BLOCKS_FILE)
    print(f"Processed {len(data['data'])} blocks, {blocks_with_variants} succeeded, {blocks_without_variants} missing.")
    print(f"Fetched {len(variants)} for association...")
    print()

    print("Associating lifted alleles...")
    variants_associated, variants, fully_assoc, partial_assoc, failed_missing_key, error_missing_data, duplicates = associate_lifted_alleles(
        BLOCKS_FILE, variants)
    print(f"Processed {len(variants_associated)} variants:\n"
          f"\tskipped {duplicates} duplicates\n"
          f"\tboth lifted alleles associated for {fully_assoc} variants\n"
          f"\tone of the lifted allele associated for {partial_assoc} variants\n"
          f"\tno alleles were lifted for {len(variants_associated) - fully_assoc - partial_assoc} variants\n")
    print()

    print("Filtering switch error-prone variants...")
    print("[WARN]", "Indels are included!" if args.indels else "All indels are being skipped")
    variants_filtered, mother_switch, grandparent_switch, both_switch, none_switch, m_snv, m_ind, gp_snv, gp_ind, b_snv, b_ind, n_snv, n_ind = select_switch_error_variants(
        variants_associated, args.indels)
    print(f"Detected {len(variants_filtered)} switch error-prone variants...")
    print(f"... {mother_switch} variants are explained by switch in mother's haplotype ({m_snv} SNVs, {m_ind} indels)")
    print(
        f"... {grandparent_switch} variants are explained by switch in grandparent's haplotype ({gp_snv} SNVs, {gp_ind} indels)")
    print(f"... {both_switch} variants are explained by switch in both haplotypes ({b_snv} SNVs, {b_ind} indels)")
    print(f"... {none_switch} variants are NOT explained by switch in any haplotype ({n_snv} SNVs, {n_ind} indels)")
    print()

    print("Saving switch error-prone variants...")
    save_sw_prone_variants(variants_filtered, PLOTS / "sw_prone_variants.vcf")
    print()

    if args.skip_clustering:
        print("[WARN] Clustering is being skipped...")
        exit(0)

    print("=== BEGIN CLUSTERING ANALYSIS ===")
    print("Clustering sw-prone variants...")
    cls = cluster_switch_variants(variants_filtered, str(REGIONS_OUT_BED), args.cluster_max_dist,
                                  args.cluster_min_members)

    sum_dist = sum_cnt = 0
    for (contig, start, end, count) in cls:
        print(f"... {contig}\t{start}\t{end}\t{count} sw-prone variants")
        sum_dist += end - start
        sum_cnt += count
    print(f"... average cluster size: {round(sum_dist / len(cls), 2)}")
    print(f"... average cluster sw-prone variant count: {round(sum_cnt / len(cls), 2)}")

    res1, res2 = summarize_variants_in_clusters(list(variants_filtered.values()), cls)
    for k, v in res1.items():
        print(k, v)

    for k, ind in res2.items():
        for k, v in ind.items():
            print(k, "indel" if k else "SNVs", v)

    from varplot import VariantSwitchErrorPlotter
    plotter = VariantSwitchErrorPlotter(variants_associated, variants_filtered, cls)
    windows_df = plotter.compute_all(chromosomes=None, haplotypes=None, workers=8)
    plotter.plot(windows_df, out_file='variants_vs_switch_errors.png')