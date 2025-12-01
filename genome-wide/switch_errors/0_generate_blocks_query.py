import glob
import json
import datetime
import os.path
from collections import defaultdict
from pathlib import Path
import argparse
from scipy import stats
from rich import print

def set_entry_field(entry: dict[str, str | int], kind: str, line: str) -> int:
    """
    Parses one part of the line from transmitted blocks definition (grandparent, mother, or granddaughter) and saves
    it to results dictionary

    Parameters
    ----------
    entry - result dictionary to be updated
    kind  - type of field to set, either grandparent, parent, or daughter
    line  - transmitted block definition, e.g. 'PAN028.chr1.haplotype2:172160018-173571062'

    Returns
    -------
    length of region defined by this transmitted block
    """

    assembly_part, region_part = line.split(":")
    region_from, region_to = region_part.split("-")

    entry[f"{kind}_assembly"] = assembly_part
    entry[f"{kind}_from"] = int(region_from)
    entry[f"{kind}_to"] = int(region_to)
    return int(region_to) - int(region_from)

def run() -> None:
    """
    Processes transmitted blocks files and generates results dictionary for further script's input

    Returns
    -------
    None (creates results.json)
    """
    query_ctr: int = 0
    sizes: list[int] = []
    for file in SRC_FILES:
        with open(file, 'r') as f:
            for line in f:
                gp_part, parent_part, daughter_part = line.strip().split("\t")
                entry: dict[str, str | int] = {}
                sz1 = set_entry_field(entry, "grandparent", gp_part)
                sz2 = set_entry_field(entry, "parent", parent_part)
                sz3 = set_entry_field(entry, "daughter", daughter_part)
                entry["id"] = query_ctr
                queries.append(entry)

                sizes.extend([sz1, sz2, sz3])
                query_ctr += 1
                debug[file] += 1

    with open(TRG_FILE, 'w') as f:
        json.dump({
            "from_step": "0_generate_blocks_query",
            "executed_at": str(datetime.datetime.now()),
            "data": queries
        }, f, indent=4)

    print(f"Transmitted blocks' size statistic: {stats.describe(sizes)}")
    print()

parser = argparse.ArgumentParser()
parser.add_argument("--input", default="./blocks", help="Path to input directory with transmitted blocks definition (.tsv files)")
parser.add_argument("--data", default="blocks/results.json", help="Path to pipeline's data file")
args = parser.parse_args()

WORKING_DIRECTORY: Path = Path(args.input)
assert os.path.exists(WORKING_DIRECTORY), f"Working directory '{WORKING_DIRECTORY}' does not exist! This script can't do anything!"

SRC_FILES: list[str] = glob.glob(str(WORKING_DIRECTORY / "*.tsv"))
TRG_FILE: Path = Path(args.data).resolve()
debug: defaultdict[str, int] = defaultdict(int)
queries: list[dict[str, str | int]] = []

print("\n" * 3)
print('[bold cyan]' + " Switch error detection pipeline ".center(50, "=") + '[/bold cyan]')
print('[bold cyan]' + "=" + "[STEP 0] Block definition generation".center(48, " ") + "=" + '[/bold cyan]')
print('[bold cyan]' + "".center(50, "=") + '[/bold cyan]')
print()

print(f"[italic]Files to process: {SRC_FILES}[/italic]")
print(f"[italic]Will save to: {TRG_FILE}[/italic]")
print()
run()

print("[green bold]Finished![/green bold]")
print(f"Processed {len(debug)} source files:")
for file, cnt in debug.items():
    print(f"\t{os.path.basename(file)}: {cnt}")
print(f"Total number of transmitted block definitions (gp -> p -> gd): {len(queries)}")
