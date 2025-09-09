#!/usr/bin/env python3

"""
make_PC-corrected_LD_command_writer.py

For each text file whose first column contains variants in the form
chrom:pos:ref:alt

this script determines a locus window (min-pos – padding, max-pos + padding)
and writes two separate command files:

1. Commands that extract those windows from one or more plink datasets as binary files
2. Commands that residualize genotypes and compute LD matrices using pysnptools

Usage
-----
python make_PC-corrected_LD_command_writer.py \
    --extract_cmd_out extract_commands.txt \
    --ld_cmd_out ld_commands.txt \
    --matrix_dir /path/to/ld_matrices \
    --plink_dirs /path/CEU,/path/YRI \
    --populations CEU,YRI \
    --variant_lists_file list_of_variant_files.txt \
    --pc_files /path/CEU/pcs.eigenvec,/path/YRI/pcs.eigenvec \
    --residualize_script /path/to/residualize_and_compute_ld.py \
    [--padding 250000] \
    [--num_pcs 10]

The file *list_of_variant_files.txt* should contain one path per line.
PC files should be in standard PLINK eigenvec format with FID, IID, PC1, PC2, ... columns.
"""

import argparse
from pathlib import Path
import re

VAR_RE = re.compile(r"^(\w+):(\d+):")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate separate command files for genotype extraction and LD matrix computation."
    )
    
    parser.add_argument("--extract_cmd_out", required=True,
                       help="Path to write the genotype extraction commands.")
    parser.add_argument("--ld_cmd_out", required=True,
                       help="Path to write the LD matrix computation commands.")
    parser.add_argument("--matrix_dir", required=True,
                       help="Directory where LD matrices should be stored.")
    parser.add_argument("--plink_dirs", required=True,
                       help="Comma-delimited directories for each population's PLINK files.")
    parser.add_argument("--populations", required=True,
                       help="Comma-delimited population names in the same order as plink_dirs.")
    parser.add_argument("--variant_lists_file", required=True,
                       help="File containing paths to variant text files, one per line.")
    parser.add_argument("--pc_files", required=True,
                       help="Comma-delimited paths to PC eigenvec files for each population.")
    parser.add_argument("--residualize_script", required=True,
                       help="Path to the residualize_and_compute_ld.py script.")
    parser.add_argument("--naming_conventions", type=str, default="chr*.plink.maf",
                       help="Comma-delimited patterns for each population's PLINK files.")
    parser.add_argument("--padding", type=int, default=250_000,
                       help="Padding (bp) added upstream & downstream of each locus (default 250 kb).")
    parser.add_argument("--num_pcs", type=int, default=10,
                       help="Number of PCs to regress out (default 10).")

    return parser.parse_args()

def find_plink_root(plink_dir, naming_convention, chrom):
    # Replace * in naming convention with the chromosome number
    pattern = naming_convention.replace("*", str(chrom))
    bed_files = list(Path(plink_dir).glob(pattern + ".bed"))
    
    if not bed_files:
        raise FileNotFoundError(f"No .bed file found for chromosome {chrom} in {plink_dir} with pattern {pattern}")
    
    # Remove .bed extension to get the root
    return str(bed_files[0])[:-4]

def get_window(variant_file, padding):
    chrom = None
    min_pos, max_pos = None, None
    
    with open(variant_file) as fh:
        for line in fh:
            fields = line.strip().split()
            match = VAR_RE.match(fields[0])
            if not match:
                continue
            
            v_chrom, v_pos = match.group(1), int(match.group(2))
            
            if chrom is None:
                chrom = v_chrom
            if chrom != v_chrom:
                raise ValueError(f"Multiple chromosomes in {variant_file}")
            
            min_pos = v_pos if min_pos is None else min(min_pos, v_pos)
            max_pos = v_pos if max_pos is None else max(max_pos, v_pos)
    
    if chrom is None:
        raise ValueError(f"No variants found in {variant_file}")
    
    return chrom, max(min_pos - padding, 0), max_pos + padding

def main():
    # Set tool paths
    plink2 = "/lustre/orion/bif154/proj-shared/mconery/tools/bin/plink2"
    
    # Parse arguments
    args = parse_args()
    
    plink_dirs = [p.strip() for p in args.plink_dirs.split(",") if p.strip()]
    populations = [p.strip() for p in args.populations.split(",") if p.strip()]
    pc_files = [p.strip() for p in args.pc_files.split(",") if p.strip()]
    
    if len(plink_dirs) != len(populations):
        raise ValueError("plink_dirs and populations must have the same length")
    if len(pc_files) != len(populations):
        raise ValueError("pc_files and populations must have the same length")
    
    # Parse naming conventions
    temp = [p.strip() for p in args.naming_conventions.split(",") if p.strip()]
    if len(populations) != len(temp) and len(temp) == 1:
        naming_conventions = [temp[0] for x in populations]
    elif len(populations) != len(temp) and len(temp) != 1:
        raise ValueError("There must be either 1 naming convention or 1 per population.")
    else:
        naming_conventions = temp[:]
    
    # Ensure output dirs exist
    matrix_root = Path(args.matrix_dir)
    for pop in populations:
        (matrix_root / pop).mkdir(parents=True, exist_ok=True)
        (matrix_root / pop / "temp").mkdir(parents=True, exist_ok=True)
    
    extract_commands = []
    ld_commands = []
    
    with open(args.variant_lists_file) as list_fh:
        for variant_file in map(str.strip, list_fh):
            if not variant_file:
                continue
            
            variant_path = Path(variant_file)
            chrom, start, end = get_window(variant_path, args.padding)
            locus_name = variant_path.stem  # base filename without extension
            
            for plink_dir, pop, naming_conv, pc_file in zip(plink_dirs, populations, naming_conventions, pc_files):
                plink_root = find_plink_root(plink_dir, naming_conv, chrom)
                out_prefix = matrix_root / pop / locus_name
                temp_prefix = matrix_root / pop / "temp" / locus_name
                
                # Extraction command - NOW USING --make-bed FOR BINARY FORMAT
                extract_cmd = f"{plink2} --bfile {plink_root} --chr {chrom} --from-bp {start} --to-bp {end} --make-bed --out {temp_prefix}_extracted"
                extract_commands.append(extract_cmd)
                
                # LD computation command - UPDATED FOR BINARY INPUT
                ld_cmd = f"python3 {args.residualize_script} "
                ld_cmd = ld_cmd + f"--geno_file {temp_prefix}_extracted "  # NO EXTENSION - pysnptools handles .bed/.bim/.fam
                ld_cmd = ld_cmd + f"--pc_file {pc_file} "
                ld_cmd = ld_cmd + f"--out_prefix {out_prefix} "
                ld_cmd = ld_cmd + f"--num_pcs {args.num_pcs} "
                ld_cmd = ld_cmd + f"--locus_name {locus_name} "
                ld_cmd = ld_cmd + f"--population {pop} "
                ld_cmd = ld_cmd + f"--chrom {chrom} "
                ld_cmd = ld_cmd + f"--start {start} "
                ld_cmd = ld_cmd + f"--end {end} "
                ld_cmd = ld_cmd + "--cleanup_temp"
                
                ld_commands.append(ld_cmd)
    
    # Write extraction commands to file
    with open(args.extract_cmd_out, "w") as out_fh:
        for cmd in extract_commands:
            out_fh.write(cmd + "\n")
    
    # Write LD computation commands to file
    with open(args.ld_cmd_out, "w") as out_fh:
        for cmd in ld_commands:
            out_fh.write(cmd + "\n")
    
    print(f"Wrote {len(extract_commands)} extraction commands to {args.extract_cmd_out}")
    print(f"Wrote {len(ld_commands)} LD computation commands to {args.ld_cmd_out}")
    print(f"Each LD command will residualize {args.num_pcs} PCs before computing LD matrices")
    print("NOTE: Now using binary PLINK format (.bed/.bim/.fam) for faster processing with pysnptools")

if __name__ == "__main__":
    main()