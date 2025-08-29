'''
select_causal_variants_ld.py

Select causal variants for GWAS simulation by tiling the genome into loci
and selecting variants based on LD constraints across multiple populations.

***Run in Conda simulate_pheno environment***
'''

import argparse
import glob
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple, Dict
import pandas as pd
import numpy as np

# -----------------------------------------------------------------------------#
# Helper functions
# -----------------------------------------------------------------------------#

def parse_bed_file(bed_path: Path) -> List[Tuple[str, int, int]]:
    """
    Parse BED file to get allowed genomic regions.
    Returns list of (chromosome, start, end) tuples.
    """
    regions = []
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                regions.append((chrom, start, end))
    return regions

def tile_genome(regions: List[Tuple[str, int, int]], locus_size: int) -> List[Tuple[str, int, int]]:
    """
    Tile the allowed genomic regions into loci of specified size.
    Returns list of (chromosome, start, end) for each locus.
    """
    loci = []
    for chrom, region_start, region_end in regions:
        pos = region_start
        while pos + locus_size <= region_end:
            loci.append((chrom, pos, pos + locus_size))
            pos += locus_size
        # Add remaining region if it's substantial (>= 50% of locus_size)
        if region_end - pos >= locus_size * 0.5:
            loci.append((chrom, pos, region_end))
    return loci

def load_bim_files(input_dir: Path, pattern: str) -> Dict[str, pd.DataFrame]:
    """
    Read BIM files matching the pattern and return a dict keyed by chromosome.
    """
    files = glob.glob(str(os.path.join(input_dir, pattern)))
    if not files:
        sys.exit(f"No BIM files matching pattern '{pattern}' found in {input_dir}")
    
    data: Dict[str, pd.DataFrame] = {}
    for f in files:
        chrom = infer_chromosome_from_filename(f)
        df = pd.read_csv(
            f,
            sep=r"\s+",
            header=None,
            names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
            dtype={"CHR": str, "SNP": str, "BP": int},
        )
        data[chrom] = df.sort_values("BP").reset_index(drop=True)
    
    return data

def infer_chromosome_from_filename(path: str) -> str:
    """
    Extract chromosome from filename (assumes 'chr' in name).
    """
    import re
    basename = os.path.basename(path)
    # Try to find chr followed by number or X/Y
    m = re.search(r"chr(\d+|X|Y)", basename, re.IGNORECASE)
    if m:
        return f"chr{m.group(1)}"
    
    # Try to find just number
    m = re.search(r"(\d+)\.plink", basename)
    if m:
        return f"chr{m.group(1)}"
    
    sys.exit(f"Cannot infer chromosome number from file name '{path}'")

def get_locus_variants(locus: Tuple[str, int, int], bim_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Get all variants within a locus from the quality-filtered BIM data.
    """
    chrom, start, end = locus
    if chrom not in bim_data:
        return pd.DataFrame(columns=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    
    df = bim_data[chrom]
    mask = (df["BP"] >= start) & (df["BP"] <= end)
    return df[mask].copy()

def calculate_ld_matrix(variants: pd.DataFrame, pop_dir: Path, naming_convention: str, temp_dir: str) -> np.ndarray:
    """
    Calculate LD matrix for variants using PLINK2.
    Returns correlation matrix as numpy array.
    """
    if len(variants) == 0:
        return np.array([])
    
    chrom = variants.iloc[0]["CHR"]
    
    # Find the PLINK files for this chromosome in this population
    plink_pattern = naming_convention.replace("*", chrom.replace("chr", ""))
    
    # Try different extensions for PLINK files
    plink_files = []
    file_type = "bfile"
    for ext in [".bed", ".pgen"]:
        files = glob.glob(str(os.path.join(pop_dir, plink_pattern + ext)))
        if files:
            plink_files = files
            file_type = "pfile" if ext == ".pgen" else "bfile"
            break
    
    if not plink_files:
        print(f"Warning: No PLINK files found for {chrom} in {pop_dir} with pattern {plink_pattern}", file=sys.stderr)
        # Return identity matrix as fallback (no LD constraints)
        return np.eye(len(variants))
    
    plink_prefix = plink_files[0].replace(".bed", "").replace(".pgen", "")
    
    # Create variant list file
    variant_file = os.path.join(temp_dir, f"variants_{chrom}_{random.randint(1000, 9999)}.txt")
    with open(variant_file, 'w') as f:
        for snp in variants["SNP"]:
            f.write(f"{snp}\n")
    
    # Calculate LD using PLINK2
    output_prefix = os.path.join(temp_dir, f"ld_{chrom}_{random.randint(1000, 9999)}")
    
    try:
        cmd = [
            "plink2",
            f"--{file_type}", plink_prefix,
            "--extract", variant_file,
            "--r-unphased", "square",
            "--out", output_prefix
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"PLINK2 error: {result.stderr}", file=sys.stderr)
            return np.eye(len(variants))
        
        # PLINK2 outputs correlation matrix as .unphased.vcor1 file
        ld_file = output_prefix + ".unphased.vcor1"
        if not os.path.exists(ld_file):
            return np.eye(len(variants))
            
        ld_matrix = np.loadtxt(ld_file, delimiter='\t')
        if ld_matrix.ndim == 1:
            ld_matrix = ld_matrix.reshape(1, 1)
        
        return ld_matrix
    
    except Exception as e:
        print(f"Error calculating LD: {e}", file=sys.stderr)
        return np.eye(len(variants))
    
    finally:
        # Clean up temporary files
        cleanup_files = [variant_file]
        for ext in [".unphased.vcor1", ".unphased.vcor1.vars", ".log"]:
            cleanup_files.append(output_prefix + ext)
        
        for file_path in cleanup_files:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
            except:
                pass

def select_variants_by_ld(variants: pd.DataFrame, ld_matrices: List[np.ndarray], 
                         num_variants: int, max_ld: float, rng: random.Random) -> List[str]:
    """
    Select variants ensuring LD constraints are met across all populations.
    """
    if len(variants) == 0 or num_variants == 0:
        return []
    
    if len(variants) < num_variants:
        return variants["SNP"].tolist()
    
    n_variants = len(variants)
    selected_indices = []
    
    # Select first variant randomly
    first_idx = rng.randrange(n_variants)
    selected_indices.append(first_idx)
    
    # Select remaining variants
    for _ in range(num_variants - 1):
        candidate_indices = []
        
        for i in range(n_variants):
            if i in selected_indices:
                continue
                
            # Check LD with all previously selected variants across all populations
            valid = True
            for selected_idx in selected_indices:
                for ld_matrix in ld_matrices:
                    if ld_matrix.size == 0:
                        continue
                    if i < ld_matrix.shape[0] and selected_idx < ld_matrix.shape[1]:
                        if abs(ld_matrix[i, selected_idx]) > max_ld:
                            valid = False
                            break
                if not valid:
                    break
            
            if valid:
                candidate_indices.append(i)
        
        if not candidate_indices:
            # If no valid candidates, we're done
            break
            
        # Randomly select from valid candidates
        selected_idx = rng.choice(candidate_indices)
        selected_indices.append(selected_idx)
    
    return variants.iloc[selected_indices]["SNP"].tolist()

def process_locus_with_retries(locus: Tuple[str, int, int], bim_data: Dict[str, pd.DataFrame], 
                              population_dirs: List[Path], plink_naming_convention: str, 
                              variants_per_set: int, max_ld: float, temp_dir: str, 
                              rng: random.Random, max_attempts: int = 10) -> Tuple[bool, List[str]]:
    """
    Process a single locus with retry logic for variant selection.
    Returns (success, selected_variants).
    """
    # Get variants in this locus
    locus_variants = get_locus_variants(locus, bim_data)
    
    if len(locus_variants) == 0:
        return False, []
    
    if len(locus_variants) < variants_per_set:
        print(f"Warning: Locus {locus[0]}:{locus[1]}-{locus[2]} has only {len(locus_variants)} variants, need {variants_per_set}", file=sys.stderr)
        # Still try to select what we can
    
    # Calculate LD matrices for each population
    ld_matrices = []
    for pop_dir in population_dirs:
        ld_matrix = calculate_ld_matrix(locus_variants, pop_dir, plink_naming_convention, temp_dir)
        ld_matrices.append(ld_matrix)
    
    # Try to select variants with retries
    for attempt in range(max_attempts):
        selected_variants = select_variants_by_ld(
            locus_variants, ld_matrices, variants_per_set, max_ld, rng
        )
        
        if len(selected_variants) == variants_per_set:
            return True, selected_variants
        elif len(selected_variants) > 0:
            # Got some variants but not the full set
            if attempt == max_attempts - 1:  # Last attempt
                print(f"Warning: Only found {len(selected_variants)}/{variants_per_set} variants after {max_attempts} attempts", file=sys.stderr)
                return True, selected_variants
        # If no variants selected, try again
    
    return False, []

def write_variant_file(out_dir: Path, set_idx: int, variant_ids: List[str]) -> None:
    """
    Write variant IDs to output file.
    """
    filename = f"Pheno{set_idx}.{len(variant_ids)}_causal.variant_set.txt"
    out_path = os.path.join(out_dir, filename)
    
    with open(out_path, "w") as fout:
        for vid in variant_ids:
            fout.write(f"{vid}\n")

# -----------------------------------------------------------------------------#
# Main
# -----------------------------------------------------------------------------#

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select causal variant sets based on genomic tiling and LD constraints."
    )
    
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Directory containing quality-filtered BIM files.",
    )
    
    parser.add_argument(
        "--naming-convention",
        default="chr*.plink.maf_0.05.shared.bim",
        help="Filename pattern with '*' as chromosome placeholder for BIM files.",
    )
    
    parser.add_argument(
        "--bed-file",
        required=True,
        type=Path,
        help="BED file specifying allowed genomic regions for tiling.",
    )
    
    parser.add_argument(
        "--population-dirs",
        required=True,
        nargs='+',
        type=Path,
        help="List of directories containing population-specific PLINK files.",
    )
    
    parser.add_argument(
        "--plink-naming-convention",
        required=True,
        help="Naming convention for chromosome-level PLINK files (use * for chromosome).",
    )
    
    parser.add_argument(
        "--max-ld",
        type=float,
        required=True,
        help="Maximum absolute LD correlation between causal variants.",
    )
    
    parser.add_argument(
        "--num-sets",
        type=int,
        default=1000,
        help="Number of loci to select for causal variant sets.",
    )
    
    parser.add_argument(
        "--variants-per-set",
        type=int,
        default=1,
        help="Number of causal variants per locus.",
    )
    
    parser.add_argument(
        "--max-distance",
        type=int,
        default=500_000,
        help="Size of each genomic locus (in bp) for tiling.",
    )
    
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory to write output files.",
    )
    
    parser.add_argument(
        "--seed",
        type=int,
        default=5,
        help="Random seed for reproducibility.",
    )
    
    parser.add_argument(
        "--max-locus-attempts",
        type=int,
        default=10,
        help="Maximum attempts to find variants in a locus before selecting a new locus.",
    )
    
    args = parser.parse_args()
    
    # Validation
    if not args.bed_file.exists():
        sys.exit(f"BED file not found: {args.bed_file}")
    
    for pop_dir in args.population_dirs:
        if not pop_dir.exists():
            sys.exit(f"Population directory not found: {pop_dir}")
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    
    # Parse BED file and tile genome
    print("Parsing BED file and tiling genome...", file=sys.stderr)
    allowed_regions = parse_bed_file(args.bed_file)
    all_loci = tile_genome(allowed_regions, args.max_distance)
    print(f"Created {len(all_loci)} genomic loci", file=sys.stderr)
    
    # Load quality-filtered BIM data
    print("Loading BIM files...", file=sys.stderr)
    bim_data = load_bim_files(args.input_dir, args.naming_convention)
    
    # Check if we have enough loci
    if len(all_loci) < args.num_sets:
        print(f"Warning: Only {len(all_loci)} loci available, but {args.num_sets} requested", file=sys.stderr)
        args.num_sets = len(all_loci)
    
    # Process each set with locus replacement logic
    successful_sets = 0
    used_loci = set()
    
    with tempfile.TemporaryDirectory() as temp_dir:
        for set_idx in range(1, args.num_sets + 1):
            max_locus_replacements = min(100, len(all_loci) - len(used_loci))  # Prevent infinite loops
            
            for replacement_attempt in range(max_locus_replacements):
                # Select a new locus that hasn't been used
                available_loci = [locus for locus in all_loci if locus not in used_loci]
                
                if not available_loci:
                    print(f"Error: No more available loci for set {set_idx}", file=sys.stderr)
                    break
                
                current_locus = rng.choice(available_loci)
                used_loci.add(current_locus)
                
                print(f"Processing set {set_idx}/{args.num_sets} (attempt {replacement_attempt + 1}): {current_locus[0]}:{current_locus[1]}-{current_locus[2]}", file=sys.stderr)
                
                # Try to process this locus
                success, selected_variants = process_locus_with_retries(
                    current_locus, bim_data, args.population_dirs, 
                    args.plink_naming_convention, args.variants_per_set, 
                    args.max_ld, temp_dir, rng, args.max_locus_attempts
                )
                
                if success and selected_variants:
                    write_variant_file(args.output_dir, set_idx, selected_variants)
                    successful_sets += 1
                    print(f"Successfully selected {len(selected_variants)} variants: {selected_variants}", file=sys.stderr)
                    break
                else:
                    print(f"Failed to select variants from locus {current_locus[0]}:{current_locus[1]}-{current_locus[2]}, trying new locus", file=sys.stderr)
            else:
                print(f"Warning: Failed to find suitable locus for set {set_idx} after {max_locus_replacements} attempts", file=sys.stderr)
            
            if set_idx % 10 == 0:
                print(f"Completed {set_idx}/{args.num_sets} sets ({successful_sets} successful)", file=sys.stderr)
    
    print(f"Completed processing. Successfully created {successful_sets}/{args.num_sets} variant sets", file=sys.stderr)

if __name__ == "__main__":
    main()