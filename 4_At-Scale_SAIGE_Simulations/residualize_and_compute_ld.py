#!/usr/bin/env python3

"""
residualize_and_compute_ld_gpu.py

PyTorch GPU-accelerated version for OLCF Frontier (AMD MI250X).

This script uses PyTorch with ROCm for GPU acceleration of genotype residualization 
and correlation matrix computation. PyTorch has much better AMD GPU support on 
Frontier than CuPy.

Installation on Frontier:
-------------------------
# Load modules
module load PrgEnv-gnu cray-python rocm/6.2.4

# Use the pre-installed PyTorch module (recommended)
module load pytorch/2.0.1

# OR create your own environment:
# conda create -n pytorch-frontier python=3.10
# conda activate pytorch-frontier
# pip install torch torchvision --index-url https://download.pytorch.org/whl/rocm6.0

Usage:
------
# Set GPU device (0-7 on Frontier nodes)
export HIP_VISIBLE_DEVICES=0

python3 residualize_and_compute_ld_pytorch.py [arguments...]
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import sys
import os

try:
    from pysnptools.snpreader import Bed
except ImportError:
    print("Error: pysnptools package is required. Install with: pip install pysnptools", file=sys.stderr)
    sys.exit(1)

# Import PyTorch for GPU acceleration
try:
    import torch
    print(f"PyTorch version: {torch.__version__}")
    # Check for GPU availability
    if torch.cuda.is_available():
        gpu_count = torch.cuda.device_count()
        current_device = torch.cuda.current_device()
        gpu_name = torch.cuda.get_device_name(current_device)
        print(f"GPU available: {gpu_name}")
        print(f"GPU count: {gpu_count}")
        print(f"Current device: {current_device}")
        # Set device
        device = torch.device('cuda')
        # Test basic operation
        test_tensor = torch.tensor([1.0, 2.0, 3.0]).to(device)
        print(f"GPU test successful: {test_tensor.sum().item()}")
        GPU_AVAILABLE = True
    else:
        print("CUDA not available, using CPU")
        device = torch.device('cpu')
        GPU_AVAILABLE = False
except ImportError:
    print("PyTorch not available, using CPU only")
    device = torch.device('cpu') if 'torch' in locals() else None
    GPU_AVAILABLE = False

def parse_args():
    parser = argparse.ArgumentParser(
        description="Residualize genotypes using PCs and compute LD correlation matrix with correct allele encoding."
    )
    parser.add_argument("--geno_file", required=True,
                       help="Path prefix to PLINK binary files (.bed/.bim/.fam) without extension.")
    parser.add_argument("--pc_file", required=True,
                       help="Path to PLINK eigenvec format PC file.")
    parser.add_argument("--out_prefix", required=True,
                       help="Output file prefix for LD matrix and log.")
    parser.add_argument("--num_pcs", type=int, default=10,
                       help="Number of PCs to regress out (default: 10).")
    parser.add_argument("--locus_name", required=True,
                       help="Name of the locus being analyzed.")
    parser.add_argument("--population", required=True,
                       help="Population name.")
    parser.add_argument("--chrom", required=True,
                       help="Chromosome number.")
    parser.add_argument("--start", type=int, required=True,
                       help="Start position of the region.")
    parser.add_argument("--end", type=int, required=True,
                       help="End position of the region.")
    parser.add_argument("--gpu_batch_size", type=int, default=400000,
                       help="Number of SNPs to process per GPU batch (default: 2000).")
    parser.add_argument("--cleanup_temp", action="store_true",
                       help="Remove temporary genotype files after processing.")
    parser.add_argument("--count_a1", action="store_true",
                       help="Count A1 alleles (like plink2 ALT). Default: count A2 alleles (like plink1.9).")
    return parser.parse_args()

def load_genotype_data(geno_file_prefix, count_a1=False):
    """Load PLINK binary format genotype data using pysnptools with correct allele encoding."""
    bed_file = f"{geno_file_prefix}.bed"
    if not Path(bed_file).exists():
        raise FileNotFoundError(f'Genotype file {bed_file} not found')
    
    print(f'Loading genotypes from {geno_file_prefix}')
    print(f'Allele counting mode: {"A1 (minor/ALT)" if count_a1 else "A2 (major/REF)"}')
    # CRITICAL FIX: Use count_A1=True to match plink2 --r-unphased behavior
    # This ensures consistent allele encoding and prevents sign flips
    bed_reader = Bed(geno_file_prefix, count_A1=count_a1)
    snp_data = bed_reader.read()
    # Debug: print the structure of the loaded data
    print(f"Debug - snp_data.iid shape: {snp_data.iid.shape}")
    print(f"Debug - snp_data.sid shape: {snp_data.sid.shape}")
    print(f"Debug - snp_data.pos shape: {snp_data.pos.shape}")
    # Extract sample IDs more robustly
    if snp_data.iid.ndim == 2:
        # Standard format with FID and IID columns
        sample_ids = pd.DataFrame({
            'FID': snp_data.iid[:, 0],
            'IID': snp_data.iid[:, 1]
        })
    elif snp_data.iid.ndim == 1:
        # Single column format - use same value for both FID and IID
        sample_ids = pd.DataFrame({
            'FID': snp_data.iid,
            'IID': snp_data.iid
        })
    else:
        raise ValueError(f"Unexpected iid structure: shape {snp_data.iid.shape}")
    # Extract SNP information more robustly
    n_snps = len(snp_data.sid)
    # Handle sid (SNP IDs) - can be 1D or 2D
    if snp_data.sid.ndim == 1:
        snp_ids = snp_data.sid
        chrom_info = [''] * n_snps  # No chromosome info available
    elif snp_data.sid.ndim == 2:
        if snp_data.sid.shape[1] >= 2:
            chrom_info = snp_data.sid[:, 0]
            snp_ids = snp_data.sid[:, 1]
        else:
            snp_ids = snp_data.sid[:, 0]
            chrom_info = [''] * n_snps
    else:
        raise ValueError(f"Unexpected sid structure: shape {snp_data.sid.shape}")
    
    # Handle pos (positions) - can be 1D or 2D
    if snp_data.pos.ndim == 1:
        positions = snp_data.pos
    elif snp_data.pos.ndim == 2:
        # Usually pos has columns: [chr, pos] or [chr, cM, pos]
        if snp_data.pos.shape[1] >= 2:
            positions = snp_data.pos[:, -1]  # Take the last column as BP position
        else:
            positions = snp_data.pos[:, 0]
    else:
        positions = np.zeros(n_snps)  # Fallback
    
    snp_info = pd.DataFrame({
        'CHR': chrom_info,
        'SNP': snp_ids,
        'BP': positions
    })
    
    genotypes = snp_data.val  # This is samples x SNPs
    n_samples, n_snps = genotypes.shape
    print(f'Loaded {n_samples} samples, {n_snps} SNPs')
    
    # Check for and handle missing values
    n_missing = np.isnan(genotypes).sum()
    if n_missing > 0:
        print(f'Found {n_missing} missing genotype calls')
    
    return sample_ids, genotypes, snp_info

def load_pc_data(pc_file, num_pcs):
    """Load PC data from PLINK eigenvec format."""
    print(f'Loading PCs from {pc_file}')
    pc_df = pd.read_csv(pc_file, sep=r'\s+', header=0,
                       names=['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)])
    
    pc_cols = ['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)]
    pc_df = pc_df[pc_cols]
    
    return pc_df

def residualize_genotypes_pytorch(sample_ids, genotypes, pc_df, num_pcs, batch_size=2000):
    """PyTorch GPU-accelerated genotype residualization."""
    print('Performing PyTorch GPU-accelerated residualization...')
    
    # Merge genotype and PC data
    merged_df = sample_ids.merge(pc_df, on=['FID', 'IID'], how='left')
    if merged_df.isnull().any().any():
        raise ValueError('Missing PC data for some samples')
    
    # Extract PC matrix and convert to float32 for GPU efficiency
    pc_matrix = merged_df.iloc[:, 2:2+num_pcs].values.astype(np.float32)
    print(f'Using {pc_matrix.shape[1]} PCs for residualization')
    
    n_samples, n_snps = genotypes.shape
    
    if GPU_AVAILABLE:
        print(f"Using GPU: {device}")
        
        # Transfer PC matrix to GPU
        pc_tensor = torch.tensor(pc_matrix, dtype=torch.float32, device=device)
        
        # Create design matrix with intercept on GPU
        ones = torch.ones((n_samples, 1), dtype=torch.float32, device=device)
        X = torch.cat([ones, pc_tensor], dim=1)
        
        print(f"Design matrix shape on GPU: {X.shape}")
        
        # Precompute (X'X)^-1 X' for efficiency
        try:
            XtX = torch.matmul(X.t(), X)
            XtX_inv = torch.inverse(XtX)
            XtX_inv_Xt = torch.matmul(XtX_inv, X.t())
            print("Precomputed (X'X)^-1 X' matrix on GPU")
        except RuntimeError as e:
            print(f"Matrix inversion failed, using pseudoinverse: {e}")
            XtX_inv_Xt = torch.pinverse(X).t()
        
        # Process genotypes in batches
        residual_genotypes = np.zeros_like(genotypes)
        n_processed = 0
        
        for batch_start in range(0, n_snps, batch_size):
            batch_end = min(batch_start + batch_size, n_snps)
            batch_size_actual = batch_end - batch_start
            
            if batch_start % (batch_size * 5) == 0:
                print(f'  Processing SNPs {batch_start}-{batch_end} / {n_snps}')
            
            # Transfer batch to GPU
            genotype_batch = torch.tensor(
                genotypes[:, batch_start:batch_end], 
                dtype=torch.float32, 
                device=device
            )
            
            residual_batch = torch.zeros_like(genotype_batch)
            
            for i in range(batch_size_actual):
                y = genotype_batch[:, i]
                
                # Check for valid data
                valid_mask = ~torch.isnan(y)
                valid_count = torch.sum(valid_mask)
                
                if valid_count < 10:
                    residual_batch[:, i] = y
                    continue
                
                # Check for invariant SNPs
                y_valid = y[valid_mask]
                if torch.var(y_valid) < 1e-8:
                    residual_batch[:, i] = y
                    continue
                
                # Compute residuals: y - X * beta
                try:
                    beta = torch.matmul(XtX_inv_Xt, y)
                    y_pred = torch.matmul(X, beta)
                    residuals = y - y_pred
                    
                    # Preserve missing values
                    residuals = torch.where(valid_mask, residuals, y)
                    residual_batch[:, i] = residuals
                    n_processed += 1
                    
                except RuntimeError:
                    residual_batch[:, i] = y
            
            # Transfer results back to CPU
            residual_genotypes[:, batch_start:batch_end] = residual_batch.cpu().numpy()
            
            # Clear GPU memory
            del genotype_batch, residual_batch
            if GPU_AVAILABLE:
                torch.cuda.empty_cache()
    
    else:
        # CPU fallback using NumPy (same as original)
        print("Using CPU fallback for residualization")
        pc_matrix = merged_df.iloc[:, 2:2+num_pcs].values
        X = np.column_stack([np.ones(n_samples), pc_matrix])
        
        try:
            XtX_inv_Xt = np.linalg.inv(X.T @ X) @ X.T
        except np.linalg.LinAlgError:
            XtX_inv_Xt = np.linalg.pinv(X).T
        
        residual_genotypes = np.zeros_like(genotypes)
        n_processed = 0
        
        for i in range(n_snps):
            if i % 1000 == 0:
                print(f'  Processing SNP {i}/{n_snps}')
            
            y = genotypes[:, i]
            valid_mask = ~np.isnan(y)
            
            if valid_mask.sum() < 10 or np.var(y[valid_mask]) < 1e-8:
                residual_genotypes[:, i] = y
                continue
            
            try:
                beta = XtX_inv_Xt @ y
                residuals = y - X @ beta
                residuals[~valid_mask] = np.nan
                residual_genotypes[:, i] = residuals
                n_processed += 1
            except:
                residual_genotypes[:, i] = y
    
    print(f'Residualization complete: {n_processed} SNPs processed')
    return residual_genotypes

def compute_correlation_matrix_pytorch(residual_genotypes, batch_size=2000):
    """PyTorch GPU-accelerated correlation matrix computation."""
    print('Computing correlation matrix with PyTorch...')
    
    n_samples, n_snps = residual_genotypes.shape
    
    if GPU_AVAILABLE and n_snps <= 8000:  # For manageable matrix sizes
        print(f'Using GPU for {n_snps}x{n_snps} correlation matrix')
        
        # Transfer data to GPU
        genotypes_tensor = torch.tensor(
            residual_genotypes.T,  # SNPs x samples
            dtype=torch.float32, 
            device=device
        )
        
        # Initialize correlation matrix on GPU
        corr_matrix = torch.full(
            (n_snps, n_snps), 
            float('nan'), 
            dtype=torch.float32, 
            device=device
        )
        
        # Compute correlations in batches to manage memory
        for i in range(0, n_snps, batch_size):
            end_i = min(i + batch_size, n_snps)
            
            if i % (batch_size * 2) == 0:
                print(f'  Computing correlations for SNPs {i}-{end_i}/{n_snps}')
            
            for j in range(i, n_snps, batch_size):
                end_j = min(j + batch_size, n_snps)
                
                # Get batch of SNPs
                snps_i = genotypes_tensor[i:end_i]  # batch_size_i x samples
                snps_j = genotypes_tensor[j:end_j]  # batch_size_j x samples
                
                # Compute correlations for this block
                for ii in range(snps_i.shape[0]):
                    for jj in range(snps_j.shape[0]):
                        global_i = i + ii
                        global_j = j + jj
                        
                        if global_j < global_i:
                            continue
                        
                        x, y = snps_i[ii], snps_j[jj]
                        
                        # Find valid observations
                        valid = ~(torch.isnan(x) | torch.isnan(y))
                        valid_count = torch.sum(valid)
                        
                        if valid_count > 10:
                            if global_i == global_j:
                                corr_matrix[global_i, global_j] = 1.0
                            else:
                                x_valid, y_valid = x[valid], y[valid]
                                
                                # Compute correlation
                                x_mean, y_mean = torch.mean(x_valid), torch.mean(y_valid)
                                x_centered = x_valid - x_mean
                                y_centered = y_valid - y_mean
                                
                                numerator = torch.sum(x_centered * y_centered)
                                x_ss = torch.sum(x_centered * x_centered)
                                y_ss = torch.sum(y_centered * y_centered)
                                
                                if x_ss > 0 and y_ss > 0:
                                    corr = numerator / torch.sqrt(x_ss * y_ss)
                                    if torch.isnan(corr):
                                        corr = torch.tensor(0.0)
                                    
                                    corr_matrix[global_i, global_j] = corr
                                    corr_matrix[global_j, global_i] = corr
                                else:
                                    corr_matrix[global_i, global_j] = 0.0
                                    corr_matrix[global_j, global_i] = 0.0
                        else:
                            corr_matrix[global_i, global_j] = 0.0
                            corr_matrix[global_j, global_i] = 0.0
        
        # Transfer result back to CPU
        corr_matrix_np = corr_matrix.cpu().numpy()
        
        # Clean up GPU memory
        del genotypes_tensor, corr_matrix
        if GPU_AVAILABLE:
            torch.cuda.empty_cache()
        
        return corr_matrix_np
    
    else:
        # CPU fallback for large matrices
        print(f'Using CPU for correlation matrix ({n_snps} SNPs)')
        corr_matrix = np.full((n_snps, n_snps), np.nan)
        
        for i in range(n_snps):
            if i % 500 == 0:
                print(f'  Computing correlations for SNP {i+1}/{n_snps}')
            
            for j in range(i, n_snps):
                x, y = residual_genotypes[:, i], residual_genotypes[:, j]
                valid = ~(np.isnan(x) | np.isnan(y))
                
                if valid.sum() > 10:
                    if i == j:
                        corr_matrix[i, j] = 1.0
                    else:
                        corr = np.corrcoef(x[valid], y[valid])[0, 1]
                        corr_matrix[i, j] = corr_matrix[j, i] = corr if not np.isnan(corr) else 0.0
                else:
                    corr_matrix[i, j] = corr_matrix[j, i] = 0.0
        
        return corr_matrix

def write_output_files(corr_matrix, out_prefix, args, snp_info=None):
    """Write LD matrix and variant files in plink2 --r-unphased square format."""
    
    # File 1: Correlation matrix in plink2 format (.unphased.vcor1)
    vcor_file = f"{out_prefix}.unphased.vcor1"
    print(f'Saving LD correlation matrix to {vcor_file}')
    np.savetxt(vcor_file, corr_matrix, fmt='%.6f', delimiter='\t')
    
    # File 2: Variant list (.unphased.vcor1.vars)
    vars_file = f"{out_prefix}.unphased.vcor1.vars"
    print(f'Saving variant list to {vars_file}')
    
    if snp_info is not None and 'SNP' in snp_info.columns:
        # Use SNP IDs from snp_info DataFrame
        variant_ids = snp_info['SNP'].values
    else:
        # Fallback: create generic variant names if snp_info not available
        n_snps = corr_matrix.shape[0]
        variant_ids = [f"SNP_{i+1}" for i in range(n_snps)]
        print("Warning: SNP info not available, using generic variant names")
    
    # Write variant list (one per line)
    with open(vars_file, 'w') as f:
        for variant_id in variant_ids:
            f.write(f"{variant_id}\n")
    
    print(f'Successfully created plink2-format LD files:')
    print(f'  - Correlation matrix: {vcor_file}')
    print(f'  - Variant list: {vars_file}')

def cleanup_temp_files(geno_file_prefix):
    """Remove temporary genotype files."""
    temp_files = [
        f"{geno_file_prefix}.bed",
        f"{geno_file_prefix}.bim", 
        f"{geno_file_prefix}.fam",
        f"{geno_file_prefix}.log",
        f"{geno_file_prefix}.nosex"
    ]
    
    for temp_file in temp_files:
        try:
            Path(temp_file).unlink(missing_ok=True)
            print(f'Removed temporary file: {temp_file}')
        except Exception as e:
            print(f'Warning: Could not remove {temp_file}: {e}')

def main():
    args = parse_args()
    
    try:
        print("=" * 60)
        print("PYTORCH GPU-ACCELERATED LD MATRIX COMPUTATION")
        print("=" * 60)
        
        if GPU_AVAILABLE:
            print(f"GPU acceleration enabled with PyTorch on {device}")
        else:
            print("Running in CPU mode with PyTorch")
        
        # Load data
        sample_ids, genotypes, snp_info = load_genotype_data(args.geno_file, args.count_a1)
        pc_df = load_pc_data(args.pc_file, args.num_pcs)
        
        # PyTorch-accelerated residualization
        residual_genotypes = residualize_genotypes_pytorch(
            sample_ids, genotypes, pc_df, args.num_pcs, args.gpu_batch_size
        )
        
        # PyTorch-accelerated correlation computation
        corr_matrix = compute_correlation_matrix_pytorch(residual_genotypes, args.gpu_batch_size)
        
        # Write output files
        write_output_files(corr_matrix, args.out_prefix, args)
        
        # Cleanup if requested
        if args.cleanup_temp:
            cleanup_temp_files(args.geno_file)
        
        print("PyTorch GPU-accelerated analysis completed successfully!")
        
        # Print memory usage if GPU was used
        if GPU_AVAILABLE:
            print(f"GPU memory allocated: {torch.cuda.memory_allocated()/1e9:.2f} GB")
            print(f"GPU memory cached: {torch.cuda.memory_reserved()/1e9:.2f} GB")
        
    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()