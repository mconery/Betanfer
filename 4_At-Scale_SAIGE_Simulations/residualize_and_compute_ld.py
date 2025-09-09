#!/usr/bin/env python3

"""
residualize_and_compute_ld_gpu.py

GPU-accelerated version for OLCF Frontier (AMD MI250X) using CuPy with ROCm.

This script uses a single GPU to accelerate the genotype residualization and 
correlation matrix computation using CuPy, which provides NumPy-like GPU operations
for AMD ROCm on Frontier.

Installation on Frontier:
-------------------------
# Load modules
module load PrgEnv-gnu cray-python rocm/6.2.4

# Create conda environment
conda create -n cupy-frontier python=3.10 numpy=1.26.4 scipy -c conda-forge
source activate cupy-frontier

# Install CuPy for ROCm
export CUPY_INSTALL_USE_HIP=1
export ROCM_HOME=${ROCM_PATH}
export HCC_AMDGPU_TARGET=gfx90a
CC=gcc CXX=g++ pip install --no-cache-dir --no-binary=cupy cupy==12.3.0

Usage:
------
# Set GPU device (0-7 on Frontier nodes)
export ROCR_VISIBLE_DEVICES=0

python3 residualize_and_compute_ld_gpu.py [arguments...]
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

# Import CuPy for GPU acceleration
try:
    import cupy as cp
    print(f"CuPy version: {cp.__version__}")
    print(f"ROCm available: {cp.cuda.runtime.getDeviceCount() > 0}")
    
    # Set memory pool to avoid memory fragmentation
    mempool = cp.get_default_memory_pool()
    pinned_mempool = cp.get_default_pinned_memory_pool()
    
    # Get GPU device info
    if cp.cuda.runtime.getDeviceCount() > 0:
        device_id = int(os.environ.get('ROCR_VISIBLE_DEVICES', '0'))
        cp.cuda.Device(0).use()  # Use first visible device
        props = cp.cuda.runtime.getDeviceProperties(0)
        print(f"Using GPU: {props['name'].decode()}")
        print(f"GPU memory: {props['totalGlobalMem'] / 1e9:.1f} GB")
    
    GPU_AVAILABLE = True
except ImportError:
    print("Warning: CuPy not available. Install with ROCm support for GPU acceleration.")
    print("Falling back to CPU-only mode.")
    GPU_AVAILABLE = False
    cp = np  # Fallback to NumPy

def parse_args():
    parser = argparse.ArgumentParser(
        description="GPU-accelerated genotype residualization and LD computation using CuPy."
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
    parser.add_argument("--cleanup_temp", action="store_true",
                       help="Remove temporary genotype files after processing.")
    parser.add_argument("--gpu_batch_size", type=int, default=1000,
                       help="Number of SNPs to process per GPU batch (default: 1000).")
    
    return parser.parse_args()

def load_genotype_data(geno_file_prefix):
    """Load PLINK binary format genotype data using pysnptools."""
    bed_file = f"{geno_file_prefix}.bed"
    if not Path(bed_file).exists():
        raise FileNotFoundError(f'Genotype file {bed_file} not found')
    
    print(f'Loading genotypes from {geno_file_prefix}')
    
    # Load using pysnptools
    bed_reader = Bed(geno_file_prefix, count_A1=False)
    snp_data = bed_reader.read()
    
    # Debug info
    print(f"Debug - snp_data.iid shape: {snp_data.iid.shape}")
    print(f"Debug - snp_data.sid shape: {snp_data.sid.shape}")
    print(f"Debug - snp_data.pos shape: {snp_data.pos.shape}")
    
    # Extract sample IDs
    if snp_data.iid.ndim == 2:
        sample_ids = pd.DataFrame({
            'FID': snp_data.iid[:, 0],
            'IID': snp_data.iid[:, 1]
        })
    else:
        sample_ids = pd.DataFrame({
            'FID': snp_data.iid,
            'IID': snp_data.iid
        })
    
    genotypes = snp_data.val.astype(np.float32)  # Use float32 for GPU efficiency
    n_samples, n_snps = genotypes.shape
    print(f'Loaded {n_samples} samples, {n_snps} SNPs')
    
    # Check for missing values
    n_missing = np.isnan(genotypes).sum()
    if n_missing > 0:
        print(f'Found {n_missing} missing genotype calls')
    
    return sample_ids, genotypes

def load_pc_data(pc_file, num_pcs):
    """Load PC data from PLINK eigenvec format."""
    print(f'Loading PCs from {pc_file}')
    pc_df = pd.read_csv(pc_file, sep=r'\s+', header=0,
                       names=['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)])
    
    pc_cols = ['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)]
    pc_df = pc_df[pc_cols]
    
    return pc_df

def residualize_genotypes_gpu(sample_ids, genotypes, pc_df, num_pcs, batch_size=1000):
    """GPU-accelerated genotype residualization using CuPy."""
    print('Performing GPU-accelerated residualization...')
    
    # Merge genotype and PC data
    merged_df = sample_ids.merge(pc_df, on=['FID', 'IID'], how='left')
    if merged_df.isnull().any().any():
        raise ValueError('Missing PC data for some samples')
    
    # Extract PC matrix and convert to float32 for GPU efficiency
    pc_matrix = merged_df.iloc[:, 2:2+num_pcs].values.astype(np.float32)
    print(f'Using {pc_matrix.shape[1]} PCs for residualization')
    
    n_samples, n_snps = genotypes.shape
    
    if GPU_AVAILABLE:
        # Transfer PC matrix to GPU
        pc_matrix_gpu = cp.asarray(pc_matrix)
        print(f"Transferred PC matrix to GPU: {pc_matrix_gpu.shape}")
        
        # Create design matrix with intercept on GPU
        X_gpu = cp.column_stack([cp.ones((n_samples, 1), dtype=cp.float32), pc_matrix_gpu])
        
        # Precompute (X'X)^-1 X' for efficiency
        try:
            XtX_inv = cp.linalg.inv(X_gpu.T @ X_gpu)
            XtX_inv_Xt = XtX_inv @ X_gpu.T
            print("Precomputed (X'X)^-1 X' matrix on GPU")
        except cp.linalg.LinAlgError:
            print("Warning: Singular matrix, using pseudoinverse")
            XtX_inv_Xt = cp.linalg.pinv(X_gpu).T
        
        # Process genotypes in batches to manage GPU memory
        residual_genotypes = np.zeros_like(genotypes)
        n_processed = 0
        
        for batch_start in range(0, n_snps, batch_size):
            batch_end = min(batch_start + batch_size, n_snps)
            batch_size_actual = batch_end - batch_start
            
            if batch_start % (batch_size * 10) == 0:
                print(f'  Processing SNPs {batch_start}-{batch_end} / {n_snps}')
            
            # Transfer batch to GPU
            genotype_batch = cp.asarray(genotypes[:, batch_start:batch_end])
            residual_batch = cp.zeros_like(genotype_batch)
            
            for i in range(batch_size_actual):
                y_gpu = genotype_batch[:, i]
                
                # Check for valid data
                valid_mask = ~cp.isnan(y_gpu)
                if cp.sum(valid_mask) < 10:
                    residual_batch[:, i] = y_gpu
                    continue
                
                # Check for invariant SNPs
                if cp.var(y_gpu[valid_mask]) < 1e-8:
                    residual_batch[:, i] = y_gpu
                    continue
                
                # Compute residuals: y - X * beta
                try:
                    beta = XtX_inv_Xt @ y_gpu
                    y_pred = X_gpu @ beta
                    residuals = y_gpu - y_pred
                    
                    # Preserve missing values
                    residuals = cp.where(valid_mask, residuals, y_gpu)
                    residual_batch[:, i] = residuals
                    n_processed += 1
                    
                except cp.linalg.LinAlgError:
                    residual_batch[:, i] = y_gpu
            
            # Transfer results back to CPU
            residual_genotypes[:, batch_start:batch_end] = cp.asnumpy(residual_batch)
            
            # Free GPU memory for this batch
            del genotype_batch, residual_batch
            if GPU_AVAILABLE:
                mempool.free_all_blocks()
    
    else:
        # CPU fallback
        print("Using CPU fallback for residualization")
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

def compute_correlation_matrix_gpu(residual_genotypes, batch_size=1000):
    """GPU-accelerated correlation matrix computation."""
    print('Computing correlation matrix on GPU...')
    
    n_samples, n_snps = residual_genotypes.shape
    
    if GPU_AVAILABLE and n_snps <= 5000:  # For manageable matrix sizes
        print(f'Using GPU for full {n_snps}x{n_snps} correlation matrix')
        
        # Transfer data to GPU
        genotypes_gpu = cp.asarray(residual_genotypes.T, dtype=cp.float32)  # SNPs x samples
        
        # Compute correlation matrix on GPU
        # Use CuPy's efficient correlation implementation
        corr_matrix_gpu = cp.full((n_snps, n_snps), cp.nan, dtype=cp.float32)
        
        for i in range(n_snps):
            if i % 100 == 0:
                print(f'  Computing correlations for SNP {i+1}/{n_snps}')
            
            for j in range(i, n_snps):
                x, y = genotypes_gpu[i], genotypes_gpu[j]
                
                # Find valid observations
                valid = ~(cp.isnan(x) | cp.isnan(y))
                valid_count = cp.sum(valid)
                
                if valid_count > 10:
                    if i == j:
                        corr_matrix_gpu[i, j] = 1.0
                    else:
                        x_valid, y_valid = x[valid], y[valid]
                        
                        # Compute correlation manually for better control
                        x_mean, y_mean = cp.mean(x_valid), cp.mean(y_valid)
                        x_centered = x_valid - x_mean
                        y_centered = y_valid - y_mean
                        
                        numerator = cp.sum(x_centered * y_centered)
                        x_ss = cp.sum(x_centered * x_centered)
                        y_ss = cp.sum(y_centered * y_centered)
                        
                        if x_ss > 0 and y_ss > 0:
                            corr = numerator / cp.sqrt(x_ss * y_ss)
                            corr = cp.where(cp.isnan(corr), 0.0, corr)
                            corr_matrix_gpu[i, j] = corr_matrix_gpu[j, i] = corr
                        else:
                            corr_matrix_gpu[i, j] = corr_matrix_gpu[j, i] = 0.0
                else:
                    corr_matrix_gpu[i, j] = corr_matrix_gpu[j, i] = 0.0
        
        # Transfer result back to CPU
        corr_matrix = cp.asnumpy(corr_matrix_gpu)
        
        # Clean up GPU memory
        del genotypes_gpu, corr_matrix_gpu
        if GPU_AVAILABLE:
            mempool.free_all_blocks()
        
    else:
        # CPU fallback for large matrices
        print(f'Using CPU for large correlation matrix ({n_snps} SNPs)')
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

def write_output_files(corr_matrix, out_prefix, args):
    """Write LD matrix and log files."""
    ld_file = f"{out_prefix}.ld"
    print(f'Saving LD matrix to {ld_file}')
    np.savetxt(ld_file, corr_matrix, fmt='%.6f', delimiter='\t')
    
    # Create log file
    log_file = f"{out_prefix}.log"
    n_samples, n_snps = corr_matrix.shape[0], corr_matrix.shape[1]
    
    with open(log_file, 'w') as f:
        f.write(f'GPU-accelerated PC-residualized LD matrix for locus {args.locus_name}\n')
        f.write(f'Population: {args.population}\n')
        f.write(f'Chromosome: {args.chrom}\n')
        f.write(f'Region: {args.start}-{args.end}\n')
        f.write(f'Number of samples: {n_samples}\n')
        f.write(f'Number of SNPs: {n_snps}\n')
        f.write(f'Number of PCs regressed out: {args.num_pcs}\n')
        f.write(f'PC file: {args.pc_file}\n')
        f.write(f'Genotype file: {args.geno_file}\n')
        f.write(f'Method: GPU-accelerated with CuPy on ROCm\n')
        f.write(f'GPU available: {GPU_AVAILABLE}\n')
        f.write(f'Batch size: {args.gpu_batch_size}\n')
        
        # Matrix statistics
        if not np.isnan(corr_matrix).all():
            f.write(f'LD matrix mean: {np.nanmean(corr_matrix):.6f}\n')
            f.write(f'LD matrix std: {np.nanstd(corr_matrix):.6f}\n')
            off_diag_mask = ~np.eye(n_snps, dtype=bool)
            f.write(f'Max LD (off-diagonal): {np.nanmax(corr_matrix[off_diag_mask]):.6f}\n')
    
    print(f'Successfully created GPU-accelerated LD matrix: {ld_file}')

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
        print("GPU-ACCELERATED LD MATRIX COMPUTATION")
        print("=" * 60)
        
        if GPU_AVAILABLE:
            print("GPU acceleration enabled with CuPy")
        else:
            print("Running in CPU mode")
        
        # Load data
        sample_ids, genotypes = load_genotype_data(args.geno_file)
        pc_df = load_pc_data(args.pc_file, args.num_pcs)
        
        # GPU-accelerated residualization
        residual_genotypes = residualize_genotypes_gpu(
            sample_ids, genotypes, pc_df, args.num_pcs, args.gpu_batch_size
        )
        
        # GPU-accelerated correlation computation
        corr_matrix = compute_correlation_matrix_gpu(residual_genotypes, args.gpu_batch_size)
        
        # Write output files
        write_output_files(corr_matrix, args.out_prefix, args)
        
        # Cleanup if requested
        if args.cleanup_temp:
            cleanup_temp_files(args.geno_file)
        
        print("GPU-accelerated analysis completed successfully!")
        
        # Print memory usage if GPU was used
        if GPU_AVAILABLE:
            print(f"Peak GPU memory usage: {mempool.total_bytes() / 1e9:.2f} GB")
        
    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()