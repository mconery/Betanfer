'''
Onekg_Simulation_Calculate_Error.py

This script calculates the error metrics for the 1000G Simulations. 

***This needs to be run in the finemap_old conda environment ***
'''

#Set key file locations
inp_file="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/master.results.csv.gz"
out_file="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/fixed_bucket_metrics.csv"

#Import needed packages
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import argparse

#Metric calculation function
def calculate_metrics(true_values, pred_values):
    mask = ~np.isnan(pred_values)
    filtered_true = true_values[mask]
    filtered_pred = pred_values[mask]
    
    if len(filtered_true) == 0:
        return np.nan, np.nan, np.nan
    
    mae = mean_absolute_error(filtered_true, filtered_pred)
    rmse = np.sqrt(mean_squared_error(filtered_true, filtered_pred))
    
    try:
        r2 = r2_score(filtered_true, filtered_pred)
        r2 = max(0, r2) if not np.isnan(r2) else np.nan
    except:
        r2 = np.nan
    
    return mae, rmse, r2


def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description='Calculate deviation metrics between genotype and beta correlations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--inp_file', 
        default='/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/master.results.csv.gz',
        help='Input gzipped CSV file containing simulation results'
    )
    parser.add_argument(
        '--out_file',
        default='/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/genotype_beta_correlation_deviation_metrics.csv',
        help='Output CSV file for calculated metrics'
    )
    args = parser.parse_args()
    
    # Define fixed bucket ranges
    GENO_BINS = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    MAF_BINS = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    
    dtypes = {
        'var1': 'category', 'var2': 'category',
        'maf1': 'float32', 'maf2': 'float32',
        'geno_corr': 'float32', 'beta_corr': 'float32',
        'beta_traits': 'uint16'
    }
    
    df = pd.read_csv(args.inp_file, 
                     compression="gzip", dtype=dtypes)
    
    # Preprocessing
    df['min_maf'] = df[['maf1', 'maf2']].min(axis=1)
    df['abs_geno_corr'] = df['geno_corr'].abs()
    df['effect_size'] = df['effect_size'].fillna(-1)
    df['case_prob'] = df['case_prob'].fillna(-1)
    
    sim_params = ['chromosome', 'population', 'n_traits', 'n_pairs', 'trait_type',
                  'case_prob', 'n_causal', 'effect_size']
    
    results = []
    for scenario_params, scenario_df in df.groupby(sim_params):
        scenario_dict = dict(zip(sim_params, scenario_params))
        
        # Overall metrics
        valid_data = scenario_df[['geno_corr', 'beta_corr']].dropna()
        overall_metrics = calculate_metrics(valid_data['geno_corr'], valid_data['beta_corr']) if not valid_data.empty else (np.nan,)*3
        results.append({
            **scenario_dict,
            'bucket_type': 'overall',
            'bucket_range': '0.0-1.0',
            'n_pairs': len(valid_data),
            'mae': overall_metrics[0],
            'rmse': overall_metrics[1],
            'r_squared': overall_metrics[2]
        })
        
        # Binary trait n_traits==beta_traits
        if scenario_dict['trait_type'] == 'binary':
            # Filter for pairs where n_traits matches beta_traits
            filtered_df = scenario_df[scenario_df['n_traits'] == scenario_df['beta_traits']]
            filtered_valid = filtered_df[['geno_corr', 'beta_corr']].dropna()
            
            # Calculate metrics
            if not filtered_valid.empty:
                filtered_mae, filtered_rmse, filtered_r2 = calculate_metrics(
                    filtered_valid['geno_corr'], filtered_valid['beta_corr']
                )
            else:
                filtered_mae = filtered_rmse = filtered_r2 = np.nan
            
            # Append new row
            results.append({
                **scenario_dict,
                'bucket_type': 'overall_filtered',
                'bucket_range': '0.0-1.0',
                'n_pairs': len(filtered_valid),
                'mae': filtered_mae,
                'rmse': filtered_rmse,
                'r_squared': filtered_r2
            })
        
        # Genotype correlation buckets
        scenario_df['geno_bucket'] = pd.cut(scenario_df['abs_geno_corr'], bins=GENO_BINS, right=False, labels=False)
        for bucket in range(len(GENO_BINS)-1):
            bucket_df = scenario_df[scenario_df['geno_bucket'] == bucket]
            bucket_data = bucket_df[['geno_corr', 'beta_corr']].dropna()
            metrics = calculate_metrics(bucket_data['geno_corr'], bucket_data['beta_corr']) if not bucket_data.empty else (np.nan,)*3
            results.append({
                **scenario_dict,
                'bucket_type': 'geno_corr',
                'bucket_range': f"{GENO_BINS[bucket]:.1f}-{GENO_BINS[bucket+1]:.1f}",
                'n_pairs': len(bucket_data),
                'mae': metrics[0],
                'rmse': metrics[1],
                'r_squared': metrics[2]
            })
        
        # MAF buckets
        scenario_df['maf_bucket'] = pd.cut(scenario_df['min_maf'], bins=MAF_BINS, right=False, labels=False)
        for bucket in range(len(MAF_BINS)-1):
            bucket_df = scenario_df[scenario_df['maf_bucket'] == bucket]
            bucket_data = bucket_df[['geno_corr', 'beta_corr']].dropna()
            metrics = calculate_metrics(bucket_data['geno_corr'], bucket_data['beta_corr']) if not bucket_data.empty else (np.nan,)*3
            results.append({
                **scenario_dict,
                'bucket_type': 'maf',
                'bucket_range': f"{MAF_BINS[bucket]:.1f}-{MAF_BINS[bucket+1]:.1f}",
                'n_pairs': len(bucket_data),
                'mae': metrics[0],
                'rmse': metrics[1],
                'r_squared': metrics[2]
            })
    
    # Clean up and save
    results_df = pd.DataFrame(results)
    results_df['effect_size'] = results_df['effect_size'].replace(-1, np.nan)
    results_df['case_prob'] = results_df['case_prob'].replace(-1, np.nan)
    results_df.to_csv(args.out_file, index=False)
    print("SUCCESS: Analysis complete with fixed bucket ranges!")

if __name__ == "__main__":
    main()