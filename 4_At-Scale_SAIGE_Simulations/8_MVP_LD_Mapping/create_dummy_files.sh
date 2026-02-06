#!/bin/bash

################################################################################
# Generate Dummy Results Commands for Failed Fine-Mapping Runs
#
# This script identifies output files that should exist but don't,
# then generates create_dummy_results.py commands to create dummy results.
#
# It mirrors the directory structure and configurations from mapping script.
################################################################################

# Set script locations
command_generator="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/8_MVP_LD_Mapping/create_dummy_commands.py"

# Set directories (mirrored from finemap_MultiSuSiE_mvp.sh)
sumstats_dir="/sw/frontier/mde/bif154_mde/betanfer_project/step2_results"
matrix_dir="/sw/frontier/mde/bif154_mde/betanfer_project/LD_matrices"
base_finemap_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/fine-mapping_results"

# Command files directory for repairs
repair_files_dir="${base_finemap_dir}/command_files"
mkdir -p "$repair_files_dir"

# Populations
populations="AFR,EUR"

# Matrix types
matrix_types="mvp_external"

# Fine-mapping method name
methods=("multisusie" "susiex" "mesusie" "susie" "carma" "rsparsepro")

# Define heritabilities (matching finemap_MultiSuSiE_mvp.sh)
heritabilities=("5e-5" "3e-4")

echo "================================================================================"
echo "GENERATING DUMMY RESULTS COMMANDS FOR FAILED MULTISUSIE RUNS"
echo "================================================================================"
echo "Start time: $(date)"
echo "GWAS directory: $sumstats_dir"
echo "LD matrix directory: $matrix_dir"
echo "Base results directory: $base_finemap_dir"
echo "Output directory: $repair_files_dir"
echo "Populations: $populations"
echo "Matrix types: $matrix_types"
echo "Method: $method_name"
echo "Heritabilities: ${heritabilities[@]}"
echo "================================================================================"
echo ""

# Counters
total_command_files=0
total_commands=0

# Loop over heritabilities
for heritability in "${heritabilities[@]}"; do
    echo "========================================="
    echo "Processing heritability: $heritability"
    echo "========================================="
    echo ""

    # Determine output directory based on heritability
    if [[ $heritability == "5e-5" ]]; then
        h_dir="fixed_homo_5e-5"
    else
        h_dir="fixed_homo_3e-4"
    fi
    
    #Loop over the methods
    for method in "${methods[@]}"; do 

        output_dir="${base_finemap_dir}/${h_dir}/${method}_results"

        # Define repair commands file
        repair_commands_file="${repair_files_dir}/make_${h_dir}_${method}_dummy_results.sh"

        echo "Heritability directory: $h_dir"
        echo "Results directory: $output_dir"
        echo "Repair commands file: $repair_commands_file"
        echo ""

        # Check if results directory exists
        if [[ ! -d "$output_dir" ]]; then
            echo "⚠️  Results directory does not exist: $output_dir"
            echo "   Skipping heritability: $heritability"
            echo ""
            continue
        fi

        # Generate dummy results commands
        echo "Generating dummy results commands..."
    
        python "$command_generator" \
            --gwas_dir "$sumstats_dir" \
            --ld_dir "$matrix_dir" \
            --output_dir "$output_dir" \
            --commands_file "$repair_commands_file" \
            --populations "$populations" \
            --matrix_types "$matrix_types" \
            --method "$method_name"

        if [[ $? -ne 0 ]]; then
            echo "❌ ERROR: Failed to generate commands for heritability: $heritability"
            exit 1
        fi

        if [[ ! -f "$repair_commands_file" ]]; then
            echo "⚠️  Commands file not created or empty: $repair_commands_file"
            echo "   No failed runs found for heritability: $heritability"
            echo ""
            continue
        fi

        # Count commands generated
        num_commands=$(grep -c "^python create_dummy_results.py" "$repair_commands_file")
        echo "✓ Generated $num_commands dummy result commands"
    
        ((total_command_files++))
        ((total_commands += num_commands))

        # Make repair script executable
        chmod +x "$repair_commands_file"
        echo "✓ Repair script created: $repair_commands_file"
        echo ""
    done
done

echo "================================================================================"
echo "COMMAND FILE GENERATION COMPLETE"
echo "================================================================================"
echo "End time: $(date)"
echo "Total repair files created: $total_command_files"
echo "Total commands generated: $total_commands"
echo "Repair files directory: $repair_files_dir"
echo "================================================================================"
echo ""

if [[ $total_command_files -eq 0 ]]; then
    echo "⚠️  No failed runs detected across all heritabilities and methods!"
    echo "   All results appear to be present."
    exit 0
fi

echo "NEXT STEPS:"
echo "==========="
echo ""
for heritability in "${heritabilities[@]}"; do
    if [[ $heritability == "5e-5" ]]; then
        h_dir="fixed_homo_5e-5"
    else
        h_dir="fixed_homo_3e-4"
    fi
    
    for method in "${methods[@]}"; do

        repair_commands_file="${repair_files_dir}/make_${h_dir}_${method}_dummy_results.sh"
    
        if [[ -f "$repair_commands_file" ]]; then
            num_cmds=$(grep -c "^python create_dummy_results.py" "$repair_commands_file")
            echo "For heritability $heritability ($num_cmds commands):"
            echo "  bash $repair_commands_file"
            echo ""
        fi
    done
done
