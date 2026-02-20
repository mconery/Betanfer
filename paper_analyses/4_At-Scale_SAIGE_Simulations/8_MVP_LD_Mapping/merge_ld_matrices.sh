#!/bin/bash

#SBATCH -A bif154_mde
#SBATCH -J merge_ld_matrices
#SBATCH -o merge_ld_matrices.log
#SBATCH -N 60
#SBATCH -p batch-spi
#SBATCH -q debug
#SBATCH -t 1:00:00

################################################################################################################
###################################### Load Modules and Environment ############################################
################################################################################################################

module load Core/25.03
module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set script locations
command_generator_script="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/generate_merge_commands.py"
matrix_combiner_script="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/combine_ld_matrices.py"
python_splitter="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Betanfer/paper_analyses/5_LD_Inference/split_commands_parallel.py"

# Set directory locations
input_dirs=(
"/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/LD_matrices/mvp_external"
)

# Output location for commands and flux output
output_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/LD_matrices"
batch_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/LD_matrices/command_files"
commands_file="${batch_dir}/ld_merge_commands.txt"

# Batch configuration
commands_per_batch=50

# Sample sizes
afr_sample_size=123328
eur_sample_size=457720

# Set number of nodes
num_nodes=${SLURM_NNODES}

################################################################################################################
################################## Create output directories and prep commands ##################################
################################################################################################################

echo "Starting LD matrix merge workflow at $(date)"
echo "Input directories: ${input_dirs[@]}"
echo "Output directory: ${output_dir}"
echo "Using ${num_nodes} nodes"

# Ensure weighted subdirectories exist in input directories
for input_dir in "${input_dirs[@]}"; do
  if [ -d "${input_dir}" ]; then
    echo "Checking directory: ${input_dir}"
    mkdir -p "${input_dir}/weighted"
  else
    echo "WARNING: Input directory does not exist: ${input_dir}"
  fi
done

echo ""

################################################################################################################
############################ Create generation job scripts for each directory ##################################
################################################################################################################

echo "Preparing parallel command generation jobs..."

# Create a temporary directory for job scripts
temp_job_dir="${batch_dir}/.temp-gen-scripts"
mkdir -p "${temp_job_dir}"

# Create a manifest file to track all generation jobs
job_manifest="${temp_job_dir}/merge-gen-manifest.txt"
> "${job_manifest}"

job_count=0

# Create a job script for each input directory
for input_dir in "${input_dirs[@]}"; do
  if [ -d "${input_dir}" ]; then
    dir_name=$(basename "${input_dir}")
    job_script="${temp_job_dir}/gen-job-${dir_name}.sh"
    
    cat > "${job_script}" <<'JOBEOF'
#!/bin/bash
set -e

python_cmd="python PLACEHOLDER_command_generator --input-dirs PLACEHOLDER_input_dir --output-commands PLACEHOLDER_output_commands --afr-sample-size PLACEHOLDER_afr_sample_size --eur-sample-size PLACEHOLDER_eur_sample_size --combiner-script PLACEHOLDER_matrix_combiner"

eval "$python_cmd"

if [ $? -ne 0 ]; then
  echo "ERROR: Failed to generate commands for directory PLACEHOLDER_input_dir" >&2
  exit 1
fi

if [ ! -f "PLACEHOLDER_output_commands" ] || [ ! -s "PLACEHOLDER_output_commands" ]; then
  echo "ERROR: Commands file is missing or empty: PLACEHOLDER_output_commands" >&2
  exit 1
fi

num_commands=$(grep -v '^#' "PLACEHOLDER_output_commands" | grep -v '^$' | wc -l)
echo "Generated ${num_commands} commands for $(basename PLACEHOLDER_input_dir)"
JOBEOF

    # Replace placeholders
    sed -i "s|PLACEHOLDER_command_generator|${command_generator_script}|g" "${job_script}"
    sed -i "s|PLACEHOLDER_input_dir|${input_dir}|g" "${job_script}"
    sed -i "s|PLACEHOLDER_output_commands|${batch_dir}/${dir_name}.tmp|g" "${job_script}"
    sed -i "s|PLACEHOLDER_afr_sample_size|${afr_sample_size}|g" "${job_script}"
    sed -i "s|PLACEHOLDER_eur_sample_size|${eur_sample_size}|g" "${job_script}"
    sed -i "s|PLACEHOLDER_matrix_combiner|${matrix_combiner_script}|g" "${job_script}"
    
    chmod +x "${job_script}"
    
    # Add to manifest for tracking
    echo "${job_script}|${batch_dir}/${dir_name}.txt|${input_dir}" >> "${job_manifest}"
    
    ((job_count++))
  fi
done

echo "Prepared ${job_count} command generation jobs"
echo ""

################################################################################################################
############################# Execute command generation jobs in parallel via Flux #############################
################################################################################################################

echo "Starting Flux-based command generation..."
echo ""

# Create output directory for generation logs
gen_log_dir="${batch_dir}/.gen-logs"
mkdir -p "${gen_log_dir}"
rm $gen_log_dir/*.log
mkdir -p "${output_dir}/jobstate/merge_gen"
rm -f "${output_dir}/jobstate/merge_gen"/*

# Start Flux instance and submit all generation jobs
echo "Submitting ${job_count} command generation jobs to Flux..."

srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start --broker-opts="-Sstatedir=${output_dir}/jobstate/merge_gen" bash -c "

echo 'Flux generation instance started at:' \$(date)
echo 'Available resources:'
flux resource list
echo ''

# Submit all generation jobs
job_count=0

while IFS='|' read -r job_script output_commands input_dir; do
  dir_name=\$(basename \${input_dir})
  flux submit -N 1 -n 1 -c 56 --output=${gen_log_dir}/\${dir_name}.log bash \"\${job_script}\"
  ((job_count++))
  
  if [ \$((job_count % 5)) -eq 0 ]; then
    echo \" Submitted \${job_count}/${job_count} generation jobs\"
  fi
done < \"${job_manifest}\"

echo \"All ${job_count} generation jobs submitted to Flux queue at: \$(date)\"
echo ''

# Wait for all generation jobs to complete
echo 'Waiting for all command generation jobs to complete...'
flux queue drain

echo 'All command generation jobs completed at:' \$(date)
echo ''
"

if [ $? -ne 0 ]; then
  echo "ERROR: Flux command generation instance failed"
  exit 1
fi

echo ""

################################################################################################################
################################## Merge Individual Command Files to Master File ##############################
################################################################################################################

echo "Merging individual command files into master file..."
echo "Start time: $(date)"

# Clear the master commands file
> "${commands_file}"

# Loop through all temporary command files and append to master file sequentially
# This avoids race conditions by processing files AFTER all parallel jobs complete
for manifest_line in $(cat "${job_manifest}"); do
  IFS='|' read -r job_script output_commands input_dir <<< "${manifest_line}"
  echo $output_commands
  dir_name=$(basename "${input_dir}")
  
  if [[ -f "${output_commands}" ]] && [[ -s "${output_commands}" ]]; then
    num_lines=$(grep -v '^#' "${output_commands}" | grep -v '^$' | wc -l)
    grep -v '^#' "${output_commands}" | grep -v '^$' >> "${commands_file}"
    echo "✓ Merged ${num_lines} commands from ${dir_name}"
  else
    echo "WARNING: Commands file missing or empty for ${dir_name}: ${output_commands}"
  fi
done

echo ""
echo "Total commands merged: $(wc -l < ${commands_file})"
echo "Merge completed at $(date)"
echo ""

# Verify master commands file
if [ ! -f "${commands_file}" ] || [ ! -s "${commands_file}" ]; then
  echo "ERROR: Master commands file is missing or empty!"
  exit 1
fi

# Count the number of commands
num_commands=$(grep -v '^#' "${commands_file}" | grep -v '^$' | wc -l)
echo "Verified: Master commands file contains ${num_commands} commands"

if [ ${num_commands} -eq 0 ]; then
  echo "WARNING: No commands were generated. Check input directories and file structure."
  exit 1
fi

# Clean up temporary command files
for manifest_line in $(cat "${job_manifest}"); do
  IFS='|' read -r job_script output_commands input_dir <<< "${manifest_line}"
  rm -f "${output_commands}"
done

echo ""

################################################################################################################
################################## Create batched command files ################################################
################################################################################################################

echo "Creating batched command files..."

# Create batch directory
mkdir -p "${batch_dir}"

# Clean up any existing batch files
rm -f "${batch_dir}"/merge_batch_*.sh

# Use Python splitter for fast parallel batch creation
python "${python_splitter}" \
  -i "${commands_file}" \
  -o "${batch_dir}" \
  -b ${commands_per_batch} \
  -p "merge_batch_" \
  -s ".sh"

# Verify batch creation
num_batch_files=$(ls -1 "${batch_dir}"/merge_batch_*.sh 2>/dev/null | wc -l)
echo "Created ${num_batch_files} batch files in ${batch_dir}"

if [ ${num_batch_files} -eq 0 ]; then
  echo "ERROR: No batch files were created"
  exit 1
fi

echo ""

################################################################################################################
################################## Execute merge commands using Flux ##########################################
################################################################################################################

echo "Starting flux-based execution of batched merge commands at $(date)"

# Create flux output directory
mkdir -p "${output_dir}/flux_output"
mkdir -p "${output_dir}/jobstate/merge_exec"
rm -f "${output_dir}/jobstate/merge_exec"/*

# Start Flux and submit batch jobs
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start --broker-opts="-Sstatedir=${output_dir}/jobstate/merge_exec" bash -c "

# Display Flux resource allocation
echo \"Flux instance started at \$(date)\"
echo \"Available resources:\"
flux resource list
echo \"\"

# Submit all batch files to Flux
echo \"Submitting ${num_batch_files} batches to Flux queue...\"
batch_count=0

for batch in ${batch_dir}/merge_batch_*.sh; do
  batch_name=\$(basename \${batch} .sh)
  flux submit -N1 -n1 -c56 \
    --output=${output_dir}/flux_output/\${batch_name}.log \
    bash \${batch}
  
  ((batch_count++))
  
  # Print progress every 10 batches
  if [ \$((batch_count % 10)) -eq 0 ]; then
    echo \" Submitted \${batch_count}/${num_batch_files} batches\"
  fi
done

echo \"All ${num_batch_files} batches submitted to Flux queue at \$(date)\"
echo \"\"

# Wait for all jobs to complete
echo \"Waiting for all batches to complete...\"
flux queue drain
echo \"\"

echo \"All batches completed at \$(date)\"

# Display final resource status
echo \"Final Flux resource status:\"
flux resource list
"

echo ""

################################################################################################################
################################## Post-processing and cleanup ##################################################
################################################################################################################

echo "Post-processing results..."

# Count successful batch completions by checking for output files
successful_batches=$(ls -1 "${output_dir}"/flux_output/batch_*.log 2>/dev/null | wc -l)
echo "Successfully completed batches: ${successful_batches}/${num_batch_files}"

# Check for any failed batches
if [ ${successful_batches} -ne ${num_batch_files} ]; then
  echo "WARNING: Some batches may have failed. Check logs in ${output_dir}/flux_output/"
fi

echo ""

# Clean up temporary job scripts
rm -rf "${temp_job_dir}"

echo "LD matrix merge workflow completed at $(date)"
echo "Batch files retained in: ${batch_dir}"
echo "Flux logs available in: ${output_dir}/flux_output"
