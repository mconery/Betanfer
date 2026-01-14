#!/bin/bash

################################################################################
# flux_node_worker.sh
#
# Node-local Flux execution script for CARMA fine-mapping
# 
# This script runs on EACH compute node and:
# 1. Starts a local Flux broker (independent scheduler)
# 2. Submits pre-assigned CARMA commands to local queue
# 3. Waits for all jobs to complete
# 4. Cleans up
#
# Usage:
#   bash flux_node_worker.sh \
#     --batch-file /path/to/batch_node_0.txt \
#     --carma-script /path/to/finemap_CARMA.R \
#     --env /path/to/carma_conda_env \
#     --node-id 0
#
# Can be tested standalone on single node for validation
################################################################################

set -e

################################################################################
# FUNCTIONS
################################################################################

print_header() {
  echo ""
  echo "=================================="
  echo "$1"
  echo "=================================="
  echo ""
}

print_status() {
  echo "[$(date +%H:%M:%S)] $1"
}

cleanup() {
  print_status "Cleaning up..."
  
  # Kill any remaining Flux processes
  killall -9 flux 2>/dev/null || true
  
  # Remove broker state directory
  if [[ -n "${BROKER_STATEDIR}" && -d "${BROKER_STATEDIR}" ]]; then
    rm -rf "${BROKER_STATEDIR}"
    print_status "Removed broker state directory"
  fi
}

trap cleanup EXIT

################################################################################
# PARSE ARGUMENTS
################################################################################

BATCH_FILE=""
CARMA_SCRIPT=""
MAP_ENVIRO=""
NODE_ID=${SLURM_NODEID}
DRY_RUN=false
TEST_MODE=false

while [[ $# -gt 0 ]]; do
  case $1 in
    --batch-file-prefix|-p)
      BATCH_FILE_PREFIX="$2"
      shift 2
      ;;
    --batch-file-suffix|-s)
      BATCH_FILE_SUFFIX="$2"
      shift 2
      ;;
    --carma-script|-c)
      CARMA_SCRIPT="$2"
      shift 2
      ;;
    --env|-e)
      MAP_ENVIRO="$2"
      shift 2
      ;;
    --dry-run|-d)
      DRY_RUN=true
      shift
      ;;
    --test)
      TEST_MODE=true
      shift
      ;;
    --help|-h)
      cat << EOF
Usage: bash flux_node_worker.sh [OPTIONS]

Options:
  --batch-file, -b FILE      Path to command batch file for this node (REQUIRED)
  --carma-script, -c FILE    Path to finemap_CARMA.R script (REQUIRED)
  --env, -e PATH             Path to conda environment (REQUIRED)
  --node-id, -n ID           Node rank/ID (default: SLURM_NODEID)
  --dry-run, -d              Print first 5 commands without executing
  --test                     Run in test mode (single command)
  --help, -h                 Show this help message

Examples:
  # Normal execution
  bash flux_node_worker.sh \\
    --batch-file ./batch_node_0.txt \\
    --carma-script ./finemap_CARMA.R \\
    --env /path/to/carma

  # Dry run (see what would execute)
  bash flux_node_worker.sh \\
    --batch-file ./batch_node_0.txt \\
    --carma-script ./finemap_CARMA.R \\
    --env /path/to/carma \\
    --dry-run

  # Test mode (execute 1 command)
  bash flux_node_worker.sh \\
    --batch-file ./batch_node_0.txt \\
    --carma-script ./finemap_CARMA.R \\
    --env /path/to/carma \\
    --test
EOF
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

BATCH_FILE="${BATCH_FILE_PREFIX}${NODE_ID}${BATCH_FILE_SUFFIX}"

################################################################################
# VALIDATE ARGUMENTS
################################################################################

if [[ -z "${BATCH_FILE}" ]]; then
    echo "ERROR: --batch-file is required"
    exit 1
fi

if [[ -z "${CARMA_SCRIPT}" ]]; then
    echo "ERROR: --carma-script is required"
    exit 1
fi

if [[ -z "${MAP_ENVIRO}" ]]; then
    echo "ERROR: --env is required"
    exit 1
fi

if [[ ! -f "${BATCH_FILE}" ]]; then
    echo "ERROR: Batch file not found: ${BATCH_FILE}"
    exit 1
fi

if [[ ! -f "${CARMA_SCRIPT}" ]]; then
    echo "ERROR: CARMA script not found: ${CARMA_SCRIPT}"
    exit 1
fi

if [[ ! -d "${MAP_ENVIRO}" ]]; then
    echo "ERROR: Conda environment not found: ${MAP_ENVIRO}"
    exit 1
fi

################################################################################
# CONFIGURATION
################################################################################

HOSTNAME=$(hostname)
NUM_CORES=$(nproc)
CORES_PER_JOB=2
MAX_CONCURRENT=$((NUM_CORES / CORES_PER_JOB))

print_header "FLUX NODE WORKER INITIALIZATION"

print_status "Hostname: ${HOSTNAME}"
print_status "Node ID: ${NODE_ID}"
print_status "Available cores: ${NUM_CORES}"
print_status "Cores per CARMA job: ${CORES_PER_JOB}"
print_status "Max concurrent jobs: ${MAX_CONCURRENT}"
print_status "Batch file: ${BATCH_FILE}"
print_status "CARMA script: ${CARMA_SCRIPT}"
print_status "Conda environment: ${MAP_ENVIRO}"

################################################################################
# COUNT COMMANDS
################################################################################

NUM_COMMANDS=$(grep -v "^#" "${BATCH_FILE}" | grep -v "^$" | wc -l)
print_status "Commands in batch: ${NUM_COMMANDS}"

if [[ ${NUM_COMMANDS} -eq 0 ]]; then
    print_status "ERROR: No commands found in batch file!"
    exit 1
fi

if [[ "${DRY_RUN}" == "true" ]]; then
    print_header "DRY RUN - FIRST 5 COMMANDS"
    echo "Commands that would be executed:"
    echo ""
    grep -v "^#" "${BATCH_FILE}" | grep -v "^$" | head -5 | nl
    echo ""
    print_status "Dry run complete (use --test to execute one command)"
    exit 0
fi

if [[ "${TEST_MODE}" == "true" ]]; then
    print_header "TEST MODE - SINGLE COMMAND"
    NUM_COMMANDS=1
    print_status "Running first command only for testing"
fi

################################################################################
# SETUP CONDA ENVIRONMENT
################################################################################

print_status "Activating conda environment..."
source activate "${MAP_ENVIRO}" 2>&1 | grep -v "^$" | head -3

# Verify R is available
if ! command -v Rscript &> /dev/null; then
    print_status "ERROR: Rscript not found in conda environment"
    exit 1
fi

print_status "R available: $(Rscript --version 2>&1 | head -1)"

################################################################################
# SETUP BROKER STATE DIRECTORY
################################################################################

# Use /tmp for local-only ephemeral state (key to avoiding locks!)
BROKER_STATEDIR="/tmp/flux-broker-${NODE_ID}-$$"
mkdir -p "${BROKER_STATEDIR}"
print_status "Broker state directory: ${BROKER_STATEDIR}"

# Verify it's writable
if ! touch "${BROKER_STATEDIR}/test" 2>/dev/null; then
    print_status "ERROR: Cannot write to ${BROKER_STATEDIR}"
    exit 1
fi

rm -f "${BROKER_STATEDIR}/test"

################################################################################
# SET FLUX ENVIRONMENT
################################################################################

# CRITICAL: Each node's broker is independent with local state
export FLUX_BROKER_OPTS="-Sstatedir=${BROKER_STATEDIR}"

# CRITICAL: No inter-node coordination (no RPC calls between brokers)
export FLUX_PMI=off

# Optional: Increase verbosity for debugging
# export FLUX_LOG_LEVEL=debug

print_status "Flux environment configured:"
print_status "  FLUX_BROKER_OPTS=${FLUX_BROKER_OPTS}"
print_status "  FLUX_PMI=${FLUX_PMI}"

################################################################################
# CREATE FLUX EXECUTION SCRIPT (FIXED: Using temp file instead of heredoc)
################################################################################

# Create a temporary script file instead of using heredoc with stdin
FLUX_SCRIPT=$(mktemp)

cat > "${FLUX_SCRIPT}" << 'FLUX_SCRIPT_EOF'
#!/bin/bash

set -e

# These variables are passed from parent shell via environment
BATCH_FILE="${BATCH_FILE_VAR}"
NODE_ID=${NODE_ID_VAR}
NUM_COMMANDS=${NUM_COMMANDS_VAR}
CORES_PER_JOB=${CORES_PER_JOB_VAR}
TEST_MODE=${TEST_MODE_VAR}

# Get current timestamp for logging
START_TIME=$(date +%s)

echo "[$(date +%H:%M:%S)] Flux broker started on $(hostname)"

# Show available resources
echo ""
echo "Available resources:"
flux resource list
echo ""

echo "Submitting ${NUM_COMMANDS} CARMA commands..."
echo ""

# Submit all commands
JOB_COUNT=0
COMMAND_NUM=0

while IFS= read -r carma_command; do

    # Skip empty lines and comments
    [[ -z "${carma_command}" ]] && continue
    [[ "${carma_command}" =~ ^# ]] && continue

    COMMAND_NUM=$((COMMAND_NUM + 1))

    # Stop after specified number of commands (for testing)
    if (( COMMAND_NUM > NUM_COMMANDS )); then
        break
    fi

    # Submit to local Flux queue
    # Each job uses CORES_PER_JOB cores, allowing multiple concurrent jobs
    JOB_ID=$(flux submit -N 1 -n 1 -c ${CORES_PER_JOB} \
        bash -c "${carma_command}" 2>/dev/null)

    JOB_COUNT=$((JOB_COUNT + 1))

    # Print progress every 25 jobs
    if (( JOB_COUNT % 25 == 0 )); then
        echo "[$(date +%H:%M:%S)] Submitted ${JOB_COUNT}/${NUM_COMMANDS} commands"
    fi

done < "${BATCH_FILE}"

echo ""
echo "[$(date +%H:%M:%S)] All ${JOB_COUNT} commands submitted"
echo ""

# Show queue status
echo "Queue status after submission:"
flux queue stats
echo ""

echo "[$(date +%H:%M:%S)] Waiting for all jobs to complete..."
echo ""

# Wait for queue to drain (all jobs complete)
# This blocks until every submitted job finishes
flux queue drain

echo ""
echo "[$(date +%H:%M:%S)] All jobs completed"

# Show final stats
FINAL_COUNT=$(flux job list --state=run --no-header 2>/dev/null | wc -l || echo "0")
echo "Running jobs after drain: ${FINAL_COUNT}"

# Calculate execution time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "[$(date +%H:%M:%S)] Flux execution took $((ELAPSED / 60)) minutes $((ELAPSED % 60)) seconds"

FLUX_SCRIPT_EOF

chmod +x "${FLUX_SCRIPT}"

################################################################################
# START FLUX BROKER AND SUBMIT JOBS
################################################################################

print_header "STARTING FLUX BROKER AND SUBMITTING JOBS"

# Export variables needed by the flux script
export BATCH_FILE_VAR="${BATCH_FILE}"
export NODE_ID_VAR="${NODE_ID}"
export NUM_COMMANDS_VAR="${NUM_COMMANDS}"
export CORES_PER_JOB_VAR="${CORES_PER_JOB}"
export TEST_MODE_VAR="${TEST_MODE}"

# FIXED: Use 'flux start <script>' instead of 'flux start bash << heredoc'
# This avoids the "stdin is not a tty" error
flux start -opmi=off "${FLUX_SCRIPT}"

FLUX_EXIT_CODE=$?

# Clean up temp script
rm -f "${FLUX_SCRIPT}"

if [[ ${FLUX_EXIT_CODE} -ne 0 ]]; then
    print_status "ERROR: Flux execution failed with code ${FLUX_EXIT_CODE}"
    exit ${FLUX_EXIT_CODE}
fi

################################################################################
# COMPLETION
################################################################################

print_header "NODE WORKER COMPLETE"

print_status "Node ${NODE_ID} finished successfully"
print_status "Processed ${NUM_COMMANDS} CARMA commands"
print_status "Results should be in output directory"

echo ""
echo "==================================="
echo "NODE WORKER FINISHED"
echo "==================================="
echo ""

exit 0
