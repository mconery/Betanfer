#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_vep_annotation.sh
#
# PBS submission script for VEP annotation of the fine-mapped variant set
# (step 12, gene nomination).
#
# Annotates the reference-validated, normalized GRCh37 VCF written by
# prep_vep_input.py, running offline VEP (cache 110) twice with --check_ref as
# a second reference-consistency check: --most_severe (feeds
# compute_vep_enrichment.py and the Figure 5 tables) and --per_gene (feeds
# nominate_coding_genes.py; the normalized VCF supplies the GRCh37 positions
# used by nominate_genes_comprehensive.py). Each run must annotate every
# variant_id of the input VCF.
#
# Input  (application_results/gene_nomination/):
#   vep_input_all_hg19.normalized.vcf  -- prep_vep_input.py
# Outputs (same directory):
#   vep_all_most_severe.tsv  -- VEP --most_severe output
#   vep_output.tsv           -- VEP --per_gene output
#
# Submit with:   qsub run_vep_annotation.sh
# Monitor with:  qstat -u $USER
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A gpu_hack
#PBS -l walltime=03:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q capacity
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/vep_annotation.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/vep_annotation.err

set -euo pipefail

# ── Environment ────────────────────────────────────────────────────────────────
# ensembl-vep is only available in this conda environment; point VEP_ENV at
# your own VEP install if you have one. Its bin/ goes first on PATH so vep runs
# under the environment's perl rather than the system perl.
VEP_ENV=/home/jianhuawang/miniforge3/envs/vep
export PATH=${VEP_ENV}/bin:${PATH}
VEP=${VEP_ENV}/bin/vep

# ── Key paths ───────────────────────────────────────────────────────────────────
out_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination
norm_vcf=${out_dir}/vep_input_all_hg19.normalized.vcf
# GRCh37 primary-assembly fasta (Ensembl release 110), only available in our tree.
grch37_fasta=/lus/flare/projects/GeomicVar/jianhuawang/betainfer/data/interim/2026-09-23_annotation-repair/reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa
vep_cache=/lus/flare/projects/GeomicVar/jianhuawang/betainfer/data/raw/vep_cache
cache_version=110
forks=16

mkdir -p "${out_dir}"

# Every input variant_id must appear in a VEP output, and nothing else.
check_coverage() {
    if ! cmp -s <(grep -v '^#' "${norm_vcf}" | cut -f3 | sort -u) \
                <(grep -v '^#' "$1" | cut -f1 | sort -u); then
        echo "ERROR: $1 does not annotate exactly the variants of ${norm_vcf}."
        exit 1
    fi
}

################################################
## Step 1: VEP --most_severe (consequence-enrichment table)
################################################
echo "NOTE: Running VEP --most_severe"
"${VEP}" --offline --cache --dir_cache "${vep_cache}" --cache_version "${cache_version}" \
    --assembly GRCh37 --fasta "${grch37_fasta}" --check_ref \
    --format vcf --tab --most_severe --force_overwrite --no_stats --fork "${forks}" \
    -i "${norm_vcf}" -o "${out_dir}/vep_all_most_severe.tsv" \
    --warning_file "${out_dir}/vep_all_most_severe.tsv.warnings"

if [ $? -ne 0 ]; then
    echo "ERROR: VEP --most_severe failed."
    exit 1
fi
check_coverage "${out_dir}/vep_all_most_severe.tsv"

################################################
## Step 2: VEP --per_gene (coding + ABC gene nomination)
################################################
echo "NOTE: Running VEP --per_gene"
"${VEP}" --offline --cache --dir_cache "${vep_cache}" --cache_version "${cache_version}" \
    --assembly GRCh37 --fasta "${grch37_fasta}" --check_ref \
    --format vcf --tab --per_gene --symbol --biotype --force_overwrite --no_stats --fork "${forks}" \
    --fields Uploaded_variation,Location,Consequence,SYMBOL,Gene,BIOTYPE,Feature,CDS_position,Protein_position,Amino_acids,Codons \
    -i "${norm_vcf}" -o "${out_dir}/vep_output.tsv" --warning_file "${out_dir}/vep_output.tsv.warnings"

if [ $? -ne 0 ]; then
    echo "ERROR: VEP --per_gene failed."
    exit 1
fi
check_coverage "${out_dir}/vep_output.tsv"

echo "NOTE: VEP annotation complete. Outputs written to:"
echo "  ${out_dir}/vep_all_most_severe.tsv"
echo "  ${out_dir}/vep_output.tsv"
