#!/bin/bash

#########################Locus-specific and Simulation parameters####################

#Set directories and locations
locus_prefix=$1
raw_geno_dir=$2
out_dir=$3
image_dir=$4
biogwas_script=$5
dependency_file=$6
gencode_gtf=$7

####################################################################################

#Activate conda directory
source ~/.bashrc
conda activate biogwas

#Extract components of the file prefix needed for processing files
if [[ $locus_prefix =~ ^([^_]+)_([^_]+)_causal_per_locus_chr([^_]+)_([^_]+)_([^_]+)$ ]]; then
    pop="${BASH_REMATCH[1]}"
    num_var="${BASH_REMATCH[2]}"
    chromosome="${BASH_REMATCH[3]}"
    start="${BASH_REMATCH[4]}"
    end="${BASH_REMATCH[5]}"
fi

################################################################################################################
####################################### Make Locus-Specific Plink Files ########################################
################################################################################################################
plink2 --bfile $raw_geno_dir/$pop/$pop"_chr"$chromosome"_biogwas" --chr $chromosome --from-bp $start --to-bp $end --make-bed --out $out_dir/$locus_prefix
echo $out_dir/$locus_prefix",chr"$chromosome > $out_dir/$locus_prefix.genotypes.list


################################################################################################################
###################################### Simulate the GWAS using bioGWAS #########################################
################################################################################################################
file_prefix=$out_dir/$locus_prefix
file_pattern=$locus_prefix
geno_file=$file_prefix".genotypes.list"
variant_file=$file_prefix".variants.txt"
out_file_pattern=$locus_prefix"_out"
gen_var=$(echo $num_var*0.0001 | bc)
mkdir -p $out_dir
cd $out_dir
python $biogwas_script -gv $gen_var --data_dir $out_dir --img_dir $image_dir -bfile --input_list $geno_file --anno_file $gencode_gtf --use_causal_snps --causal_snps $variant_file --skip_simulation -d $dependency_file --pattern $out_file_pattern