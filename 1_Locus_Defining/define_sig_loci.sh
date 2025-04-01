#!/bin/bash
#PBS -A GeomicVar
#PBS -N define_sig_loci
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:grand
#PBS -l select=1

#Set mapping script location
map_script=/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/1_Locus_Defining/define_sig_loci_for_trait.py

#Set key directories
INPUT_DIR="/grand/projects/GeomicVar/proj-mvp-ukb-meta/output_gia/metaanalysis_for_export/ALL"
OUTPUT_DIR="/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci"
SCRIPT_PATH="/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/1_Locus_Defining/define_sig_loci.py"
TEMP_SCRIPT_DIR="/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/temp_scripts"

#Change directory to temp script directory so output files accumulate there
cd $TEMP_SCRIPT_DIR

#Loop over files and process each one
ls $INPUT_DIR/*.ALL.mrmega.sumstats.txt.gz | while read file; do
	#Extract the trait name and create output file name
	base_file=`basename $file`
	trait=$(echo $base_file | sed 's/.ALL.mrmega.sumstats.txt.gz//g')
	out_file=$OUTPUT_DIR/$trait.sig_loci.csv
	#Check whether the individual trait already exists
	if [ -f "$out_file" ]; then
		echo "Skipping trait $trait since output file already exists."
		continue
	fi
	#Send info to new shell script
	echo ' 
		source ~/.bashrc  
		conda activate finemap
		python '$map_script' -g '$file' -o '$OUTPUT_DIR' -p 5e-8 -s 1e-6 -t 250000 -b '$trait' --chr-col CHR --pos-col POS --p-col P' > $TEMP_SCRIPT_DIR/$trait.define_sig_loci.sh
        sed -i '1i#!/bin/bash'  $TEMP_SCRIPT_DIR/$trait.define_sig_loci.sh
	#Submit the shell script
        qsub -l filesystems=home:grand -N define_$trait -q preemptable -l select=1:ncpus=1 -A GeomicVar -l walltime=1:00:00 $TEMP_SCRIPT_DIR/$trait.define_sig_loci.sh
done

