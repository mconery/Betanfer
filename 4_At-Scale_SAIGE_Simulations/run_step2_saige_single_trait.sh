#!/bin/bash
#Extract command line parameters
chunk=$1
chromo=$2
pop_dir=$3
bgen_file=$4
manifest_file=$5
sample_file=$6
var_list=$7
saige_step2_script=$8
minMAF=$9
minMAC=$10
LOCO=$11
is_fastTest=$12
is_Firth_beta=$13
pCutoffforFirth=$14
nThreads=$15

# Start timing
start_time=$(date +%s.%N)

#Strip chr from chromo
chrom=${chromo#chr}

#Run step 2
Rscript $saige_step2_script --bgenFile=$bgen_file --bgenFileIndex=$bgen_file.bgi --manifestFile=$manifest_file --sampleFile=$sample_file --chrom=$chrom --AlleleOrder=ref-first --minMAF=$minMAF --minMAC=$minMAC --LOCO=$LOCO --is_fastTest=$is_fastTest --is_Firth_beta=$is_Firth_beta --idstoIncludeFile=$var_list --pCutoffforFirth=$pCutoffforFirth --nThreads=$nThreads
# End timing and calculate duration
end_time=$(date +%s.%N)
duration=$(echo "$end_time - $start_time" | bc -l)

#Move Files
mv $pop_dir/$chunk/*.$chromo.step2.txt $pop_dir

# Save timing results
echo "$chunk,$chromo,$duration" >> $pop_dir/timing2.csv