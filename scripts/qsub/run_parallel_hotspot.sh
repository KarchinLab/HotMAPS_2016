#!/bin/bash
#$ -S /bin/bash
#$ -o sge_log
#$ -e sge_log
#$ -cwd
#$ -t 1-1000
#$ -tc 100
hostname
date
split_dir=$1
num_sim=$2
radius=$3
out_dir=$4
echo "python hotspot.py --log-level=INFO \
    -m $split_dir/mut_info_split_$((SGE_TASK_ID-1)).txt \
    -a $split_dir/pdb_info_split_$((SGE_TASK_ID-1)).txt \
    -t EVERY \
    -n $num_sim \
    -r $radius \
    -o $out_dir/data/hotspot/full_output/output_${signif_lvl}_$((SGE_TASK_ID-1)).txt \
    -e $out_dir/error/error_pdb_${signif_lvl}_$((SGE_TASK_ID-1)).txt \
    --log=stdout "
python hotspot.py --log-level=INFO \
    -m $split_dir/mut_info_split_$((SGE_TASK_ID-1)).txt \
    -a $split_dir/pdb_info_split_$((SGE_TASK_ID-1)).txt \
    -t EVERY \
    -n $num_sim \
    -r $radius \
    -o $out_dir/data/hotspot/full_output/output_${signif_lvl}_$((SGE_TASK_ID-1)).txt \
    -e $out_dir/error/error_pdb_${signif_lvl}_$((SGE_TASK_ID-1)).txt \
    --log=stdout
date
