#!/bin/bash
#$ -S /bin/bash
#$ -o sge_log
#$ -e sge_log
#$ -cwd
#$ -t 1-1000
#$ -tc 100
hostname
date
out_dir=$1
alpha_lvl=$2
python scripts/get_hotspot_residues.py \
    -i $out_dir/data/hotspot/full_output/output_`echo $alpha_lvl`_$((SGE_TASK_ID-1)).txt \
    -s $alpha_lvl \
    -o $out_dir/data/hotspot/residues/hotspot_residues_$((SGE_TASK_ID-1)).txt 
date
