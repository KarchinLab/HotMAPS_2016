#########################################
# Parameters for running hotspot pipeline
#########################################
# mysql information
MYSQL_HOST=karchin-db01.icm.jhu.edu 
MYSQL_USER=collin
MYSQL_DB=mupit_modbase

# directory containing output files
OUTPUT_DIR=output/all_pdb_run/

# q-value threshold for significance
Q_VALUE=.01
# number of simulations 
NUM_SIM=10000
# sphere radius for residue neighbors
RADIUS=10.0
# window size to look in sequence
# for 1D algorithm
FLANK_1D=3

# annotation input for pdb files
PDB_INFO=data/fully_described_pdb_info.txt
# mutation file from mupit
MUT_FILE=data/mutations.txt
# directory containing PDB_INFO and MUT_FILE
# split into pieces for parallelization
SPLIT_DIR=data/split_pdbs/

# temp data files
pdb_info_init=data/pdb_info.txt
TEMP_DIR=tmp/

GROUP_FUNC=min

###################################
# Paths to external tools/libraries
###################################

##################################################
# Directories containing mutations and their 
# annotations
##################################################
# Directory for merged annotation info
MUPIT_ANNOTATION_DIR=data/mupit/mupit_annotations_10_27_2015/

##################################
# Prepare input files for hot spot
# detetection
##################################
# Get info about PDB, chain, and gene names
# important for running on known structures
getPDBInfo:
	mkdir -p data
	mysql -u ${MYSQL_USER} -A -p -h ${MYSQL_HOST} ${MYSQL_DB} < scripts/sql/get_pdb_info.sql > ${pdb_info_init}

# get mutations from mupit database
getMutations:
	mysql -u ${MYSQL_USER} -A -p -h ${MYSQL_HOST} ${MYSQL_DB} < scripts/sql/get_mutations.sql > ${MUT_FILE}

# add file path information for pdb files
getPDBPath:
	python scripts/add_path_info.py -p ${pdb_info_init} -o data/pdb_info.path.txt

# get the chain desciption for the PDB files
getPDBDescription:
	python scripts/chain_description.py -i data/pdb_info.path.txt -o ${PDB_INFO}

# split input files for parallelization
splitInputFiles:
	python scripts/divide_pdb_info.py \
		-f ${PDB_INFO} \
		-m ${MUT_FILE} \
		--split-dir ${SPLIT_DIR}

# Run all commands for preparing input for hot spot detection code
prepareHotspotInput: getPDBInfo getMutations getPDBPath getPDBDescription splitInputFiles
annotateStructures: getPDBPath getPDBDescription splitInputFiles

#####################################
# Run hotspot code
#####################################
# run the 3D hotspot code in parallel on the cluster
runParallelHotspot:
	# create output directories if needed
	mkdir -p ${OUTPUT_DIR}
	mkdir -p ${OUTPUT_DIR}/data/hotspot/full_output
	mkdir -p ${OUTPUT_DIR}/data/hotspot/residues
	mkdir -p ${OUTPUT_DIR}/error
	# run hotspot code in parallel
	qsub -N PDB2HOTSPOT -v PATH=$$PATH scripts/qsub/run_parallel_hotspot.sh ${SPLIT_DIR} ${NUM_SIM} ${RADIUS} ${OUTPUT_DIR}

# run hotspot without the parallel aspect
runNormalHotspot:
	mkdir -p ${OUTPUT_DIR}
	mkdir -p ${OUTPUT_DIR}/data/hotspot/full_output
	mkdir -p ${OUTPUT_DIR}/data/hotspot/residues
	mkdir -p ${OUTPUT_DIR}/error
	python hotspot.py --log-level=INFO \
		-m ${MUT_FILE} \
		-a ${PDB_INFO} \
		-t EVERY \
		-n ${NUM_SIM} \
		-r ${RADIUS} \
		-o ${OUTPUT_DIR}/output_merged.txt \
		-e ${OUTPUT_DIR}/error/error_pdb_${PDB_INFO}.txt \
		--log=stdout

# merge all files about hotspots together
mergeHotspotFiles:
	rm -f ${OUTPUT_DIR}/output_merged.txt
	cat ${OUTPUT_DIR}/data/hotspot/full_output/output_* | awk -F"\t" 'NR==1 || $$1!="Structure"' > ${OUTPUT_DIR}/output_merged.txt

# Multiple testing correction
#
# NOTE: the annotation results from CRAVAT
# are needed. Thus, please run commands in the next
# section before doing multiple testing correction.
multipleTestCorrect:
	python multiple_testing_correction.py \
		-i ${OUTPUT_DIR}/output_merged.txt \
		-f ${GROUP_FUNC} \
		-m ${MUPIT_ANNOTATION_DIR} \
		-q ${Q_VALUE} \
		-o ${OUTPUT_DIR}/mtc_output_${GROUP_FUNC}_${Q_VALUE}.txt \
		-s ${OUTPUT_DIR}/significance_level_${Q_VALUE}.txt

# find hotspot regions (i.e. collection of residues) in structures
findHotregionStruct:
	python find_hotspot_regions_struct.py \
		-i ${OUTPUT_DIR}/output_merged.txt \
		-a ${MUPIT_ANNOTATION_DIR} \
		-p ${PDB_INFO} \
		-r ${RADIUS} \
		-o ${OUTPUT_DIR}/hotspot_regions_structure_${Q_VALUE}.txt \
		-s ${OUTPUT_DIR}/significance_level_${Q_VALUE}.txt \
	    --log=stdout

# find hotspot regions (i.e. collection of residues) for gene 
# rather then for a single structure
findHotregionGene:
	python find_hotspot_regions_gene.py \
		-m ${OUTPUT_DIR}/mtc_output_min_${Q_VALUE}.txt \
		-a ${MUPIT_ANNOTATION_DIR} \
		-p ${PDB_INFO} \
		-r ${RADIUS} \
		-q ${Q_VALUE} \
		-o ${OUTPUT_DIR}/hotspot_regions_gene_${Q_VALUE}.txt \
		--log=stdout
