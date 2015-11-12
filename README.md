# HotMAPS

Hotspot Missense mutation Areas in Protein Structures (HotMAPS) detects 3D somatic mutation hotspots in PDB structures.
Identical chains are accounted for in the model to prevent errors induced by assuming chains
are independent. For example, an annotated mutation in a gene forming a homodimer
will always contain the same mutation in both chains. PDB biological assemblies are
used when available to best reveal biologically meaningful 3D hotspots.

Running the is generally done through using `make` commands that invoke scripts
found in the scripts directory or top-level scripts.

## Install

### Python

The python code depends on the following python libraries (listed in requirements.txt):

* numpy
* scipy
* biopython

The necessary python libraries can be installed through pip.

```bash
$ pip install -r requirements.txt
```

The exact version numbers may not be necessary but are what it was tested on.

### Java

For consistency of processing protein chain descriptions found in the PDB file to that performed in MuPit, we use BioJava to annotate each chain and add that label as a column in the input file to our pipeline. Download the BioJava 3.1 [core](http://biojava.org/download/maven/org/biojava/biojava3-core/3.1.0/biojava3-core-3.1.0.jar) and [structure](http://biojava.org/download/maven/org/biojava/biojava3-structure/3.1.0/biojava3-structure-3.1.0.jar) jar files and place them in a directory called "lib".

## 3D Hotspot Detection Pipeline

### Initial Setup

First, download the Protein Data Bank (PDB) structures from [ftp://ftp.wwpdb.org/pub/pdb/](ftp://ftp.wwpdb.org/pub/pdb/) and the theoretical protein structure models from [here](ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013/). Then update the config.txt to point toward the directories that you save the structure files at. A MySQL dump of the MuPIT database containing mutation counts in our study and associated tables that map genome coordinates to PDB structures is available [here](http://karchinlab.org/data/HotMAPS/mupit_modbase.sql.gz). The MuPIT database has a fairly large file size, you may want to directly download and uploade to MYSQL.

```bash
$ wget http://karchinlab.org/data/HotMAPS/mupit_modbase.sql.gz
$ gunzip mupit_modbase.sql.gz
$ mysql [options] < mupit_modbase.sql
```

This will create a database named `mupit_modbase`. Additionally, download the mutation annotations for the CRAVAT reference transcript [here](http://karchinlab.org/data/HotMAPS/mupit_annotations.tar.gz), and then place in a sub-directory called "data".

### Running 3D HotMAPS

First, the input files need to be generated. The initial input information
is retrieved from the Mupit MySQL database. To prepare the input files
simply invoke the following make command.

```bash
$ make MYSQL_USER=myuser MYSQL_DB=mydb prepareHotspotInput
```

Where `myuser` is your MySQL user name and `mydb` is the database name
for Mupit (Default: mupit_modbase).

To run the code in parallel using Sun Grid Engine (SGE) execute the following make command:

```bash
$ make OUTPUT_DIR=myoutput_dir runParallelHotspot
```

To run the code normally (no parallelization) execute:

```bash
$ make OUTPUT_DIR=myoutput_dir runNormalHotspot
```

`myoutput_dir` is the output directory (Default: output/all_pdb_run).

To merge the output from the parallel runs use the following make command:
Note if you ran the normal version instead of parallel, you need not run this step
the merged file will already be produced

```bash
$ make OUTPUT_DIR=myoutput mergeHotspotFiles
```

Next, the p-values need to be adjusted for multiple hypotheses testing. 
This needs the mutation annotation files noted in the `Initital Setup`
section that was saved in the "data" sub-directory (parameter `MUPIT_ANNOTATION_DIR` in the make command).

```bash
$ make multipleTestCorrect OUTPUT_DIR=myoutput MUPIT_ANNOTATION_DIR=annotation_dir Q_VALUE=myqvalue 
```

`myqvalue` is the q-value for the False Discovery Rate (FDR) correction (.01 by default). The next step is group significant residues
into regions. If you are interested in regions on the actual PDB protein structure,
then you will utilize the `find_hotspot_regions_struct.py` script use the following
commmand:

```bash
$ make findHotregionStruct OUTPUT_DIR=myoutput_dir Q_VALUE=myqvalue MUPIT_ANNOTATION_DIR=annotation_dir
```

Where like before `myoutput_dir` is the output directory and `myqvalue` is the
q-value (Default: .01). Similarly the regions can be constructed for each gene
using the reference transcript selected by CRAVAT for each mutation. 

```bash
$ make findHotregionGene OUTPUT_DIR=myoutput_dir Q_VALUE=myqvalue MUPIT_ANNOTATION_DIR=annotation_dir
```

### 1D Hotspot Detection Pipeline

This is very similar to the 3D detection one but instead of the 3D versions
of the make commands, the 1D versions are used instead. If you have 
a SGE cluster, the code can be executed in parallel:

```bash
$ make runParallelHotspot1D OUTPUT_DIR=myoutput_dir
```

If not, then run the non-parallel version:

```bash
$ make runNormalHotspot1D OUTPUT_DIR=myoutput_dir
```

Only for the parallel execution do you need to execute the
next merging step. Otherwise continue without running.

```bash
$ make mergeHotspotFiles1D OUTPUT_DIR=myoutput_dir
```

Perform multiple testing correction on the p-values.

```bash
$ make multipleTestCorrect1D OUTPUT_DIR=myoutput MUPIT_ANNOTATION_DIR=annotation_dir Q_VALUE=myqvalue 
```

Construct hotspot regions for each structure.

```bash
$ make findHotregionStruct1D OUTPUT_DIR=myoutput_dir Q_VALUE=myqvalue MUPIT_ANNOTATION_DIR=annotation_dir
```

Construct merged hotspot regions evidenced from all PDBs to form the gene-level
regions.

```bash
$ make findHotregionGene1D OUTPUT_DIR=myoutput_dir Q_VALUE=myqvalue MUPIT_ANNOTATION_DIR=annotation_dir
```
