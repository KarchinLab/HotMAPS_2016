# HotMAPS

Hotspot Missense mutation Areas in Protein Structures (HotMAPS) detects 3D somatic mutation hotspots in PDB structures.
Identical chains are accounted for in the model to prevent errors induced by assuming chains
are independent. For example, an annotated mutation in a gene forming a homodimer
will always contain the same mutation in both chains. PDB biological assemblies are
used when available to best reveal biologically meaningful 3D hotspots.

Running the is generally done through using `make` commands that invoke scripts
found in the scripts directory or top-level scripts.

## Install

HotMAPS is intended to be ran on **linux** operating systems.

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

For consistency of processing protein chain descriptions found in the PDB file to that performed in MuPIT, we use BioJava to annotate each chain and add that label as a column in the input file to our pipeline. Download the BioJava 3.1 [core](http://biojava.org/download/maven/org/biojava/biojava3-core/3.1.0/biojava3-core-3.1.0.jar) and [structure](http://biojava.org/download/maven/org/biojava/biojava3-structure/3.1.0/biojava3-structure-3.1.0.jar) jar files and place them in a directory called "lib".

## HotMAPS Pipeline

### Initial Setup

First, download the Protein Data Bank (PDB) structures from [ftp://ftp.wwpdb.org/pub/pdb/](ftp://ftp.wwpdb.org/pub/pdb/) and the theoretical protein structure models ([ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013/](ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013/)). You will need both the RefSeq and Ensembl theoretical protein structure models (H_sapiens_2013.tar.xz and ModBase_H_sapiens_2013_refseq.tar.xz, respectively). Then update the config.txt to point toward the directories that you save the structure files at after extracting from compressed format.  Additionally, download the [mutations file](http://karchinlab.org/data/HotMAPS/mutations.txt.gz), [protein structure annotation file](http://karchinlab.org/data/HotMAPS/pdb_info.txt.gz), and  annotations for the CRAVAT reference transcript available [here](http://karchinlab.org/data/HotMAPS/mupit_annotations.tar.gz). Place all three files in a sub-directory called "data". Assuming you are already in the HotMAPS directory:

```bash
$ mkdir -p data
$ cd data
$ wget http://karchinlab.org/data/HotMAPS/mutations.txt.gz
$ gunzip mutations.txt.gz
$ wget http://karchinlab.org/data/HotMAPS/pdb_info.txt.gz
$ gunzip pdb_info.txt.gz
$ wget http://karchinlab.org/data/HotMAPS/mupit_annotations.tar.gz
$ tar xvzf mupit_annotations.tar.gz
$ cd ..
```

Assuming you have changed the `config.txt` file to point towards where you downloaded the protein structure files, an additional step is needed to annotate those protein structures.

```bash
$ make annotateStructures
```

### Running 3D HotMAPS

To run the code in parallel using Sun Grid Engine (SGE) execute the following make command:

```bash
$ make OUTPUT_DIR=myoutput_dir runParallelHotspot
```

To run the code normally (no parallelization) execute:

```bash
$ make OUTPUT_DIR=myoutput_dir runNormalHotspot
```

`myoutput_dir` is the output directory (Default: output/all_pdb_run).

Note if you ran the normal version instead of parallel, you need not run this next step
as the merged file will already be produced. To merge the output from the parallel 
runs use the following make command:

```bash
$ make OUTPUT_DIR=myoutput mergeHotspotFiles
```

Next, the p-values need to be adjusted for multiple hypotheses testing. 
This needs the CRAVAT reference transcript files noted in the `Initital Setup`
section that was saved in the "data" sub-directory (parameter `MUPIT_ANNOTATION_DIR` in the make command).

```bash
$ make multipleTestCorrect OUTPUT_DIR=myoutput MUPIT_ANNOTATION_DIR=annotation_dir Q_VALUE=myqvalue 
```

`myqvalue` is the q-value for the False Discovery Rate (FDR) correction (.01 by default). The next step is group significant residues
into regions. If you are interested in regions on the actual PDB protein structure,
script use the following command:

```bash
$ make findHotregionStruct OUTPUT_DIR=myoutput_dir Q_VALUE=myqvalue MUPIT_ANNOTATION_DIR=annotation_dir
```

Where like before `myoutput_dir` is the output directory and `myqvalue` is the
q-value (Default: .01). Similarly the regions can be constructed for each gene
using the reference transcript selected by CRAVAT for each mutation. 

```bash
$ make findHotregionGene OUTPUT_DIR=myoutput_dir Q_VALUE=myqvalue MUPIT_ANNOTATION_DIR=annotation_dir
```

### 1D Version

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

### Advanced (optional)

The input data for HotMAPS can also be prepared by directly using the MuPIT MySQL database. A MySQL dump of the MuPIT database containing mutation counts in our study and associated tables that map genome coordinates to PDB structures is available [here](http://karchinlab.org/data/HotMAPS/mupit_modbase.sql.gz). The MuPIT database has a fairly large file size, you may want to directly download and upload to MYSQL.

```bash
$ wget http://karchinlab.org/data/HotMAPS/mupit_modbase.sql.gz
$ gunzip mupit_modbase.sql.gz
$ mysql [options] < mupit_modbase.sql
```

This will create a database named `mupit_modbase`, where `[options]` is the necessary MySQL parameters to login.

Next, the input files need to be generated before starting the `Running 3D HotMAPS` section. The initial input information
is retrieved from the MuPIT MySQL database, in contrast to downloading already made files (done in the `Initial Setup` section). To prepare the input files simply invoke the following make command.

```bash
$ make MYSQL_USER=myuser MYSQL_HOST=myhost MYSQL_DB=mydb prepareHotspotInput
```

Where `myuser` is your MySQL user name, `myhost` is the host name for MySQL, and `mydb` is the database name
for Mupit (Default: mupit_modbase).
