# RNAseq Workflow
#### Thomas W. Battaglia
#### tb1280@nyu.edu
#### 5/17/16
  
  
### Introduction  
RNAseq is an emerging technology for deep RNA sequencing. Not only does RNAseq technology have the ability to analyze differences in gene exporession between samples, but can discover new isoforms of genes and analyze SNP variation of particular genes. This tutorial will cover the basic workflow for analysing RNAseq data on host tissue. It will not cover the workflow for analyzing RNAseq on microbial community samples. For analysis of the meta-transcriptome using RNAseq, see the workflow for Whole Genome Sequencing Metagenomics (WGS Metagenomics) in a differnt tutorial.

### 1. Setup new directory for analysis
Many files are generated during the RNAseq workflow, so it is best to keep them stored in an organized fashion. The first step before running/installing any tools, is to generate a directory structure for easier organization. 

**input** : per sample sequencing files (.fastq)  
**genome** : directory for storing genome of interest (.fasta/.fna)  
**annotation** : directory for storing the annotation file assocaited with host genome (.gtf/.gff3)  
**index** : directory thatwill store STAR aligner's index  
**tools** : directory to store depdencies required for analysis (SortMeRNA)  

```
# Make new root directory
mkdir new_analysis

# Make sub-directories to store:
mkdir -p new_analysis/input
mkdir -p new_analysis/genome
mkdir -p new_analysis/annotation
mkdir -p new_analysis/index
mkdir -p new_analysis/tools
```


### 2. Installing required dependencies
This workflow requires many different tools, but many to all of them are available on the phoenix cluster. The tutorial below assumes you are using the NYULMC phoneix cluster. The only package that may not be available directly on the cluster is SortMeRNA. Follow the steps below to install a local copy of SortMeRNA on the cluster environment.

Other tools that are required for processing are:
FASTQC (link)
Trim Galore! (link)
STAR-aligner (link)
samstat (link)
Subread (link)

#### Install SortMeRNA and generate index's

```bash
# Download SortMeRNA to cluster
wget

# Unzip

# 
```


### 3. Download host genome and annotation file 

### 4. Generate STAR aligner index

### 5. Place sequence files in 'input' folder and change parameters
Once the index has been generated and all of the tools are installed, it is time to start the alignment workflow. The workflow is stored as a shell file which has ```for``` loops to iterate over each input ```fastq``` file. This setup allows the shell file to be submitted to the cluster as a job. The workflow will place files generated during processing into the correct folders.

Before submitting the command, make sure you change the variables within the file to reflect your dataset.  
**sequenceLength** :  
**adapterSequence** : The nucloetide sequence associated with  
**adapterSequence** : The nucloetide sequence associated with  

The workflow shell file can be downloaded here: (TODO)

### 6. Run the workflow 

### 7. Import count table to DESeq2 and run differential expression analysis
### 8. Pathway analysis on resulting DE genes

