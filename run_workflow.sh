#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N RNAseq_Workflow
#$ -e worflow_error.txt
#$ -o workflow_stdout.txt
#$ -pe threaded 16

# - - - - - - - - - - - - - - - - - - - - - - - -
# @ Thomas W. Battaglia
# @ tb1280@nyu.edu

# Description:

# - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - -
# Variables to set
# - - - - - - - - - - - 
sequenceLength=50
qualityCutoff=20
trimLength=20
trimAdapterSeq="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


# - - - - - - - - - - -
# Load required modules
# - - - - - - - - - - - 
module load STAR


# - - - - - - - - - - - -
# Verify workflow can run
# - - - - - - - - - - - -
# No Trim Galore! in PATH
# No Cutadapt in PATH
# No FASTQ in PATH
# No STAR in PATH
# No SortMeRNA in PATH
# No index files in index folder
# More than one file in genome or annotation folder


# - - - - - - - - - - - -
# Set default variables
# - - - - - - - - - - - -

# Input/Output
inputFiles=input/*.fastq.gz
outputFolder="output.$(date +%F_%R)/"
mkdir -p outputFolder

# Genome/Annotation
genomeFile=genome/genome/*
indexFolder="index/"
annotationFile="annotation/*"

# SortMeRNA Location
sortMeRNALoc="tools/sortmerna/sortmerna"
sortMeRNADB="tools/sortmerna"


# - - - - - - - - - - -
# Variables to set
# - - - - - - - - - - - 


# - - - - - - - - - - -
# Run FASTQC
# - - - - - - - - - - - 


# - - - - - - - - - - -
# Run Trim Galore!
# - - - - - - - - - - - 


# - - - - - - - - - - -
# Run SortMeRNA
# - - - - - - - - - - - 


# - - - - - - - - - - -
# Run STAR-aligner
# - - - - - - - - - - - 


# - - - - - - - - - - -
# Run samtools
# - - - - - - - - - - - 


# - - - - - - - - - - - - - - 
# Run Subread (featureCounts)
# - - - - - - - - - - - - - -







