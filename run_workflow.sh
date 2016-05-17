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
qualityCutoff=20
trimLength=20


# - - - - - - - - - - -
# Load required modules
# - - - - - - - - - - - 
module load samstat/1.09
module load star/2.5.0c
module load subread/1.4.6-p3
module load fastqc/0.11.4


# - - - - - - - - - - - -
# Verify workflow can run
# - - - - - - - - - - - -

# No Trim Galore! in PATH
if [[ -z $(STAR -version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1

# No Cutadapt in PATH
elif [[ -z $(cutadapt -version) ]]; then
    echo "No Cutadapt installation found! Exiting now."
    exit 1

# No FASTQ in PATH
elif [[ -z $(fastqc -version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1
    
# No STAR in PATH
elif [[ -z $(STAR -version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1
    
# No SortMeRNA in PATH
elif [[ -z $(tools/sortmerna/sortmerna -version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1
    
# No index files in index folder
elif [[ -z $(STAR -version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1
    
# More than one file in genome or annotation folder
elif [[ -z $(STAR -version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1
fi


# - - - - - - - - - - - -
# Set default variables
# - - - - - - - - - - - -

# Input/Output
inputFiles=input/*.fastq.gz
outputFolder="output.$(date +%F_%R)"
mkdir -p outputFolder

# Genome/Annotation
genomeFile=genome/*
indexFolder="index/"
annotationFile=annotation/*

# SortMeRNA Location
sortMeRNALoc="tools/sortmerna/sortmerna"
sortMeRNADB="tools/sortmerna"


# - - - - - - - - - - -
# Variables to set
# - - - - - - - - - - - 
export $outputQcFolder="${outputFolder}/initial_qc"
export $outputTrimFolder="${outputFolder}/trimmed_output"


# - - - - - - - - - - -
# Run FASTQC
# - - - - - - - - - - - 
echo "Quality analysis of raw reads..."
mkdir -p $outputQcFolder

# Run Quality Analysis on Raw Reads
for seq in $inputFiles; do
  fastqc -o $outputQcFolder --noextract -t 12 $seq
done


# - - - - - - - - - - -
# Run Trim Galore!
# - - - - - - - - - - - 
echo "Trimming raw reads..."
mkdir -p $outputQcFolder

# Run Trim Galore to remove adapters and low base quality scores
for trim in $inputFiles*; do
  tools/trim_galore_zip/trim_galore \
  --quality $qualityCutoff \
  --fastqc \
  --gzip \
  --length $trimLength \
  --output_dir $outputTrimFolder \
  $trim
done


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

