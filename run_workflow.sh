#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N RNAseq_Workflow
#$ -e workflow_error.txt
#$ -o workflow_stdout.txt
#$ -pe threaded 4
#$ -l mem_free=8G

# - - - - - - - - - - - - - - - - - - - - - - - -
# @ Thomas W. Battaglia
# @ tb1280@nyu.edu

# Description:
# This is the workflow script for processing and
# analyzing RNAseq data from raw FASTQ sequences.
# There must be an index made for both the STAR
# aligner as well as the SortMeRNA removal tool
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
if [[ -z $(STAR --version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1

# No Cutadapt in PATH
elif [[ -z $(cutadapt --version) ]]; then
    echo "No CutAdapt installation found! Exiting now."
    exit 1

# No FASTQ in PATH
elif [[ -z $(fastqc --version) ]]; then
    echo "No FASTQC installation found! Exiting now."
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
outputQcFolder="${outputFolder}/initial_qc"
outputTrimFolder="${outputFolder}/trimmed_output"

sortMeRnaAligned="${outputFolder}/rRNA/aligned/"
sortMeRnaFiltered="${outputFolder}/rRNA/filtered/"
sortMeRnaLogs="${outputFolder}/rRNA/logs/"
sortmernaDB="tools/sortmerna-2.1-linux-64"

alignedSequences="${outputFolder}/aligned_sequences"
alignedBAM="${outputFolder}/aligned_sequences/aligned_bam/"
alignedLog="${outputFolder}/aligned_sequences/aligned_logs/"
alignedStat="${outputFolder}/aligned_sequences/aligned_stats/"

finalCounts="${outputFolder}/final_counts"


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
mkdir -p $outputTrimFolder

# Run Trim Galore to remove adapters and low base quality scores
for trim in $inputFiles*; do
  tools/trim_galore_zip/trim_galore \
  --quality $qualityCutoff \
  --fastqc \
  --length $trimLength \
  --output_dir $outputTrimFolder \
  $trim
done

# unzip all sequences
gunzip $outputTrimFolder/*.fq.gz


# - - - - - - - - - - -
# Run SortMeRNA
# - - - - - - - - - - -
echo "Removing rRNA sequences..."
mkdir -p $sortMeRnaAligned
mkdir -p $sortMeRnaFiltered
mkdir -p $sortMeRnaLogs

# Initialize Database
sortmernaREF=${sortmernaDB}/rRNA_databases/silva-arc-16s-id95.fasta,${sortmernaDB}/index/silva-arc-16s-id95:\
${sortmernaDB}/rRNA_databases/silva-arc-23s-id98.fasta,${sortmernaDB}/index/silva-arc-23s-id98:\
${sortmernaDB}/rRNA_databases/silva-bac-16s-id90.fasta,${sortmernaDB}/index/silva-bac-16s-id95:\
${sortmernaDB}/rRNA_databases/silva-bac-23s-id98.fasta,${sortmernaDB}/index/silva-bac-23s-id98:\
${sortmernaDB}/rRNA_databases/silva-euk-18s-id95.fasta,${sortmernaDB}/index/silva-euk-18s-id95:\
${sortmernaDB}/rRNA_databases/silva-euk-28s-id98.fasta,${sortmernaDB}/index/silva-euk-28s-id98


# Align to rRNA databases
for trim in $outputTrimFolder/*.fq; do

    # Get basename
    FQ=`basename $trim .fq`

    # run sortmerna
    tools/sortmerna-2.1-linux-64/sortmerna \
    --ref $sortmernaREF \
    --reads $trim \
    --aligned ${sortMeRnaAligned}/${FQ}_aligned \
    --other ${sortMeRnaFiltered}/${FQ}_filtered \
    --fastx \
    --log \
    -a 12 \
    -v

    # Move Log Files into correct order
    mv ${sortMeRnaAligned}/${FQ}_aligned.log $sortMeRnaLogs

    # Removed Cutadapt Trimmed Files to save space
    rm -rf $trim

    # Gzip the sequences that aligned to rRNA databases to save space
    gzip ${sortMeRnaAligned}${FQ}_aligned.fq
done


# - - - - - - - - - - -
# Run STAR-aligner
# - - - - - - - - - - -
echo "Aligning sequences..."
mkdir -p $alignedSequences
mkdir -p $alignedBAM
mkdir -p $alignedLog
mkdir -p $alignedStat

# Run STAR-aligner


# Move Log Files into correct folder
for log in $alignedSequences/*Log.final.out; do
    mv $log $alignedLog/
done

# Move Samstat to correct folder
for html in $alignedSequences/*.html; do
    mv $html $alignedStat/
done

# Remove *Log.out files
for log in $alignedSequences/*Log.out; do
    rm -rf $log
done

# Remove *Log.out files
for log in $alignedSequences/*Log.progress.out; do
    rm -rf $log
done

# Remove (.Tab) Files
for seq in $alignedSequences/*out.tab; do
    rm -rf $seq
done


# - - - - - - - - - - - - -
# Run samtools + move BAM
# - - - - - - - - - - - - -
for seq in $alignedSequences/*.bam; do
    samstat $seq
    mv $seq $alignedBAM/
done







