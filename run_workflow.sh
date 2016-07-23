#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N RNAseq_Workflow
#$ -e workflow_error.txt
#$ -o workflow_stdout.txt
#$ -pe threaded 24
#$ -l mem_free=64G

# - - - - - - - - - - - - - - - - - - - - - - - -
# @ Thomas W. Battaglia
# @ tb1280@nyu.edu

# Description:
# This is the workflow script for processing and
# analyzing RNAseq data from raw FASTQ sequences.
# There must be an index made for both the STAR
# aligner as well as the SortMeRNA removal tool.
# - - - - - - - - - - - - - - - - - - - - - - - -

# Source virtual environment that contains all
# the required tools.
source venv/bin/activate

# - - - - - - - - - - -
# Variables to set
# - - - - - - - - - - -
qualityCutoff=20
trimLength=50


# - - - - - - - - - - - - - - - - -
# Print cluster info
# - - - - - - - - - - - - - - - - -
echo "- - - - Diagnostics - - - - -"
echo "Number of slots: $NSLOTS"
echo "Number of hosts: $NHOSTS"
echo "Number in Queue: $QUEUE"
echo "OS type: $SGE_ARCH"
echo ""
echo "Run parameters:"
echo "Host name: "$(hostname -f)" "
echo "Date started: $(date)"
echo "Currently in" $(pwd)
echo "Python version: "$(python -V 2>&1)" "
echo "- - - - - - - - - - - - - - -"


# - - - - - - - - - - - -
# Verify workflow can run
# - - - - - - - - - - - -

# No Trim Galore! in PATH
if [[ -z $(type STAR) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1

# No Trim Galore! in PATH
elif [[ -z $(type trim_galore) ]]; then
    echo "No TrimeGaore! installation found! Exiting now."
    exit 1

# No cutadapt in PATH
elif [[ -z $(type cutadapt) ]]; then
    echo "No cutadapt installation found! Exiting now."
    exit 1

# No sortmerna in PATH
elif [[ -z $(type sortmerna) ]]; then
    echo "No SortMeRNA installation found! Exiting now."
    exit 1

# No featureCounts in PATH
elif [[ -z $(type featureCounts) ]]; then
    echo "No featureCounts installation found! Exiting now."
    exit 1

# No FASTQ in PATH
elif [[ -z $(type fastqc) ]]; then
    echo "No FASTQC installation found! Exiting now."
    exit 1

elif [[ -z "$(ls -A index/)" ]]; then
    echo "Error genome index is empty. You must first generate an index before alignment."
    exit 1
fi


# - - - - - - - - - - - -
# Set default variables
# - - - - - - - - - - - -

# Input/Output
inputFiles=input/*
#outputFolder="output.$(date +%F_%R)"
outputFolder="output"
mkdir -p $outputFolder

# Genome/Annotation
genomeIndexDir="index/"
annotationFile=$(ls -d -1 annotation/*)

# QC data
outputQcFolder="${outputFolder}/1_initial_qc"
outputTrimFolder="${outputFolder}/2_trimmed_output"

# SortMeRNA
sortMeRnaAligned="${outputFolder}/3_rRNA/aligned/"
sortMeRnaFiltered="${outputFolder}/3_rRNA/filtered/"
sortMeRnaLogs="${outputFolder}/3_rRNA/logs/"
sortmernaDB="tools/sortmerna-2.1-linux-64"

# Alignment
alignedSequences="${outputFolder}/4_aligned_sequences"
alignedBAM="${outputFolder}/4_aligned_sequences/aligned_bam/"
alignedLog="${outputFolder}/4_aligned_sequences/aligned_logs/"
alignedStat="${outputFolder}/4_aligned_sequences/aligned_stats/"

# featureCounts
finalCounts="${outputFolder}/5_final_counts"


# - - - - - - - - - - -
# Run FASTQC
# - - - - - - - - - - -
echo "Quality analysis of raw reads..."
mkdir -p $outputQcFolder

# Run Quality Analysis on Raw Reads
for seq in $inputFiles; do
    fastqc \
    -o $outputQcFolder \
    --noextract \
    -t $NSLOTS \
    $seq
done


# - - - - - - - - - - -
# Run Trim Galore!
# - - - - - - - - - - -
echo "Trimming raw reads..."
mkdir -p $outputTrimFolder

# Run Trim Galore to remove adapters and low base quality scores
for trim in $inputFiles*; do
    trim_galore \
    --quality $qualityCutoff \
    --fastqc \
    --length $trimLength \
    --output_dir $outputTrimFolder \
    $trim
done

# unzip all sequences (if needed)
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
    sortmerna \
    --ref $sortmernaREF \
    --reads $trim \
    --aligned ${sortMeRnaAligned}/${FQ}_aligned \
    --other ${sortMeRnaFiltered}/${FQ}_filtered \
    --fastx \
    --log \
    -a $NSLOTS \
    -v

    # Move Log Files into correct order
    mv ${sortMeRnaAligned}/${FQ}_aligned.log $sortMeRnaLogs

    # Removed Cutadapt trimmed files to save space
    rm -rf $trim

    # Gzip the aligned rRNA sequences to save space
    gzip ${sortMeRnaAligned}${FQ}_aligned.fq

    # Gzip the sequences that didn't aligned to save space
    gzip ${sortMeRnaFiltered}${FQ}_filtered.fq
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
for seq in $sortMeRnaFiltered*; do

    # Get basename and directory
    seq_dir=$(dirname ${seq})/
    baseName=`basename $seq .fq.gz`

    # Remove the last '_trimmed_filtered' to make the naming cleaner
    baseNameClean=${baseName%???????????????????}

    # STAR
    STAR \
    --genomeDir $genomeIndexDir \
    --readFilesIn $seq \
    --runThreadN $NSLOTS \
    --outFileNamePrefix $alignedSequences/$baseNameClean \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts
done

# Move Log Files into correct folder
for log in $alignedSequences/*Log.final.out; do
    mv $log $alignedLog/
done

# Move *Log.out files
for log in $alignedSequences/*Log.out; do
    mv $log $alignedLog/
done

# Remove *Log.progress.out files
for log in $alignedSequences/*Log.progress.out; do
    rm -rf $log
done

# Remove (.Tab) Files
for seq in $alignedSequences/*out.tab; do
    rm -rf $seq
done

# Remove leftover TEMP folders
for tmp in $alignedSequences/*_STARtmp; do
    rm -rf $tmp
done

# Move Sorted BAM
for seq in $alignedSequences/*.bam; do
    mv $seq $alignedBAM/
done


# - - - - - - - - - - - - - -
# Run Subread (featureCounts)
# - - - - - - - - - - - - - -
echo "Summarizing gene counts..."
mkdir -p $finalCounts

# Store list of files as a variable
dirlist=$(ls -t $alignedBAM/* | tr '\n' ' ')

# Run featureCounts
featureCounts \
-a $annotationFile \
-o $finalCounts/final_counts.txt \
-T $NSLOTS \
$dirlist


# - - - - - - - - - - - - - -
# Run MultiQC 
# - - - - - - - - - - - - - -
multiqc $outputFolder \
--o $outputFolder


# End of script
echo "RNAseq completed!"
echo "Date finished: $(date)"
deactivate
