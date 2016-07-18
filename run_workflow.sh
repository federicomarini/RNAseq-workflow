#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N RNAseq_Workflow
#$ -e workflow_error.txt
#$ -o workflow_stdout.txt
#$ -pe threaded 16
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

# Source virtual environment that contains all the required tools
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
hostname -f
date
pwd
python --version
echo "- - - - - - - - - - - - - - -"


# - - - - - - - - - - - -
# Verify workflow can run
# - - - - - - - - - - - -

# No Trim Galore! in PATH
if [[ -z $(STAR --version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1

# No Trim Galore! in PATH
elif [[ -z $(trim_galore --version) ]]; then
    echo "No TrimeGaore! installation found! Exiting now."
    exit 1
    
# No sortmerna in PATH
elif [[ -z $(sortmerna --version) ]]; then
    echo "No SortMeRNA installation found! Exiting now."
    exit 1
    
# No Cutadapt in PATH
elif [[ -z $(featureCounts --version) ]]; then
    echo "No featureCounts installation found! Exiting now."
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
inputFiles=input/*
#outputFolder="output.$(date +%F_%R)"
outputFolder="outputFolder"
mkdir -p $outputFolder

# Genome/Annotation
genomeFile=genome/*
indexFolder=index/
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
  tools/trim_galore_zip/trim_galore \
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
    tools/sortmerna-2.1-linux-64/sortmerna \
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
for seq in $sortMeRnaFiltered*; do

    # Get basename and directory
    seq_dir=$(dirname ${seq})/
    baseName=`basename $seq .fq.gz`

    # Remove the last '_trimmed_filtered' to make the naming cleaner
    baseNameClean=${baseName%?????????????????}

    # STAR
    STAR \
    --genomeDir $indexFolder \
    --readFilesIn $seq \
    --runThreadN $NSLOTS \
    --outFileNamePrefix $alignedSequences/$baseNameClean \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --readFilesCommand zcat

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


# - - - - - - - - - - - - -
# Run samtools + move BAM
# (TODO)
# - - - - - - - - - - - - -
for seq in $alignedSequences/*.bam; do
    samstat $seq
    mv $seq $alignedBAM/
done

# Move Samstat to correct folder
for html in $alignedSequences/*.html; do
    mv $html $alignedStat/
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


# End of script
echo "RNAseq completed!"
deactivate
