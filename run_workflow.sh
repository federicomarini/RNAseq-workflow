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









# Create folder locations for outputs
#-----------------------------------------------------------#
InitialQCFolder="${OutputFolder}initial_qc/"

TrimmedFolder="${OutputFolder}trimmed_files/reads/"
TrimmedReports="${OutputFolder}trimmed_files/post_trim_qc/"

SortmernaAligned="${OutputFolder}rRNA/aligned/"
SortmernaFiltered="${OutputFolder}rRNA/filtered/"
SortmernaLogs="${OutputFolder}rRNA/logs/"

AlignedSequences="${OutputFolder}aligned_sequences"
AlignedBAM="${OutputFolder}/aligned_sequences/aligned_bam/"
AlignedLog="${OutputFolder}/aligned_sequences/aligned_logs/"
AlignedStat="${OutputFolder}/aligned_sequences/aligned_stats/"




#-----------------------------------------------------------#
# Start Analysis!
#-----------------------------------------------------------#

# FastQC
#-----------------------------
if [ "$InitalQC" -eq 1 ]; then

    echo "Processing Raw Reads for Quality Analysis..."
    mkdir -p $InitialQCFolder

    # Run Quality Analysis on Raw Reads
    for seq in $InputFiles*; do
        fastqc -o $InitialQCFolder --noextract -t 12 $seq
    done


fi



# Cutadapt
#-----------------------------------
if [ "$RemoveAdapters" -eq 1 ]; then

    echo "Filtering Low Quality Bases and Adapters..."
    mkdir -p $TrimmedFolder
    mkdir -p $TrimmedReports

    # Cutadapt Trimming
    for seq in $InputFiles*; do
        FQ=`basename $seq .fastq.gz`
        cutadapt -q $QualityCutoff -m $TrimLength -a $TrimAdapterSeq -o ${TrimmedFolder}${FQ}_trimmed.fq $seq
    done


    # FastQC Trimmed Data
    for seq in $TrimmedFolder*.fq; do
        fastqc --noextract -o $TrimmedReports -t 12 $seq
    done

fi



# SortMeRNA
#-----------------------------
if [ "$RemoveRNA" -eq 1 ]; then
    echo "Removing rRNA Sequences..."
    chmod +x packages/sortmerna/sortmerna
    mkdir -p $SortmernaAligned
    mkdir -p $SortmernaFiltered
    mkdir -p $SortmernaLogs

    # Initialize Database
    SortmernaREF=${SortmernaDB}/rRNA_databases/silva-arc-16s-id95.fasta,${SortmernaDB}/index/silva-arc-16s-id95:${SortmernaDB}/rRNA_databases/silva-arc-23s-id98.fasta,${SortmernaDB}/index/silva-arc-23s-id98:${SortmernaDB}/rRNA_databases/silva-bac-16s-id90.fasta,${SortmernaDB}/index/silva-bac-16s-id95:${SortmernaDB}/rRNA_databases/silva-bac-23s-id98.fasta,${SortmernaDB}/index/silva-bac-23s-id98:${SortmernaDB}/rRNA_databases/silva-euk-18s-id95.fasta,${SortmernaDB}/index/silva-euk-18s-id95:${SortmernaDB}/rRNA_databases/silva-euk-28s-id98.fasta,${SortmernaDB}/index/silva-euk-28s-id98


    # Align to rRNA databases
    for trim in $TrimmedFolder*.fq; do
        FQ=`basename $trim .fq`

        packages/sortmerna/sortmerna --ref $SortmernaREF --reads $trim --aligned ${SortmernaAligned}/${FQ}_aligned --other ${SortmernaFiltered}/${FQ}_filtered --fastx --log -a 16 -v


        # Move Log Files into correct order
        mv ${SortmernaAligned}/${FQ}_aligned.log $SortmernaLogs

        # Removed Cutadapt Trimmed Files to save space
        rm $trim

        # Gzip the sequences that aligned to rRNA databases to save space
        gzip ${SortmernaAligned}${FQ}_aligned.fq
    done

fi




# STAR Index Builder
#-----------------------------
# Build an index for alignment. Only needs to be run once.
if [ "$IndexGenome" -eq 1 ]; then
    chmod +x packages/STAR/bin/Linux_x86_64/STAR
    packages/STAR/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir $GenomeIndex --genomeFastaFiles $Genome --sjdbGTFfile $GenomeGTF --sjdbOverhang 49 --runThreadN 16
fi





# Align Sequences using STAR
#-----------------------------
if [ "$AlignSequences" -eq 1 ]; then

    # Make Directories
    mkdir -p $AlignedSequences
    mkdir -p $AlignedBAM
    mkdir -p $AlignedLog
    mkdir -p $AlignedStat

    # Not needed each time. Only after first run due to permissions error
    chmod +x packages/STAR/bin/Linux_x86_64/STAR

    # Align to STAR aligner
    for seq in $SortmernaFiltered*; do

        # Variables--------
        seq_dir=$(dirname ${seq})/
        baseName=`basename $seq .fq`

        # STAR--------------
        packages/STAR/bin/Linux_x86_64/STAR --genomeDir $GenomeIndex --readFilesIn $seq --runThreadN 16 --outFileNamePrefix $AlignedSequences/$baseName --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat

        # gzip to save space
        gzip $seq
    done



    # Move results into the correct places-----------

    # Run Samstat to get alignment stats. Move (.BAM) files into proper folder.
    for seq in $AlignedSequences/*.bam; do
        samstat $seq
        mv $seq $AlignedBAM/
    done

    # Move Log Files into correct folder
    for log in $AlignedSequences/*Log.final.out; do
        mv $log $AlignedLog/
    done

    # Move Samstat to correct folder
    for html in $AlignedSequences/*.html; do
        mv $html $AlignedStat/
    done





    # Remove unnecessary output files----------------

    # Remove *Log.out files
    for log in $AlignedSequences/*Log.out; do
        rm $log
    done

    # Remove *Log.out files
    for log in $AlignedSequences/*Log.progress.out; do
        rm $log
    done

    # Remove (.Tab) Files
    for seq in $AlignedSequences/*out.tab; do
        rm $seq
    done

fi


