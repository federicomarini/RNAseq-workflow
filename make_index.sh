#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N STAR_Index
#$ -e star_index_error.txt
#$ -o star_index_stdout.txt
#$ -pe threaded 12
#$ -l mem_free=16G

# - - - - - - - - - - - - - - - - - - - - - - - - - - -
# @ Thomas W. Battaglia
# tb1280@nyu.edu

# make_index.sh
# Argument 1 = The length of your RNAseq sequences.

# This script will generate index files for STAR-aligner
# before use with RNAseq data alignment. This file can
# can be submitted as a job if on a proper cluster
# enviroment.
# - - - - - - - - - - - - - - - - - - - - - - - - - -

# Source virtualevn
source venv/bin/activate


# - - - - - - - - - - - -
# Set variables
# - - - - - - - - - - - -
export readLength=`expr $1 - 1`
export genomeFile=$(ls genome/)
export indexFolder="index/"
export annotationFile=$(ls annotation/)


# - - - - - - - - - - - -
# Verify script can run
# - - - - - - - - - - - -

# No STAR installation
if [[ -z $(STAR --version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1
fi

# Unzip genome (if compressed)
if [ ${$genomeFile: -3} == ".gz" ]; then
    gunzip $genomeFile
fi

# Build an index for alignment. Only needs to be run once.
STAR \
--runMode genomeGenerate \
--genomeDir $indexFolder \
--genomeFastaFiles $genomeFile \
--sjdbGTFfile $annotationFile \
--sjdbOverhang $readLength \
--limitGenomeGenerateRAM 17179869184 \
--runThreadN 12

# Message completed
echo 'STAR Index generated!'
deactivate
