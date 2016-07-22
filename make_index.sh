#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N STAR_Index
#$ -e star_index_error.txt
#$ -o star_index_stdout.txt
#$ -pe threaded 24
#$ -l mem_free=64G

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
export genomeFile=$(ls -d -1 genome/*)
export indexFolder="index/"
export annotationFile=$(ls -d -1 annotation/*)


# - - - - - - - - - - - -
# Verify script can run
# - - - - - - - - - - - -

# No STAR installation
if [[ -z $(STAR --version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1
fi

# Unzip genome (if compressed)
gunzip $genomeFile


# Build an index for alignment. Only needs to be run once.
STAR \
--runMode genomeGenerate \
--genomeDir $indexFolder \
--genomeFastaFiles $(pwd)/$genomeFile \
--sjdbGTFfile $(pwd)/$annotationFile \
--sjdbOverhang $readLength \
--runThreadN $NSLOTS

# Clean up temporary folder
rm -rf _STARtmp/

# Message completed
echo 'STAR Index generated!'
deactivate