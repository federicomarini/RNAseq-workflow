RNAseq Workflow
================
Thomas W. Battaglia
14 Sept 2016

### Introduction

RNA seq is an emerging technology for deep RNA sequencing. Not only does RNA seq technology have the ability to analyze differences in gene expression between samples, but can discover new isoforms of genes and analyze SNP variation of particular genes. This tutorial will cover the basic workflow for analyzing RNA seq data on host tissue. It will not cover the workflow for analyzing RNA seq on microbial community samples. For analysis of the meta-transcriptome using RNA seq, see the workflow for Whole Genome Sequencing Metagenomics.

#### 1. Setup new bio-cookie template for a new RNAseq analysis

The easiest way to set up the RNAseq environment is by using a bio-cookie template. This template will download all the necessary tools into a new virtual environment, so that there are no disruptions with your current computer packaging versions. The template will also ask if you would like to download a human or mouse genome for transcriptome alignment. If you want to align to a genome which is neither human or murine, then you must select "None" in the option and download the genome and annotation file separately. Below are a list of the tools which are downloaded during the installation steps.

#### 1a. Load python version 2.7 or higher

``` bash
# Load python 2.7.3 on cluster (if needed)
module load python/2.7.3

# Make sure python 2.7.3 was installed into environment
python --version
>> Python 2.7.3
```

#### 1b. Install cookiecutter template engine

``` bash
pip install cookiecutter --user
```

#### 1c. Set up bio-cookie-WMS environment.

Be sure to choose the correct options related to your data. For more information see: <https://github.com/twbattaglia/bio-cookie-RNAseq>

``` bash
cookiecutter https://github.com/twbattaglia/bio-cookie-RNAseq

Tools that are downloaded during bio-cookie creation:
a. cutadapt - adapter and base quality filtering
b. FastQC - sequence quality analysis
c. Trim Galore! - more easily run cutadapt and fastqc
d. SortMeRNA - remove rRNA sequence contamination
e. STAR-aligner - quickly align sequences to genome
f. Subread/featureCounts - summarize alignments into gene counts
g. MultiQC - generate a readable report of all logs
```
