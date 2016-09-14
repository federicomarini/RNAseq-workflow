RNAseq Workflow
================
Thomas W. Battaglia
14 Sept 2016

### Introduction

RNA seq is an emerging technology for deep RNA sequencing. Not only does RNA seq technology have the ability to analyze differences in gene expression between samples, but can discover new isoforms of genes and analyze SNP variation of particular genes. This tutorial will cover the basic workflow for analyzing RNA seq data on host tissue. It will not cover the workflow for analyzing RNA seq on microbial community samples. For analysis of the meta-transcriptome using RNA seq, see the workflow for Whole Genome Sequencing Metagenomics.

#### 1. Setup new bio-cookie template for a new RNAseq analysis

The easiest way to set up the RNAseq environment is by using a bio-cookie template. This template will download all the necessary tools into a new virtual environment, so that there are no disruptions with your current computer packaging versions. It is high reccomended to run on all of these commands on a high performance cluster.

The template will also ask if you would like to download a human or mouse genome for transcriptome alignment. If you want to align to a genome which is neither human or murine, then you must select "None" in the option and download the genome and annotation file separately. Below are a list of the tools which are downloaded during the installation steps. See here for a list of availabel genomes (<https://genome.ucsc.edu/cgi-bin/hgGateway> / <http://useast.ensembl.org/info/data/ftp/index.html> )

##### 1a. Load python version 2.7 or higher

``` bash
# Load python 2.7.3 on cluster (if needed)
module load python/2.7.3

# Make sure python 2.7.3 was installed into environment
python --version
>> Python 2.7.3
```

##### 1b. Install cookiecutter template engine

``` bash
pip install cookiecutter --user
```

##### 1c. Set up bio-cookie-WMS environment.

Be sure to choose the correct options related to your data. For more information see: <https://github.com/twbattaglia/bio-cookie-RNAseq>

``` bash
# Generate a bio-cookie RNAseq template project
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

#### 1d. Change directory into the new folder

``` bash
cd new_rnaseq_analysis/
```

### 2. Generate STAR aligner index

The STAR aligner requires an index to be created before aligning any `fastq` sequences. The index command requires the host genome of interest and associated annotation file as inputs. The command is best run when submitted as a job. Download the shell file and submit it as a job. Additionally you must specific the length of the sequences in the RNAseq run in the command. This is typically 50 or 101, but will depend on what was requested.

``` bash
# Submit the make_index.sh as a job. This can take ~20 minutes to finish.
qsub make_index.sh
```

### 3. Run the workflow in parallel

The job can be submitted to the cluster once you generated a STAR-aligner index and you have placed the `.fastq.gz` files into the folder, `input`.

The command is written to submit an array job, which means you must specify how many times the job should run. In this case it is dependent upon how many individual `.fastq.gz` files that are located in the `input` folder. If there are 6 different samples within the folder you would run the command `qsub -t 1:6 run_workflow.sh`, but if you have 10 files, you would run the command `qsub -t 1:10 run_workflow.sh`. The `-t` flag corresponds to how many different commands (files) should be submitted.

``` bash
# Submit the array workflow script
qsub -t 1:4 run_workflow.sh
```

**Note:** The command is written for SGE clusters, but if you have a different cluster set up, you can change the variable `$SGE_TASK_ID` to the correct array style variable.

### 4. Process the output files and generate a report

After the main workflow has finished running, you will have generated aligned `.bam` files for each sample. These files only give the coordinates for the sequence alignments on the genome, and must be summarized to give gene counts per sample, before any statistical analysis can be performed. To do this we perform further processing the output files using the `postprocessing.sh` script There is no need for a parallel array script and the command can be submitted as is as long as no folder names/locations were not changed in main workflow step.

``` bash
qsub postprocessing.sh
```

Once completed, many files are generated after running the main workflow and the post-processing scripts. Below is a breakdown of each of the folders and their respective description of their outputs. Each step reports a log file for every sample which can be used to calculcated the number of sequences.

``` bash
── output/
  │   └── 1_initial_qc/ - FastQC quality reports for each sample
  │   └── 2_trimmed_output/ - Trimmed reports for each sample
  │   └── 3_rRNA
  │       ├── aligned/ - Sequences aligned to rRNA databases
  │       ├── filtered/ - Sequences with rRNA sequences removed
  │       ├── logs/ - Logs from running SortMeRNA
  │   └── 4_aligned_sequences
  │       ├── aligned_bam/ - Main alignment files for each sample
  │       ├── aligned_logs/ - Log from running STAR alignment step
  │       ├── aligned_counts/ - STAR alignment counts output (for comparison with featureCounts)
  │   └── 5_final_counts - Final gene count summarization output from running featureCounts
  │   └── 6_multiQC - Overall workflow report
```

### 5. Differential expression analysis using R

Once the workflow has completed, you can now use the gene count table as an input into DESeq2 for statistical analysis using the R-programming language. It is highly reccomended to use **RStudio** when writing R code and generating R-related analyses. You can download RStudio for your system here: <https://www.rstudio.com/products/rstudio/download/>

#### 5a. Install required R-libraries

Multiple libraries are required to perform analysis, plotting and enrichment.

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("ggplot2")
biocLite("clusterProfiler")
biocLite("biomaRt")
biocLite("ReactomePA")
biocLite("DOSE")
biocLite("KEGG.db")
biocLite("pathview")
biocLite("org.Mm.eg.db")
biocLite("org.Hs.eg.db")
biocLite("pheatmap")
biocLite("genefilter")
biocLite("RColorBrewer")
biocLite("topGO")
biocLite("dplyr")
```

#### 5b. Import featureCounts output

One you have an R environment appropriatley set up, you can begin to import the featureCounts output of genes counts, found within the `5_final_counts` folder. This tutorial will use DESeq2 to normalize and perform the statistical analysis between sample groups. Be sure to copy the `final_counts.txt` file generate from featureCounts step to your set working directory, or specify the full location when importing the table. <br>
If you would like to use an example featureCounts output, download the gene counts and metadata here:
[Download gene counts per sample](https://raw.githubusercontent.com/twbattaglia/RNAseq-workflow/master/example/final_counts.txt)
[Download sample metadata](https://raw.githubusercontent.com/twbattaglia/RNAseq-workflow/master/example/metadata.txt)
[Download R command script](https://raw.githubusercontent.com/twbattaglia/RNAseq-workflow/master/example/rnaseq_commands.R)

``` r
# - - - - - - - - - - - - - - -
# Import featureCounts data
# - - - - - - - - - - - - - - -
library(DESeq2) # statistical analysis 
library(ggplot2) # plotting 
library(knitr) # for better formatting

# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
countdata <- read.table("example/final_counts.txt", header = TRUE, skip = 1, row.names = 1)

# Remove .bam + '..' from column identifiers
colnames(countdata) <- gsub(".bam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("..", "", colnames(countdata), fixed = T)

# Remove length/char columns
countdata <- countdata[ ,c(-1:-5)]

# Make sure ID's are correct
head(countdata)
```

    ##               SRR1374924 SRR1374923 SRR1374921 SRR1374922
    ## 4933401J01Rik          0          0          0          0
    ## Gm26206                0          0          0          0
    ## Xkr4                 214        302        459        425
    ## Gm18956                0          0          0          0
    ## Gm37180                4          2          3          1
    ## Gm37363                1          0          0          1

<br>

##### Import metadata text file. The SampleID's must be the first column.

``` r
# Import metadata file
# - make row names the matching sampleID's from the countdata
metadata <- read.delim("example/metadata.txt", row.names = 1)

# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)

# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]

# Make sure ID's are correct
head(metadata)
```

    ##            Group Replicate   sampleid
    ## SRR1374924 HiGlu      Rep2 SRR1374924
    ## SRR1374923 HiGlu      Rep1 SRR1374923
    ## SRR1374921 LoGlu      Rep1 SRR1374921
    ## SRR1374922 LoGlu      Rep2 SRR1374922

##### Make DESeq2 object from counts and metadata

``` r
# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~Group)


# Find differential expressed genes
# Run DESEq2
ddsMat <- DESeq(ddsMat)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing
