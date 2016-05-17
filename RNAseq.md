# RNAseq Workflow
#### Thomas W. Battaglia
#### tb1280@nyu.edu
#### 5/17/16
  
  
### Introduction  
RNAseq is an emerging technology for deep RNA sequencing. Not only does RNAseq technology have the ability to analyze differences in gene expression between samples, but can discover new isoforms of genes and analyze SNP variation of particular genes. This tutorial will cover the basic workflow for analyzing RNA seq data on host tissue. It will not cover the workflow for analyzing RNA seq on microbial community samples. For analysis of the meta-transcriptome using RNA seq, see the workflow for Whole Genome Sequencing Metagenomics (WGS Metagenomics) in a different tutorial.

### 1. Setup new directory for analysis
Many files are generated during the RNAseq workflow, so it is best to keep them stored in an organized directory. The first step before running/installing any tools, is to generate a directory structure for easier organization. The rest of this tutorial assumes your working directory ```pwd``` is within the new folder created below.

**input** : per sample sequencing files (.fastq)  
**genome** : directory for storing genome of interest (.fasta/.fna)  
**annotation** : directory for storing the annotation file associated with host genome (.gtf/.gff3)  
**index** : directory that will store STAR aligner's index  
**tools** : directory to store dependencies required for analysis (SortMeRNA)  

```
# Make new root directory
mkdir new_analysis

# Change directory into new folder
cd new_analysis

# Make sub-directories to store:
mkdir -p input
mkdir -p genome
mkdir -p annotation
mkdir -p index
mkdir -p tools
```

-----

### 2. Installing required dependencies
This workflow requires many different tools, but many to all of them are available on the phoenix cluster. The tutorial below assumes you are using the NYULMC phoenix cluster. The only package that may not be available directly on the cluster is SortMeRNA. Follow the steps below to install a local copy of SortMeRNA on the cluster environment.

Other tools that are required for processing are:  
**FASTQC** (link)  
**Trim Galore!** (link)  
**STAR-aligner** (link)  
**samstat** (link)  
**Subread** (link)  
**CutAdapt** (link)  
If the python package CutAdapt is not on your cluster environment, you can install a local copy by running the command below:
```bash
pip install cutadapt --user
```

-----

#### 2A. Install Trim Galore!
Trim Galore is a tool developed to run FASTQC and Cutadapt in one step. The tool can analyze and trim in the same step, saving time and computation processing performance. It is not normally installed on any cluster environment, so it must be installed separately.

```bash
# Download Trim Galore! to the 'tools' folder
wget -P tools/ http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip 

# Unzip and remove compressed file
unzip tools/trim_galore_v0.4.1.zip -d tools/
rm -rf tools/trim_galore_v0.4.1.zip

# Check if Trim galore! was installed properly
tools/trim_galore_zip/trim_galore --version
```


#### 2B. Install SortMeRNA and generate index's

```bash
# Download SortMeRNA to the 'tools' folder
wget -P tools/ http://bioinfo.lifl.fr/RNA/sortmerna/code/sortmerna-2.1-linux-64.tar.gz 

# Unzip and remove compressed file
tar -zxvf tools/sortmerna-2.1-linux-64.tar.gz -C tools/
rm -rf tools/sortmerna-2.1-linux-64.tar.gz

# Check if SortMeRNA was installed properly
tools/sortmerna-2.1-linux-64/indexdb_rna --help
tools/sortmerna-2.1-linux-64/sortmerna --help

# Generate indexes for SortMeRNA
(TODO)

```

-----

### 3. Download host genome and annotation file
A host genome and annotation file are required for a complete RNAseq analysis. The annotation file contains the gene information associated with the coordinates of an alignment. The two files need to places in their respective folders before running the workflow.

#### 3A. Human genome
```bash
# Download human genome to the 'genome' folder
wget -p genome/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.transcripts.fa.gz

# Download genome annotation file to the 'annotation' folder
wget -P annotation/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz

# Decompress
gunzip genome/gencode.v24.transcripts.fa.gz
gunzip annotation/gencode.v24.annotation.gtf.gz
```

##### 3B. Mouse genome
```bash
# Download mouse genome to the 'genome' folder
wget -P genome/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.transcripts.fa.gz

# Download genome annotation file to the 'annotation' folder
wget -P annotation/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gtf.gz

# Decompress
gunzip genome/gencode.vM9.transcripts.fa.gz
gunzip annotation/gencode.vM9.annotation.gtf.gz
```

##### 3C. Other genome
If you are not working with a mouse or human genome, see the following websites for a complete listing.
https://genome.ucsc.edu/cgi-bin/hgGateway  
http://useast.ensembl.org/info/data/ftp/index.html  

-----

### 4. Generate STAR aligner index
The STAR aligner requires an index to be created before aligning any ```fastq``` sequences. The index command requires the host genome of interest and associated annotation file as inputs. The command is best run on the cluster, but should be submitted as a job. Download the shell file and submit it as a job. There is only one argument needed, which is the total length of the reads in your dataset.

```bash
# Download index job file
wget https://raw.githubusercontent.com/twbattaglia/RNAseq-workflow/master/make_index.sh

# Submit the make_index.sh as a job with read legnth as an argument
qsub make_index.sh 50
```

-----

### 5. Download workflow and change parameters
Once the index has been generated and all of the tools are installed, it is time to start the alignment workflow. The workflow is stored as a shell file which has ```for``` loops to iterate over each input ```fastq``` file. This setup allows the shell file to be submitted to the cluster as a job. The workflow will place files generated during processing into the correct folders.

```bash
wget https://raw.githubusercontent.com/twbattaglia/RNAseq-workflow/master/run_workflow.sh
```

Before submitting the command, make sure you change the variables within the file to reflect your data set.  
**qualityCutoff** : The minimum Phred score to retain nucleotides  
**trimLength** :  The minimum length of the sequence to keep after quality trimming. (Typically 50% or greater)  

-----

### 6. Run the workflow 
The job can be submitted to the cluster once you have changed the parameters to reflect your data and you have placed the ```fastq``` files into the folder, 'input'.
```bash
qsub run_workflow.sh
```

-----

### 7. Import count table to DESeq2 and run differential expression analysis
Once the workflow has completed, you can now use the gene count table as an input into DESeq2 for statistical analysis. 

One additional required file which is needed is a type of mapping file. This can be created in R or it can be imported as a text file. The mapping file must have sample identifiers that match the the resulting table with more columns that describe the sample (e.g Treatment).
```R
# Install required libraries

# Load required libraries
library(DESEq2)
library(ggplot2)

# Import counts table from featureCounts

# Import metadata (or create a new dataframe)

# Relevel so we know whats control group
metadatas$Group <- relevel(targets$Group, ref = "CON")

# Make Deseq2 object from featureCounts object -----------
ddsMat <- DESeqDataSetFromMatrix(countData = countsTable$counts,
                                 colData = metadata,
                                 design = ~Group)
                                 
# Run DESEq2
ddsMat <- DESeq(ddsMat)

# Get results from testing
res_out <- results(ddsMat, pAdjustMethod = "fdr")

# Generate summary of testing
summary(res_out)
```

-----

### 8. Add gene annotation information to results table
Depending upon the dataset, you may have to change the database for gene annotation.  
**Human** : ```hsapiens_gene_ensembl```   
**Mouse** : ```itridecemlineatus_gene_ensembl```   
**Squirrel** : ```itridecemlineatus_gene_ensembl```  
**More** :    
```R
# Run this command to get a list of available marts
biomaRt::listDatasets(useEnsembl(biomart="ensembl"))
```

```R
# Load library
library(biomaRt)

# Get information from ENSEMBL (change mart depending upon species)
geneIDs <- getBM(filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id", 
                                "external_gene_name", 
                                "description", 
                                "entrezgene"), 
                 values = row.names(res_out), 
                 mart = useMart("ensembl", dataset = ""), 
                 uniqueRows = T)
                 
idx <- match(row.names(res_out), geneIDs$ensembl_gene_id)
res_out$ensembl_gene_id <- row.names(res_out)
res_out$gene_name <- geneIDs$external_gene_name[idx]
res_out$description <- geneIDs$description[idx]
res_out$entrez <- geneIDs$entrezgene[idx]

# Show only significant genes
sig_res <- subset(res_out, padj < 0.05)

# Get summary of significant genes
summary(sig_res)

```


-----

### 9. Pathway analysis on resulting DE genes
```R
# Load required libraries
library(clusterProfiler)

# Convert geneID's to EntrezID or use column create from annotation step.


```


### Citations:

