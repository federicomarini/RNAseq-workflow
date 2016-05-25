---
output: html_document
---
# RNAseq Workflow
#### Thomas W. Battaglia
#### tb1280@nyu.edu
#### 5/17/16
  
  
### Introduction  
RNA seq is an emerging technology for deep RNA sequencing. Not only does RNA seq technology have the ability to analyze differences in gene expression between samples, but can discover new isoforms of genes and analyze SNP variation of particular genes. This tutorial will cover the basic workflow for analyzing RNA seq data on host tissue. It will not cover the workflow for analyzing RNA seq on microbial community samples. For analysis of the meta-transcriptome using RNA seq, see the workflow for Whole Genome Sequencing Metagenomics (WGS Metagenomics) in a different tutorial.

### 1. Setup new directory for analysis
Many files are generated during the RNA seq workflow, so it is best to keep them stored in an organized directory. The first step before running/installing any tools, is to generate a directory structure for easier organization. The rest of this tutorial assumes your working directory ```pwd``` is within the new folder created below.

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
mkdir input
mkdir genome
mkdir annotation
mkdir index
mkdir tools
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
If the python package Cutadapt is not on your cluster environment, you can install a local copy by running the command below:
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
```

Next we need to generate an index for the SortMeRNA database. This can be run on the head node as it doesn't require too much computation.
```bash
# Set variable for location of index files
sortmernaDB="tools/sortmerna-2.1-linux-64"

# Set location for database
sortmernaREF=${sortmernaDB}/rRNA_databases/silva-arc-16s-id95.fasta,${sortmernaDB}/index/silva-arc-16s-id95:\
${sortmernaDB}/rRNA_databases/silva-arc-23s-id98.fasta,${sortmernaDB}/index/silva-arc-23s-id98:\
${sortmernaDB}/rRNA_databases/silva-bac-16s-id90.fasta,${sortmernaDB}/index/silva-bac-16s-id95:\
${sortmernaDB}/rRNA_databases/silva-bac-23s-id98.fasta,${sortmernaDB}/index/silva-bac-23s-id98:\
${sortmernaDB}/rRNA_databases/silva-euk-18s-id95.fasta,${sortmernaDB}/index/silva-euk-18s-id95:\
${sortmernaDB}/rRNA_databases/silva-euk-28s-id98.fasta,${sortmernaDB}/index/silva-euk-28s-id98

# Generate indexs
tools/sortmerna-2.1-linux-64/indexdb_rna --ref $sortmernaREF 
```

-----

### 3. Download host genome and annotation file
A host genome and annotation file are required for a complete RNA seq analysis. The annotation file contains the gene information associated with the coordinates of an alignment. The two files need to places in their respective folders before running the workflow.


#### 3A. Mouse genome
There are different releases of the full mouse genome, so be aware which genome and which annotation file you are using for the sequencing alignment. ***See this paper for more information:***(TODO)  
  

##### Gencodes
http://www.gencodegenes.org/mouse_releases/current.html
```bash
# Download mouse genome to the 'genome' folder
wget -P genome/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.transcripts.fa.gz

# Download genome annotation file to the 'annotation' folder
wget -P annotation/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gtf.gz

# Decompress
gunzip genome/*
gunzip annotation/*
```

#### Ensembl
http://www.ensembl.org/info/data/ftp/index.html
```bash
# Download mouse genome to the 'genome' folder
wget -P genome/ ftp://ftp.ensembl.org/pub/release-84/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

# Download genome annotation file to the 'annotation' folder
wget -P annotation/ ftp://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.gtf.gz

# Decompress
gunzip genome/*
gunzip annotation/*
```

#### UNSC
```bash

```


#### 3B. Human genome

```bash
# Download human genome to the 'genome' folder
wget -p genome/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.transcripts.fa.gz

# Download genome annotation file to the 'annotation' folder
wget -P annotation/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz

# Decompress
gunzip genome/gencode.v24.transcripts.fa.gz
gunzip annotation/gencode.v24.annotation.gtf.gz
```

##### 3C. Other genome
If you are not working with a mouse or human genome, see the following websites for a complete listing.
https://genome.ucsc.edu/cgi-bin/hgGateway  
http://useast.ensembl.org/info/data/ftp/index.html  

-----

### 4. Generate STAR aligner index
The STAR aligner requires an index to be created before aligning any ```fastq``` sequences. The index command requires the host genome of interest and associated annotation file as inputs. The command is best run on the cluster, but should be submitted as a job. Download the shell file and submit it as a job.

```bash
# Download index job file
wget https://raw.githubusercontent.com/twbattaglia/RNAseq-workflow/master/make_index.sh

# Submit the make_index.sh as a job 
qsub make_index.sh
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
Once the workflow has completed, you can now use the gene count table as an input into DESeq2 for statistical analysis. See this tutorial for more information about using DESeq2 http://www.bioconductor.org/help/workflows/rnaseqGene/

One additional required file which is needed, is a type of mapping file. This can be created in R or it can be imported as a text file. The mapping file must have sample identifiers that match the the resulting table with more columns that describe the sample (e.g Treatment).

#####  First install the required packages for analysis. It is best to use an IDE like RSudio for the analysis: https://www.rstudio.com/products/rstudio/download/
```R
# Install required libraries
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("ggplot2")
biocLite("clusterProfiler")
biocLite("biomaRt")
biocLite("ReactomePA")
biocLite("DOSE")
biocLite("pathview")
biocLite("org.Mm.eg.db")
biocLite("pheatmap")
biocLite("genefilter")
biocLite("fdrtool")
biocLite("RColorBrewer")
```

Be sure to copy the final_counts.txt file generate from featureCounts, to your working directory, or specific it when importing the table.  
```R
# Load required libraries
library(DESeq2)
library(ggplot2)

# Set working directory. (if needed)
setwd("~/path/to/working/directory/")

# Import counts table from featureCounts
# Skip first row (or delete before import)
counts <- read.delim("~/path/to/working/directory/final_counts.txt", skip = 1, row.names = 1)

# Remove Length column
counts <- counts[,-1]

# Import metadata (or create a new dataframe)
# The sample identifers must be the row names for the dataframe and must match the names of the counts table columns.
metadata <- read.delim("sample_mapping_file.txt", row.names = 1)

# Relevel (if needed) to know which group is the reference (control)
metadata$Group <- relevel(metadata$Group, ref = "Control")

# Make DESeq2 object from featureCounts object 
# countData : count dataframe
# colData : sample metadata in the dataframe with row names as sampleID's
# design : The design of the comparisons to use. Use (~) before the name of the column variable to compare
ddsMat <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = metadata,
                                 design = ~Group)
                                 
# Run DESEq2
ddsMat <- DESeq(ddsMat)

# Get results from testing with FDR adjust pvalues
res_out <- results(ddsMat, pAdjustMethod = "fdr")

# order results by padj value 
res_out <- res_out[order(res_out$padj), ]

# Generate summary of testing. 
# q-value cutoff is 1% by default.
summary(res_out)

```

-----

### 8a. Add gene annotation information to results table
Depending upon the data set, you may have to change the database for gene annotation.  
**Human** : ```hsapiens_gene_ensembl```   
**Mouse** : ```mmusculus_gene_ensembl```   
**Squirrel** : ```itridecemlineatus_gene_ensembl```  
```R
# Run this command to get a list of available marts
biomaRt::listDatasets(useEnsembl(biomart="ensembl"))
```
##### Get GeneID's for each identifier and add it to the results table
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
                 mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), 
                 uniqueRows = T)

# Get row numbers that correspond to geneID's and annotation data            
idx <- match(row.names(res_out), geneIDs$ensembl_gene_id)

# Make a new column of gene annotation data
res_out$ensembl_gene_id <- row.names(res_out)
res_out$gene_name <- geneIDs$external_gene_name[idx]
res_out$description <- geneIDs$description[idx]
res_out$entrez <- geneIDs$entrezgene[idx]
            
# Subset for only significant genes (q<0.05)
res_out_sig <- subset(res_out, padj < 0.05)


```

##### Write all the important results to .txt files
```R
# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ddsMat), normalized = T), 
            file = 'DESeq2_normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(ddsMat[row.names(res_out_sig)], normalized = T), 
            file = 'DESeq2_normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)
            
# Write the annotated results table to a .txt file
write.table(x = as.data.frame(res_out), 
            file = "DESEq2_results_gene_annotated.txt", 
            sep = '\t', 
            quote = F)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(res_out_sig), 
            file = "DESEq2_results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)
```

-----


### 8b. Make heatmap of significant genes 
Make a heatmap of the significant genes using the pheatmap package.
```R
# Load libraries
library(pheatmap)
library(RColorBrewer)

# rlog transformation of the raw count data. 
# Necessary to show data in plots
rld <- rlog(ddsMat)

# Gather top genes and make matrix
mat <- assay(rld[row.names(res_out_sig)])

# Make column annotation.
# Choose which column variables you want to annotate
# the columns by.
annotation_col = data.frame(
  Treatment = factor(colData(rld)$Group),
  row.names = colData(rld)$SampleID
)

# Make Heatmap with pheatmap function.
# See more in documentation for customization
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", 
         annotation_col = annotation_col,
         fontsize_row = 8)
```

-----

### 9. Pathway analysis on DE genes. 
The code below assumes you have RNA seq results from the mouse genome.  
Get more information about clusterProfiler here: http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
```R
# Load required libraries
library(clusterProfiler)
library(ReactomePA)
library(KEGG.db)
library(DOSE)
library(pathview)
library(org.Mm.eg.db)

# Create a matrix of gene entrez ID's and log fold changes
gene_matrix <- res_out_sig$logFC
names(gene_matrix) <- res_out_sig$entrez


# - - - - - - - - - - - - - 
# Enrich with KEGG database
# - - - - - - - - - - - - - 
kegg_enrich <- enrichKEGG(gene = gene_matrix,
                 organism = 'mouse',
                 pvalueCutoff = 0.05, 
                 readable = TRUE)
                 
# Get table of results
summary(kegg_enrich)

# Plot results
barplot(kegg_enrich, drop=TRUE, showCategory=12)

# Plot specific KEGG pathways (w fold change) with pathview
# pathway.id : KEGG pathway
pathview(gene.data = gene_matrix, pathway.id = "04940", species = "mouse", map.symbol = T)


# - - - - - - - - - - - - - 
# Enrich with GO databse
# - - - - - - - - - - - - - 
go_enrich <- groupGO(gene = gene_matrix,
               organism = 'mouse',
               ont = "BP",
               level = 3,
               readable = TRUE)

# Get table of results
summary(go_enrich)

# Plot results
barplot(go_enrich, drop=TRUE, showCategory=12)


# - - - - - - - - - - - - - 
# Enrich with ReactomeDB
# - - - - - - - - - - - - - 
reactome_enrich <- enrichPathway(gene = gene_matrix, organism = 'mouse', pvalueCutoff = 0.05)
                 
# Get table of results
summary(reactome_enrich)
```


### Citations:

Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

Martin, Marcel. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, [S.l.], v. 17, n. 1, p. pp. 10-12, may. 2011. ISSN 2226-6089. Available at: <http://journal.embnet.org/index.php/embnetjournal/article/view/200>. doi:http://dx.doi.org/10.14806/ej.17.1.200. 

Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611

Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21. doi:10.1093/bioinformatics/bts635.

Lassmann et al. (2010) "SAMStat: monitoring biases in next generation sequencing data." Bioinformatics doi:10.1093/bioinformatics/btq614 [PMID: 21088025] 

Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30. 

Love MI, Huber W and Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, pp. 550.

Yu G, Wang L, Han Y and He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” OMICS: A Journal of Integrative Biology, 16(5), pp. 284-287. 
