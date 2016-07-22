---
output: html_document
---
# RNAseq Workflow
#### Thomas W. Battaglia
#### tb1280@nyu.edu
#### 5/17/16
  
  
### Introduction  
RNA seq is an emerging technology for deep RNA sequencing. Not only does RNA seq technology have the ability to analyze differences in gene expression between samples, but can discover new isoforms of genes and analyze SNP variation of particular genes. This tutorial will cover the basic workflow for analyzing RNA seq data on host tissue. It will not cover the workflow for analyzing RNA seq on microbial community samples. For analysis of the meta-transcriptome using RNA seq, see the workflow for Whole Genome Sequencing Metagenomics (WGS Metagenomics) in a different tutorial.

### 1A. Setup new directory for analysis
Many files are generated during the RNA seq workflow, so it is best to keep them stored in an organized directory. The first step before running/installing any tools, is to generate a directory structure for easier organization. The rest of this tutorial assumes your working directory ```pwd``` is within the ```new_analysis``` folder created during the first step.

**input** : per sample sequencing files (.fastq)  
**genome** : directory for storing genome of interest (.fasta/.fna)  
**annotation** : directory for storing the annotation file associated with host genome (.gtf/.gff3)  
**index** : directory that will store STAR aligner's index files   
**tools** : directory to store dependencies required for analysis (SortMeRNA)  
**venv** : a virtual environment of python to store all dependencies

```
# Make new root directory
mkdir rnaseq_analysis

# Change directory into new folder
cd rnaseq_analysis

# Make sub-directories to store:
mkdir input
mkdir genome
mkdir annotation
mkdir index
mkdir tools
```

### 1B. Setup new virtual environment to store tools
A new virtual environment will keep all the tools that are required for analysis in its own folder so that there are no errors intorduced into the system python installation and/or break and tools that require a particular version. Most of the tools require python 2.7+. After creating and activating the virtualenv, you should see ```(venv)``` on the command line.

```bash
# Install the virtualenv package
pip install virtualenv --user

# Make a new virtual environment folder in the 
virtualenv venv --python=python2.7

# Start the virtualenv
source venv/bin/activate

# Create a function to link tools to virtualenv
# Sourced from http://huttenhower.sph.harvard.edu/docs/anadama_workflows/install.html#id1
function link() { ln -sv $(readlink -f "$1") $(pwd)/venv/bin/; }
```

#### Note:
If there is an error during the virtualenv creation step that states ```The executable python2.7 (from --python=python2.7) does not exist```, then the cluster environment's python version is not set for Python 2.7+, so the command ```module load python/2.7.3``` must be run first before running the ```virtualenv venv --python=python2.7``` command.


### 1C. Demo dataset
If you want to use a demo dataset to practice the RNAseq alignment workflow, run one of the commands below to place a fastq file in the ```input``` folder.
```
# Download publically available mouse RNAseq fastq file. (2 biological replicates)
# https://www.encodeproject.org/experiments/ENCSR648YEP/
wget -P input/ https://www.encodeproject.org/files/ENCFF377KCE/@@download/ENCFF377KCE.fastq.gz
wget -P input/ https://www.encodeproject.org/files/ENCFF473WMT/@@download/ENCFF473WMT.fastq.gz

# Download publically available human RNAseq fastq file.
# https://www.encodeproject.org/experiments/ENCFF001RNK/
wget -P input/ (TODO)
```

-----

### 2. Installing required tools
This workflow requires many different tools. Some of these tools maybe available on your cluster environment, but to ensure the correct versioning of the tools, it is suggested to install them again inside a virtual environment. See http://docs.python-guide.org/en/latest/dev/virtualenvs/ for more information about using virtual-environments.

Other tools that are required for processing are:  
A. **cutadapt**  
B. **FastQC**  
C. **Trim Galore!**  
D. **SortMeRNA**  
E. **STAR-aligner**  
F. **Subread**  


-----
#### 2A. Install cutadapt
**Link:** http://cutadapt.readthedocs.io/en/stable/guide.html  
"Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads."  
```bash
pip install cutadapt
```

#### 2B. Install FastQC
**Link:** http://www.bioinformatics.babraham.ac.uk/projects/fastqc/  
"FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis."
```bash
# Download FastQC to the 'tools' folder  
wget -P tools/ http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip  

# Unzip and remove compressed file
unzip tools/fastqc_v0.11.5.zip -d tools/
rm -rf tools/fastqc_v0.11.5.zip

# Change permissions (only for fastqc tool)
chmod 777 tools/FastQC/fastqc

# Link package to virtual environment 
link tools/FastQC/fastqc

# Check if FastQC was installed properly
fastqc --version
```

#### 2C. Install Trim Galore!
**Link:** http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/  
"Trim Galore! is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing)."

```bash
# Download Trim Galore! to the 'tools' folder
wget -P tools/ http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip 

# Unzip and remove compressed file
unzip tools/trim_galore_v0.4.1.zip -d tools/
rm -rf tools/trim_galore_v0.4.1.zip

# Link package to virtual environment 
link tools/trim_galore_zip/trim_galore

# Check if Trim galore! was installed properly
trim_galore --version
```

#### 2D. Install SortMeRNA
**Link:** http://bioinfo.lifl.fr/RNA/sortmerna/  
"SortMeRNA is a biological sequence analysis tool for filtering, mapping and OTU-picking NGS reads. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of nucleotide sequences. The main application of SortMeRNA is filtering rRNA from metatranscriptomic data."  

```bash
# Download SortMeRNA to the 'tools' folder
wget -P tools/ http://bioinfo.lifl.fr/RNA/sortmerna/code/sortmerna-2.1-linux-64.tar.gz 

# Unzip and remove compressed file
tar -zxvf tools/sortmerna-2.1-linux-64.tar.gz -C tools/
rm -rf tools/sortmerna-2.1-linux-64.tar.gz

# Link package to virtual environment 
link tools/sortmerna-2.1-linux-64/indexdb_rna
link tools/sortmerna-2.1-linux-64/sortmerna

# Check if SortMeRNA was installed properly
indexdb_rna -h
sortmerna --version
```

Next we need to generate an index for the SortMeRNA database.
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
indexdb_rna --ref $sortmernaREF 
```

#### 2E. Install STAR-aligner
**Link:** https://github.com/alexdobin/STAR  
"Spliced Transcripts Alignment to a Reference"

```bash
# Download STAR to the 'tools' folder
wget -P tools/ https://github.com/alexdobin/STAR/archive/2.5.2a.zip

# Unzip and remove compressed file
unzip tools/2.5.2a.zip -d tools/
rm -rf tools/2.5.2a.zip

# Link package to virtual environment 
link tools/STAR-2.5.2a/bin/Linux_x86_64/STAR

# Check if STAR was installed properly
STAR --version
```

#### 2F. Install Subread
**Link:** http://subread.sourceforge.net/  
**Link:** http://bioinf.wehi.edu.au/featureCounts/  
"Subread package: high-performance read alignment, quantification and mutation discovery"  

```bash
# Download Subread to the 'tools' folder
wget -P tools/ https://sourceforge.net/projects/subread/files/subread-1.5.0-p3/subread-1.5.0-p3-Linux-x86_64.tar.gz

# Unzip and remove compressed file
tar -zxvf tools/subread-1.5.0-p3-Linux-x86_64.tar.gz -C tools/
rm -rf tools/subread-1.5.0-p3-Linux-x86_64.tar.gz

# Link package to virtual environment 
link tools/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts

# Check if Subread (featureCounts) was installed properly
featureCounts -v
```

-----

### 3. Download host genome and annotation file
A host genome and annotation file are required for a complete RNA seq analysis. The annotation file contains the gene information associated with the coordinates of an alignment. The two files need to places in their respective folders before running the workflow.

#### 3A. Mouse genome
There are different releases of the full mouse genome, so be aware which genome and which annotation file you are using for the sequencing alignment. ***See this paper for more information:***(TODO)  


##### Gencodes (recommended)
http://www.gencodegenes.org/mouse_releases/current.html  
```bash
# Download mouse genome to the 'genome' folder
wget -P genome/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M10/GRCm38.p4.genome.fa.gz

# Download genome annotation file to the 'annotation' folder
wget -P annotation/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gtf.gz

# Download genome annotation file to the 'annotation' folder
wget -P annotation/ http://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm38.85.gtf.gz

# Decompress
gunzip genome/GRCm38.p4.genome.fa.gz
gunzip annotation/gencode.vM9.annotation.gtf.gz
```

##### UNSC
```bash
# (TODO)
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

#### 3C. Other genome
If you are not working with a mouse or human genome, see the following websites for a complete listing.
https://genome.ucsc.edu/cgi-bin/hgGateway  
http://useast.ensembl.org/info/data/ftp/index.html  

-----

### 4. Generate STAR aligner index
The STAR aligner requires an index to be created before aligning any ```fastq``` sequences. The index command requires the host genome of interest and associated annotation file as inputs. The command is best run on the cluster, but should be submitted as a job. Download the shell file and submit it as a job. Additionally you must specific the length of the sequences in the RNAseq run in the command. This is typically 50 or 101, but will depend on what was requested.

```bash
# Download index job file
wget https://raw.githubusercontent.com/twbattaglia/RNAseq-workflow/master/make_index.sh

# Submit the make_index.sh as a job. Set the read length to 101bp.
qsub make_index.sh 101
```

-----

### 5. Download workflow and change parameters
Once the index has been generated and all of the tools are installed, it is time to start the alignment workflow. The workflow is stored as a shell file which has ```for``` loops to iterate over each input ```fastq``` file. This setup allows the shell file to be submitted to the cluster as a job. The workflow will place files generated during processing into the correct folders.

Before submitting the command, make sure you change the variables within the file to reflect your data set.  
- [ ] **qualityCutoff** : The minimum Phred score to retain nucleotides  
- [ ] **trimLength** :  The minimum length of the sequence to keep after quality trimming. (Typically 50% or greater)  

```bash
wget https://raw.githubusercontent.com/twbattaglia/RNAseq-workflow/master/run_workflow.sh
```

-----

### 6. Run the workflow 
The job can be submitted to the cluster once you have changed the parameters to reflect your data and you have placed the ```fastq``` files into the folder, ```input```.
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

Be sure to copy the ```final_counts.txt``` file generate from featureCounts step to your set working directory, or specify the full location when importing the table.  
```R
# Load required libraries
library(DESeq2)
library(ggplot2)

# Set working directory. (if needed)
setwd("~/path/to/working/directory/")

# Import counts table from featureCounts
# Skip first row (or delete before import)
counts <- read.delim("~/path/to/working/directory/final_counts.txt", skip = 1, row.names = 1)

# Remove 'Length' column
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
            quote = F,
            col.names = NA)

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

# Make colored column annotations.
# Choose which column variables you want to annotate the columns by.
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
