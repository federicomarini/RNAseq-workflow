# RNAseq Workflow
#### Thomas W. Battaglia
#### 2/13/16
  
  
### Introduction  
RNAseq is an emerging technology for deep RNA sequencing. Not only does RNAseq technology have the ability to analyze differences in gene exporession between samples, but can discover new isoforms of genes and analyze SNP variation of particular genes. This tutorial will cover the basic workflow for analysing RNAseq data on host tissue. It will not cover the workflow for analyzing RNAseq on microbial community samples. For analysis of the meta-transcriptome using RNAseq, see the workflow for Whole Genome Sequencing Metagenomics (WGS Metagenomics) in a differnt tutorial.


### Installation
This workflow requires many different tools, but many to all of them are available on the phoenix cluster. The tutorial below assumes you are using the NYULMC phoneix cluster. The only package that may not be available directly on the cluster is SortMeRNA. Follow the steps below to install a local copy of SortMeRNA on the cluster environment.

#### SortMeRNA


 
#### Step 1. Download and place sequences into a folder.  
  
#### Step 2. FastQC
Download the script below to use for analysis of your sequencing files. You must provide the name of the folder where your sequences are located as an input when using this script.

```
fastqc
``` 

  
#### Step 3. Cutadapt 
Cutadapt is used to remove contamination and bases with very low quality scoring. Typically a phred score of 20 (99%) is optimal for most alignment analysis, but if you are using this data to analyze SNP and variant data, it is best to use a score of 30 (99.9%)

```
cutadapt 
```

#### Step 4. SortMeRNA

```
sortmerna 
```

#### Step 5. STAR-aligner

```
STAR 
```

#### Step 6. SAMStat

```
samstat 
```

#### Step 7. Summarize alignments to gene counts using ```featureCounts```
Now that you have created the alignment files (.bam) for each sample, you can now move ahead to analysis and plotting. The .BAM file will be used to generate read counts per gene and subsequently generate a gene abundance table across all samples. This table will then be used within the statistical analysis program DESeq2.

```
featureCounts 
```



### Analysis Workflow




#### Step 1. Compare groups using DESeq2
A more detailed guide can be found HERE(LINK), but this will cover the important steps for performing the statistical analysis.
```
# Load required libraries
library(DESeq2)
library(fdrtool)

# Import gene count data txt file

# Import sample metadata txt file

# Relevel so we know whats control group
metadatas$Group <- relevel(targets$Group, ref = "CON")

# Make Deseq2 object from featureCounts object -----------
ddsMat <- DESeqDataSetFromMatrix(countData = fc$counts,
                                 colData = metadata,
                                 design = ~Group)
                                 
# Run DESEq2
ddsMat <- DESeq(ddsMat)

# Get results from testing
res_out <- results(ddsMat, pAdjustMethod = "fdr")

# Generate summary of testing
summary(res_out)

```


#### Step 2. Annotate geneID's to gene names
```
# Load libraries
library(biomaRt)

# Get information from ENSEMBL (change mart depending upon species)
geneIDs <- getBM(filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id", 
                                "external_gene_name", 
                                "description", 
                                "entrezgene"), 
                 values = row.names(res_out), 
                 mart = useMart("ensembl", dataset = "itridecemlineatus_gene_ensembl"), 
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




