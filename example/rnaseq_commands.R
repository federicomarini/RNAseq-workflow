# @ Thomas W. Battaglia
# tb1280@nyu.edu

# Install and load required libraries
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("ggplot2") 
biocLite("clusterProfiler") 
biocLite("biomaRt")
biocLite("ReactomePA") 
biocLite("DOSE")
biocLite("KEGG.db")
biocLite("pathview") 
biocLite("org.Hs.eg.db") 
biocLite("pheatmap")
biocLite("genefilter") 
biocLite("RColorBrewer") 
biocLite("topGO") 
biocLite("dplyr")


# - - - - - - - - - - - - - - -
# Import featureCounts data
# - - - - - - - - - - - - - - -
library(DESeq2)
library(ggplot2)
library(knitr)

# - - - - - - - - - - - - - 
# Import gene counts table
# - - - - - - - - - - - - - 

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


# - - - - - - - - - - - - - 
# Import metadata file
# - - - - - - - - - - - - - 

# Import and make row names the matching sampleID's from the countdata
metadata <- read.delim("example/metadata.txt", row.names = 1)

# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)

# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]

# Make sure ID's are correct
head(metadata)


# - - - - - - - - - - - - - 
# Run Deseq2
# - - - - - - - - - - - - - 
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



