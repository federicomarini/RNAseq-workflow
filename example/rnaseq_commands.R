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
biocLite("org.Mm.eg.db")
biocLite("org.Hs.eg.db")
biocLite("pheatmap")
biocLite("genefilter")
biocLite("RColorBrewer")
biocLite("GO.db")
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

# Get results from testing with FDR adjust pvalues
res <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(res)

# Check directionality of fold change
mcols(res, use.names = T)


# - - - - - - - - - - - - - 
# Gene annotation
# - - - - - - - - - - - - - 
library(AnnotationDbi)
library(org.Mm.eg.db)

# Add gene full name
res$description <- mapIds(x = org.Mm.eg.db,
                          keys = row.names(res),
                          column = "GENENAME",
                          keytype = "SYMBOL",
                          multiVals = "first")

# Add gene symbol
res$symbol <- row.names(res)

# Add ENTREZ ID
res$entrez <- mapIds(x = org.Mm.eg.db,
                     keys = row.names(res),
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

# Subset for only significant genes (q < 0.05)
res_sig <- subset(res, padj < 0.05)
head(res_sig)

# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ddsMat), normalized = T), 
            file = 'DESeq2_normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(ddsMat[row.names(res_sig)], normalized = T), 
            file = 'DESeq2_normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(res), 
            file = "DESEq2_results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(res_sig), 
            file = "DESEq2_results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)


# - - - - - - - - - - - - - 
# PCA plot
# - - - - - - - - - - - - - 
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Plot PCA by column variable
plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 500) +
  theme_bw() +
  ggsave('figures/pca_plot.png')


# - - - - - - - - - - - - - 
# Heatmap plot
# - - - - - - - - - - - - - 
# Load libraries
library(pheatmap) 
library(RColorBrewer) 

# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Gather top 30 genes and make matrix
mat <- assay(ddsMat_rlog[row.names(res_sig)])[1:30, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group),
  Replicate = factor(colData(ddsMat_rlog)$Replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c("LoGlu" = "blue", "HiGlu" = "orange"),
  Replicate = c(Rep1 = "red", Rep2 = "green")
)

# Make Heatmap with pheatmap function.
# See more in documentation for customization
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors, 
         show_colnames = F)


# - - - - - - - - - - - - - 
# Volcano plot
# - - - - - - - - - - - - - 
# Load libraries
library(ggplot2)
library(RColorBrewer)

# Gather Log-fold change and FDR-corrected pvalues from deseq2 results
data <- data.frame(pval = -log10(res$padj), 
                   lfc = res$log2FoldChange, 
                   row.names = row.names(res))

# Remove any rows that have NA as an entry
data <- na.omit(data)

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval))

# Add ggplot2 layers
vol +   
  ggtitle(label = "Volcano Plot") +
  geom_point(aes(colour = pval), size = 5, alpha = 0.7, na.rm = T) + # color the dots
  theme_bw(base_size = 18) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("Treatment" / "Control"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_vline(xintercept = c(-2, 2), colour = "darkgrey") + # Add fold change cutoff lines
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_colour_gradient(low = "red", high = "green") # Add red-green scaling by intensity


# - - - - - - - - - - - - - 
# MA plot
# - - - - - - - - - - - - - 
plotMA(res, ylim = c(-5, 5))


# - - - - - - - - - - - - - 
# Single gene plot
# - - - - - - - - - - - - - 
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Get gene with highest expression
top_gene <- rownames(res)[which.min(res$log2FoldChange)]

# Plot single gene
plotCounts(ddsMat, gene = top_gene, intgroup = "Group")


# - - - - - - - - - - - - - - -
# Pathway analysis of DE genes
# - - - - - - - - - - - - - - -

# Load required libraries
library(clusterProfiler)
library(ReactomePA)
library(KEGG.db)
library(DOSE)
library(org.Mm.eg.db)
library(pathview)

# Remove any genes that do not have any entrez identifiers
res_sig_entrez <- subset(res_sig, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
gene_matrix <- res_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- res_sig_entrez$entrez


# - - - - - - - - - - - - -
# Enrich with KEGG database
# - - - - - - - - - - - - -
kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'mouse',
                          pvalueCutoff = 0.05, 
                          readable = TRUE)

# Get table of results
head(summary(kegg_enrich))

# KEGG plot
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)


# - - - - - - - - - - - - -
# Enrich with ReactomeDB
# - - - - - - - - - - - - -
reactome_enrich <- enrichPathway(gene = names(gene_matrix), 
                                 organism = 'mouse', 
                                 readable = TRUE,
                                 pvalueCutoff = 0.05)

# Get table of results
head(summary(reactome_enrich))

# Plot results
barplot(reactome_enrich, 
        drop = TRUE, 
        showCategory=10, 
        title = "ReactomeDB Enrichment Pathways",
        font.size = 8)


# - - - - - - - - - - - - -
# Enrich with GO
# - - - - - - - - - - - - -
go_enrich <- enrichGO(gene = names(gene_matrix), 
                      organism = 'mouse', 
                      readable = TRUE,
                      ont = "BP",
                      pvalueCutoff = 0.05)

# Get table of results
head(summary(go_enrich))

# Plot results
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)


# - - - - - - - - - - - - -
# Plot specific KEGG pathways
# (with fold change) 
# - - - - - - - - - - - - -
# pathway.id : KEGG pathway identifier
pathview(gene.data = gene_matrix, 
         pathway.id = "04070", 
         species = "mouse", 
         map.symbol = T)

