# ============================================================================
# EPENDYMOMA RNA-SEQ ANALYSIS PIPELINE - PART 1
# Differential Gene Expression Analysis using DESeq2
# ============================================================================
# Project: End-to-End RNA-Seq Analysis Pipeline: Ependymoma Transcriptomics
# Purpose: Load, preprocess, and perform differential gene expression analysis
# ============================================================================

# ============================================================================
# SECTION 1: SETUP & ENVIRONMENT INITIALIZATION
# ============================================================================
# Load required libraries for RNA-Seq analysis
library(tximeta)              # For reading transcript-level quantification data
library(SummarizedExperiment) # For managing genomic experiment data
library(org.Hs.eg.db)         # Human gene annotation database
library(rtracklayer)          # For reading GTF annotation files
library(dplyr)                # For data manipulation and transformation
library(DESeq2)               # For differential gene expression analysis

# Set working directory to where your data files are located
setwd()  # TODO: Specify your working directory path

# ============================================================================
# SECTION 2: LOAD & VALIDATE SAMPLE DATA
# ============================================================================
# Description: Read sample metadata files for two tumor subtypes (RELA and PFA2)
# These files should contain sample information, file paths, and experimental conditions

# Read RELA subtype sample metadata
coldata_RELA <- read.table("samples_RELA.txt", header=T, sep = "\t")

# Read PFA2 subtype sample metadata
coldata_PFA2 <- read.table("samples_PFA2.txt", header=T, sep = "\t")

# Verify that all data files referenced in the metadata actually exist
# This prevents errors later when trying to read non-existent files
file.exists(coldata_RELA$files)
file.exists(coldata_PFA2$files)

# Display the metadata to verify correct loading
print("RELA Sample Metadata:")
coldata_RELA 

print("PFA2 Sample Metadata:")
coldata_PFA2 

# ============================================================================
# SECTION 3: IMPORT TRANSCRIPT QUANTIFICATION DATA
# ============================================================================
# Description: Import transcript-level abundance estimates using tximeta
# tximeta automatically links quantification data with genomic annotations

# Load PFA2 quantification data (transcript level)
se_PFA2 <- tximeta(coldata_PFA2)

# Display the structure of the SummarizedExperiment object
# This shows dimensions, colnames, metadata, and available assays
se_PFA2

# ============================================================================
# SECTION 4: EXPLORE DATA STRUCTURE & GENOMIC ANNOTATIONS
# ============================================================================
# Description: Inspect the imported data and its associated genomic information

# Display sequence information (chromosome lengths, genome version)
seqinfo(se_PFA2)

# Display sample metadata (conditions, treatments, batch effects)
colData(se_PFA2)

# Display available assay types (counts, TPM, abundance measures)
assayNames(se_PFA2)

# Display genomic ranges for each transcript
rowRanges(se_PFA2)

# ============================================================================
# SECTION 5: SUMMARIZE TRANSCRIPT DATA TO GENE LEVEL
# ============================================================================
# Description: Aggregate transcript-level counts to gene level
# This is necessary for downstream gene-level differential expression analysis
# assignRanges="abundant" uses the most abundant transcript per gene

gse_PFA2 <- summarizeToGene(se_PFA2, assignRanges="abundant")

# Display gene-level SummarizedExperiment object
gse_PFA2

# Display genomic ranges now at gene level
rowRanges(gse_PFA2)

# ============================================================================
# SECTION 6: ADD GENE SYMBOLS & ANNOTATIONS (Method 1: Using org.Hs.eg.db)
# ============================================================================
# Description: Annotate genes with human gene symbols using Bioconductor database
# This adds SYMBOL column to gene annotations for easier interpretation

# Display info about the organism annotation database
org.Hs.eg.db

# Add SYMBOL (gene name) column to gene information
gse_PFA2 <- addIds(gse_PFA2, "SYMBOL", gene=TRUE)

# Display updated genomic ranges with gene symbols
rowRanges(gse_PFA2)

# Check for missing gene symbols (genes without human annotation)
# Output: [1] 25988 (number of genes without annotated symbols)
missing_symbols <- sum(is.na(rowRanges(gse_PFA2)$SYMBOL))
print(paste("Missing gene symbols:", missing_symbols))

# Total number of genes in the dataset
total_genes <- length(rowRanges(gse_PFA2))
print(paste("Total genes:", total_genes))

# ============================================================================
# SECTION 7: ADD GENE ANNOTATIONS (Method 2: Using GTF File)
# ============================================================================
# Description: Read GTF annotation file and extract gene_id, gene_name, gene_type
# This method provides more complete annotations than org.Hs.eg.db alone
# The GTF file contains detailed genomic features and gene metadata

# Read the GENCODE annotation file (v45, primary assembly)
gtf.df <- readGFF("gencode.v45.primary_assembly.annotation.gtf.gz")

# Preview the GTF dataframe structure
head(gtf.df)

# Extract unique gene-level information: gene_id, gene_name, and gene_type
# This creates a mapping dataframe for gene annotation enrichment
geneId2Name.df <- gtf.df %>% 
  dplyr::select(gene_id, gene_name, gene_type) %>% 
  distinct()

# Preview the extracted gene annotations
head(geneId2Name.df)
tail(geneId2Name.df)

# Set gene_id as row names for easy matching with counts data
rownames(geneId2Name.df) <- geneId2Name.df$gene_id

# Verify row names are set correctly
head(geneId2Name.df)

# ============================================================================
# SECTION 8: MATCH ANNOTATIONS TO EXPRESSION DATA
# ============================================================================
# Description: Subset annotation dataframe to match gene order in expression data
# This ensures annotations align correctly with count data

# Subset geneId2Name.df to include only genes present in gse_PFA2
# Maintains the same gene order as the expression data (critical!)
geneId2Name.df <- geneId2Name.df[rownames(gse_PFA2), ]

# Display matched annotations
head(geneId2Name.df)

# ============================================================================
# SECTION 9: ADD MATCHED ANNOTATIONS TO EXPRESSION DATA
# ============================================================================
# Description: Incorporate gene_name and gene_type into the SummarizedExperiment

# Add gene_name column to genomic ranges
rowRanges(gse_PFA2)$gene_name <- geneId2Name.df$gene_name

# Add gene_type column (e.g., protein_coding, lncRNA, miRNA)
rowRanges(gse_PFA2)$gene_type <- geneId2Name.df$gene_type

# Display updated genomic ranges with all annotations
rowRanges(gse_PFA2)

# ============================================================================
# SECTION 10: PREPARE DATA FOR DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================
# Description: Create a DESeq2 object with appropriate experimental design
# The design formula specifies what variables explain the gene expression differences

# Create DESeq2 DataSet object with condition as the experimental design variable
# This will test for expression differences between different conditions
dds_PFA2 <- DESeqDataSet(gse_PFA2, design = ~ condition)

# Display the DESeq2 object structure and sample information
dds_PFA2

# ============================================================================
# SECTION 11: DATA QUALITY CONTROL - VARIANCE STABILIZATION & PCA
# ============================================================================
# Description: Transform data for visualization and identify potential outliers
# Variance stabilizing transformation (VST) reduces the dependence of variance on mean

# Apply variance stabilizing transformation for stable variance across expression levels
vst_PFA2 <- varianceStabilizingTransformation(dds_PFA2)

# Generate PCA plot to visualize sample clustering by condition
# returnData=TRUE returns the PCA coordinates instead of just plotting
# Used to identify outlier samples or batch effects
pca_data <- plotPCA(vst_PFA2, intgroup = "condition", returnData = TRUE)
print("Check the PCA plot to identify outliers and batch effects")
print("Samples that cluster separately from their condition group may be outliers")

# ============================================================================
# SECTION 12: FILTER DATA & PERFORM DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================
# Description: Apply filtering criteria and run DESeq2 statistical test
# Note: dds_PFA2_filtered should be defined (current code has potential bug)

# TODO: Define filtering criteria (e.g., remove low-count genes)
# Example: dds_PFA2_filtered <- dds_PFA2[rowSums(counts(dds_PFA2)) >= 10, ]

# Run DESeq2 analysis: estimates size factors, dispersions, and fits negative binomial GLM
dds_PFA2_filtered <- DESeq(dds_PFA2_filtered)

# Display the analyzed DESeq2 object with statistical results
dds_PFA2_filtered

# ============================================================================
# SECTION 13: EXTRACT DIFFERENTIAL EXPRESSION RESULTS
# ============================================================================
# Description: Get results with specific statistical thresholds
# lfcThreshold=1: log2 fold change threshold (2-fold change)
# alpha=0.05: FDR (False Discovery Rate) significance level

# Extract results with log2 fold change >= 1 and adjusted p-value <= 0.05
FDR_0.05_2x_change_PFA2 <- results(dds_PFA2_filtered, lfcThreshold=1, alpha=0.05)

# Display summary of results (number of up/down regulated genes)
summary(FDR_0.05_2x_change_PFA2)

# Display the results dataframe
FDR_0.05_2x_change_PFA2

# ============================================================================
# SECTION 14: CONVERT RESULTS TO DATAFRAME FORMAT
# ============================================================================
# Description: Convert DESeq2 results object to a standard R dataframe
# This allows for easier data manipulation and export

# Convert results to dataframe
result_PFA2.df <- as.data.frame(FDR_0.05_2x_change_PFA2)

# Display dimensions (rows = genes, columns = statistics)
dim(result_PFA2.df)

# Preview first and last few rows of results
head(result_PFA2.df)
tail(result_PFA2.df)

# ============================================================================
# SECTION 15: EXTRACT & ORGANIZE GENE ANNOTATION INFORMATION
# ============================================================================
# Description: Get gene-level metadata from the DESeq2 object
# Includes: gene_id, SYMBOL (gene symbol), gene_name, and gene_type

# Extract gene annotation columns from DESeq2 object
gene_names_PFA2.df <- as.data.frame(
  rowData(dds_PFA2_filtered)[, c("gene_id", "SYMBOL", "gene_name", "gene_type")]
)

# Display dimensions and preview of gene annotations
dim(gene_names_PFA2.df)
head(gene_names_PFA2.df)
tail(gene_names_PFA2.df)

# ============================================================================
# SECTION 16: MERGE GENE ANNOTATIONS WITH DIFFERENTIAL EXPRESSION RESULTS
# ============================================================================
# Description: Combine gene metadata with statistical test results
# First verify that row order matches between annotations and results

# Verify that gene order is identical between annotations and results
# This is critical - if FALSE, alignment is incorrect
all(rownames(gene_names_PFA2.df) == rownames(result_PFA2.df))  # Should output [1] TRUE

# Combine annotations and results into a single comprehensive dataframe
allresult_PFA2.df <- data.frame(gene_names_PFA2.df, result_PFA2.df)

# Display combined dataframe dimensions and content
dim(allresult_PFA2.df)
head(allresult_PFA2.df)
tail(allresult_PFA2.df)

# ============================================================================
# SECTION 17: QUALITY FILTERING - REMOVE GENES WITH MISSING STATISTICS
# ============================================================================
# Description: Filter out genes that don't have statistical results
# These are typically genes filtered out during DESeq2 analysis (low counts, etc.)

# Identify genes with missing adjusted p-values (NA values)
excluded_genes_PFA2 <- is.na(allresult_PFA2.df$padj)

# Count how many genes have missing values
num_excluded <- sum(excluded_genes_PFA2)
print(paste("Genes excluded (NA p-values):", num_excluded))

# Remove genes with missing adjusted p-values
allresult_PFA2.df <- allresult_PFA2.df[!excluded_genes_PFA2, ]

# Display filtered dataframe dimensions
# Output should be: 27641 genes x 10 columns
dim(allresult_PFA2.df)

# ============================================================================
# SECTION 18: FINAL RESULT FORMATTING & SORTING
# ============================================================================
# Description: Select relevant columns and sort by significance and effect size
# Sorting by padj then log2FoldChange helps identify most important genes

# Select relevant columns for final output
extracted_result_PFA2.df <- allresult_PFA2.df %>%
  dplyr::select(
    gene_id,           # Ensembl gene identifier
    gene_name,         # Gene name from GTF annotation
    gene_type,         # Biotype (protein_coding, lncRNA, etc.)
    baseMean,          # Average expression level across samples
    log2FoldChange,    # Log2 fold change between conditions
    padj               # Adjusted p-value (FDR corrected)
  ) %>%
  arrange(padj, desc(log2FoldChange))  # Sort by significance, then by effect size

# Display final result dimensions
# Output: 27641 genes x 6 columns
dim(extracted_result_PFA2.df)

# Preview top and bottom genes
head(extracted_result_PFA2.df)
tail(extracted_result_PFA2.df)

# ============================================================================
# SECTION 19: EXPORT RESULTS TO FILE
# ============================================================================
# Description: Save differential expression results to tab-separated file
# This file can be used for downstream pathway analysis and visualization

write.table(
  extracted_result_PFA2.df,
  file = "DESeq2_FDR_0.05_2x_genes_PFA2.txt",
  quote = FALSE,          # Don't quote string values
  sep = "\t",             # Tab-separated format
  row.names = FALSE,      # Don't include row numbers
  col.names = TRUE        # Include column headers
)

print("Analysis complete! Results saved to: DESeq2_FDR_0.05_2x_genes_PFA2.txt")

# ============================================================================
# NEXT STEPS (For Part 2 - Suggested)
# ============================================================================
# 1. Repeat analysis for RELA subtype using coldata_RELA
# 2. Compare results between RELA and PFA2 subtypes
# 3. Perform pathway enrichment analysis using ClusterProfiler
# 4. Generate publication-ready visualizations:
#    - Volcano plot (log2FoldChange vs -log10(padj))
#    - Heatmaps of top differentially expressed genes
#    - PCA plots with better formatting
# 5. Create summary report with key findings
# ============================================================================
