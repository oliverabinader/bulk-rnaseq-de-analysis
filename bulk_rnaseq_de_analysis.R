############################################################
# Bulk RNA-seq Differential Expression and Visualization
#
# Author: Oliver Abinader
#
# Description:
# This script performs a bulk RNA-seq workflow for comparing control and treated samples using count-based input data.
# The workflow includes:
#   1. Data preparation and sample selection
#   2. Differential expression analysis using edgeR
#   3. Filtering of significantly differentially expressed genes
#   4. Volcano plot generation
#   5. Correlation analysis across log fold-change comparisons
#
# Notes:
# - Update file paths before running the script.
############################################################

############################
# 1. Load libraries
############################
library(edgeR)
library(ggplot2)
library(ggrepel)

############################
# 2. Read input count data
############################
# Provide full path from where you need to read the data 
raw_counts <- read.delim("data/raw_counts.tsv", check.names = FALSE) # Read count matrix (rows = genes, columns = samples)

raw_counts <- as.data.frame(raw_counts)
rownames(raw_counts) <- raw_counts$EnsemblId
raw_counts$EnsemblId <- NULL

############################
# 3. Subset samples of interest
############################
# Example: Select columns corresponding to a specific cell line (e.g., LS174T)
cell_line <- grep("LS174T", names(raw_counts), value = TRUE)
cell_line_df <- raw_counts[, cell_line]

# Define sample groups based on naming pattern
control_columns <- grep("CD19", names(cell_line_df), value = TRUE)
treated_columns <- grep("5F9", names(cell_line_df), value = TRUE)

# Prepare count matrix for DE analysis
df1 <- cell_line_df[, c(control_columns, treated_columns)]

# Define group labels
group <- factor(c(
  rep("control", length(control_columns)),
  rep("treated", length(treated_columns))
))

############################
# 4. Run differential expression analysis with edgeR
############################
# Create DGEList object (counts + group information)
y <- DGEList(counts = df1, group = group)

# Filter out lowly expressed genes to improve statistical power
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes across samples
y <- calcNormFactors(y)

# Estimate common dispersion across genes
y <- estimateCommonDisp(y, verbose = TRUE)
# Estimate gene-specific (tagwise) dispersion
y <- estimateTagwiseDisp(y, verbose = TRUE)

# Perform differential expression testing 
et <- exactTest(y, pair = c("control", "treated"))
et$table$FDR <- p.adjust(et$table$PValue, method = "fdr")
# Convert results to data frame
et_df <- as.data.frame(et)

############################
# 5. Export DE results
############################
# First, Add gene symbols to DE results
# Read gene annotation file (Ensembl ID → Gene Symbol)
annotation <- read.delim("data/gene_id_to_symbol.tsv", check.names = FALSE)
# Map gene symbols to DE results using Ensembl IDs
et_df$gene <- annotation$GeneSymbol[match(rownames(et_df), annotation$Geneid)]

# Save full DE results table
write.table(
  et_df,
  file = "results/DE_results.tsv",
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

############################
# 6. Filter significant genes
############################
# Make sure the columns are numeric and not character for filtering values (str(object))
# Identify significantly upregulated genes (FDR < 0.05, log2FC >= 1)
significant_up_genes <- et_df[et_df$FDR < 0.05 & et_df$logFC >= 1, ]
# Identify significantly downregulated genes (FDR < 0.05, log2FC < -1)
significant_dn_genes <- et_df[et_df$FDR < 0.05 & et_df$logFC < -1, ]

write.table(
  significant_up_genes,
  file = "results/significant_up_genes.tsv", # provide full path to where you want to create the file 
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  significant_dn_genes,
  file = "results/significant_down_genes.tsv",
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

############################
# 7. Generate volcano plot
############################
de <- et_df
# Initialize category for differential expression
de$diffexpressed <- "NO"
# Assign categories based on thresholds
de$diffexpressed[de$logFC >= 1 & de$FDR < 0.05] <- "UP"
de$diffexpressed[de$logFC < -1 & de$FDR < 0.05] <- "DOWN"

# Define colors for plot
mycolors <- c("DOWN" = "blue", "UP" = "red", "NO" = "grey")

p <- ggplot(de, aes(x = logFC, y = -log10(FDR), col = diffexpressed)) +
  geom_point() +
  geom_text_repel(
    data = de[de$diffexpressed != "NO", ],
    aes(label = gene),
    max.overlaps = 10,
    box.padding = 1,
    point.padding = 0.2,
    size = 4
  ) +
  theme_minimal() +
  scale_colour_manual(values = mycolors, name = NULL) +
  labs(
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

ggsave(
   filename = "results/volcano_plot.png",
   plot = p,
   device = "tiff",
   width = 8,
   height = 6,
   units = "in",
   dpi = 350,
   compression = "lzw"
 )

############################
# 9. Correlation analysis
############################
# Example correlation analysis between two logFC vectors
# Update column names and input object as needed

# correlation <- cor(RNAseq_results$H5_logFC, RNAseq_results$D7_logFC)
# r_squared <- correlation^2
#
# plot(
#   RNAseq_results$H5_logFC,
#   RNAseq_results$D7_logFC,
#   xlab = "H5_logFC",
#   ylab = "D7_logFC",
#   main = "Correlation between H5_logFC and D7_logFC",
#   pch = 16,
#   col = "blue"
# )
#
# abline(lm(RNAseq_results$D7_logFC ~ RNAseq_results$H5_logFC), col = "red")
#
# text(
#   x = min(RNAseq_results$H5_logFC),
#   y = max(RNAseq_results$D7_logFC),
#   labels = paste(
#     "Correlation =", round(correlation, 2),
#     "\nR squared =", round(r_squared, 2)
#   ),
#   pos = 4,
#   col = "black"
# )

############################################################
# End of script
############################################################
