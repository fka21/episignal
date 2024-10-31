# Installing and loading required packages

# Differential Gene Expression Analysis with limma-voom and Wilcoxon Rank-Sum Test
# Author: Ferenc Kagan
# Description: This script performs differential gene expression analysis on RNA-seq data 
# using both limma-voom and Wilcoxon rank-sum test methods, comparing results with a Venn diagram. 
# It visualizes significant genes via a heatmap and performs pathway enrichment analysis using 
# human GO and KEGG databases. The gene IDs are Ensembl IDs.

# Load required libraries (install if necessary)
packages <- c("tidyverse", "ggfortify", "limma", "here", "VennDiagram", "clusterProfiler", "org.Hs.eg.db", "pheatmap")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# Create "output" directory if it doesn't exist
if(!dir.exists(here("output"))) dir.create(here("output"))

# Load data
data <- read.table(here("input", "6_samples_tpm_human.csv"), header = TRUE)

# Define normal sample names
normal_samples <- c("CG.in_S31", "Nav.in_S36", "QB.in_S33")

# Create sample_info data frame
sample_info <- data.frame(
  sampleId = colnames(data)[2:length(colnames(data))],
  status = ifelse(colnames(data)[2:length(colnames(data))] %in% normal_samples, "normal", "inversed")
)

# Data Visualization with Boxplot and PCA
data %>%
  pivot_longer(2:length(.)) %>%
  mutate(value = log(value)) %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  theme_bw() +
  xlab(NULL) +
  ylab('log(TPM + 1)')
# Save boxplot
ggsave(here("output", "boxplot_log_transformed_values.png"))


# Prepare PCA input
pca_input <- t(data[, 2:length(data)])
colnames(pca_input) <- data$Geneid
pca_input <- cbind(data.frame(sample = sample_info$status[match(sample_info$sampleId, rownames(pca_input))], pca_input))
data_pca <- prcomp(pca_input[, 2:length(pca_input)])
autoplot(data_pca, data = pca_input, color = "sample", size = 5) +
theme_bw()
# Save boxplot
ggsave(here("output", "pca.png"))

# Log-transform the data for limma analysis
data_transformed <- log(data[, 2:length(data)] + 1)
rownames(data_transformed) <- data$Geneid

# Prepare design matrix for limma analysis
rownames(sample_info) <- sample_info$sampleId
sample_info <- sample_info[, 2, drop = F]
design <- model.matrix(~ status, data = sample_info)

# Perform limma-voom differential expression analysis
fit <- lmFit(data_transformed, design)
fit <- eBayes(fit, trend = TRUE)
limma_results <- topTable(fit, coef = ncol(design), number = Inf)
limma_sig <- limma_results[limma_results$adj.P.Val < 0.05, ]

# Wilcoxon Rank-Sum Test
wilcoxon_results <- apply(data_transformed, 1, function(x) {
  group1 <- x[sample_info$status == "normal"]
  group2 <- x[sample_info$status == "inversed"]
  wilcox.test(group1, group2, paired = TRUE)$p.value
})

# Adjust p-values and filter significant genes
wilcoxon_adj_p <- p.adjust(wilcoxon_results, method = "BH")
wilcoxon_sig <- names(wilcoxon_adj_p[wilcoxon_adj_p < 0.05])

# Inspect if Wilcoxon Rank-Sum Test has any significant hits
length(na.omit(wilcoxon_sig))

# Heatmap of Significant Genes from Limma-Voom
pheatmap(
  data_transformed[rownames(limma_sig), ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = sample_info,
  show_rownames = FALSE,
  filename=here("output", "heatmap.png"),
  main = "Heatmap of Significant Genes"
)

# Pathway Enrichment Analysis with ClusterProfiler
gene_ids <- rownames(limma_sig)
gene_list <- bitr(gene_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO Enrichment Analysis
go_enrichment <- enrichGO(
  gene = gene_list$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# KEGG Pathway Enrichment Analysis
kegg_enrichment <- enrichKEGG(
  gene = gene_list$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Visualize GO and KEGG enrichment results
dotplot(go_enrichment, showCategory = 20, title = "GO Enrichment")
ggsave(here("output", "GO_dotplot.png"))
dotplot(kegg_enrichment, showCategory = 20, title = "KEGG Enrichment")
ggsave(here("output", "KEGG_dotplot.png"))

