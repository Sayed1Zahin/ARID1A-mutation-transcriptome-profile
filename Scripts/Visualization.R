install.packages("ggfortify")
library(ggplot2)
library(limma)
library(ggrepel)
library(ggplot2)
library(ggrepel)



# Volcano plot ARID1A vs Control
res_arid1a_annot$threshold <- as.factor(res_arid1a_annot$adj.P.Val < 0.05 & abs(res_arid1a_annot$logFC) > 1)

# Filter top 10 significant hits
top_hits_1 <- res_arid1a_annot[res_arid1a_annot$threshold == TRUE, ][1:10, ]

# Plot
ggplot(res_arid1a_annot, aes(x = logFC, y = -log10(P.Value), color = threshold)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_hits_1,
    aes(label = SYMBOL),  # Use SYMBOL column for gene names
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("gray", "red")) +
  labs(
    title = "Volcano Plot: ARID1A vs Control",
    x = "Log2 Fold Change",
    y = "-log10(P-value)"
  ) +
  theme_minimal()



# Volcano Plot GSK126 vs Control

res_gsk126_annot$threshold <- as.factor(res_gsk126_annot$adj.P.Val < 0.05 & abs(res_gsk126_annot$logFC) > 1)
# Filter top 10 significant hits
top_hits_2 <- res_gsk126_annot[res_gsk126_annot$threshold == TRUE, ][1:10, ]

ggplot(res_gsk126_annot, aes(x = logFC, y = -log10(P.Value), color = threshold)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_hits_2,
    aes(label = SYMBOL),  # Use SYMBOL column for gene names
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("gray", "red")) +
  labs(
    title = "Volcano Plot: GSK126 vs Control",
    x = "Log2 Fold Change",
    y = "-log10(P-value)"
  ) +
  theme_minimal()



# Volcano Plot GSSSK126 vs ARID1A 
res_gsk126_arid1a_annot$threshold <- as.factor(res_gsk126_arid1a_annot$adj.P.Val < 0.05 & abs(res_gsk126_arid1a_annot$logFC) > 1)
# Filter top 10 significant hits
top_hits_3 <- res_gsk126_arid1a_annot[res_gsk126_arid1a_annot$threshold == TRUE, ][1:10, ]

ggplot(res_gsk126_arid1a_annot, aes(x = logFC, y = -log10(P.Value), color = threshold)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_hits_3,
    aes(label = SYMBOL),  # Use SYMBOL column for gene names
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("gray", "red")) +
  labs(
    title = "Volcano Plot: GSK126 vs ARID1A",
    x = "Log2 Fold Change",
    y = "-log10(P-value)"
  ) +
  theme_minimal()


ggplot(res_gsk126_arid1a_annot, aes(x = logFC, y = -log10(P.Value), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Volcano Plot: GSK126 vs ARID1A", x = "Log2 Fold Change", y = "-log10(P-value)") +
  theme_minimal()




# 2. PCA Plot
library(ggplot2)
library(ggfortify)  # for autoplot

# Suppose your expression matrix is `exprs_data` (genes as rows, samples as columns)
expr_t <- t(expr)  # transpose for PCA
group <- factor(c(rep("Control", 3), rep("ARID1A", 3), rep("GSK126", 3)))  # adjust accordingly

pca <- prcomp(expr_t, scale. = TRUE)
autoplot(pca, data = data.frame(group), colour = 'group') +
  ggtitle("PCA Plot") + theme_minimal()




# MA plot
limma::plotMA(fit2, coef = "ARID1A_vs_Control", main = "MA Plot: ARID1A vs Control")
limma::plotMA(fit2, coef = "GSK126_vs_Control", main = "MA Plot: GSK126 vs Control")
limma::plotMA(fit2, coef = "GSK126_vs_ARID1A", main = "MA Plot: GSK126_vs_ARID1A")


# histogram
hist(res_arid1a$adj.P.Val, breaks = 50, col = "skyblue",
     main = "Adjusted p-value Distribution", xlab = "Adjusted p-value (FDR)")




# heatmap
library(pheatmap)

# Select top 50 genes by adjusted p-value
top_genes <- rownames(res_arid1a)[order(res_arid1a$adj.P.Val)][1:50]
heat_data <- expr[top_genes, ]

# Optional: Z-score scaling per gene
heat_data_scaled <- t(scale(t(heat_data)))

annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(heat_data)

pheatmap(heat_data_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Heatmap of Top 50 DEGs (ARID1A vs Control)")




