# Install if not already installed

if (!requireNamespace("illuminaHumanv4.db", quietly = TRUE))
  BiocManager::install("illuminaHumanv4.db")
install.packages("limma")
BiocManager::install("limma")

library(limma)
library(illuminaHumanv4.db)
# Load expression data (assuming it's in a CSV)
data <- read.csv("Data/rawcounts.csv", row.names = 1, header = TRUE)
expr <- data[, seq(1, ncol(data), by = 2)]
# Optionally, add probe IDs if you have them separately
# data$ID_REF <- id_vector  # Or merge with an annotation file

# Extract expression matrix only (remove detection p-values)
expr <- data[, c("SAMPLE.1", "SAMPLE.2", "SAMPLE.3", 
                 "SAMPLE.4", "SAMPLE.5", "SAMPLE.6", 
                 "SAMPLE.7", "SAMPLE.8", "SAMPLE.9")]  # Use actual column names

# Add probe IDs as rownames if available
# rownames(expr) <- data$ID_REF

# Create group factor
group <- factor(c("Control", "Control", "Control", 
                  "ARID1A", "ARID1A", "ARID1A", 
                  "GSK126", "GSK126", "GSK126"))

# Create design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Fit linear model
fit <- lmFit(expr, design)

# Define contrasts of interest
contrast.matrix <- makeContrasts(
  ARID1A_vs_Control = ARID1A - Control,
  GSK126_vs_Control = GSK126 - Control,
  GSK126_vs_ARID1A  = GSK126 - ARID1A,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


results <- topTable(fit2, coef = "ARID1A_vs_Control", number = Inf, adjust.method = "BH")

# Create a data frame with annotation
annot <- select(illuminaHumanv4.db,
                keys = rownames(results),
                columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                keytype = "PROBEID")

# Merge with limma results
results_annotated <- merge(results, annot, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)
rownames(results_annotated) <- results_annotated$Row.names
results_annotated$Row.names <- NULL


# Get full results for all comparisons
res_arid1a <- topTable(fit2, coef = "ARID1A_vs_Control", number = Inf, adjust.method = "BH")
res_arid1a_annot <- merge(res_arid1a, annot, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)
rownames(res_arid1a_annot) <- res_arid1a_annot$Row.names
res_arid1a_annot$Row.names <- NULL
write.csv(res_arid1a_annot, "Data/ARID1A_vs_Control_limma_results.csv")




res_gsk126 <- topTable(fit2, coef = "GSK126_vs_Control", number = Inf, adjust.method = "BH")
res_gsk126_annot <- merge(res_gsk126, annot, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)
rownames(res_gsk126_annot) <- res_gsk126_annot$Row.names
res_gsk126_annot$Row.names <- NULL
write.csv(res_gsk126_annot, "Data/GSK126_vs_Control_limma_results.csv")



res_gsk126_arid1a <- topTable(fit2, coef = "GSK126_vs_ARID1A", number = Inf, adjust.method = "BH")
res_gsk126_arid1a_annot <- merge(res_gsk126_arid1a, annot, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)
rownames(res_gsk126_arid1a_annot) <- res_gsk126_arid1a_annot$Row.names
res_gsk126_arid1a_annot$Row.names <- NULL
write.csv(res_gsk126_arid1a_annot, "Data/GSK126_vs_ARID1A_limma_results.csv")


# Top 50 DEGs for each comparison
top50_ARID1A_vs_Control <- topTable(fit2, coef = "ARID1A_vs_Control", adjust.method = "BH", number = 50)
top50_ARID1A_vs_Control_annot <- merge(top50_ARID1A_vs_Control, annot, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)
rownames(top50_ARID1A_vs_Control_annot) <- top50_ARID1A_vs_Control_annot$Row.names
top50_ARID1A_vs_Control_annot$Row.names <- NULL
write.csv(top50_ARID1A_vs_Control_annot, "Data/ARID1A_vs_Control_top50_results.csv")


top50_GSK126_vs_Control <- topTable(fit2, coef = "GSK126_vs_Control", adjust.method = "BH", number = 50)
top50_GSK126_vs_Control_annot <- merge(top50_GSK126_vs_Control, annot, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)
rownames(top50_GSK126_vs_Control_annot) <- top50_GSK126_vs_Control_annot$Row.names
top50_GSK126_vs_Control_annot$Row.names <- NULL
write.csv(top50_GSK126_vs_Control_annot, "Data/GSK126_vs_Control_top50_results.csv")


top50_GSK126_vs_ARID1A <- topTable(fit2, coef = "GSK126_vs_ARID1A", adjust.method = "BH", number = 50)
top50_GSK126_vs_ARID1A_annot <- merge(top50_GSK126_vs_ARID1A, annot, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)
rownames(top50_GSK126_vs_ARID1A_annot) <- top50_GSK126_vs_ARID1A_annot$Row.names
top50_GSK126_vs_ARID1A_annot$Row.names <- NULL
write.csv(top50_GSK126_vs_ARID1A_annot, "Data/GSK126_vs_ARID1A_top50_results.csv")


annotate_and_save <- function(result, name) {
  result_annot <- merge(result, annot, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)
  rownames(result_annot) <- result_annot$Row.names
  result_annot$Row.names <- NULL
  write.csv(result_annot, paste0("Data/", name, "_annotated.csv"))
}


annotate_and_save(top50_ARID1A_vs_Control_annot, "ARID1A_vs_Control_top50")



