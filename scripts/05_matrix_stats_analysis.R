#Load in the filtered matrix
project_root <- getwd()
filtered_matrix <- read.csv(file.path(project_root, "results", "filtered_matrix.csv"))

# Set feature IDs as rownames
rownames(filtered_matrix) <- filtered_matrix$X
filtered_matrix <- filtered_matrix[, -1]  # remove the X column

# Separate QC and biological sample columns
qc <- grep("QCPool", colnames(filtered_matrix), value = TRUE)
bio <- grep("STDMIX", colnames(filtered_matrix), value = TRUE)

# Extract group labels from biological sample names
groups <- gsub(".*STDMIX_([A-D]1)_.*", "\\1", bio)

# cat("QC samples:", length(qc), "\n")
# cat("Biological samples:", length(bio), "\n")
# cat("Groups:", unique(groups), "\n")

# Extract just the intensity values as a numeric matrix
qc_matrix <- as.matrix(filtered_matrix[, qc])
bio_matrix <- as.matrix(filtered_matrix[, bio])

# Log transform (standard in metabolomics to normalise the data)
qc_log <- log2(qc_matrix + 1)
bio_log <- log2(bio_matrix + 1)

# Quick check
cat("QC matrix dimensions:", dim(qc_log), "\n")
cat("Bio matrix dimensions:", dim(bio_log), "\n")

# Calculate CV for each feature across QC samples
qc_cv <- apply(qc_log, 1, function(x) {
  (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
})

# Summarise
cat("Median CV:", round(median(qc_cv, na.rm = TRUE), 1), "%\n")
cat("Features with CV < 30%:", sum(qc_cv < 30, na.rm = TRUE), "/", length(qc_cv), "\n")

# Plot CV distribution
png(file.path(project_root, "results", "statistical_analysis", "qc_cv_distribution.png"),
    width = 800,
    height = 600)

hist(qc_cv, 
     breaks = 50,
     xlim = c(0, 35),
     ylim = c(0, 200),
     main = "QC Coefficient of Variation Distribution",
     xlab = "CV (%)",
     col = "steelblue")
abline(v = 30, col = "red", lty = 2, lwd = 2)
legend("topright", legend = "30% threshold", col = "red", lty = 2)

dev.off()

# Keep only features with CV < 30%
reliable_features <- qc_cv < 30
reliable_features[is.na(reliable_features)] <- FALSE

# Filter both matrices to reliable features only
qc_log_filtered <- qc_log[reliable_features, ]
bio_log_filtered <- bio_log[reliable_features, ]

# Impute missing values with half the minimum value per feature
impute_min_half <- function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE) / 2
  return(x)
}

bio_log_filtered <- t(apply(bio_log_filtered, 1, impute_min_half))
qc_log_filtered <- t(apply(qc_log_filtered, 1, impute_min_half))

# Check no NAs remain
cat("NAs remaining in bio:", sum(is.na(bio_log_filtered)), "\n")
cat("NAs remaining in QC:", sum(is.na(qc_log_filtered)), "\n")

library(ggplot2)

# Combine QC and biological samples for PCA
all_log_filtered <- cbind(bio_log_filtered, qc_log_filtered)

# PCA requires samples as rows and features as columns so we transpose
pca_result <- prcomp(t(all_log_filtered), scale. = TRUE, center = TRUE)

# Calculate variance explained by each PC
var_explained <- summary(pca_result)$importance[2, ] * 100

#Save the PC scree plot
png(file.path(project_root, "results", "statistical_analysis", "pc_scree_plot.png"),
    width = 800,
    height = 600)

barplot(var_explained[1:20],
        names.arg = paste0("PC", 1:20),
        ylim = c(0, 20),
        xlab = "Principal Component",
        ylab = "Variance Explained (%)",
        main = "Scree Plot",
        col = "steelblue",
        las = 2)

dev.off()

# Create a dataframe for plotting
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  sample = colnames(all_log_filtered),
  group = c(groups, rep("QC", ncol(qc_log_filtered)))
)

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, colour = group, label = sample)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Score Plot",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
  ) +
  theme_bw()

#Save the first pca plot
ggsave(file.path(project_root, "results", "statistical_analysis", "pca_plot1.png"), 
       width = 8, 
       height = 6, 
       dpi = 300)

# Identify and remove outliers
outlier_threshold <- 10
outliers <- rownames(pca_df)[abs(pca_df$PC1) > outlier_threshold]
cat("Outlier samples:", outliers, "\n")
pca_df_clean <- pca_df[!rownames(pca_df) %in% outliers, ]

# Replot without outliers and add convex hulls
hull_data <- do.call(rbind, lapply(unique(pca_df_clean$group), function(g) {
  sub_df <- pca_df_clean[pca_df_clean$group == g, ]
  sub_df[chull(sub_df$PC1, sub_df$PC2), ]
}))

ggplot(pca_df_clean, aes(x = PC1, y = PC2, colour = group, fill = group)) +
  geom_point(size = 3) +
  geom_polygon(data = hull_data, aes(group = group), alpha = 0.1) +
  labs(title = "PCA Score Plot (outliers removed)",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  theme_bw()

#Save the updated pca plot
ggsave(file.path(project_root, "results", "statistical_analysis", "pca_plot2.png"), 
       width = 8, 
       height = 6, 
       dpi = 300)

library(ropls)

# Prepare matrix - samples as rows, features as columns
X <- t(bio_log_filtered)

# Group labels as factor
y <- as.factor(groups)

# Run PLS-DA
plsda_result <- opls(X, y, 
                     predI = 2,
                     fig.pdfC = "none")

# Extract scores
scores <- getScoreMN(plsda_result)

# Build plot dataframe
plsda_df <- data.frame(
  LV1 = scores[, 1],
  LV2 = scores[, 2],
  group = y
)

# Plot
hull_data_plsda <- do.call(rbind, lapply(unique(plsda_df$group), function(g) {
  sub_df <- plsda_df[plsda_df$group == g, ]
  sub_df[chull(sub_df$LV1, sub_df$LV2), ]
}))

ggplot(plsda_df, aes(x = LV1, y = LV2, colour = group, fill = group)) +
  geom_point(size = 3) +
  geom_polygon(data = hull_data_plsda, aes(group = group), alpha = 0.1) +
  labs(title = "PLS-DA Score Plot",
       x = "LV1",
       y = "LV2") +
  theme_bw()

# Save
ggsave(file.path(project_root, "results", "statistical_analysis", "plsda_plot.png"),
       width = 8,
       height = 6,
       dpi = 300)


library(stats)

# Run one-way ANOVA for each feature
p_values <- apply(bio_log_filtered, 1, function(x) {
  df <- data.frame(intensity = x, group = groups)
  fit <- aov(intensity ~ group, data = df)
  summary(fit)[[1]][["Pr(>F)"]][1]
})

# Apply FDR correction
p_adjusted <- p.adjust(p_values, method = "BH")

# Summarise results
cat("Significant features (FDR < 0.05):", sum(p_adjusted < 0.05, na.rm = TRUE), "\n")
cat("Significant features (FDR < 0.01):", sum(p_adjusted < 0.01, na.rm = TRUE), "\n")

#Save ANOVA results
anova_results <- data.frame(
  feature = rownames(bio_log_filtered),
  p_value = p_values,
  p_adjusted = p_adjusted,
  significant = p_adjusted < 0.05
)

write.csv(anova_results, 
          file.path(project_root, "results", "statistical_analysis", "anova_results.csv"),
          row.names = FALSE)

# Calculate mean intensity for A1 vs rest
a1_samples <- groups == "A1"
rest_samples <- groups != "A1"

mean_a1 <- rowMeans(bio_log_filtered[, a1_samples], na.rm = TRUE)
mean_rest <- rowMeans(bio_log_filtered[, rest_samples], na.rm = TRUE)

# Log2 fold change
log2fc <- mean_a1 - mean_rest  # difference in log2 space = log2 fold change

# Build plot dataframe
volcano_df <- data.frame(
  feature = rownames(bio_log_filtered),
  log2fc = log2fc,
  p_adjusted = p_adjusted,
  neg_log10_p = -log10(p_adjusted),
  significant = p_adjusted < 0.05
)

# Plot
ggplot(volcano_df, aes(x = log2fc, y = neg_log10_p, colour = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(title = "Volcano Plot: A1 vs Rest",
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  theme_bw()

library(pheatmap)

# Extract significant features
sig_features <- rownames(bio_log_filtered)[p_adjusted < 0.05]
sig_matrix <- bio_log_filtered[sig_features, ]

# Scale each feature (row) to have mean 0 and sd 1
# This makes patterns visible across features with different intensity scales
sig_matrix_scaled <- t(scale(t(sig_matrix)))

# Create annotation for samples showing which group they belong to
sample_annotation <- data.frame(
  Group = groups,
  row.names = colnames(sig_matrix)
)

# Plot heatmap
pheatmap(sig_matrix_scaled,
         annotation_col = sample_annotation,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Significant Features Heatmap (FDR < 0.05)",
         clustering_method = "ward.D2",
         filename = file.path(project_root, "results", "statistical_analysis", "heatmap.png"),
         width = 8,
         height = 8) 