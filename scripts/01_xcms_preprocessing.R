#Project Setup

#Establish project_root and data_path
project_root <- getwd()
data_path <- file.path(project_root,"data","raw","raw_data_POS")

#Load packages
library(xcms)
library(MsExperiment)

########################################################

#XCMS preprocessing

#Find mzML files
mzml_files <- list.files(
  path = data_path,
  pattern = "\\.mzML$",
  recursive = TRUE,
  full.names = TRUE
)

#Build metadata
sample_info <- data.frame(
  sample_name = basename(dirname(mzml_files)),
  file_path = mzml_files
)

#Build MS experiment object
ms_data <- readMsExperiment(
  spectraFiles = mzml_files,
  sampleData = sample_info
)

#Define peak detection settings
cwp <- CentWaveParam()

#Run peak detection
add_peaks_data <- findChromPeaks(ms_data, param = cwp)

#Define QC samples
sample_type <- ifelse(grepl("QCPool", sample_info$sample_name), "QC", "Study")
qc_indices <- which(sample_type == "QC")

#Define parameters for initial peak grouping
pdp <- PeakDensityParam(
  sampleGroups = sample_type
)

#Perform initial peak grouping
add_peaks_data <- groupChromPeaks(add_peaks_data, param = pdp)

#Define parameters for alignment
pgp <- PeakGroupsParam(
  minFraction = 0.8,
  subset = qc_indices,
  subsetAdjust = "average",
  smooth = "loess",
  span = 0.4
)

#Align retention times
rt_aligned_data <- adjustRtime(add_peaks_data, param = pgp)

#Correspondence for peak regrouping
rt_aligned_data <- groupChromPeaks(rt_aligned_data, param = pdp)

#Create matrix
feature_matrix <- featureValues(rt_aligned_data)

#Define features separately
feature_info <- featureDefinitions(rt_aligned_data)

########################################################

#Reduce to 1000 features

#First identify features that are most often missing or 0
remove_missing <- !is.na(feature_matrix) & feature_matrix > 0
feature_detection_rate <- rowMeans(remove_missing)

#Remove features with a mean detection rate <=0.5
keep_features <- feature_detection_rate > 0.5

filtered_matrix <- feature_matrix[keep_features, ]
feature_info <- feature_info[keep_features, ]

#Update feature_detection_rate to apply to filtered_matrix
remove_missing <- !is.na(filtered_matrix) & filtered_matrix > 0
feature_detection_rate <- rowMeans(remove_missing)

#Filter further by detection rate and intensity
#Find mean intensity
mean_intensity <- rowMeans(filtered_matrix, na.rm = TRUE)

#Log transform intensity
log_intensity <- log10(mean_intensity + 1)

#Make a combined score
score <- feature_detection_rate * log_intensity

# Sort scores descending
sorted_scores <- sort(score, decreasing = TRUE)

# Plot
png(file.path(project_root, "results", "stats", "feature_selection_score.png"),
    width = 800, height = 600)

plot(sorted_scores,
     type = "l",
     xlab = "Feature rank",
     ylab = "Score (detection rate × log intensity)",
     main = "Feature Selection Score Distribution",
     col = "#1B9AAA",
     lwd = 2)

# Add vertical line at rank 1000
abline(v = 1000, col = "#F4A261", lty = 2, lwd = 2)

# Add threshold score value
threshold_score <- sorted_scores[1000]
abline(h = threshold_score, col = "#F4A261", lty = 3, lwd = 1.5)

legend("topright", 
       legend = c("Feature scores", "Top 1000 cutoff"),
       col = c("#1B9AAA", "#F4A261"),
       lty = c(1, 2), lwd = 2)

dev.off()

#Choose top 1000 features
top_features <- order(score, decreasing = TRUE)[1:min(1000, length(score))]

#Filter matrix further
filtered_matrix <- filtered_matrix[top_features, ]
feature_info <- feature_info[top_features, ]

########################################################

#Prepare output for ipa processing

#Create feature_info_reduced
feature_info_reduced <- data.frame(
  ids = rownames(feature_info),
  mzs = feature_info$mzmed,
  RTs = feature_info$rtmed,
  Int = rowMeans(filtered_matrix, na.rm = TRUE)
)

#Combine filtered_matrix and feature_info into an output ready for ipaPy2
rownames(feature_info_reduced) <- feature_info_reduced$ids
xcms_output <- cbind(feature_info_reduced[, c("ids", "mzs", "RTs")], filtered_matrix)

#Convert to .csv
write.csv(xcms_output,
          file.path("results", "xcms_output.csv"),
          row.names = FALSE)

#Also output filtered matrix for stats analysis
write.csv(filtered_matrix, 
          file.path("results", "filtered_matrix"))

