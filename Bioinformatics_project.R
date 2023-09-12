if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("curatedTCGAData", force = TRUE) 
BiocManager::install("TCGAutils") 
BiocManager::install("TCGAbiolinks") 

install.packages("SNFtool") 
install.packages("caret") 
install.packages("cluster")  
install.packages("mclustcomp") 
install.packages("ggplot2") 
install.packages("mclust") 

library("curatedTCGAData") 
library("TCGAutils") 
library("TCGAbiolinks") 
library("SNFtool") 
library("caret") 
library("cluster") 
library("mclustcomp") 
library("NetPreProc") 
library("ggplot2") 
library("mclust")
library("ggplot2")

# This code follows the steps proposed by the project guidelines
# for this reason (since I chose to do the project individually and not in a group) steps 7, 9 and 10 aren't mentioned

# Step 1) Download PRAD dataset from TCGA

assays <- c("miRNASeqGene", "RNASeq2Gene", "RPPAArray")
mo <- curatedTCGAData(diseaseCode = "PRAD", 
                      assays = assays, 
                      version = "2.0.1", dry.run = FALSE)

mo <- mo[, , paste0("PRAD", "_", assays, "-20160128")]

mo


# Step 2) Pre-processing steps

# 2.1) Consider only primary solid tumors
# Select primary solid tumors based on the barcode pattern
primary <- TCGAutils::TCGAsampleSelect(colnames(mo), sampleCodes = "01")
mo <- mo[, primary, ]

# 2.2) Check for replicates
check_rep <- anyReplicated(mo)
print(check_rep)

# 2.3) Remove FFPE samples:
no_ffpe <- which(as.data.frame(colData(mo))$patient.samples.sample.is_ffpe == "no")
mo <- mo[, no_ffpe, ]

# Obtain samples having all the considered omics:
complete <- intersectColumns(mo)

# Extract assays in list:
complete <- assays(complete) 

# Obtain matrices samples x features:
complete <- lapply(complete, FUN=t)

# 2.4) Remove features having NAs:
complete[[3]] <- complete[[3]][, colSums(is.na(complete[[3]])) == 0] 

# 2.5) Remove features with near zero variance and retain top 100 features having higher variance:
nf <- 100
for (i in 1:length(complete)) {
  idx <- caret::nearZeroVar(complete[[i]])
  message(paste("Removed", length(idx), "features from", names(complete)[i]))
  if (length(idx) != 0) {
    complete[[i]] <- complete[[i]][, -idx]
  }
  
  if (ncol(complete[[i]]) <= nf) next
  
  vars <- apply(complete[[i]], 2, var)
  idx <- sort(vars, index.return = TRUE, decreasing = TRUE)$ix
  
  complete[[i]] <- complete[[i]][, idx[1:nf]]
}

# 2.6) Perform features standardization using z-score:
zscore <- function(data){
  
  zscore_vec <- function(x) { return ((x - mean(x)) / sd(x))}
  data <- apply(data, 2, zscore_vec)
  
  
  return(data)
}

complete <- lapply(complete, zscore) 

# 2.7) Clean barcodes retaining only "Project-TSS-Participant":
for (v in 1:length(complete)) {
  # Access a specific row by its name
  selected_row <- complete[[v]][rownames(complete[[v]]) %in% "Project-TSS-Participant", ]
}


# Step 3) Download disease subtypes from TCGAbiolinks
subtypes <- as.data.frame(TCGAbiolinks::PanCancerAtlas_subtypes())
subtypes <- subtypes[subtypes$cancer.type == "PRAD", ]

# Retain only samples with an associated subtype
subtypes <- subtypes[!is.na(subtypes$Subtype_Integrative), ]
#table(subtypes$Subtype_Integrative)

#head(subtypes)

# Get the patient IDs from the complete dataset
patients_complete <- substr(colnames(complete[[1]]), 1, 12)

# Get the patient IDs from the subtypes dataset
patients_subtypes <- substr(subtypes$pan.samplesID, 1, 12)

# Find the common patients in both datasets
common_patients <- intersect(patients_complete, patients_subtypes)

# Get the indices of the common patients in the complete dataset
idx_complete <- match(common_patients, patients_complete)

# Get the indices of the common patients in the subtypes dataset
idx_subtypes <- match(common_patients, patients_subtypes)

# Reorder the rows of the complete dataset and the subtypes dataset
complete_ordered <- complete
subtypes_ordered <- subtypes[idx_subtypes, ]
for (i in 1:length(complete)) {
  complete_ordered[[i]] <- complete[[i]][, idx_complete]
}


# Step 4) Check if the patients in the multi-omics dataset and the subtypes dataset are in the same order
same_order <- identical(common_patients, substr(subtypes_ordered$pan.samplesID, 1, 12))

# Print the result
if (same_order) {
  message("The patients in the multi-omics dataset and the subtypes dataset are in the same order.")
} else {
  warning("The patients in the multi-omics dataset and the subtypes dataset are not in the same order.")
}

# Step 5) Similarity Network Fusion with the scaled exponential euclidean distance
W_list <- list();
for(i in 1:length(complete)){
  
  # Calculate similarity matrix with k=20
  K = 20 
  # using scaled exponential euclidean distance, to be returned ad a matrix for using it in SNF
  Dist <- as.matrix(dist(scale(as.matrix(complete[[i]])), method="euclidean"))
  Dist <- Dist/max(Dist)
  W = exp(-Dist^2/K)
  
  W_list[[i]] <- W
  
}

# Integration using SNF 
W_int <- SNF(W_list, K=20, t=20)


# Step 6) Integrate the similarity matrices from each data source using an average of the matrices
W_avg <- Reduce("+", W_list) / length(W_list)

# Convert similarity matrices into distance matrices:
dist_int <- 1 - W_int
dist_avg <- 1 - W_avg


# Step 8) Perform disease subtype discovery using PAM algorithm
k <- length(unique(subtypes$Subtype_Selected)) #number of clusters = number of iCluster disease subtypes
pam_res_int <- pam(as.dist(dist_int), k = k)
pam_res_avg <- pam(as.dist(dist_avg), k = k)

# Convert disease subtypes to numeric vectors:
labels <- as.numeric(factor(subtypes$Subtype_Selected, levels = unique(subtypes$Subtype_Selected)))

#labels, pam_res_int$clustering and pam_res_avg$clustering must have the same length
length <- min(length(labels), length(pam_res_int$clustering), length(pam_res_avg$clustering))
pam_res_int_clustering <- pam_res_int$clustering[1:length]
pam_res_avg_clustering <- pam_res_avg$clustering[1:length]

# Compute measures for SNF integration:
types <- c("rand", "adjrand", "nmi1")

# Step 8a) Compute similarity matrices for each data source
similarity_matrices <- list()
for (i in 1:length(complete)) {
  similarity_matrix <- SNFtool::affinityMatrix(1 - NetPreProc::Prob.norm(W_list[[i]]))
  similarity_matrices[[i]] <- similarity_matrix
}

# Perform disease subtype discovery using PAM algorithm on each similarity matrix
pam_results <- list()
for (i in 1:length(complete)) {
  distance_matrix <- as.dist(1 - NetPreProc::Prob.norm(similarity_matrices[[i]]))
  pam_result <- pam(distance_matrix, k = length(unique(subtypes$Subtype_Integrative)))
  pam_results[[i]] <- pam_result
}

# Run PAM clustering on each individual omics similarity matrix and save the cluster assignments 
pam_clusters <- list()
for(i in 1:length(complete)) {
  
  pam_result <- pam_results[[i]]
  
  pam_clusters[[i]] <- pam_result$clustering
  
}

# Compute evaluation metrics comparing individual omics clusterings to known subtypes
source_metrics <- list()
for(i in 1:length(complete)) {
  
  labels_sub <- labels[1:length(pam_clusters[[i]])]
  
  source_metrics[[i]] <- mclustcomp(pam_clusters[[i]],
                                    labels_sub,
                                    types = types)
  
}

# Step 8b) Compute the average of the similarity matrices
average_similarity_matrix <- Reduce("+", similarity_matrices) / length(similarity_matrices)
# Perform disease subtype discovery using PAM algorithm on the average similarity matrix
average_distance_matrix <- as.dist(1 - NetPreProc::Prob.norm(average_similarity_matrix))
average_pam_result <- pam(average_distance_matrix, k = length(unique(subtypes$Subtype_Integrative)))

# Compute measures for average integration
metrics_pam_avg <- mclustcomp(pam_res_avg_clustering, labels, types = types)

# Step 8c) Perform disease subtype discovery using PAM algorithm on the integrated similarity matrix (W_int)
integrated_distance_matrix <- as.dist(1 - NetPreProc::Prob.norm(W_int))
integrated_pam_result <- pam(integrated_distance_matrix, k = length(unique(subtypes$Subtype_Integrative)))

# Compute measures for SNF integration
metrics_pam_int <- mclustcomp(pam_res_int_clustering, labels, types = types)


# Step 11) Compare the clusterings obtained by each considered approach and show the results
# Create a data frame to store the evaluation metrics for each approach
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("curatedTCGAData", force = TRUE) 
BiocManager::install("TCGAutils") 
BiocManager::install("TCGAbiolinks") 

install.packages("SNFtool") 
install.packages("caret") 
install.packages("cluster")  
install.packages("mclustcomp") 
install.packages("ggplot2") 
install.packages("mclust") 

library("curatedTCGAData") 
library("TCGAutils") 
library("TCGAbiolinks") 
library("SNFtool") 
library("caret") 
library("cluster") 
library("mclustcomp") 
library("NetPreProc") 
library("ggplot2") 
library("mclust")
library("ggplot2")


# Step 1) Download PRAD dataset from TCGA

assays <- c("miRNASeqGene", "RNASeq2Gene", "RPPAArray")
mo <- curatedTCGAData(diseaseCode = "PRAD", 
                      assays = assays, 
                      version = "2.0.1", dry.run = FALSE)

mo <- mo[, , paste0("PRAD", "_", assays, "-20160128")]

mo


# Step 2) Pre-processing steps

# 2.1) Consider only primary solid tumors
# Select primary solid tumors based on the barcode pattern
primary <- TCGAutils::TCGAsampleSelect(colnames(mo), sampleCodes = "01")
mo <- mo[, primary, ]

# 2.2) Check for replicates
check_rep <- anyReplicated(mo)
print(check_rep)

# 2.3) Remove FFPE samples:
no_ffpe <- which(as.data.frame(colData(mo))$patient.samples.sample.is_ffpe == "no")
mo <- mo[, no_ffpe, ]

# Obtain samples having all the considered omics:
complete <- intersectColumns(mo)

# Extract assays in list:
complete <- assays(complete) 

# Obtain matrices samples x features:
complete <- lapply(complete, FUN=t)

# 2.4) Remove features having NAs:
complete[[3]] <- complete[[3]][, colSums(is.na(complete[[3]])) == 0] 

# 2.5) Remove features with near zero variance and retain top 100 features having higher variance:
nf <- 100
for (i in 1:length(complete)) {
  idx <- caret::nearZeroVar(complete[[i]])
  message(paste("Removed", length(idx), "features from", names(complete)[i]))
  if (length(idx) != 0) {
    complete[[i]] <- complete[[i]][, -idx]
  }
  
  if (ncol(complete[[i]]) <= nf) next
  
  vars <- apply(complete[[i]], 2, var)
  idx <- sort(vars, index.return = TRUE, decreasing = TRUE)$ix
  
  complete[[i]] <- complete[[i]][, idx[1:nf]]
}

# 2.6) Perform features standardization using z-score:
zscore <- function(data){
  
  zscore_vec <- function(x) { return ((x - mean(x)) / sd(x))}
  data <- apply(data, 2, zscore_vec)
  
  
  return(data)
}

complete <- lapply(complete, zscore) 

# 2.7) Clean barcodes retaining only "Project-TSS-Participant":
for (v in 1:length(complete)) {
  # Access a specific row by its name
  selected_row <- complete[[v]][rownames(complete[[v]]) %in% "Project-TSS-Participant", ]
}


# Step 3) Download disease subtypes from TCGAbiolinks
subtypes <- as.data.frame(TCGAbiolinks::PanCancerAtlas_subtypes())
subtypes <- subtypes[subtypes$cancer.type == "PRAD", ]

# Retain only samples with an associated subtype
subtypes <- subtypes[!is.na(subtypes$Subtype_Integrative), ]
#table(subtypes$Subtype_Integrative)

#head(subtypes)

# Get the patient IDs from the complete dataset
patients_complete <- substr(colnames(complete[[1]]), 1, 12)

# Get the patient IDs from the subtypes dataset
patients_subtypes <- substr(subtypes$pan.samplesID, 1, 12)

# Find the common patients in both datasets
common_patients <- intersect(patients_complete, patients_subtypes)

# Get the indices of the common patients in the complete dataset
idx_complete <- match(common_patients, patients_complete)

# Get the indices of the common patients in the subtypes dataset
idx_subtypes <- match(common_patients, patients_subtypes)

# Reorder the rows of the complete dataset and the subtypes dataset
complete_ordered <- complete
subtypes_ordered <- subtypes[idx_subtypes, ]
for (i in 1:length(complete)) {
  complete_ordered[[i]] <- complete[[i]][, idx_complete]
}


# Step 4) Check if the patients in the multi-omics dataset and the subtypes dataset are in the same order
same_order <- identical(common_patients, substr(subtypes_ordered$pan.samplesID, 1, 12))

# Print the result
if (same_order) {
  message("The patients in the multi-omics dataset and the subtypes dataset are in the same order.")
} else {
  warning("The patients in the multi-omics dataset and the subtypes dataset are not in the same order.")
}

# Step 5) Similarity Network Fusion with the scaled exponential euclidean distance
W_list <- list();
for(i in 1:length(complete)){
  
  # Calculate similarity matrix with k=20
  K = 20 
  # using scaled exponential euclidean distance, to be returned ad a matrix for using it in SNF
  Dist <- as.matrix(dist(scale(as.matrix(complete[[i]])), method="euclidean"))
  Dist <- Dist/max(Dist)
  W = exp(-Dist^2/K)
  
  W_list[[i]] <- W
  
}

# Integration using SNF 
W_int <- SNF(W_list, K=20, t=20)


# Step 6) Integrate the similarity matrices from each data source using an average of the matrices
W_avg <- Reduce("+", W_list) / length(W_list)

# Convert similarity matrices into distance matrices:
dist_int <- 1 - W_int
dist_avg <- 1 - W_avg


# Step 8) Perform disease subtype discovery using PAM algorithm
k <- length(unique(subtypes$Subtype_Selected)) #number of clusters = number of iCluster disease subtypes
pam_res_int <- pam(as.dist(dist_int), k = k)
pam_res_avg <- pam(as.dist(dist_avg), k = k)

# Convert disease subtypes to numeric vectors:
labels <- as.numeric(factor(subtypes$Subtype_Selected, levels = unique(subtypes$Subtype_Selected)))

#labels, pam_res_int$clustering and pam_res_avg$clustering must have the same length
length <- min(length(labels), length(pam_res_int$clustering), length(pam_res_avg$clustering))
pam_res_int_clustering <- pam_res_int$clustering[1:length]
pam_res_avg_clustering <- pam_res_avg$clustering[1:length]

# Compute measures for SNF integration:
types <- c("rand", "adjrand", "nmi1")

# Step 8a) Compute similarity matrices for each data source
similarity_matrices <- list()
for (i in 1:length(complete)) {
  similarity_matrix <- SNFtool::affinityMatrix(1 - NetPreProc::Prob.norm(W_list[[i]]))
  similarity_matrices[[i]] <- similarity_matrix
}

# Perform disease subtype discovery using PAM algorithm on each similarity matrix
pam_results <- list()
for (i in 1:length(complete)) {
  distance_matrix <- as.dist(1 - NetPreProc::Prob.norm(similarity_matrices[[i]]))
  pam_result <- pam(distance_matrix, k = length(unique(subtypes$Subtype_Integrative)))
  pam_results[[i]] <- pam_result
}

# Run PAM clustering on each individual omics similarity matrix and save the cluster assignments 
pam_clusters <- list()
for(i in 1:length(complete)) {
  
  pam_result <- pam_results[[i]]
  
  pam_clusters[[i]] <- pam_result$clustering
  
}

# Compute evaluation metrics comparing individual omics clusterings to known subtypes
source_metrics <- list()
for(i in 1:length(complete)) {
  
  labels_sub <- labels[1:length(pam_clusters[[i]])]
  
  source_metrics[[i]] <- mclustcomp(pam_clusters[[i]],
                                    labels_sub,
                                    types = types)
  
}

# Step 8b) Compute the average of the similarity matrices
average_similarity_matrix <- Reduce("+", similarity_matrices) / length(similarity_matrices)
# Perform disease subtype discovery using PAM algorithm on the average similarity matrix
average_distance_matrix <- as.dist(1 - NetPreProc::Prob.norm(average_similarity_matrix))
average_pam_result <- pam(average_distance_matrix, k = length(unique(subtypes$Subtype_Integrative)))

# Compute measures for average integration
metrics_pam_avg <- mclustcomp(pam_res_avg_clustering, labels, types = types)

# Step 8c) Perform disease subtype discovery using PAM algorithm on the integrated similarity matrix (W_int)
integrated_distance_matrix <- as.dist(1 - NetPreProc::Prob.norm(W_int))
integrated_pam_result <- pam(integrated_distance_matrix, k = length(unique(subtypes$Subtype_Integrative)))

# Compute measures for SNF integration
metrics_pam_int <- mclustcomp(pam_res_int_clustering, labels, types = types)


# Step 11) Compare the clusterings obtained by each considered approach and show the results
# Create a data frame to store the evaluation metrics for each approach
metrics_df <- data.frame(
  Approach = c("SNF integration", "Avg integration", "miRNA", "mRNA", "proteins"),
  AdjRand = c(
    metrics_pam_int$scores[which(metrics_pam_int$types == "adjrand")],
    metrics_pam_avg$scores[which(metrics_pam_avg$types == "adjrand")],
    source_metrics[[1]]$scores[which(source_metrics[[1]]$types == "adjrand")],
    source_metrics[[2]]$scores[which(source_metrics[[2]]$types == "adjrand")],
    source_metrics[[3]]$scores[which(source_metrics[[3]]$types == "adjrand")]
  ),
  NMI = c(
    metrics_pam_int$scores[which(metrics_pam_int$types == "nmi1")],
    metrics_pam_avg$scores[which(metrics_pam_avg$types == "nmi1")],
    source_metrics[[1]]$scores[which(source_metrics[[1]]$types == "nmi1")],
    source_metrics[[2]]$scores[which(source_metrics[[2]]$types == "nmi1")],
    source_metrics[[3]]$scores[which(source_metrics[[3]]$types == "nmi1")]
  ),
  Rand = c(
    metrics_pam_int$scores[which(metrics_pam_int$types == "rand")],
    metrics_pam_avg$scores[which(metrics_pam_avg$types == "rand")],
    source_metrics[[1]]$scores[which(source_metrics[[1]]$types == "rand")],
    source_metrics[[2]]$scores[which(source_metrics[[2]]$types == "rand")],
    source_metrics[[3]]$scores[which(source_metrics[[3]]$types == "rand")]
  )
)

# Print the table
print(metrics_df)

# Convert the data frame to long format
metrics_df_long <- tidyr::gather(metrics_df, "Metric", "Score", -Approach)

# Create the bar plot
ggplot(metrics_df_long, aes(x = Approach, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Approach", y = "Score", fill = "Metric") +
  ggtitle("Comparison of evaluation metrics for each approach")

# Print the table
print(metrics_df)

# Convert the data frame to long format
metrics_df_long <- tidyr::gather(metrics_df, "Metric", "Score", -Approach)

# Create the bar plot
ggplot(metrics_df_long, aes(x = Approach, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Approach", y = "Score", fill = "Metric") +
  ggtitle("Comparison of evaluation metrics for each approach")
