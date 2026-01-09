#!/usr/bin/env Rscript
# Script to process kallisto tabular files and generate pre-compiled *_all_values.csv files
# 
# Usage: 
#   Rscript create_precompiled_data.R <input_directory> [output_directory]
#
# Example:
#   Rscript create_precompiled_data.R /path/to/kallisto_files MSI_data
#
# The script will:
#   1. Read all .tabular files from the input directory
#   2. Process them using tximport (same as Custom workflow)
#   3. Run differential expression analysis for each contrast
#   4. Generate *_all_values.csv files (D_O_all_values.csv, L_O_all_values.csv, etc.)
#
# Expected file naming: B1D_S12.tabular, B2L_S13.tabular, G1_S4.tabular, etc.
#   - B = Brown, G = Green
#   - D = Low, L = High, no suffix = Medium
#   - _S## = sample number (ignored)

library(tximport)
library(edgeR)
library(dplyr)
library(tidyr)
library(stringr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript create_precompiled_data.R <input_directory> [output_directory]\n",
       "  input_directory: directory containing kallisto .tabular files\n",
       "  output_directory: directory to save *_all_values.csv files (default: MSI_data)")
}

input_dir <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], "MSI_data")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Reading kallisto files from:", input_dir, "\n")
cat("Output directory:", output_dir, "\n")

# ==============================================================================
# Helper functions (from utils.R)
# ==============================================================================

rename_msi_df_columns <- function(colnames_all, old_name, new_name) {
  colnames_all[colnames_all %in% c(old_name)] <- new_name
  return(colnames_all)
}

calculate_y <- function(filtered_counts, group, normMat_subset = NULL) {
  y <- DGEList(counts = filtered_counts, group = group)
  keep <- rowSums(cpm(y) > 1) >= 2
  y <- y[keep, keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "none")  # tximport norm before
  
  # Apply normalization offset if provided (from full dataset normalization)
  if (!is.null(normMat_subset)) {
    # Subset normMat to match the kept genes
    normMat_subset <- normMat_subset[rownames(y), , drop = FALSE]
    y <- scaleOffset(y, normMat_subset)
  }
  
  return(y)
}

make_fit_y <- function(single_replicates, y, design) {
  if (dim(single_replicates)[1] > 1) {
    y$common.dispersion <- 0.1
    fit <- glmFit(y, design, robust = TRUE)
    warning_message <- "WARNING: too few replicates per treatment/subgroup."
    cat(warning_message, "\n")
  } else {
    y.Disp <- estimateDisp(y, design, robust = TRUE)
    warning_message <- ""
    fit <- glmFit(y.Disp, design, robust = TRUE)
  }
  return(list(y, fit))
}

# ==============================================================================
# Step 1: Read kallisto files using tximport
# ==============================================================================

cat("\n=== Step 1: Reading kallisto files ===\n")

# Get all .tabular files
tabular_files <- list.files(input_dir, pattern = "\\.tabular$", full.names = TRUE)
if (length(tabular_files) == 0) {
  stop("No .tabular files found in ", input_dir)
}

cat("Found", length(tabular_files), "tabular files\n")

# Import using tximport
txi <- tximport(tabular_files, type = "kallisto", ignoreAfterBar = TRUE, txOut = TRUE)

cts <- txi$counts
normMat <- txi$length

# Extract sample names from filenames
sample_names <- basename(tabular_files) %>% 
  str_remove("\\.tabular$")

colnames(cts) <- sample_names
colnames(normMat) <- sample_names

# Normalization (same as server_combine_galaxy_output.R)
normMat <- normMat / exp(rowMeans(log(normMat)))
normCts <- cts / normMat

# Computing effective library sizes
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Create DGEList object
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

# Filtering
keep <- filterByExpr(y)
y <- y[keep, ]

# Store the normalized y object and offset for later use
y_full <- y
normMat_full <- normMat[keep, , drop = FALSE]  # Subset offset matrix to match filtered genes

# Get CPM data for display/filtering
data <- as.data.frame(cpm(y))
data$Geneid <- rownames(data)
rownames(data) <- NULL
data <- data %>% select(Geneid, everything())

cat("Data loaded:", nrow(data), "genes,", ncol(data) - 1, "samples\n")

# ==============================================================================
# Step 2: Process and filter data
# ==============================================================================

cat("\n=== Step 2: Processing and filtering data ===\n")

# Remove duplicates
c <- data[!duplicated(data$Geneid), ]
rownames(c) <- c$Geneid

# Create count matrix from CPM (for filtering purposes)
countData.all <- as.matrix(c[2:ncol(c)])
rownames(countData.all) <- c$Geneid
countData <- countData.all

# Filter: keep genes with >10 counts in at least 1/5 of samples
filter <- apply(countData, 1, function(x) length(x[x > 10]) >= ncol(countData) / 5)
filtered <- countData[filter, ]

# IMPORTANT: Use raw counts from y_full, not CPM, for DE analysis
# The normalization offset is already applied to y_full
# Ensure rownames match between filtered and y_full$counts
common_genes <- intersect(rownames(filtered), rownames(y_full$counts))
if (length(common_genes) < nrow(filtered)) {
  cat("WARNING: Some genes in filtered are not in y_full$counts\n")
  cat("  filtered genes:", nrow(filtered), ", common genes:", length(common_genes), "\n")
}
filtered_counts <- y_full$counts[common_genes, , drop = FALSE]
filtered <- filtered[common_genes, , drop = FALSE]  # Also subset filtered to match

cat("After filtering:", nrow(filtered), "genes\n")

# ==============================================================================
# Step 3: Extract metadata from sample names
# ==============================================================================

cat("\n=== Step 3: Extracting metadata from sample names ===\n")

# Parse sample names: B1D_S12 -> Brown1.Low
# B = Brown, G = Green
# D = Low, L = High, no suffix = Medium/O

parse_sample_name <- function(name) {
  # Remove _S## suffix
  base_name <- str_remove(name, "_S\\d+$")
  
  # Extract group (B1, B2, G1, G2, etc.)
  group_match <- str_extract(base_name, "^([BG]\\d+)")
  group <- ifelse(is.na(group_match), base_name, group_match)
  
  # Extract treatment (D, L, or empty for Medium)
  if (str_detect(base_name, "D$")) {
    treatment <- "Low"
  } else if (str_detect(base_name, "L$")) {
    treatment <- "High"
  } else if (str_detect(base_name, "Dayo$")) {
    # Special case: B2Dayo -> Brown2, Medium
    treatment <- "Medium"
    group <- str_remove(base_name, "Dayo$")
  } else {
    treatment <- "Medium"
  }
  
  # Extract subgroup class (Brown or Green)
  subgroup_class <- ifelse(str_detect(group, "^B"), "Brown", "Green")
  
  # Convert B1 -> Brown1, G1 -> Green1 (for internal use)
  group_label <- ifelse(str_detect(group, "^B"), 
                        str_replace(group, "^B", "Brown"),
                        str_replace(group, "^G", "Green"))
  
  # Create column name for internal processing: Brown1.Low, Green2.Medium, etc.
  colname <- paste0(group_label, ".", treatment)
  
  # Create short format for output: B1.D, G2.L, etc.
  # Convert treatment to short code: Low->D, High->L, Medium->O
  treatment_code <- ifelse(treatment == "Low", "D",
                          ifelse(treatment == "High", "L", "O"))
  colname_short <- paste0(group, ".", treatment_code)
  
  return(list(
    subgroup = group,
    subgroup_class = subgroup_class,
    treatment = treatment,
    colname = colname,
    colname_short = colname_short
  ))
}

# Create metadata dataframe
metadata_list <- lapply(sample_names, parse_sample_name)
metadata_df <- do.call(rbind, lapply(metadata_list, function(x) {
  data.frame(
    sample = x$colname,
    sample_short = x$colname_short,
    subgroup = x$subgroup,
    subgroup_class = x$subgroup_class,
    treatment = x$treatment,
    stringsAsFactors = FALSE
  )
}))

# Rename columns in filtered data (use long format for internal processing)
# Apply to both CPM (filtered) and counts (filtered_counts)
colnames(filtered) <- metadata_df$sample[match(colnames(filtered), sample_names)]
colnames(filtered_counts) <- metadata_df$sample[match(colnames(filtered_counts), sample_names)]

cat("Sample metadata:\n")
print(metadata_df)

# ==============================================================================
# Step 4: Get contrasts
# ==============================================================================

cat("\n=== Step 4: Identifying contrasts ===\n")

subgroup <- factor(metadata_df$subgroup)
subgroup_class <- factor(metadata_df$subgroup_class)
treatment <- factor(metadata_df$treatment)

# Get inter contrasts (between treatments)
if (length(unique(treatment)) > 1) {
  inter_choices <- combn(levels(treatment), 2, simplify = FALSE) %>%
    lapply(function(x) paste0(x[1], "-", x[2])) %>%
    unlist()
} else {
  inter_choices <- c()
}

# Get intra contrasts (between groups within treatment)
intra_choices <- c()
for (treat in unique(metadata_df$treatment)) {
  groups_in_treat <- unique(metadata_df$subgroup_class[metadata_df$treatment == treat])
  if (length(groups_in_treat) >= 2) {
    contrasts <- combn(groups_in_treat, 2, simplify = FALSE) %>%
      lapply(function(x) paste0(x[1], "-", x[2], "_in_", treat))
    intra_choices <- c(intra_choices, unlist(contrasts))
  }
}

# Combine all contrasts
all_contrasts <- c(inter_choices, intra_choices)

# Map to expected condition names
contrast_mapping <- list(
  "Low-Medium" = "D_O",
  "Medium-Low" = "D_O",  # Reverse
  "High-Medium" = "L_O",
  "Medium-High" = "L_O",  # Reverse
  "Low-High" = "D_L",
  "High-Low" = "D_L",  # Reverse
  "Brown-Green_in_Medium" = "B_G_in_O",
  "Green-Brown_in_Medium" = "B_G_in_O",  # Reverse
  "Brown-Green_in_High" = "B_G_in_L",
  "Green-Brown_in_High" = "B_G_in_L",  # Reverse
  "Brown-Green_in_Low" = "B_G_in_D",
  "Green-Brown_in_Low" = "B_G_in_D"  # Reverse
)

cat("Available contrasts:\n")
print(all_contrasts)

# ==============================================================================
# Step 5: Run DE analysis for each contrast
# ==============================================================================

cat("\n=== Step 5: Running DE analysis ===\n")

lfcT <- 1.5
fdrT <- 0.05

# Process each contrast
for (contrast_str in all_contrasts) {
  cat("\nProcessing contrast:", contrast_str, "\n")
  
  # Parse contrast
  if (str_detect(contrast_str, "_in_")) {
    # Intra contrast: Brown-Green_in_Medium
    parts <- str_split(contrast_str, "_in_", simplify = TRUE)
    group_contrast <- parts[1]  # "Brown-Green"
    conditional_treatment <- parts[2]  # "Medium"
    
    group_1 <- str_split(group_contrast, "-", simplify = TRUE)[1]
    group_2 <- str_split(group_contrast, "-", simplify = TRUE)[2]
    
    # Filter columns for this treatment (use counts, not CPM)
    filtered_to_contrast <- filtered_counts[, grep(paste0("\\.", conditional_treatment, "$"), colnames(filtered_counts)), drop = FALSE]
    filtered_to_contrast <- filtered_to_contrast[, c(
      grep(paste0("^", group_1, "\\d*\\."), colnames(filtered_to_contrast), perl = TRUE),
      grep(paste0("^", group_2, "\\d*\\."), colnames(filtered_to_contrast), perl = TRUE)
    ), drop = FALSE]
    
    # Get group_by from subgroup_class
    sst <- metadata_df %>%
      filter(sample %in% colnames(filtered_to_contrast)) %>%
      select(subgroup_class)
    group_by <- factor(sst$subgroup_class)
    
  } else {
    # Inter contrast: Low-Medium
    parts <- str_split(contrast_str, "-", simplify = TRUE)
    group_1 <- parts[1]
    group_2 <- parts[2]
    conditional_treatment <- ""
    
    # Filter columns (use counts, not CPM)
    filtered_to_contrast <- filtered_counts[, c(
      grep(paste0("\\.", group_1, "$"), colnames(filtered_counts)),
      grep(paste0("\\.", group_2, "$"), colnames(filtered_counts))
    ), drop = FALSE]
    
    # Get group_by from treatment
    sst <- metadata_df %>%
      filter(sample %in% colnames(filtered_to_contrast)) %>%
      select(treatment)
    group_by <- factor(sst$treatment)
  }
  
  if (length(unique(group_by)) < 2) {
    cat("  Skipping: insufficient groups\n")
    next
  }
  
  # Create design matrix
  design <- model.matrix(~0 + group_by)
  rownames(design) <- colnames(filtered_to_contrast)
  colnames(design) <- levels(factor(group_by))
  
  # Count replicates per group
  single_replicates <- as.data.frame(design) %>%
    gather(key = "k", value = "v") %>%
    group_by(k) %>%
    summarise(n = sum(v == 1), .groups = "drop") %>%
    filter(n == 1)
  
  # Calculate y and fit
  # Get the subset of normalization offset matrix for this contrast
  normMat_subset <- normMat_full[, colnames(filtered_to_contrast), drop = FALSE]
  y <- calculate_y(filtered_to_contrast, group_by, normMat_subset)
  y_fit_list <- make_fit_y(single_replicates, y, design)
  y <- y_fit_list[[1]]
  fit <- y_fit_list[[2]]
  
  # Create contrast
  cond_clean <- paste0(group_1, "-", group_2)
  cmd <- paste0("my.contrasts <- makeContrasts(", cond_clean, ", levels = design)")
  eval(parse(text = cmd))
  
  # Run GLM likelihood ratio test
  lrt <- glmLRT(fit, contrast = my.contrasts)
  
  # Get results
  v0 <- lrt$table
  v0$FDR <- p.adjust(v0$PValue, method = "BH")
  
  # Add gene labels and change
  deg.edger0 <- v0[abs(v0$logFC) >= lfcT & v0$FDR <= fdrT, ]
  
  if (dim(deg.edger0)[1] == 0) {
    lfcT_tmp <- 0
    fdrT_tmp <- 1
    deg.edger0 <- v0[abs(v0$logFC) >= lfcT_tmp & v0$FDR <= fdrT_tmp, ]
    cat("  WARNING: Not enough DEGs, using relaxed thresholds\n")
  }
  
  deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC) * (-log10(deg.edger0$FDR)), 
                                 decreasing = TRUE)[seq_len(min(150, nrow(deg.edger0)))], ]
  degnames1 <- rownames(deg.edger1)
  degnames1 <- grep("^Gm\\d+|Rik$|^RP\\d+", degnames1, invert = TRUE, value = TRUE)
  
  v0$genelabels <- NA
  if (dim(deg.edger1 %>% filter(!is.na(PValue)))[1] > 0) {
    v0[rownames(v0) %in% degnames1, ]$genelabels <- TRUE
  }
  
  v0$change <- as.factor(ifelse(v0$FDR <= fdrT & abs(v0$logFC) >= lfcT,
                                ifelse(v0$logFC >= lfcT, "Up", "Down"), "NotSig"))
  
  # Prepare output dataframe - keep genelabels and change, keep rownames
  # Column order: logFC, logCPM, LR, PValue, FDR, genelabels, change
  v0 <- v0 %>%
    select(logFC, logCPM, LR, PValue, FDR, genelabels, change, everything())
  
  # Get count data from the full filtered dataset (all samples)
  # Use the original y_full object which has all samples with proper normalization
  counts <- as.data.frame(y_full$counts)
  
  # Map long column names to short format for output (B1.D, G2.L, etc.)
  colname_mapping <- metadata_df %>%
    filter(sample %in% colnames(counts)) %>%
    select(sample, sample_short)
  
  # Rename count columns to short format
  for (i in seq_len(nrow(colname_mapping))) {
    old_name <- colname_mapping$sample[i]
    new_name <- colname_mapping$sample_short[i]
    if (old_name %in% colnames(counts)) {
      colnames(counts)[colnames(counts) == old_name] <- new_name
    }
  }
  
  # Select only relevant columns for this contrast (matching MSI script logic)
  # For inter contrasts: select by treatment suffix (e.g., .D and .O for D_O)
  # For intra contrasts: select by group prefix and treatment suffix (e.g., B.*.D and G.*.D for B_G_in_D)
  if (conditional_treatment != "") {
    # Intra contrast: select by group prefix (B or G) and treatment suffix
    # e.g., for B_G_in_D: select B1.D, B2.D, ... G1.D, G2.D, ... (only .D columns)
    # Convert "Brown" -> "B", "Green" -> "G" for matching short format column names
    group_1_short <- ifelse(group_1 == "Brown", "B", "G")
    group_2_short <- ifelse(group_2 == "Brown", "B", "G")
    treatment_code <- ifelse(conditional_treatment == "Low", "D",
                            ifelse(conditional_treatment == "High", "L", "O"))
    selectedCols <- c(
      grep(paste0("^", group_1_short, "[1-9,a-z,A-Z]*\\.", treatment_code, "$"), colnames(counts), perl = TRUE),
      grep(paste0("^", group_2_short, "[1-9,a-z,A-Z]*\\.", treatment_code, "$"), colnames(counts), perl = TRUE)
    )
  } else {
    # Inter contrast: select by treatment suffix
    # e.g., for D_O: select all .D and .O columns
    treatment_code_1 <- ifelse(group_1 == "Low", "D",
                              ifelse(group_1 == "High", "L", "O"))
    treatment_code_2 <- ifelse(group_2 == "Low", "D",
                              ifelse(group_2 == "High", "L", "O"))
    selectedCols <- c(
      grep(paste0("\\.", treatment_code_1, "$"), colnames(counts)),
      grep(paste0("\\.", treatment_code_2, "$"), colnames(counts))
    )
  }
  
  # Subset counts to only selected columns
  counts_selected <- counts[, selectedCols, drop = FALSE]
  
  # Combine v0 (with rownames) and selected counts
  # Use cbind to preserve rownames (matching MSI script format)
  v0 <- cbind(v0, counts_selected[rownames(v0), , drop = FALSE])
  
  # Map to condition name
  # Build the mapping key
  if (conditional_treatment != "") {
    mapping_key <- paste0(group_1, "-", group_2, "_in_", conditional_treatment)
  } else {
    mapping_key <- paste0(group_1, "-", group_2)
  }
  
  cond_name <- contrast_mapping[[mapping_key]]
  if (is.null(cond_name)) {
    # Try reverse
    if (conditional_treatment != "") {
      reverse_key <- paste0(group_2, "-", group_1, "_in_", conditional_treatment)
    } else {
      reverse_key <- paste0(group_2, "-", group_1)
    }
    cond_name <- contrast_mapping[[reverse_key]]
  }
  
  if (is.null(cond_name)) {
    cat("  WARNING: Could not map contrast '", mapping_key, "' to condition name, skipping\n", sep = "")
    next
  }
  
  cat("  Mapped to condition:", cond_name, "\n")
  
  # Save to file
  # Format: rownames as first column (empty header), then logFC, logCPM, LR, PValue, FDR, genelabels, change, then sample counts
  output_file <- file.path(output_dir, paste0(cond_name, "_all_values.csv"))
  write.table(v0, output_file, sep = "\t", row.names = TRUE, quote = FALSE, col.names = NA)
  cat("  Saved to:", output_file, "\n")
  cat("  Genes:", nrow(v0), "\n")
}

cat("\n=== Done! ===\n")
cat("Output files saved to:", output_dir, "\n")

