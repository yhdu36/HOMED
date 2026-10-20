
suppressWarnings(suppressMessages({
    library(limma)
    library(nlme)
    library(parallel)
    library(doParallel)
    library(foreach)
    library(EpiDISH)
    library(lumi)
}))

source("R/HOMED_Estimate_multilevel.r")
source("R/HOMED_Plot_multilevel.r")
source("R/HOMED_Reference_multilevel.r")

# Run from the HOMED repository root.
# Download input data from:
# https://huggingface.co/datasets/Durianhead/HOMED_data
# Place all RDS files directly in benchmarking/PBMC/.

data_dir <- file.path("benchmarking", "PBMC")

required_files <- c(
    "PBMC_pure_methyl.rds",
    "PBMC_pure_methyl_pdata.rds",
    "PBMC_bulk_methyl_train.rds",
    "PBMC_bulk_pdata_train.rds",
    "DWLS_bulk_prop.rds"
)

missing_files <- required_files[
    !file.exists(file.path(data_dir, required_files))
]

if (length(missing_files) > 0) {
    stop(
        "Missing required PBMC files:\n",
        paste(missing_files, collapse = "\n"),
        "\nDownload input datasets from: ",
        "https://huggingface.co/datasets/Durianhead/HOMED_data"
    )
}

pure_methyl <- readRDS(
    file.path(data_dir, "PBMC_pure_methyl.rds")
)

pure_methyl_pdata <- readRDS(
    file.path(data_dir, "PBMC_pure_methyl_pdata.rds")
)

training_bulk_methyl <- readRDS(
    file.path(data_dir, "PBMC_bulk_methyl_train.rds")
)

bulk_pdata <- readRDS(
    file.path(data_dir, "PBMC_bulk_pdata_train.rds")
)

DWLS_bulk_prop <- readRDS(
    file.path(data_dir, "DWLS_bulk_prop.rds")
)

# Remove cell types with zero proportions across all samples
keep_samples <- colSums(DWLS_bulk_prop == 0) == 0
DWLS_bulk_prop <- DWLS_bulk_prop[, keep_samples, drop = FALSE]

# Names are Expression_id, values are Methyl_id
expr_to_methyl <- setNames(
    bulk_pdata$Methyl_id,
    bulk_pdata$Expression_id
)

# Rename rows from Expression_id to Methyl_id
rownames(DWLS_bulk_prop) <- expr_to_methyl[
    rownames(DWLS_bulk_prop)
]

# Match training samples
stopifnot(
    all(colnames(training_bulk_methyl) %in%
        rownames(DWLS_bulk_prop))
)

DWLS_training_bulk_prop <- DWLS_bulk_prop[
    colnames(training_bulk_methyl), ,
    drop = FALSE
]

# Keep CpGs shared by pure and training methylation matrices
common_cpgs <- intersect(
    rownames(pure_methyl),
    rownames(training_bulk_methyl)
)

pure_methyl <- pure_methyl[
    common_cpgs, ,
    drop = FALSE
]

training_bulk_methyl <- training_bulk_methyl[
    common_cpgs, ,
    drop = FALSE
]

two_level_celltypes <- list(
    mono = NULL,
    t = c("cd4t", "cd8t"),
    nk = NULL,
    b = NULL
)

# Generate HOMED reference
start_wall <- Sys.time()

PBMC_ref_DWLS <- HOMED_Reference(
    FACS_beta = pure_methyl,
    FACS_pd = pure_methyl_pdata,
    training_beta = training_bulk_methyl,
    training_prop = DWLS_training_bulk_prop,
    celltypes = two_level_celltypes
)

end_wall <- Sys.time()

elapsed_wall_sec <- as.numeric(
    difftime(end_wall, start_wall, units = "secs")
)

cat("Elapsed wall time (seconds):", elapsed_wall_sec, "\n")

# Save generated reference
saveRDS(
    PBMC_ref_DWLS,
    file.path(data_dir, "nonlog_PBMC_ref_DWLS.rds")
)
