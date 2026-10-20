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
source("R/HOMED_Reference_multilevel.r")

# Run from the HOMED repository root.
# Download input data from:
# https://huggingface.co/datasets/Durianhead/HOMED_data

data_dir <- file.path("benchmarking", "Placenta")

required_files <- c(
    "Placenta_methyl_pure.rds",
    "Placenta_meta_pure.rds",
    "Placenta_methyl_train.rds",
    "Placenta_methyl_test.rds",
    "Placenta_meta_train.rds",
    "Placenta_meta_test.rds"
)

missing_files <- required_files[
    !file.exists(file.path(data_dir, required_files))
]

if (length(missing_files) > 0) {
    stop(
        "Missing required Placenta datasets:\n",
        paste(missing_files, collapse = "\n"),
        "\n\nDownload from: ",
        "https://huggingface.co/datasets/Durianhead/HOMED_data"
    )
}

# Load Placenta datasets
pure_methyl <- readRDS(
    file.path(data_dir, "Placenta_methyl_pure.rds")
)

pure_metadata <- readRDS(
    file.path(data_dir, "Placenta_meta_pure.rds")
)

mm_train <- readRDS(
    file.path(data_dir, "Placenta_methyl_train.rds")
)

mm_test <- readRDS(
    file.path(data_dir, "Placenta_methyl_test.rds")
)

meta_train <- readRDS(
    file.path(data_dir, "Placenta_meta_train.rds")
)

meta_test <- readRDS(
    file.path(data_dir, "Placenta_meta_test.rds")
)

# Specify cell types
two_level_celltypes <- list(
    Hofbauer = NULL,
    Endothelial = NULL,
    Trophoblast = c(
        "Villoustrophoblast",
        "Syncytiotrophoblast",
        "Extravilloustrophoblast"
    ),
    nRBC = NULL,
    Stromal = NULL
)

# Start timing
start_wall <- Sys.time()
start_cpu <- proc.time()

cat(sprintf(
    "Run started: %s\n",
    format(start_wall, "%Y-%m-%d %H:%M:%S %Z")
))

# Generate HOMED reference
Placenta_ref <- HOMED_Reference(
    FACS_beta = pure_methyl,
    FACS_pd = pure_metadata,
    training_beta = mm_train,
    training_prop = meta_train,
    celltypes = two_level_celltypes
)

# End timing
end_wall <- Sys.time()
end_cpu <- proc.time()

elapsed_wall_sec <- as.numeric(
    difftime(end_wall, start_wall, units = "secs")
)

elapsed_cpu_sec <- (end_cpu - start_cpu)[["elapsed"]]

cat(sprintf(
    "\nTotal wall time: %.1f sec (%.2f min)\n",
    elapsed_wall_sec,
    elapsed_wall_sec / 60
))

cat(sprintf(
    "Total CPU time: %.1f sec (%.2f min)\n",
    elapsed_cpu_sec,
    elapsed_cpu_sec / 60
))

cat(sprintf(
    "Run ended: %s\n",
    format(end_wall, "%Y-%m-%d %H:%M:%S %Z")
))

# Save generated reference
saveRDS(
    Placenta_ref,
    file.path(data_dir, "Placenta_ref_run1.rds")
)

cat("HOMED analysis completed successfully\n")
cat("Results saved to benchmarking/Placenta/Placenta_ref_run1.rds\n")
