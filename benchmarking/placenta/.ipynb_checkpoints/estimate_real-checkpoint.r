
suppressWarnings(suppressMessages({
    library(EpiDISH)
    library(ggplot2)
    library(dplyr)
    library(Seurat)
    library(Matrix)
}))

source("R/HOMED_Estimate_multilevel.r")
source("R/HOMED_Plot_multilevel.r")
source("R/HOMED_Reference_multilevel.r")

# Run from the HOMED repository root.
# Input data downloaded from:
# https://huggingface.co/datasets/Durianhead/HOMED_data

data_dir <- file.path("benchmarking", "Placenta")

required_files <- c(
    "Placenta_methyl_test.rds",
    "Placenta_meta_test.rds",
    "Placenta_ref_run1.rds"
)

missing_files <- required_files[
    !file.exists(file.path(data_dir, required_files))
]

if (length(missing_files) > 0) {
    stop(
        "Missing required Placenta files:\n",
        paste(missing_files, collapse = "\n"),
        "\nRun reference generation and download ",
        "the required Hugging Face datasets first."
    )
}

# Load testing data and generated reference
mm_test <- readRDS(
    file.path(data_dir, "Placenta_methyl_test.rds")
)

meta_test <- readRDS(
    file.path(data_dir, "Placenta_meta_test.rds")
)

ref_matrix_1 <- readRDS(
    file.path(data_dir, "Placenta_ref_run1.rds")
)

two_level_celltypes <- list(
    nRBC = NULL,
    Trophoblast = c(
        "Extravilloustrophoblast",
        "Villoustrophoblast",
        "Syncytiotrophoblast"
    ),
    Endothelial = NULL,
    Hofbauer = NULL,
    Stromal = NULL
)

# HOMED deconvolution on testing Placenta data
estimated_prop_test <- HOMED_Estimate(
    celltypes = two_level_celltypes,
    deconv_beta = mm_test,
    reference_obj = ref_matrix_1
)$deconvolution_proportion

# Save estimated proportions
saveRDS(
    estimated_prop_test,
    file.path(data_dir, "estimated_prop_test.rds")
)

cat("Placenta deconvolution completed successfully\n")
cat("Results saved to benchmarking/Placenta/estimated_prop_test.rds\n")
