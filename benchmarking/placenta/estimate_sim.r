
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
# Input datasets downloaded from:
# https://huggingface.co/datasets/Durianhead/HOMED_data

data_dir <- file.path("benchmarking", "Placenta")

required_files <- c(
    "Sim100_beta_random_newdata.rds",
    "Sim100_prop_random_newdata.rds",
    "Placenta_ref_run1.rds"
)

missing_files <- required_files[
    !file.exists(file.path(data_dir, required_files))
]

if (length(missing_files) > 0) {
    stop(
        "Missing required Placenta files:\n",
        paste(missing_files, collapse = "\n"),
        "\nDownload the required Hugging Face datasets ",
        "and run reference generation first."
    )
}

# Load simulated methylation data and ground-truth proportions
sim_mm_test <- readRDS(
    file.path(data_dir, "Sim100_beta_random_newdata.rds")
)

sim_meta_test <- readRDS(
    file.path(data_dir, "Sim100_prop_random_newdata.rds")
)

rownames(sim_meta_test) <- gsub(
    " ", "", rownames(sim_meta_test)
)

# Load generated HOMED reference
ref_matrix <- readRDS(
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

# HOMED deconvolution on simulated Placenta data
estimated_prop_test <- HOMED_Estimate(
    celltypes = two_level_celltypes,
    deconv_beta = sim_mm_test,
    reference_obj = ref_matrix
)$deconvolution_proportion

# Save estimated proportions
saveRDS(
    estimated_prop_test,
    file.path(data_dir, "sim_estimated_prop_test.rds")
)

cat("Simulated Placenta deconvolution completed successfully\n")
cat("Results saved to benchmarking/Placenta/sim_estimated_prop_test.rds\n")
