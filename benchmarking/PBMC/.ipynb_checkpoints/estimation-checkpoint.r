
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
# Download input data from:
# https://huggingface.co/datasets/Durianhead/HOMED_data
# All inputs and outputs are stored in benchmarking/PBMC/.

data_dir <- file.path("benchmarking", "PBMC")

required_files <- c(
    "PBMC_bulk_methyl_test.rds",
    "PBMC_bulk_pdata_test.rds",
    "PBMC_sim100_beta.rds",
    "PBMC_sim100_prop.rds",
    "DWLS_bulk_prop.rds",
    "nonlog_PBMC_ref_DWLS.rds"
)

missing_files <- required_files[
    !file.exists(file.path(data_dir, required_files))
]

if (length(missing_files) > 0) {
    stop(
        "Missing required PBMC files:\n",
        paste(missing_files, collapse = "\n"),
        "\nRun reference generation and download ",
        "the required Hugging Face datasets first."
    )
}

# Load testing methylation data and generated reference
testing_bulk_methyl <- readRDS(
    file.path(data_dir, "PBMC_bulk_methyl_test.rds")
)

DWLS_ref_matrix <- readRDS(
    file.path(data_dir, "nonlog_PBMC_ref_DWLS.rds")
)

two_level_celltypes <- list(
    mono = "mono",
    t = c("cd4t", "cd8t"),
    nk = "nk",
    b = "b"
)

# HOMED deconvolution on real PBMC data
PBMC_HOMED_DWLS_estimated_prop <- HOMED_Estimate(
    celltypes = two_level_celltypes,
    deconv_beta = testing_bulk_methyl,
    reference_obj = DWLS_ref_matrix
)$deconvolution_proportion

saveRDS(
    PBMC_HOMED_DWLS_estimated_prop,
    file.path(data_dir, "est_prop.rds")
)

# Load DWLS pseudo-ground-truth proportions
DWLS_testing_bulk_prop <- readRDS(
    file.path(data_dir, "DWLS_bulk_prop.rds")
)

bulk_pdata <- readRDS(
    file.path(data_dir, "PBMC_bulk_pdata_test.rds")
)

# Names are Expression_id, values are Methyl_id
expr_to_methyl <- setNames(
    bulk_pdata$Methyl_id,
    bulk_pdata$Expression_id
)

# Assumes DWLS_bulk_prop.rds has samples as rows
# and cell types as columns.
rownames(DWLS_testing_bulk_prop) <- expr_to_methyl[
    rownames(DWLS_testing_bulk_prop)
]

stopifnot(
    all(rownames(PBMC_HOMED_DWLS_estimated_prop) %in%
        rownames(DWLS_testing_bulk_prop))
)

DWLS_testing_bulk_prop <- DWLS_testing_bulk_prop[
    rownames(PBMC_HOMED_DWLS_estimated_prop), ,
    drop = FALSE
]

meta_test <- DWLS_testing_bulk_prop
est_prop <- PBMC_HOMED_DWLS_estimated_prop

ct_order <- c("mono", "b", "cd8t", "nk", "cd4t")

cts <- ct_order[
    ct_order %in% colnames(meta_test) &
    ct_order %in% colnames(est_prop)
]

ct_names <- c(
    "cd4t" = "CD4+ T cells",
    "b" = "B cells",
    "mono" = "Monocytes",
    "nk" = "NK cells",
    "cd8t" = "CD8+ T cells"
)

# Plot real PBMC performance
png(
    file.path(data_dir, "PBMC_scatter_HOMED.png"),
    width = 5000,
    height = 1000,
    res = 300
)

par(
    mfrow = c(1, 5),
    mar = c(4, 4, 2, 1)
)

ylim_all <- range(
    c(
        as.matrix(meta_test[, cts, drop = FALSE]),
        as.matrix(est_prop[, cts, drop = FALSE])
    ),
    na.rm = TRUE
)

xlim_all <- ylim_all

for (ct in cts) {

    truth <- meta_test[, ct]
    est <- est_prop[, ct]

    pcc <- cor(truth, est)
    rmse <- sqrt(mean((truth - est)^2))

    plot(
        truth,
        est,
        main = paste("HOMED real", ct_names[ct]),
        cex.main = 1.5,
        font.main = 2,
        cex.lab = 1.5,
        cex.axis = 1.3,
        xlab = "True",
        ylab = "Estimated",
        pch = 19,
        xlim = xlim_all,
        ylim = ylim_all,
        asp = 1
    )

    abline(0, 1, col = "red", lwd = 2)

    legend(
        "topleft",
        legend = c(
            paste0("PCC = ", round(pcc, 2)),
            paste0("RMSE = ", round(rmse, 2))
        ),
        bty = "n",
        cex = 1.4
    )
}

dev.off()

# HOMED deconvolution on simulated PBMC data
sim_methyl <- readRDS(
    file.path(data_dir, "PBMC_sim100_beta.rds")
)

simPBMC_HOMED_DWLS_estimated_prop <- HOMED_Estimate(
    celltypes = two_level_celltypes,
    deconv_beta = sim_methyl,
    reference_obj = DWLS_ref_matrix
)$deconvolution_proportion

saveRDS(
    simPBMC_HOMED_DWLS_estimated_prop,
    file.path(data_dir, "sim_est_prop.rds")
)

# Load simulated ground-truth proportions
sim_prop <- readRDS(
    file.path(data_dir, "PBMC_sim100_prop.rds")
)
