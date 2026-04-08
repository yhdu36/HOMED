#' Validate Cell Types Using Linear/Mixed Models
#'
#' Fits linear or mixed-effects models for each probe to estimate cell type-specific
#' methylation coefficients. Used by pickProbes_limma for 3+ cell types.
#'
#' @param Y Matrix of methylation values (probes x samples).
#' @param pheno Data frame with sample information (design matrix).
#' @param modelFix Formula for fixed effects model.
#' @param modelBatch Optional formula for random effects (batch correction). Default is NULL.
#' @param L.forFstat Optional contrast matrix for F-statistics. Default is NULL.
#' @param verbose Logical indicating whether to print progress. Default is FALSE.
#'
#' @return A list containing:
#'   \item{coefEsts}{Matrix of coefficient estimates (probes x cell types)}
#'   \item{coefVcovs}{List of variance-covariance matrices for each probe}
#'   \item{Pval}{P-values from F-tests}
#'   \item{Fstat}{F-statistics}
#'   And other diagnostic information.
#'
#' @export
validationCellType <- function(Y, pheno, modelFix, modelBatch = NULL, L.forFstat = NULL, verbose = FALSE) {
    
    N <- dim(pheno)[1]
    pheno$y <- rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel <- dim(xTest)[2]
    M <- dim(Y)[1]
    
    if (is.null(L.forFstat)) {
        # All non-intercept coefficients
        L.forFstat <- diag(sizeModel)[-1, ]
        colnames(L.forFstat) <- colnames(xTest)
        rownames(L.forFstat) <- colnames(xTest)[-1]
    }
    
    # Initialize containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()
    
    if (verbose) cat("[validationCellType] ")
    
    # Loop over each probe
    for (j in seq_len(M)) {
        # Remove missing methylation values
        ii <- !is.na(Y[j, ])
        nObserved[j] <- sum(ii)
        pheno$y <- Y[j, ]
        
        if (j %% round(M / 10) == 0 && verbose) cat(".")
        
        # Try to fit model
        try({
            if (!is.null(modelBatch)) {
                fit <- try(
                    lme(modelFix, random = modelBatch, data = pheno[ii, ])
                )
                # If LME can't be fit, use OLS
                OLS <- inherits(fit, "try-error")
            } else {
                OLS <- TRUE
            }
            
            if (OLS) {
                fit <- lm(modelFix, data = pheno[ii, ])
                fitCoef <- fit$coef
                sigmaResid[j] <- summary(fit)$sigma
                sigmaIcept[j] <- 0
                nClusters[j] <- 0
            } else {
                fitCoef <- fit$coef$fixed
                sigmaResid[j] <- fit$sigma
                sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
                nClusters[j] <- length(fit$coef$random[[1]])
            }
            
            coefEsts[j, ] <- fitCoef
            coefVcovs[[j]] <- vcov(fit)
            
            useCoef <- L.forFstat %*% fitCoef
            useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
            Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef)) / sizeModel
        })
    }
    
    if (verbose) cat(" done\n")
    
    # Name the rows for matching to target data
    rownames(coefEsts) <- rownames(Y)
    colnames(coefEsts) <- names(fitCoef)
    degFree <- nObserved - nClusters - sizeModel + 1
    
    # Get P-values from F-statistics
    Pval <- 1 - pf(Fstat, sizeModel, degFree)
    
    list(
        coefEsts = coefEsts,
        coefVcovs = coefVcovs,
        modelFix = modelFix,
        modelBatch = modelBatch,
        sigmaIcept = sigmaIcept,
        sigmaResid = sigmaResid,
        L.forFstat = L.forFstat,
        Pval = Pval,
        orderFstat = order(-Fstat),
        Fstat = Fstat,
        nClusters = nClusters,
        nObserved = nObserved,
        degFree = degFree
    )
}


#' Optimize Probe Selection Using IDOL Algorithm
#'
#' Iteratively refines probe selection to minimize RMSE when predicting cell type
#' proportions. Uses leave-one-out testing and probability weighting to identify
#' the most informative probes.
#'
#' @param candDMRFinderObject Candidate probe coefficient matrix or IDOL object.
#' @param trainingBetas Matrix of methylation beta values from training samples.
#' @param trainingCovariates Matrix of true cell type proportions for training.
#' @param libSize Number of probes to sample each iteration. Default is 500.
#' @param maxIt Maximum number of iterations. Default is 200.
#' @param numCores Number of CPU cores for parallel processing. Default is 4.
#' @param IDOLobj Logical indicating if input is an IDOL object. Default is TRUE.
#' @param rmse_improve_thresh Minimum RMSE improvement to continue iterating.
#'   Default is 0.0001.
#' @param rmse_for Character vector of cell types to optimize for.
#' @param batch_size Batch size for leave-one-out testing. Default is 5.
#'
#' @return A list containing:
#'   \item{IDOL Optimized Library}{Character vector of optimized probe names}
#'   \item{IDOL Optimized CoefEsts}{Matrix of coefficients for optimized probes}
#'   \item{RMSE_CelltypeLevel}{Vector of RMSE values per iteration}
#'   \item{PCC_CelltypeLevel}{Vector of correlation values per iteration}
#'   \item{Number of Iterations}{Total iterations run}
#'   \item{LibrarySize}{Number of probes in library}
#'
#' @export
IDOL_opt_rmse <- function(candDMRFinderObject, trainingBetas, trainingCovariates, libSize = 500, maxIt = 200, numCores = 4, 
                          IDOLobj = TRUE, rmse_improve_thresh = 0.0001, rmse_for = c("Hofbauer"), batch_size = 5) {
    
    # Helper functions
    expit <- function(w) exp(w) / (1 + exp(w))
    
    R2_celltype_level <- function(obs, pred) {
        pcc <- numeric(ncol(obs))
        for (i in 1:ncol(obs)) {
            pcc[i] <- cor(obs[, i], pred[, i], method = "pearson")
        }
        mean(pcc, na.rm = TRUE)
    }
    
    RMSE_celltype_level <- function(obs, pred) {
        cols_use <- intersect(rmse_for, colnames(obs))
        if (length(cols_use) < 1) stop("Missing cell type in the data.")
        obs <- obs[, cols_use, drop = FALSE]
        pred <- pred[, cols_use, drop = FALSE]
        rmse <- numeric(ncol(obs))
        for (i in 1:ncol(obs)) {
            rmse[i] <- sqrt(mean((obs[, i] - pred[, i])^2, na.rm = TRUE))
        }
        mean(rmse, na.rm = TRUE)
    }
    
    polar <- function(x, y, scale = 1) {
        r <- sqrt(x^2 + y^2)
        theta <- atan2(y, x)
        r * cos(theta - (scale * pi / 4))
    }
    
    # Extract probes and coefficients
    if (IDOLobj) {
        trainingProbes1 <- rownames(candDMRFinderObject$coefEsts)
        coefEsts <- candDMRFinderObject$coefEsts
    } else {
        trainingProbes1 <- rownames(candDMRFinderObject)
        coefEsts <- candDMRFinderObject
    }
    
    P <- length(trainingProbes1)
    ProbVector <- rep(1 / P, P)
    V <- libSize
    B <- maxIt
    
    R2_best <- 0
    RMSE_best <- 10000
    
    R2Vals <- numeric(B)
    RMSEVals <- numeric(B)
    
    cellTypes <- colnames(coefEsts)
    K <- length(cellTypes)
    
    if (!all(cellTypes %in% colnames(trainingCovariates))) {
        stop("Cell type names in covariates do not match coefEsts columns")
    }
    
    # Ensure numeric matrices
    trainingBetas <- as.matrix(trainingBetas)
    coefEsts <- as.matrix(coefEsts)
    omega.obs <- as.matrix(trainingCovariates[, cellTypes, drop = FALSE])
    
    # Set up parallel processing
    cl <- makeCluster(numCores)
    
    clusterEvalQ(cl, {
        Sys.setenv(
            OMP_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1"
        )
        library(EpiDISH)
    })
    
    registerDoParallel(cl)
    use_par <- numCores > 1
    
    # Export stable objects to workers
    clusterExport(cl, c("trainingBetas", "omega.obs"), envir = environment())
    
    # Main iteration loop
    for (i in 1:B) {
        # Sample library
        Probes <- sample(1:P, V, prob = ProbVector)
        CpGNames <- trainingProbes1[Probes]
        Beta <- coefEsts[CpGNames, ]
        
        # Evaluate full library
        estF <- suppressMessages(
            epidish(
                trainingBetas[CpGNames, , drop = FALSE],
                Beta,
                method = "CP",
                maxit = 50,
                nu.v = c(0.25, 0.5, 0.75),
                constraint = "inequality"
            )$estF
        )
        omega.tilde <- as.matrix(estF)
        
        R2_ct <- R2_celltype_level(omega.obs, omega.tilde)
        RMSE_ct <- RMSE_celltype_level(omega.obs, omega.tilde)
        
        # Leave-one-out testing
        CpG_idx <- match(CpGNames, rownames(trainingBetas))
        X_full <- trainingBetas[CpG_idx, , drop = FALSE]
        Beta_full <- Beta
        V_iter <- nrow(X_full)
        
        n_batches <- ceiling(V_iter / batch_size)
        chunks <- parallel::splitIndices(n_batches, min(numCores, n_batches))
        
        if (use_par) {
            clusterExport(cl, c("X_full", "Beta_full", "V_iter", "batch_size"), 
                         envir = environment())
        }
        
        Perform.q <- foreach(ch = chunks, .packages = "EpiDISH") %dopar% {
            chunk_len <- sum(pmin(batch_size, V_iter - (ch - 1L) * batch_size))
            R2_vec <- numeric(chunk_len)
            RMSE_vec <- numeric(chunk_len)
            
            write_pos <- 1L
            
            for (batch in ch) {
                start_idx <- (batch - 1L) * batch_size + 1L
                end_idx <- min(batch * batch_size, V_iter)
                
                if (start_idx == 1L) {
                    keep_idx <- seq(end_idx + 1L, V_iter)
                } else if (end_idx == V_iter) {
                    keep_idx <- seq_len(start_idx - 1L)
                } else {
                    keep_idx <- c(seq_len(start_idx - 1L), seq(end_idx + 1L, V_iter))
                }
                
                X.q <- X_full[keep_idx, , drop = FALSE]
                Beta.q <- Beta_full[keep_idx, , drop = FALSE]
                
                estF.q <- suppressMessages(
                    epidish(
                        X.q, Beta.q,
                        method = "CP",
                        maxit = 50,
                        nu.v = c(0.25, 0.5, 0.75),
                        constraint = "inequality"
                    )$estF
                )
                omega.tilde.q <- as.matrix(estF.q)
                
                R2.q <- R2_celltype_level(omega.obs, omega.tilde.q)
                RMSE.q <- RMSE_celltype_level(omega.obs, omega.tilde.q)
                len <- end_idx - start_idx + 1L
                
                idx <- write_pos:(write_pos + len - 1L)
                R2_vec[idx] <- R2.q
                RMSE_vec[idx] <- RMSE.q
                write_pos <- write_pos + len
            }
            
            list(R2 = R2_vec, RMSE = RMSE_vec)
        }
        
        # Flatten results
        R2.q <- unlist(lapply(Perform.q, `[[`, "R2"), use.names = FALSE)
        RMSE.q <- unlist(lapply(Perform.q, `[[`, "RMSE"), use.names = FALSE)
        
        # Update probability vector
        rmse.dq <- (RMSE_ct - RMSE.q) * (-1)
        norm.rmse <- rmse.dq / sd(rmse.dq)
        r2.dq <- R2_ct - R2.q
        norm.r2 <- r2.dq / sd(r2.dq)
        
        p1 <- polar(norm.rmse, norm.r2)
        for (j in 1:length(Probes)) {
            p0 <- ProbVector[[Probes[j]]]
            ProbVector[[Probes[j]]] <- expit(p1[j]) * p0 + p0 / 2
        }
        ProbVector <- ProbVector / sum(ProbVector)
        
        # Track best library
        if ((RMSE_best - RMSE_ct > rmse_improve_thresh)) {
            RMSE_best <- RMSE_ct
            R2_best <- R2_ct
            IDOL.optim.DMRs <- CpGNames
            IDOL.optim.coefEsts <- coefEsts[CpGNames, ]
        }
        
        RMSEVals[i] <- RMSE_ct
        R2Vals[i] <- R2_ct
    }
    
    stopCluster(cl)
    
    IDOLObjects <- list(
        IDOL.optim.DMRs,
        IDOL.optim.coefEsts,
        RMSEVals,
        R2Vals,
        B,
        V
    )
    names(IDOLObjects) <- c(
        "IDOL Optimized Library",
        "IDOL Optimized CoefEsts",
        "RMSE_CelltypeLevel",
        "PCC_CelltypeLevel",
        "Number of Iterations",
        "LibrarySize"
    )
    
    print(paste("Final Library -", rmse_for, ": Avg Cell-RMSE =", round(RMSE_best, 4),
                "; Cell-PCC =", round(R2_best, 4), sep = " "))
    
    return(IDOLObjects)
}

# ---- Helper: normalize hierarchy so old 2-level format still works ----
normalize_celltypes <- function(celltypes) {
    if (is.null(celltypes)) return(NULL)
    if (!is.list(celltypes)) stop("celltypes must be a named list")

    out <- vector("list", length(celltypes))
    names(out) <- names(celltypes)

    for (nm in names(celltypes)) {
        node <- celltypes[[nm]]

        if (is.null(node)) {
            out[[nm]] <- NULL
        } else if (is.list(node)) {
            out[[nm]] <- normalize_celltypes(node)
        } else {
            # old format: parent = c("child1", "child2")
            child_names <- as.character(node)
            child_list <- as.list(rep(list(NULL), length(child_names)))
            names(child_list) <- child_names
            out[[nm]] <- child_list
        }
    }

    out
}

# ---- Helper: depth of hierarchy ----
get_tree_depth <- function(celltypes) {
    celltypes <- normalize_celltypes(celltypes)

    node_depth <- function(node) {
        if (is.null(node)) return(1)
        1 + max(vapply(node, node_depth, numeric(1)))
    }

    max(vapply(celltypes, node_depth, numeric(1)))
}

# ---- Helper: get all node names at a given level ----
get_nodes_at_level <- function(celltypes, target_level, current_level = 1) {
    celltypes <- normalize_celltypes(celltypes)

    if (current_level == target_level) {
        return(names(celltypes))
    }

    out <- character(0)
    for (nm in names(celltypes)) {
        node <- celltypes[[nm]]
        if (!is.null(node)) {
            out <- c(out, get_nodes_at_level(node, target_level, current_level + 1))
        }
    }
    unique(out)
}

# ---- Helper: get subtrees at a given level ----
get_subtrees_at_level <- function(celltypes, target_level, current_level = 1) {
    celltypes <- normalize_celltypes(celltypes)

    if (current_level == target_level) {
        return(celltypes)
    }

    out <- list()
    for (nm in names(celltypes)) {
        node <- celltypes[[nm]]
        if (!is.null(node)) {
            out <- c(out, get_subtrees_at_level(node, target_level, current_level + 1))
        }
    }
    out
}

# ---- Helper: get leaf descendants under one node ----
get_leaf_names <- function(node, node_name = NULL) {
    if (is.null(node)) return(node_name)

    out <- character(0)
    for (nm in names(node)) {
        out <- c(out, get_leaf_names(node[[nm]], nm))
    }
    unique(out)
}

# ---- Helper: map each leaf to the level where it lives ----
get_leaf_levels <- function(celltypes, current_level = 1) {
    celltypes <- normalize_celltypes(celltypes)

    out <- c()
    for (nm in names(celltypes)) {
        node <- celltypes[[nm]]
        if (is.null(node)) {
            out[nm] <- current_level
        } else {
            out <- c(out, get_leaf_levels(node, current_level + 1))
        }
    }
    out
}

# ---- Generic summarizer: aggregate leaf proportions to any level ----
summarize_to_level <- function(cell_types_list, prop_data, target_level) {
    cell_types_list <- normalize_celltypes(cell_types_list)
    subtrees <- get_subtrees_at_level(cell_types_list, target_level)

    summarized_data <- data.frame(
        matrix(
            NA,
            nrow = ncol(prop_data),
            ncol = length(subtrees)
        )
    )
    colnames(summarized_data) <- names(subtrees)
    rownames(summarized_data) <- colnames(prop_data)

    for (node_name in names(subtrees)) {
        subtree <- subtrees[[node_name]]

        if (is.null(subtree)) {
            valid_leaves <- intersect(node_name, rownames(prop_data))
        } else {
            leaf_names <- get_leaf_names(subtree)
            valid_leaves <- intersect(leaf_names, rownames(prop_data))
        }

        if (length(valid_leaves) > 0) {
            summarized_data[, node_name] <- rowSums(
                t(prop_data)[, valid_leaves, drop = FALSE]
            )
        } else {
            summarized_data[, node_name] <- NA
            warning("No valid leaves found for ", node_name)
        }
    }

    summarized_data
}

#' Summarize Level 2 Proportions to Level 1
#'
#' Aggregates subtype proportions to their parent cell type proportions
#' by summing subtypes within each major cell type.
#'
#' @param cell_types_list Named list defining hierarchy (same as celltypes).
#' @param prop_data Matrix of cell type proportions (cell types x samples).
#'
#' @return Data frame of summarized proportions (samples x level 1 cell types).
#'
#' @export
summarize_level1 <- function(cell_types_list, prop_data) {
    summarize_to_level(cell_types_list, prop_data, target_level = 1)
}

#' Select Differentially Methylated Probes and Estimate Reference Coefficients
#'
#' Identifies probes that distinguish between cell types using limma differential
#' methylation analysis, then calculates reference methylation coefficients for
#' deconvolution at a specified hierarchical level.
#'
#' @param p Matrix of methylation beta values (probes x samples) from purified 
#'   cell populations.
#' @param pd Data frame of sample phenotype information. Must include 
#'   'cellType_level1' and 'cellType_level2' columns.
#' @param level Integer (1 or 2) indicating hierarchical level.
#' @param numProbes Number of candidate probes to select per cell type. 
#'   Default is 500.
#' @param p_val Adjusted p-value threshold for probe significance. 
#'   Default is 0.05.
#' @param by Character string for selection method: 'logfc' (default) selects
#'   by log fold-change (equal numbers of hyper/hypo-methylated), otherwise 
#'   selects by smallest adjusted p-values.
#'
#' @return A list with three elements:
#'   \item{coefEsts}{Numeric matrix (probes x cell types) of reference 
#'     coefficients representing expected methylation for each probe in each 
#'     pure cell type. This is used for deconvolution.}
#'   \item{limmaStats}{Named list where each element contains the complete limma
#'     results (logFC, p-values, etc.) for one cell type vs all others.}
#'   \item{selectedProbes}{Character vector of probe names selected across all
#'     cell types.}
#'
#' @details
#' The function performs the following steps:
#' 1. For each cell type, runs limma to compare that cell type vs all others
#' 2. Selects top differentially methylated probes per cell type
#' 3. Takes the union of probes across all cell types
#' 4. Estimates reference coefficients via linear regression (direct matrix
#'    solution for 2 cell types, validationCellType for 3+ cell types)
#'
#' @export
pickProbes_limma <- function(p, pd, level, numProbes = 500, p_val = 0.05, by = 'logfc') {
    
    level_col <- paste0("cellType_level", level)

    if (!level_col %in% colnames(pd)) {
        stop("Missing column in pd: ", level_col)
    }

    pd$CellType <- pd[[level_col]]
    cellTypes <- unique(pd[[level_col]])
    cellTypes <- cellTypes[!is.na(cellTypes)]
    
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    
    trainingProbes <- list()
    limmaStats <- list()
    
    for (celltype in cellTypes) {
        message("Running limma for ", celltype, " vs others...")
        
        pd$status <- ifelse(pd$CellType == celltype, "interest", "other")
        pd$status <- relevel(factor(pd$status), ref = "interest")
        
        design <- model.matrix(~ pd$status)
        
        fit <- lmFit(beta2m(p), design)
        fit <- eBayes(fit)
        
        table_limma <- topTable(fit, coef = 2, n = Inf)
        
        table_limma_complete <- table_limma[is.finite(rowSums(table_limma)), ]
        
        limmaStats[[celltype]] <- table_limma_complete
        
        sig_limma <- table_limma_complete[table_limma_complete$adj.P.Val < p_val, ]
        
        if (by == 'logfc') {
            hyper_limma <- sig_limma[order(sig_limma$logFC, decreasing = TRUE), ]
            hypo_limma <- sig_limma[order(sig_limma$logFC, decreasing = FALSE), ]
            n_per_direction <- floor(numProbes / 2)
            probes <- c(
                rownames(hyper_limma)[seq_len(min(n_per_direction, nrow(hyper_limma)))],
                rownames(hypo_limma)[seq_len(min(n_per_direction, nrow(hypo_limma)))]
            )
        } else {
            sig_limma_sorted <- sig_limma[order(sig_limma$adj.P.Val, decreasing = FALSE), ]
            probes <- rownames(sig_limma_sorted)[seq_len(min(numProbes, nrow(sig_limma_sorted)))]
        }
        
        trainingProbes <- append(trainingProbes, list(probes))
    }
    
    trainingProbes <- unique(unlist(trainingProbes))
    trainingProbes <- trainingProbes[!is.na(trainingProbes)]
    
    message("Selected ", length(trainingProbes), " unique probes across all cell types")
    p <- p[trainingProbes, ]
    
    form <- as.formula(
        sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+"))
    )
    
    phenoDF <- as.data.frame(model.matrix(~ pd$CellType - 1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    
    if (ncol(phenoDF) == 2) {
        message("Estimating coefficients using direct matrix solution (2 cell types)")
        X <- as.matrix(phenoDF)
        coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
    } else {
        message("Estimating coefficients using validationCellType (", 
                ncol(phenoDF), " cell types)")
        tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
        coefEsts <- tmp$coefEsts
    }
    
    message("Reference building complete!")
    
    return(list(
        coefEsts = coefEsts,
        limmaStats = limmaStats,
        selectedProbes = trainingProbes
    ))
}

#' Build Hierarchical Reference Profiles for Cell Type Deconvolution
#'
#' Creates reference methylation profiles for two-level hierarchical 
#' deconvolution using IDOL optimization. Builds references for both major
#' cell types (level 1) and subtypes (level 2).
#'
#' @param FACS_beta Matrix of methylation beta values (probes x samples) from
#'   FACS-sorted purified cell populations.
#' @param FACS_pd Data frame of phenotype information for FACS samples with
#'   'cellType_level1' and 'cellType_level2' columns.
#' @param training_beta Matrix of methylation beta values (probes x samples)
#'   from training samples with known cell type proportions.
#' @param training_prop Matrix of true cell type proportions for training
#'   samples (samples x cell types).
#' @param celltypes Named list defining cell type hierarchy. Each element
#'   represents a level 1 cell type; NULL values indicate no subtypes,
#'   character vectors indicate subtypes. Example:
#'   list(Stromal = c("Fibroblast", "Endothelial"), Hofbauer = NULL)
#' @param p_val Adjusted p-value threshold for probe selection. Default is 0.05.
#' @param numProbes Number of candidate probes to select per cell type via
#'   differential methylation analysis. Default is 500.
#' @param maxIt Vector of length 2 specifying maximum iterations for IDOL
#'   optimization [level1, level2]. Default is c(100, 100).
#' @param libSize Library size for IDOL optimization (number of probes to
#'   sample each iteration). Default is 500.
#' @param rmse_improve_thresh RMSE improvement threshold for IDOL convergence.
#'   Algorithm stops if improvement is less than this value. Default is 0.001.
#' @param batch_size Batch size for parallel leave-one-out testing in IDOL.
#'   Default is 5.
#' @param numCores Number of CPU cores for parallel processing. Default is 4.
#'
#' @return A list with two elements:
#'   \item{layer1_reference_list}{Named list where each element is a coefficient
#'     matrix (probes x cell types) for one level 1 cell type, optimized by IDOL}
#'   \item{layer2_reference_list}{Named list where each element is a coefficient
#'     matrix (probes x cell types) for one level 2 subtype, optimized by IDOL}
#'
#' @details
#' The function performs hierarchical reference building:
#' 
#' Level 1 (Major cell types):
#' 1. Selects candidate probes that distinguish major cell types
#' 2. Aggregates subtype proportions to major cell type level
#' 3. Runs IDOL optimization for each major cell type
#' 
#' Level 2 (Subtypes):
#' 1. Selects candidate probes that distinguish subtypes
#' 2. Runs IDOL optimization for each subtype within parent cell types
#' 
#' IDOL iteratively refines probe selection to minimize RMSE when predicting
#' cell type proportions in the training data.
#'
#' @export
HOMED_Reference <- function(FACS_beta, FACS_pd, training_beta, training_prop, celltypes, p_val = 0.05, numProbes = 500,
                            maxIt = c(100, 100), libSize = 500, rmse_improve_thresh = 0.001, batch_size = 5, numCores = 4) {

    celltypes <- normalize_celltypes(celltypes)
    max_level <- get_tree_depth(celltypes)

    if (length(maxIt) < max_level) {
        maxIt <- c(maxIt, rep(tail(maxIt, 1), max_level - length(maxIt)))
    }

    reference_out <- list()

    for (lvl in seq_len(max_level)) {
        message("=== Building Level ", lvl, " References ===")

        probes_limma_lvl <- pickProbes_limma(
            p = FACS_beta,
            pd = FACS_pd,
            level = lvl,
            p_val = p_val,
            numProbes = numProbes
        )

        summarized_props_lvl <- summarize_to_level(
            cell_types_list = celltypes,
            prop_data = t(training_prop),
            target_level = lvl
        )

        layer_reference_list <- list()

        for (ct in colnames(probes_limma_lvl$coefEsts)) {
            message("Optimizing reference for level ", lvl, " cell type: ", ct)

            idol_res <- IDOL_opt_rmse(
                candDMRFinderObject = probes_limma_lvl$coefEsts,
                trainingBetas = training_beta,
                trainingCovariates = summarized_props_lvl,
                maxIt = maxIt[lvl],
                libSize = libSize,
                IDOLobj = FALSE,
                rmse_improve_thresh = rmse_improve_thresh,
                rmse_for = ct,
                batch_size = batch_size,
                numCores = numCores
            )

            layer_reference_list[[ct]] <- idol_res$`IDOL Optimized CoefEsts`
        }

        reference_out[[paste0("layer", lvl, "_reference_list")]] <- layer_reference_list
    }

    message("\n=== Reference building complete! ===")
    return(reference_out)
}