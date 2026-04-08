#' Calculate Cell Type Proportions Using Reference-Based Deconvolution
#'
#' Uses constrained projection (CP) method via EpiDISH to estimate the
#' proportion of a specific cell type.
#'
#' @param beta Matrix of methylation beta values (probes x samples).
#' @param reference Reference coefficient matrix (probes x cell types).
#' @param celltype Character string specifying which cell type to extract
#'   from the deconvolution results.
#'
#' @return Numeric vector of estimated proportions for the specified cell type
#'   (one value per sample).
#'
#' @export
calculate_cell_proportions <- function(beta, reference, celltype) {

    common_probes <- intersect(rownames(beta), rownames(reference))
    
    cell_prop_res <- suppressMessages(
        epidish(
            beta.m = beta[common_probes, ],
            ref.m = reference[common_probes, ],
            method = "CP"
        )$estF[, celltype]
    )
    
    return(cell_prop_res)
}


#' Merge Hierarchical Deconvolution Results
#'
#' Combines level 1 and level 2 deconvolution results into a single
#' proportion matrix with all cell types, then normalizes to sum to 1.
#'
#' @param layer1_df Data frame of level 1 proportions (samples x major cell types).
#' @param layer2_df Data frame of level 2 proportions (samples x subtypes).
#' @param cell_types Named list defining hierarchy.
#'
#' @return Matrix of final cell type proportions (samples x all cell types),
#'   normalized to sum to 1 per sample.
#'
#' @export
merge_layers <- function(layer1_df, layer2_df, cell_types) {
    
    group_1 <- names(cell_types)[sapply(cell_types, is.null)]
    group_2 <- names(cell_types)[sapply(cell_types, function(x) !is.null(x))]
    
    normalized_level2 <- NULL
    
    for (cell_type in group_2) {
        subtypes <- cell_types[[cell_type]]
        
        level2_sub <- layer2_df[, subtypes, drop = FALSE]
        
        normalized_level2_temp <- level2_sub
        
        if (is.null(normalized_level2)) {
            normalized_level2 <- normalized_level2_temp
        } else {
            normalized_level2 <- cbind(normalized_level2, normalized_level2_temp)
        }
    }
    
    merged_result <- cbind(
        layer1_df[, group_1, drop = FALSE],
        normalized_level2
    )
    
    merged_result <- merged_result / rowSums(merged_result)
    
    return(merged_result)
}

merge_layers_general <- function(layer_results, cell_types) {
    cell_types <- normalize_celltypes(cell_types)
    leaf_levels <- get_leaf_levels(cell_types)

    final_cols <- list()

    for (leaf_name in names(leaf_levels)) {
        lvl <- unname(leaf_levels[[leaf_name]])
        layer_name <- paste0("layer", lvl)

        if (!layer_name %in% names(layer_results)) {
            stop("Missing deconvolution results for ", layer_name)
        }

        if (!leaf_name %in% colnames(layer_results[[layer_name]])) {
            stop("Missing leaf ", leaf_name, " in ", layer_name)
        }

        final_cols[[leaf_name]] <- layer_results[[layer_name]][, leaf_name]
    }

    final_prop <- as.data.frame(final_cols)
    rownames(final_prop) <- rownames(layer_results[[1]])

    final_prop <- as.matrix(final_prop)
    final_prop <- final_prop / rowSums(final_prop)

    return(final_prop)
}
                                        
#' Estimate Cell Type Proportions Using Hierarchical Deconvolution
#'
#' Performs two-level hierarchical deconvolution to estimate cell type
#' proportions in mixed samples using pre-built reference profiles.
#'
#' @param celltypes Named list defining cell type hierarchy (same structure
#'   as used in HOMED_Reference).
#' @param deconv_beta Matrix of methylation beta values (probes x samples)
#'   for samples to be deconvolved.
#' @param reference_obj List containing reference profiles from HOMED_Reference,
#'   with elements 'layer1_reference_list' and 'layer2_reference_list'.
#'
#' @return A list with one element:
#'   \item{deconvolution_proportion}{Matrix of estimated cell type proportions
#'     (samples x cell types). Proportions are normalized to sum to 1 for each
#'     sample.}
#'
#' @details
#' The function performs hierarchical deconvolution in two steps:
#' 
#' Level 1: Estimates proportions of major cell types
#' Level 2: Estimates proportions of subtypes within each major cell type
#' 
#' Final proportions combine both levels, with subtype proportions scaled
#' by their parent cell type proportion (optional, depending on merge_layers
#' implementation).
#'
#' @export
HOMED_Estimate <- function(celltypes, deconv_beta, reference_obj) {

    celltypes <- normalize_celltypes(celltypes)
    max_level <- get_tree_depth(celltypes)

    deconv_by_level <- list()

    for (lvl in seq_len(max_level)) {
        message("=== Level ", lvl, " Deconvolution ===")

        ref_name <- paste0("layer", lvl, "_reference_list")
        if (!ref_name %in% names(reference_obj)) {
            stop("Missing reference object entry: ", ref_name)
        }

        level_nodes <- get_nodes_at_level(celltypes, lvl)
        level_proportion_list <- list()

        for (ct in level_nodes) {
            message("Deconvolving level ", lvl, " cell type: ", ct)

            pred_prop <- calculate_cell_proportions(
                beta = deconv_beta,
                reference = reference_obj[[ref_name]][[ct]],
                celltype = ct
            )

            level_proportion_list[[ct]] <- data.frame(pred_prop = pred_prop)
        }

        wide_df_level <- as.data.frame(
            bind_cols(lapply(level_proportion_list, function(x) x$pred_prop))
        )
        rownames(wide_df_level) <- colnames(deconv_beta)
        colnames(wide_df_level) <- level_nodes

        deconv_by_level[[paste0("layer", lvl)]] <- wide_df_level
    }

    message("\n=== Merging hierarchical results ===")
    final_prop <- merge_layers_general(deconv_by_level, celltypes)

    message("Deconvolution complete!")

    return(list(
        deconvolution_proportion = final_prop,
        deconvolution_by_level = deconv_by_level
    ))
}