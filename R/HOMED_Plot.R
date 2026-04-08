#' Visualize Deconvolution Results with Correlation Metrics
#'
#' Creates scatter plots comparing predicted vs true cell type proportions
#' with Pearson correlation (PCC) and RMSE metrics displayed for each cell type.
#'
#' @param beta Data frame of methylation beta values (probes x samples).
#' @param ref Reference object containing coefficient estimates. Can be from
#'   IDOL optimization or pickProbes_limma.
#' @param true_prop Data frame of true cell type proportions (cell types x samples).
#' @param title Character string for plot title. Default is 'Deconvolution Results'.
#' @param idol Logical indicating if ref is an IDOL object (TRUE) or regular
#'   coefficient matrix (FALSE). Default is TRUE.
#' @param method Deconvolution method to use: 'RPC', 'CBS', or 'CP'. 
#'   Default is 'RPC'.
#' @param con Constraint type: 'equality' or 'inequality'. Default is 'equality'.
#' @param pred_prop Optional pre-computed predicted proportions. If provided,
#'   skips deconvolution step. Default is NULL.
#' @param xylim Axis limits for plots. Default is 0.5.
#'
#' @return A list with four elements:
#'   \item{cell_type_pcc_df}{Data frame with Pearson correlation per cell type}
#'   \item{cell_type_rmse_df}{Data frame with RMSE per cell type}
#'   \item{plot}{ggplot2 object with faceted scatter plots}
#'   \item{predictions}{Matrix of predicted proportions}
#'
#' @export
HOMED_plot <- function(
    beta,
    ref,
    true_prop,
    title = 'Deconvolution Results',
    idol = TRUE,
    method = 'RPC',
    con = 'equality',
    pred_prop = NULL,
    xylim = 0.5
) {
    
    y_true <- as.data.frame(true_prop)
    
    if (is.null(pred_prop)) {
        if (idol) {
            res_epidish <- epidish(
                beta.m = beta,
                ref.m = ref$`IDOL Optimized CoefEsts`,
                method = method,
                maxit = 50,
                nu.v = c(0.25, 0.5, 0.75),
                constraint = con
            )
        } else {
            res_epidish <- epidish(
                beta.m = beta,
                ref.m = ref$coefEsts,
                method = method,
                maxit = 50,
                nu.v = c(0.25, 0.5, 0.75),
                constraint = con
            )
        }
        y_pred <- as.data.frame(res_epidish$estF)[, colnames(y_true)]
    } else {
        common_ct <- intersect(colnames(y_true), colnames(pred_prop))
        y_pred <- pred_prop[, common_ct]
        y_true <- y_true[, common_ct]
    }
    
    y_true$sample <- rownames(y_true)
    y_pred$sample <- rownames(y_pred)
    
    y_true_long <- reshape(
        y_true,
        varying = names(y_true)[-ncol(y_true)],
        v.names = "y_true",
        timevar = "cell_type",
        times = names(y_true)[-ncol(y_true)],
        direction = "long"
    )
    
    y_pred_long <- reshape(
        y_pred,
        varying = names(y_pred)[-ncol(y_pred)],
        v.names = "y_pred",
        timevar = "cell_type",
        times = names(y_pred)[-ncol(y_pred)],
        direction = "long"
    )
    
    df <- merge(
        y_true_long[, c("sample", "cell_type", "y_true")],
        y_pred_long[, c("sample", "cell_type", "y_pred")],
        by = c("sample", "cell_type")
    )
    
    cell_type_pcc_df <- df %>%
        group_by(cell_type) %>%
        summarise(pcc = cor(y_true, y_pred, method = "pearson", use = "complete.obs"))
    
    cell_type_rmse_df <- df %>%
        group_by(cell_type) %>%
        summarise(rmse = sqrt(mean((y_true - y_pred)^2, na.rm = TRUE)))
    
    cell_colors <- rainbow(length(unique(df$cell_type)))
    
    df <- df %>%
        left_join(cell_type_pcc_df, by = "cell_type") %>%
        left_join(cell_type_rmse_df, by = "cell_type")
    
    p1 <- ggplot(df, aes(x = y_true, y = y_pred, color = cell_type)) +
        geom_point(alpha = 0.8, size = 1.8, show.legend = FALSE) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
        facet_wrap(~cell_type, scales = "free", ncol = 4) +
        theme_minimal() +
        ylim(0, xylim) + 
        xlim(0, xylim) +
        scale_color_manual(values = cell_colors) +
        labs(
            x = "True Cell Proportion",
            y = "Predicted Cell Proportion",
            title = title
        ) +
        geom_text(
            data = distinct(df, cell_type, pcc),
            aes(x = 0.05, y = xylim - 0.1, label = paste0("PCC = ", round(pcc, 2))),
            color = "black", 
            size = 3, 
            hjust = 0
        ) +
        geom_text(
            data = distinct(df, cell_type, rmse),
            aes(x = 0.05, y = xylim - 0.2, label = paste0("RMSE = ", round(rmse, 2))),
            color = "black", 
            size = 3, 
            hjust = 0
        )
    
    return(list(
        cell_type_pcc_df = cell_type_pcc_df,
        cell_type_rmse_df = cell_type_rmse_df,
        plot = p1,
        predictions = y_pred
    ))
}