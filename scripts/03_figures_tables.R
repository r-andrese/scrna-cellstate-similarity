## ---------------------------------------------------------------------------
## 03_figures_tables.R
##
## Draws the similarity heatmaps and SingleR diagnostics, and writes all
## results to an .xlsx workbook.
##
## Usage (from the repository root):
##   Rscript scripts/03_figures_tables.R [config/config.R]
##
## Input:  <OUT_DIR>/rds/02_results.rds
## Output: <OUT_DIR>/figures/*.pdf, <OUT_DIR>/tables/similarity_results.xlsx
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleR)
  library(ComplexHeatmap)
  library(circlize)
  library(openxlsx)
  library(grid)
})

source(file.path("R", "functions.R"))
cfg <- load_config()

res_file <- file.path(cfg$dir_rds, "02_results.rds")
if (!file.exists(res_file)) stop("Not found: ", res_file, ". Run 02_similarity.R first.")
res <- readRDS(res_file)

## 01 rerun after 02: the results no longer correspond to the prepared data.
id_file <- file.path(cfg$dir_rds, "01_run_id.txt")
if (file.exists(id_file) && !identical(readLines(id_file), res$prep_run_id)) {
  stop("01_prepare.R was rerun after 02_similarity.R. Rerun 02_similarity.R first.",
       call. = FALSE)
}
check_params(res$params, cfg, PARAMS_01, "01_prepare.R")
check_params(res$params, cfg, PARAMS_02, "02_similarity.R")
check_inputs(res$inputs, cfg)

## All reported analysis parameters come from the stored results, not the
## current config.
p <- res$params

pred        <- res$pred
glm_res     <- res$glm_res
test_labels <- res$test_labels
cor_res <- res$cor_res


## ---------------------------------------------------------------------------
## SingleR diagnostics
## ---------------------------------------------------------------------------

n_ref_states <- nlevels(res$ref_labels)

pdf(file.path(cfg$dir_fig, "singler_score_heatmap.pdf"),
    width = 15, height = max(6, 2 + 0.25 * n_ref_states))
plotScoreHeatmap(pred, order.by = "clusters", clusters = test_labels,
                 fontsize = cfg$BASE_SIZE)
invisible(dev.off())

pdf(file.path(cfg$dir_fig, "singler_delta_distribution.pdf"),
    width = 10, height = max(8, 2.5 * ceiling(n_ref_states / 3)))
print(plotDeltaDistribution(pred, ncol = 3) +
        ggplot2::theme_classic(base_size = cfg$BASE_SIZE))
invisible(dev.off())


## ---------------------------------------------------------------------------
## Row and column order
## ---------------------------------------------------------------------------
## Derived once from the Spearman matrix and applied to every heatmap, so the
## panels can be read against each other.

cor_sp <- cor_res$spearman

if (cfg$CLUSTER_HM) {
  hc_row  <- hclust(dist(cor_sp),    method = cfg$CLUSTER_METHOD)
  hc_col  <- hclust(dist(t(cor_sp)), method = cfg$CLUSTER_METHOD)
  ord_row <- rownames(cor_sp)[hc_row$order]
  ord_col <- colnames(cor_sp)[hc_col$order]
} else {
  hc_row  <- FALSE
  hc_col  <- FALSE
  ord_row <- natural_order(rownames(cor_sp))
  ord_col <- natural_order(colnames(cor_sp))
  cor_sp  <- cor_sp[ord_row, ord_col, drop = FALSE]
}


## ---------------------------------------------------------------------------
## Logistic regression heatmap
## ---------------------------------------------------------------------------
## A reference state whose regularisation path did not converge returns a
## logit of zero for every test cell, i.e. a similarity of exactly 0.5. Those
## rows are masked to NA so they are not read as intermediate similarity.

lr_mat <- t(res$lr_similarity)   # rows = reference, columns = test

failed <- names(glm_res$converged)[!glm_res$converged]
if (length(failed)) {
  lr_mat[rownames(lr_mat) %in% failed, ] <- NA_real_
  message("Masking non-converged reference states: ", paste(failed, collapse = ", "))
}

lr_mat <- lr_mat[intersect(ord_row, rownames(lr_mat)),
                 intersect(ord_col, colnames(lr_mat)), drop = FALSE]

draw_hm(lr_mat, "Predicted\nsimilarity",
        file.path(cfg$dir_fig, "lr_similarity_heatmap.pdf"),
        col_fun = colorRamp2(c(0, 0.5, 1), c(cfg$LR_LOW, cfg$LR_MID, cfg$LR_HIGH)),
        cfg = cfg)


## ---------------------------------------------------------------------------
## Correlation heatmap
## ---------------------------------------------------------------------------

cor_col_fun <- function(m) {
  lim <- max(abs(m), na.rm = TRUE)
  colorRamp2(c(-lim, 0, lim), c(cfg$COR_LOW, cfg$COR_MID, cfg$COR_HIGH))
}

draw_hm(cor_sp, "Spearman\nrho",
        file.path(cfg$dir_fig, "pseudobulk_spearman_heatmap.pdf"),
        col_fun = cor_col_fun(cor_sp), cfg = cfg,
        cluster_rows = hc_row, cluster_columns = hc_col)


## ---------------------------------------------------------------------------
## Tables (.xlsx)
## ---------------------------------------------------------------------------

converge_tab <- data.frame(
  reference_state = names(glm_res$converged),
  n_train_cells   = as.integer(glm_res$n_train[names(glm_res$converged)]),
  converged       = as.logical(glm_res$converged),
  n_nonzero_coef  = as.integer(glm_res$n_nonzero[names(glm_res$converged)]),
  stringsAsFactors = FALSE
)

tv <- as.character(test_labels)

singler_tab <- as.data.frame(table(
  test_state    = tv,
  singler_label = as.character(pred$labels)
), stringsAsFactors = FALSE)
singler_tab <- singler_tab[singler_tab$Freq > 0, ]
singler_tab <- singler_tab[order(match(singler_tab$test_state, natural_order(unique(tv))),
                                 -singler_tab$Freq), ]
rownames(singler_tab) <- NULL
colnames(singler_tab)[3] <- "n_cells"

n_cells_by_state  <- table(tv)
n_pruned_by_state <- tapply(is.na(pred$pruned.labels), tv, sum)
pruned_tab <- data.frame(
  test_state = names(n_pruned_by_state),
  n_cells    = as.integer(n_cells_by_state[names(n_pruned_by_state)]),
  n_pruned   = as.integer(n_pruned_by_state),
  stringsAsFactors = FALSE
)
pruned_tab$prop_pruned <- pruned_tab$n_pruned / pruned_tab$n_cells
pruned_tab <- pruned_tab[match(natural_order(pruned_tab$test_state), pruned_tab$test_state), ]
rownames(pruned_tab) <- NULL

ps <- res$prep_summary
fmt <- function(x) paste(format(x, scientific = FALSE, trim = TRUE), collapse = ", ")

meta_tab <- data.frame(
  item = c("reference_file", "test_file",
           "reference_species", "test_species", "assay",
           "reference_label_column", "test_label_column",
           "ortholog_min_support", "ortholog_pairs_total", "ortholog_pairs_1to1",
           "shared_genes", "cells_dropped_reference", "cells_dropped_test",
           "n_reference_states", "n_test_states",
           "seed", "glm_alpha", "glm_nfolds", "glm_thresh", "glm_maxit",
           "n_states_nonconverged", "singler_de_method",
           "hvg_n", "n_features_correlation",
           "states_removed_no_cells",
           "clustering", "cluster_method"),
  value = c(ps$ref_file, ps$test_file,
            p$REF_SPECIES, p$TEST_SPECIES, p$ASSAY,
            p$REF_LABEL_COL, p$TEST_LABEL_COL,
            fmt(p$ORTHOLOG_MIN_SUPPORT), fmt(ps$n_pairs_total), fmt(ps$n_pairs_1to1),
            fmt(nrow(ps$gene_map)), fmt(ps$n_dropped_ref), fmt(ps$n_dropped_test),
            fmt(nlevels(res$ref_labels)), fmt(nlevels(res$test_labels)),
            fmt(p$SEED), fmt(p$GLM_ALPHA), fmt(p$GLM_NFOLDS),
            fmt(p$GLM_THRESH), fmt(p$GLM_MAXIT),
            fmt(sum(!glm_res$converged)), p$SINGLER_DE_METHOD,
            fmt(p$HVG_N), fmt(cor_res$n_features),
            if (length(ps$states_lost)) paste(ps$states_lost, collapse = ", ") else "none",
            fmt(cfg$CLUSTER_HM), if (cfg$CLUSTER_HM) cfg$CLUSTER_METHOD else "none"),
  stringsAsFactors = FALSE
)

features_tab <- data.frame(gene = cor_res$features, stringsAsFactors = FALSE)
if ("human_symbol" %in% colnames(ps$gene_map)) {
  features_tab$human_symbol <- ps$gene_map$human_symbol[
    match(features_tab$gene, ps$gene_map$mouse_symbol)
  ]
}

write.xlsx(
  list(
    metadata              = meta_tab,
    LR_similarity         = as_df(lr_mat, "reference_state"),
    LR_convergence        = converge_tab,
    cor_spearman          = as_df(cor_sp[ord_row, ord_col, drop = FALSE], "reference_state"),
    features_used         = features_tab,
    SingleR_assignments   = singler_tab,
    SingleR_pruning       = pruned_tab
  ),
  file = file.path(cfg$dir_tab, "similarity_results.xlsx"),
  overwrite = TRUE
)
write_session_info(cfg, "03")

message("Done. Output written to: ", cfg$OUT_DIR)
