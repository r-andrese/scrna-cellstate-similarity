## ---------------------------------------------------------------------------
## 02_similarity.R
##
## Computes three readouts of similarity between reference and test states:
##   1. SingleR label transfer (per cell)
##   2. Elastic-net logistic regression trained on the reference, applied to
##      the test cells, summarised per test state
##   3. Correlation of z-scored pseudobulk profiles in a shared HVG space
##
## Usage (from the repository root):
##   Rscript scripts/02_similarity.R [config/config.R]
##
## Input:  <OUT_DIR>/rds/01_prepared.rds
## Output: <OUT_DIR>/rds/02_results.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(SingleR)
  library(glmnet)
})

source(file.path("R", "functions.R"))
cfg <- load_config()

prep_file <- file.path(cfg$dir_rds, "01_prepared.rds")
if (!file.exists(prep_file)) stop("Not found: ", prep_file, ". Run 01_prepare.R first.")
prep <- readRDS(prep_file)
check_params(prep$params, cfg, PARAMS_01, "01_prepare.R")
check_inputs(prep$inputs, cfg)

set.seed(cfg$SEED)


## ---------------------------------------------------------------------------
## 1. SingleR
## ---------------------------------------------------------------------------

message("Running SingleR.")
pred <- SingleR(
  test      = prep$test_data,
  ref       = prep$ref_data,
  labels    = prep$ref_labels,
  de.method = cfg$SINGLER_DE_METHOD
)


## ---------------------------------------------------------------------------
## 2. Logistic regression
## ---------------------------------------------------------------------------

message("Running elastic-net logistic regression.")
glm_res <- glm_predict(
  train.data  = prep$ref_data,  train.group = prep$ref_labels,
  test.data   = prep$test_data, test.group  = prep$test_labels,
  alpha  = cfg$GLM_ALPHA,  nfolds = cfg$GLM_NFOLDS, seed = cfg$SEED,
  thresh = cfg$GLM_THRESH, maxit  = cfg$GLM_MAXIT
)

if (any(!glm_res$converged)) {
  message("  states with non-converged paths: ",
          paste(names(glm_res$converged)[!glm_res$converged], collapse = ", "))
}

## Mean logit per test state, transformed to the probability scale.
## Rows = test states, columns = reference states.
lr_mean_logit <- apply(
  glm_res$logits, 2, function(e) vapply(split(e, glm_res$test.group), mean, numeric(1))
)
lr_similarity <- (1 + exp(-lr_mean_logit)) ** -1


## ---------------------------------------------------------------------------
## 3. Pseudobulk correlation
## ---------------------------------------------------------------------------

message("Building pseudobulk profiles.")
ref_pb  <- pseudobulk_logcpm(pseudobulk_counts(prep$ref_counts,  prep$ref_labels))
test_pb <- pseudobulk_logcpm(pseudobulk_counts(prep$test_counts, prep$test_labels))

message("Ranking HVGs independently within each dataset.")
hvg_ref  <- rank_hvg(prep$ref_counts,  cfg$HVG_N)
hvg_test <- rank_hvg(prep$test_counts, cfg$HVG_N)

cor_res <- cor_hvg(cfg$HVG_N, ref_pb = ref_pb, test_pb = test_pb,
                   hvg_ref = hvg_ref, hvg_test = hvg_test)


## ---------------------------------------------------------------------------
## Save
## ---------------------------------------------------------------------------

saveRDS(
  list(
    params         = cfg[c(PARAMS_01, PARAMS_02)],
    prep_run_id    = prep$run_id,
    inputs         = prep$inputs,
    prep_summary   = prep[c("ref_file", "test_file", "gene_map",
                            "n_pairs_total", "n_pairs_1to1",
                            "n_dropped_ref", "n_dropped_test", "states_lost")],
    ref_labels     = prep$ref_labels,
    test_labels    = prep$test_labels,
    pred           = pred,
    glm_res        = glm_res,
    lr_similarity  = lr_similarity,
    ref_pb         = ref_pb,
    test_pb        = test_pb,
    hvg_ref        = hvg_ref,
    hvg_test       = hvg_test,
    cor_res        = cor_res
  ),
  file.path(cfg$dir_rds, "02_results.rds")
)
write_session_info(cfg, "02")

message("Done: ", file.path(cfg$dir_rds, "02_results.rds"))
