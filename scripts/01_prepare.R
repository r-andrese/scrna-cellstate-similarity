## ---------------------------------------------------------------------------
## 01_prepare.R
##
## Loads the reference and test Seurat objects, maps genes to a shared gene
## universe (1:1 human-mouse orthologs across species, or the intersection of
## gene names within a species), and renormalises both datasets from raw
## counts on that universe.
##
## Usage (from the repository root):
##   Rscript scripts/01_prepare.R [config/config.R]
##
## Output: <OUT_DIR>/rds/01_prepared.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(babelgene)
})

source(file.path("R", "functions.R"))
cfg <- load_config()


## ---------------------------------------------------------------------------
## Load
## ---------------------------------------------------------------------------

inputs <- input_fingerprint(cfg)

message("Loading objects.")
ref_obj  <- readRDS(cfg$REF_PATH)
test_obj <- readRDS(cfg$TEST_PATH)

ref_counts  <- get_layer(ref_obj,  cfg$ASSAY, "counts")
test_counts <- get_layer(test_obj, cfg$ASSAY, "counts")

ref_labels  <- get_labels(ref_obj,  cfg$REF_LABEL_COL,  cfg$REF_NAME)
test_labels <- get_labels(test_obj, cfg$TEST_LABEL_COL, cfg$TEST_NAME)
ref_labels  <- ref_labels[colnames(ref_counts)]
test_labels <- test_labels[colnames(test_counts)]
if (anyNA(ref_labels) || anyNA(test_labels)) {
  stop("Cell names in the counts matrix do not match the object metadata.")
}

rm(ref_obj, test_obj)
invisible(gc())


## ---------------------------------------------------------------------------
## Shared gene universe
## ---------------------------------------------------------------------------

if (cfg$REF_SPECIES != cfg$TEST_SPECIES) {

  message("Mapping 1:1 human-mouse orthologs (babelgene).")
  human_counts <- if (cfg$REF_SPECIES == "human") ref_counts else test_counts
  mouse_counts <- if (cfg$REF_SPECIES == "mouse") ref_counts else test_counts

  orth <- orthologs_1to1(rownames(human_counts), cfg$ORTHOLOG_MIN_SUPPORT)
  n_pairs_total <- orth$n_pairs_total
  n_pairs_1to1  <- nrow(orth$table)
  message(sprintf(
    "  ortholog pairs: %d total, %d retained as 1:1, %d dropped as ambiguous.",
    n_pairs_total, n_pairs_1to1, n_pairs_total - n_pairs_1to1
  ))

  gene_map <- orth$table[
    orth$table$human_symbol %in% rownames(human_counts) &
      orth$table$mouse_symbol %in% rownames(mouse_counts), , drop = FALSE
  ]
  if (anyDuplicated(gene_map$human_symbol) || anyDuplicated(gene_map$mouse_symbol)) {
    stop("Shared gene universe contains duplicates after 1:1 filtering.")
  }

  human_counts <- human_counts[gene_map$human_symbol, , drop = FALSE]
  rownames(human_counts) <- gene_map$mouse_symbol
  mouse_counts <- mouse_counts[gene_map$mouse_symbol, , drop = FALSE]

  if (cfg$REF_SPECIES == "human") {
    ref_counts  <- human_counts
    test_counts <- mouse_counts
  } else {
    ref_counts  <- mouse_counts
    test_counts <- human_counts
  }
  rm(human_counts, mouse_counts)

} else {

  message("Same species: using the intersection of gene names.")
  shared <- intersect(rownames(ref_counts), rownames(test_counts))
  gene_map <- data.frame(gene = shared, stringsAsFactors = FALSE)
  n_pairs_total <- NA_integer_
  n_pairs_1to1  <- NA_integer_
  ref_counts  <- ref_counts[shared, , drop = FALSE]
  test_counts <- test_counts[shared, , drop = FALSE]
}

message(sprintf("  shared gene universe: %d genes.", nrow(ref_counts)))
if (nrow(ref_counts) == 0) stop("No shared genes. Check species and gene naming.")


## ---------------------------------------------------------------------------
## Cells, labels, normalisation
## ---------------------------------------------------------------------------

## Drop cells left with zero counts on the shared universe.
ref_keep  <- Matrix::colSums(ref_counts)  > 0
test_keep <- Matrix::colSums(test_counts) > 0
if (any(!ref_keep) || any(!test_keep)) {
  message(sprintf(
    "  dropping %d reference and %d test cells with zero shared-universe counts.",
    sum(!ref_keep), sum(!test_keep)
  ))
}
ref_counts  <- ref_counts[,  ref_keep,  drop = FALSE]
test_counts <- test_counts[, test_keep, drop = FALSE]
states_before <- c(levels(droplevels(ref_labels)), levels(droplevels(test_labels)))
ref_labels  <- droplevels(ref_labels[colnames(ref_counts)])
test_labels <- droplevels(test_labels[colnames(test_counts)])

## States all of whose cells were dropped disappear from every output.
states_lost <- setdiff(states_before, c(levels(ref_labels), levels(test_labels)))
if (length(states_lost)) {
  message("  states removed (no cells left): ", paste(states_lost, collapse = ", "))
}

if (nlevels(ref_labels) < 2 || nlevels(test_labels) < 2) {
  stop("At least two states are needed in each dataset.")
}
## Checked here as well as in 02, so the run fails before the slow steps.
n_ref <- table(ref_labels)
if (any(n_ref < 2)) {
  stop("Reference states with fewer than 2 cells cannot be modelled: ",
       paste(names(n_ref)[n_ref < 2], collapse = ", "), ". Remove or merge them.")
}

message("Renormalising both datasets on the shared gene universe.")
ref_data  <- normalise_shared(ref_counts)
test_data <- normalise_shared(test_counts)


## ---------------------------------------------------------------------------
## Save
## ---------------------------------------------------------------------------

run_id <- format(Sys.time(), "%Y%m%d-%H%M%OS3")

saveRDS(
  list(
    params        = cfg[PARAMS_01],
    run_id        = run_id,
    inputs        = inputs,
    ref_file      = basename(cfg$REF_PATH),
    test_file     = basename(cfg$TEST_PATH),
    gene_map      = gene_map,
    n_pairs_total = n_pairs_total,
    n_pairs_1to1  = n_pairs_1to1,
    n_dropped_ref  = sum(!ref_keep),
    n_dropped_test = sum(!test_keep),
    states_lost    = states_lost,
    ref_counts    = ref_counts,
    test_counts   = test_counts,
    ref_data      = ref_data,
    test_data     = test_data,
    ref_labels    = ref_labels,
    test_labels   = test_labels
  ),
  file.path(cfg$dir_rds, "01_prepared.rds")
)
writeLines(run_id, file.path(cfg$dir_rds, "01_run_id.txt"))
write_session_info(cfg, "01")

message("Done: ", file.path(cfg$dir_rds, "01_prepared.rds"))
