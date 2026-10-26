## ---------------------------------------------------------------------------
## make_synthetic_data.R
##
## Builds a small synthetic human reference and mouse test object, used only
## to generate the example figures in docs/img/ and to check that the
## pipeline runs end to end. 
##
## Design: four human states (A_Macro, B_Mono, C_DC, D_rare) and three mouse
## states (mac1, mac2, mac10). Each state up-regulates 60 random genes; mouse
## mac1, mac2 and mac10 share their programmes with human A_Macro, B_Mono and
## C_DC respectively, so the expected matches are known. D_rare has no mouse
## counterpart.
##
## Usage (from the repository root):
##   Rscript examples/make_synthetic_data.R <output_folder>
##
## Then set REF_PATH and TEST_PATH in the config to
## <output_folder>/synthetic_ref.rds and <output_folder>/synthetic_test.rds,
## with REF_LABEL_COL and TEST_LABEL_COL set to "cell_state".
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(babelgene)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript examples/make_synthetic_data.R <output_folder>")
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

set.seed(1)

N_GENES       <- 1500
N_PROG_GENES  <- 60
FOLD_CHANGE   <- 8
N_CELLS_REF   <- 80   # per human state
N_CELLS_TEST  <- 60   # per mouse state

## Gene universe: 1:1 human-mouse ortholog pairs taken from babelgene's
## internal table (the package exposes no exported gene list).
db <- babelgene:::orthologs_df
human_pool <- unique(db$human_symbol[db$taxon_id == 10090])[1:4000]
pairs <- unique(babelgene::orthologs(
  genes = human_pool, species = "mouse", human = TRUE, min_support = 3
)[, c("human_symbol", "symbol")])
pairs <- pairs[!duplicated(pairs$human_symbol) & !duplicated(pairs$symbol), ]
pairs <- pairs[seq_len(N_GENES), ]

## Expression programmes.
base <- rgamma(N_GENES, shape = 0.5, rate = 0.5)
prog <- lapply(1:4, function(k) {
  v <- base
  idx <- sample(N_GENES, N_PROG_GENES)
  v[idx] <- v[idx] * FOLD_CHANGE
  v
})

simulate <- function(n_per_state, programmes, genes, prefix) {
  mats <- lapply(programmes, function(p) {
    depth <- runif(1, 0.8, 1.2)
    matrix(rpois(length(p) * n_per_state, rep(p, n_per_state) * depth),
           nrow = length(p))
  })
  m <- do.call(cbind, mats)
  rownames(m) <- genes
  colnames(m) <- paste0(prefix, seq_len(ncol(m)))
  obj <- CreateSeuratObject(counts = as(m, "dgCMatrix"))
  obj$cell_state <- rep(names(programmes), each = n_per_state)
  obj
}

ref <- simulate(
  N_CELLS_REF,
  setNames(prog[1:4], c("A_Macro", "B_Mono", "C_DC", "D_rare")),
  pairs$human_symbol, "h"
)
test <- simulate(
  N_CELLS_TEST,
  setNames(prog[1:3], c("mac1", "mac2", "mac10")),
  pairs$symbol, "m"
)

saveRDS(ref,  file.path(out, "synthetic_ref.rds"))
saveRDS(test, file.path(out, "synthetic_test.rds"))
message("Written: ", file.path(out, "synthetic_ref.rds"), ", ",
        file.path(out, "synthetic_test.rds"))
