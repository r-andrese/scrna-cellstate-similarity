## ---------------------------------------------------------------------------
## Shared functions for scrna-cellstate-similarity.
## Sourced by every script in scripts/.
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------

## Read the config file given on the command line (default config/config.R),
## validate it, and create the output folders.
load_config <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  cfg_path <- if (length(args) >= 1) args[1] else file.path("config", "config.R")
  if (!file.exists(cfg_path)) {
    stop("Config file not found: ", cfg_path,
         "\nCopy config/config.example.R to config/config.R and edit it, ",
         "or pass the config path as the first argument.")
  }
  cfg <- new.env()
  sys.source(cfg_path, envir = cfg)
  cfg <- as.list(cfg)

  required <- c(
    "REF_PATH", "TEST_PATH", "OUT_DIR",
    "REF_SPECIES", "TEST_SPECIES", "ASSAY",
    "REF_LABEL_COL", "TEST_LABEL_COL", "REF_NAME", "TEST_NAME",
    "ORTHOLOG_MIN_SUPPORT",
    "SEED", "GLM_ALPHA", "GLM_NFOLDS", "GLM_THRESH", "GLM_MAXIT",
    "SINGLER_DE_METHOD", "HVG_N",
    "BASE_SIZE", "CLUSTER_HM", "CLUSTER_METHOD", "CELL_MM",
    "COR_LOW", "COR_MID", "COR_HIGH", "LR_LOW", "LR_MID", "LR_HIGH", "NA_COL"
  )
  missing <- setdiff(required, names(cfg))
  if (length(missing)) {
    stop("Config is missing: ", paste(missing, collapse = ", "))
  }

  species_ok <- c("human", "mouse")
  if (!cfg$REF_SPECIES %in% species_ok || !cfg$TEST_SPECIES %in% species_ok) {
    stop("REF_SPECIES and TEST_SPECIES must each be \"human\" or \"mouse\".")
  }
  if (identical(cfg$REF_NAME, cfg$TEST_NAME)) {
    stop("REF_NAME and TEST_NAME must differ.")
  }
  if (!is.numeric(cfg$HVG_N) || length(cfg$HVG_N) != 1 || cfg$HVG_N < 1) {
    stop("HVG_N must be a single positive number.")
  }
  if (!is.numeric(cfg$GLM_NFOLDS) || cfg$GLM_NFOLDS < 3) {
    stop("GLM_NFOLDS must be at least 3.")
  }
  if (!is.logical(cfg$CLUSTER_HM) || length(cfg$CLUSTER_HM) != 1 || is.na(cfg$CLUSTER_HM)) {
    stop("CLUSTER_HM must be TRUE or FALSE.")
  }
  hc_methods <- c("ward.D", "ward.D2", "single", "complete", "average",
                  "mcquitty", "median", "centroid")
  if (!cfg$CLUSTER_METHOD %in% hc_methods) {
    stop("CLUSTER_METHOD must be one of: ", paste(hc_methods, collapse = ", "))
  }

  cfg$config_path <- cfg_path
  cfg$dir_rds     <- file.path(cfg$OUT_DIR, "rds")
  cfg$dir_fig     <- file.path(cfg$OUT_DIR, "figures")
  cfg$dir_tab     <- file.path(cfg$OUT_DIR, "tables")
  for (d in c(cfg$dir_rds, cfg$dir_fig, cfg$dir_tab)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  cfg
}

## Parameters that determine the output of each step. Input files are tracked
## separately by input_fingerprint().
PARAMS_01 <- c("REF_SPECIES", "TEST_SPECIES", "ASSAY",
               "REF_LABEL_COL", "TEST_LABEL_COL", "REF_NAME", "TEST_NAME",
               "ORTHOLOG_MIN_SUPPORT")
PARAMS_02 <- c("SEED", "GLM_ALPHA", "GLM_NFOLDS", "GLM_THRESH", "GLM_MAXIT",
               "SINGLER_DE_METHOD", "HVG_N")

## Exact comparison; integers and doubles of equal value are treated as equal.
same_value <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) {
    a <- as.double(a)
    b <- as.double(b)
  }
  identical(unname(a), unname(b))
}

## Stop if the stored parameters of an upstream step differ from the config.
check_params <- function(stored, cfg, keys, rerun) {
  diff <- keys[!vapply(keys, function(k) same_value(stored[[k]], cfg[[k]]),
                       logical(1))]
  if (length(diff)) {
    stop("Config changed since ", rerun, " was run (",
         paste(diff, collapse = ", "), "). Rerun ", rerun, " first.",
         call. = FALSE)
  }
  invisible(TRUE)
}

## Path, size and modification time of the two input objects.
input_fingerprint <- function(cfg) {
  f <- c(ref = cfg$REF_PATH, test = cfg$TEST_PATH)
  if (!all(file.exists(f))) {
    stop("Input object not found: ", paste(f[!file.exists(f)], collapse = ", "))
  }
  info <- file.info(f)
  data.frame(role  = names(f),
             path  = normalizePath(f),
             size  = info$size,
             mtime = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
             stringsAsFactors = FALSE, row.names = NULL)
}

## Stop if the input objects named in the config are not those used by 01.
check_inputs <- function(stored, cfg) {
  now <- input_fingerprint(cfg)
  if (!identical(stored, now)) {
    stop("REF_PATH or TEST_PATH, or the files they point to, changed since ",
         "01_prepare.R was run. Rerun 01_prepare.R first.", call. = FALSE)
  }
  invisible(TRUE)
}

write_session_info <- function(cfg, step) {
  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(cfg$OUT_DIR, paste0("session_info_", step, ".txt")))
}


## ---------------------------------------------------------------------------
## Data access and normalisation
## ---------------------------------------------------------------------------

## Retrieve a layer from an assay across Seurat v3/v4/v5 object versions.
get_layer <- function(obj, assay = "RNA", layer = "counts") {
  if (!assay %in% names(obj@assays)) {
    stop("Assay not found in object: ", assay)
  }
  a <- obj[[assay]]
  if (inherits(a, "Assay5")) {
    lyrs <- SeuratObject::Layers(obj, assay = assay)
    if (any(grepl(paste0("^", layer, "\\."), lyrs))) {
      obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]])
    }
    return(SeuratObject::LayerData(obj, assay = assay, layer = layer))
  }
  m <- try(SeuratObject::LayerData(obj, assay = assay, layer = layer), silent = TRUE)
  if (inherits(m, "try-error")) {
    m <- SeuratObject::GetAssayData(obj, assay = assay, slot = layer)
  }
  m
}

## Cell-state labels from a metadata column, prefixed with the dataset name.
get_labels <- function(obj, col, prefix) {
  if (!col %in% colnames(obj@meta.data)) {
    stop("Label column not found in object metadata: ", col)
  }
  lab <- as.character(obj@meta.data[[col]])
  if (anyNA(lab)) {
    stop(sum(is.na(lab)), " cells have a missing label in column ", col,
         ". Subset the object to labelled cells before running.")
  }
  lab <- factor(paste(prefix, lab, sep = "."))
  names(lab) <- colnames(obj)
  lab
}

## log1p CP10K over the gene universe the matrix is restricted to.
normalise_shared <- function(counts) {
  Seurat::NormalizeData(
    counts,
    normalization.method = "LogNormalize",
    scale.factor = 1e4,
    verbose = FALSE
  )
}


## ---------------------------------------------------------------------------
## Orthologs
## ---------------------------------------------------------------------------

## 1:1 human-mouse orthologs (babelgene), queried with the human gene list.
## Pairs are retained only if unambiguous in both directions within the
## returned table. Returns a data.frame with columns human_symbol and
## mouse_symbol, plus the number of pairs before filtering.
orthologs_1to1 <- function(human_genes, min_support = 3) {
  orth <- babelgene::orthologs(
    genes = human_genes, species = "mouse",
    human = TRUE, min_support = min_support
  )
  req <- c("human_symbol", "symbol")
  if (!all(req %in% colnames(orth))) {
    stop("babelgene::orthologs did not return the expected columns: ",
         paste(setdiff(req, colnames(orth)), collapse = ", "))
  }
  orth <- unique(orth[, req])

  h_tab <- table(orth$human_symbol)
  m_tab <- table(orth$symbol)
  keep  <- orth$human_symbol %in% names(h_tab)[h_tab == 1] &
           orth$symbol       %in% names(m_tab)[m_tab == 1]

  out <- data.frame(human_symbol = orth$human_symbol[keep],
                    mouse_symbol = orth$symbol[keep],
                    stringsAsFactors = FALSE)
  list(table = out, n_pairs_total = nrow(orth))
}


## ---------------------------------------------------------------------------
## Pseudobulk correlation
## ---------------------------------------------------------------------------

## Sum raw counts within each label. Returns genes x labels.
pseudobulk_counts <- function(counts, labels) {
  labels <- droplevels(factor(labels))
  ind <- Matrix::sparse.model.matrix(~ 0 + labels)
  colnames(ind) <- levels(labels)
  as.matrix(counts %*% ind)
}

## Counts-per-million then log1p.
pseudobulk_logcpm <- function(pb) {
  ls <- colSums(pb)
  if (any(ls == 0)) stop("Pseudobulk column with zero total counts.")
  log1p(t(t(pb) / ls) * 1e6)
}

## Rank genes by vst standardised variance, independently within a dataset.
rank_hvg <- function(counts, n_max) {
  ## Seurat replaces underscores in feature names, which would make the
  ## returned names disagree with the pseudobulk row names.
  n_us <- sum(grepl("_", rownames(counts)))
  if (n_us > 0) {
    warning(n_us, " gene names contain '_' and will not match after HVG ",
            "ranking; they are effectively excluded from the correlation.")
  }
  o <- Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  o <- Seurat::FindVariableFeatures(
    o, selection.method = "vst", nfeatures = n_max, verbose = FALSE
  )
  hv <- Seurat::HVFInfo(o)
  vcol <- grep("variance.standardized", colnames(hv), value = TRUE)[1]
  if (is.na(vcol)) stop("Could not locate vst.variance.standardized in HVFInfo output.")
  rownames(hv)[order(hv[[vcol]], decreasing = TRUE)]
}

## z-score each gene across the labels of one dataset.
zscore_rows <- function(m) {
  mu <- rowMeans(m)
  sd <- apply(m, 1, stats::sd)
  list(z = (m - mu) / sd, sd = sd)
}

## Spearman correlation of z-scored pseudobulk profiles, in the
## feature space of the top-n HVGs shared by both datasets.
cor_hvg <- function(n, ref_pb, test_pb, hvg_ref, hvg_test) {
  feats <- intersect(head(hvg_ref, n), head(hvg_test, n))
  feats <- intersect(feats, intersect(rownames(ref_pb), rownames(test_pb)))
  if (length(feats) < 50) {
    warning(sprintf("Only %d shared HVGs at n = %d.", length(feats), n))
  }
  zr <- zscore_rows(ref_pb[feats, , drop = FALSE])
  zt <- zscore_rows(test_pb[feats, , drop = FALSE])
  ## Genes invariant across labels in either dataset carry no information
  ## and would produce NaN on z-scoring; removed from both.
  ok <- zr$sd > 0 & zt$sd > 0
  zr <- zr$z[ok, , drop = FALSE]
  zt <- zt$z[ok, , drop = FALSE]
  list(
    hvg_n       = n,
    n_features  = nrow(zr),
    features    = rownames(zr),
    spearman    = cor(zr, zt, method = "spearman")
  )
}


## ---------------------------------------------------------------------------
## Elastic-net logistic regression
## ---------------------------------------------------------------------------

## Offset equal to the log-odds of the positive class, so that the fitted
## coefficients describe deviation from the class prior.
getPopulationOffset <- function(y) {
  if (!is.factor(y)) y <- factor(y)
  if (length(levels(y)) != 2) stop("y must be a two-level factor")
  off <- sum(y == levels(y)[2]) / length(y)
  off <- log(off / (1 - off))
  rep(off, length(y))
}

## One-vs-rest elastic-net logistic regression. A model is trained for each
## reference state on all reference cells and applied to every test cell with
## a zero offset. One fold partition is shared by every fit.
glm_predict <- function(train.data, train.group, test.data, test.group,
                        alpha, nfolds, seed, thresh, maxit) {
  if (nrow(train.data) < 50) {
    message("Fewer than 50 features in the training data.")
  }

  set.seed(seed)
  foldid <- sample(rep(seq_len(nfolds), length.out = ncol(train.data)))

  train.group <- droplevels(factor(train.group))
  groups      <- levels(train.group)
  converged   <- setNames(rep(TRUE, length(groups)), groups)
  n_train     <- setNames(as.integer(table(train.group)[groups]), groups)

  ## glmnet cannot fit a class with fewer than two cells. With fewer cells
  ## than folds, some folds contain no cell of the state.
  too_few <- names(n_train)[n_train < 2]
  if (length(too_few)) {
    stop("Reference states with fewer than 2 cells cannot be modelled: ",
         paste(too_few, collapse = ", "), ". Remove or merge them.")
  }
  few <- names(n_train)[n_train < nfolds]
  if (length(few)) {
    warning("Reference states with fewer cells than CV folds: ",
            paste(few, collapse = ", "), ". Their lambda selection is unreliable.")
  }
  nonzero <- setNames(rep(NA_integer_, length(groups)), groups)

  x_train <- Matrix::t(train.data)
  x_test  <- Matrix::t(test.data)
  logits  <- vector("list", length(groups))
  names(logits) <- groups

  for (g in groups) {
    fac <- factor(train.group == g)
    ok  <- TRUE
    fit <- withCallingHandlers(
      glmnet::cv.glmnet(
        x = x_train, y = fac,
        offset = getPopulationOffset(fac),
        family = "binomial", intercept = FALSE,
        alpha = alpha, foldid = foldid, type.measure = "class",
        thresh = thresh, maxit = maxit
      ),
      warning = function(w) {
        if (grepl("Convergence for", conditionMessage(w))) {
          ok <<- FALSE
          invokeRestart("muffleWarning")
        }
      }
    )
    converged[g] <- ok
    if (!ok) {
      message(sprintf("  glmnet path did not converge for: %s (n = %d)",
                      g, n_train[[g]]))
    }
    ## Coefficients retained at lambda.min. Zero means the null model was
    ## selected, which predicts a similarity of exactly 0.5 for every cell.
    nonzero[g] <- as.integer(sum(as.numeric(
      stats::coef(fit, s = "lambda.min"))[-1] != 0))
    logits[[g]] <- as.numeric(predict(
      object    = fit,
      newx      = x_test,
      newoffset = rep(0, ncol(test.data)),
      s         = "lambda.min"
    ))
  }

  logits <- as.data.frame(logits, check.names = FALSE)
  rownames(logits) <- colnames(test.data)
  prob <- (1 + exp(-logits)) ** -1

  list(
    test.group  = test.group,
    logits      = logits,
    probability = prob,
    cluster     = colnames(prob)[apply(prob, 1, which.max)],
    foldid      = foldid,
    converged   = converged,
    n_train     = n_train,
    n_nonzero   = nonzero
  )
}


## ---------------------------------------------------------------------------
## Plotting
## ---------------------------------------------------------------------------

## Natural sort
natural_order <- function(x) {
  num  <- suppressWarnings(as.numeric(sub("^.*?([0-9]+)$", "\\1", x)))
  stem <- sub("([0-9]+)$", "", x)
  x[order(stem, ifelse(is.na(num), 0, num), x)]
}

## Heatmap with fixed cell size.
draw_hm <- function(m, title, file, col_fun, cfg,
                    cluster_rows = FALSE, cluster_columns = FALSE) {
  bs <- cfg$BASE_SIZE

  hm <- ComplexHeatmap::Heatmap(
    m,
    col                  = col_fun,
    na_col               = cfg$NA_COL,
    name                 = title,
    column_title         = "Test",
    row_title            = "Training",
    cluster_rows         = cluster_rows,
    cluster_columns      = cluster_columns,
    show_row_names       = TRUE,
    row_names_side       = "right",
    row_title_gp         = grid::gpar(fontsize = bs),
    column_title_gp      = grid::gpar(fontsize = bs),
    row_names_gp         = grid::gpar(fontsize = bs),
    column_names_gp      = grid::gpar(fontsize = bs),
    width                = ncol(m) * grid::unit(cfg$CELL_MM, "mm"),
    height               = nrow(m) * grid::unit(cfg$CELL_MM, "mm"),
    heatmap_legend_param = list(
      direction = "horizontal",
      title_gp  = grid::gpar(fontsize = bs),
      labels_gp = grid::gpar(fontsize = bs - 2)
    )
  )

  pad <- grid::unit(c(4, 4, 4, 4), "mm")

  ## Measure on a null device, then draw at the measured size.
  grDevices::pdf(NULL)
  d <- ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom", padding = pad)
  w <- grid::convertWidth(ComplexHeatmap::width.HeatmapList(d), "inch", valueOnly = TRUE)
  h <- grid::convertHeight(ComplexHeatmap::height.HeatmapList(d), "inch", valueOnly = TRUE)
  grDevices::dev.off()

  grDevices::pdf(file, width = w + 0.5, height = h + 0.5)
  ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom", padding = pad)
  grDevices::dev.off()
  invisible(file)
}

## Matrix to data.frame with the row names as the first column.
as_df <- function(m, rowname) {
  d <- data.frame(rownames(m), as.data.frame(m, check.names = FALSE),
                  check.names = FALSE, stringsAsFactors = FALSE)
  colnames(d)[1] <- rowname
  d
}
