## ---------------------------------------------------------------------------
## Configuration for scrna-cellstate-similarity
##
## Copy this file to config/config.R and edit the values as needed.
## Every script takes the config path as its only argument:
##
##   Rscript scripts/01_prepare.R config/config.R
##
## If no argument is given, config/config.R is used.
##
## Parameters are grouped by the script that uses them. Scripts 02 and 03
## check that the parameters of upstream steps still match the results they
## read, and stop with a message naming the script to rerun if they do not.
## Plotting parameters (section 03) can be changed freely and only require
## rerunning 03.
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## Inputs and outputs
## ---------------------------------------------------------------------------

## Seurat objects. Raw counts must be present in the
## counts layer (Seurat v5) or counts slot (Seurat v3/v4) of ASSAY.
REF_PATH  <- "/path/to/reference_object.rds"
TEST_PATH <- "/path/to/test_object.rds"

## All outputs are written here
OUT_DIR   <- "/path/to/output_folder"

## ---------------------------------------------------------------------------
## 01_prepare.R
## ---------------------------------------------------------------------------

## Species of each object: "human" or "mouse".
## If both are the same species, no ortholog mapping is performed and the
## shared gene universe is the intersection of gene names.
REF_SPECIES  <- "human"
TEST_SPECIES <- "mouse"

## Assay holding raw counts in both objects.
ASSAY <- "RNA"

## Metadata columns holding the cell-state labels to compare. The object should
## be subset beforehand so there are no missing labels.
REF_LABEL_COL  <- "cell_state"
TEST_LABEL_COL <- "cell_state"

## Short prefixes added to the labels of each dataset.
REF_NAME  <- "ref"
TEST_NAME <- "test"

## Minimum number of databases supporting an ortholog pair (babelgene).
ORTHOLOG_MIN_SUPPORT <- 3


## ---------------------------------------------------------------------------
## 02_similarity.R
## ---------------------------------------------------------------------------

SEED <- 1234

## Elastic-net logistic regression (glmnet), one-vs-rest per reference state.
GLM_ALPHA  <- 0.99   # mixing parameter; 1 = lasso, 0 = ridge
GLM_NFOLDS <- 10     # cross-validation folds for lambda selection
GLM_THRESH <- 1e-6   # convergence tolerance (glmnet default 1e-7)
GLM_MAXIT  <- 1e6    # maximum passes over the data (glmnet default 1e5)

## SingleR marker detection method.
SINGLER_DE_METHOD <- "wilcox"

## Pseudobulk correlation. Highly variable genes are ranked independently in
## each dataset; the feature space is the intersection of the top HVG_N of each.
HVG_N <- 2000


## ---------------------------------------------------------------------------
## 03_figures_tables.R (plotting only; can be changed without rerunning 01/02)
## ---------------------------------------------------------------------------

BASE_SIZE <- 14

## Hierarchical clustering of the Spearman heatmap. The resulting order is
## also applied to the logistic regression heatmap.
CLUSTER_HM     <- TRUE
CLUSTER_METHOD <- "ward.D2"

## Heatmap cell size in mm.
CELL_MM <- 8

## Correlation heatmap diverging colour scale.
COR_LOW  <- "#1965B0"
COR_MID  <- "white"
COR_HIGH <- "#DC050C"

## Logistic regression heatmap colour scale.
LR_LOW  <- "#e9e9e9"
LR_MID  <- "white"
LR_HIGH <- "red"

## Colour for masked (non-converged) cells.
NA_COL <- "#e9e9e9"
