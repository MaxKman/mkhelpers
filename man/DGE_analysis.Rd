% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DGE_analysis.R
\name{DGE_analysis}
\alias{DGE_analysis}
\title{DEseq2 wrapper for pseudobulk differential gene expression analysis in single cell data}
\usage{
DGE_analysis(
  m,
  md,
  cluster_col,
  sample_col,
  group_col,
  batch_col = NULL,
  balance_batches = FALSE,
  add_var = NULL,
  title,
  group1,
  group2,
  design = ~group,
  n_cells_normalize,
  n_cells_min,
  min_n_samples_group,
  cell_name_col = cell_name,
  exp_percentage = 5,
  save_results = FALSE,
  savepath = "../../RDS/DGE"
)
}
\arguments{
\item{m}{A single cell matrix with genes as rows and cells as columns}

\item{md}{A metadata dataframe (must contain cluster_col, sample_col and group_col)}

\item{cluster_col}{The metadata column identifying the cluster of each cell}

\item{sample_col}{The metadata column identifying the sample to which each cell belongs (e.g. the donor)}

\item{group_col}{The metadata column identifying the group / condition to which each cell belongs (e.g. treatment vs. control)}

\item{batch_col}{The metadata column identifying the batch to which each cell belongs. Note, that to correct for a batch effect in the analysis the design formula needs to be adapted accordingly}

\item{balance_batches}{If set to TRUE, will ensure that only samples are included in the comparison for which a corresponding sample from the same batch is available in the other group. Requires batch_col to be defined.}

\item{add_var}{An additional column that should be preserved, for example because it is used in the design formula}

\item{title}{The title of comparison (will be used for file naming)}

\item{group1}{The name of group 1 in group_col}

\item{group2}{The name of group 2 in group_col}

\item{design}{The DEseq2 design formula to be used}

\item{n_cells_normalize}{The number of cells to which the pseudobulking will be standardized. E.g. if set to 10000, the counts aggregated from a cluster/sample combination will be normalized as follows: count_matrix / n_cells * n_cells_normalize.}

\item{n_cells_min}{The minimum number of cells for a cluster/sample combination to be included in the analysis. E.g. if set to 20, samples will be excluded for clusters, for which they contain less than 20 cells}

\item{min_n_samples_group}{Minimum number of samples in each group to perform DGE analysis. E.g. if set to 3 only 3 versus 3 comparisons will be carried out.}

\item{cell_name_col}{The name of the column that contains the cell identifier}

\item{exp_percentage}{Minimum percentage of cells in a cluster from one group, which need to express a gene for it to be included in the DGE analysis}

\item{save_results}{Whether to save results to an RDS file (default: FALSE)}

\item{savepath}{Path where the results of the analysis should be saved}
}
\value{
A dataframe with the DEseq2 results with an additional column indicating the cluster
}
\description{
This function can be run in parallel by setting up a future environment. See examples.
}
\examples{
\dontrun{
library(Seurat)
library(SeuratData)
InstallData("ifnb")
pbmc_ifnb <- LoadData("ifnb")

# Simulate some batches and donors
add_donors <- pbmc_ifnb@meta.data \%>\%
  rownames_to_column("cell_name") \%>\%
  group_by(stim) \%>\%
  mutate(simulated_batches = sample(4, n(), replace = TRUE)) \%>\%
  ungroup \%>\%
  mutate(simulated_donors = simulated_batches) \%>\%
  mutate(simulated_donors = str_c(simulated_donors, "_", stim)) \%>\%
  column_to_rownames("cell_name") \%>\%
  select(simulated_donors, simulated_batches)
add_donors$cell_name <- rownames(add_donors)
pbmc_ifnb <- AddMetaData(pbmc_ifnb, add_donors)

# Remove a donor from one batch to test whether the entire batch is
# correctly removed from the analysis (if balance_batches is set to true)
donor_remove <- pbmc_ifnb@meta.data$simulated_donors \%>\% unique \%>\% .[[1]]
pbmc_ifnb <- subset(pbmc_ifnb, simulated_donors != donor_remove)

# run sequential:
options(future.globals.maxSize = 20*1024^3)
future::plan(future::sequential())
tictoc::tic()
test_dge <- DGE_analysis(m = pbmc_ifnb@assays$RNA@counts,
                         md = pbmc_ifnb@meta.data,
                         cluster_col = seurat_annotations,
                         sample_col = simulated_donors,
                         group_col = stim,
                         group1 = "CTRL",
                         group2 = "STIM",
                         design = ~simulated_batches + stim,
                         balance_batches = TRUE,
                         batch_col = simulated_batches,
                         n_cells_normalize = 10000,
                         n_cells_min = 10,
                         min_n_samples_group = 3,
                         cell_name_col = cell_name,
                         exp_percentage = 1)
tictoc::toc()

# run in parallel:
options(future.globals.maxSize = 20*1024^3)
future::plan(future::multisession(workers = 10, gc = TRUE))
tictoc::tic()
test_dge <- DGE_analysis(m = pbmc_ifnb@assays$RNA@counts,
                         md = pbmc_ifnb@meta.data,
                         cluster_col = seurat_annotations,
                         sample_col = simulated_donors,
                         group_col = stim,
                         group1 = "CTRL",
                         group2 = "STIM",
                         design = ~simulated_batches + stim,
                         balance_batches = TRUE,
                         batch_col = simulated_batches,
                         n_cells_normalize = 10000,
                         n_cells_min = 10,
                         min_n_samples_group = 3,
                         cell_name_col = cell_name,
                         exp_percentage = 1)
tictoc::toc()
future:::ClusterRegistry("stop")
 }
}
