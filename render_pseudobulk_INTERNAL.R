# This script generates pseudobulk SummarizedExperiment data for cellNexus hca_2024 and hca_2025 samples.
# cellNexus metadata version: metadata.v2024.2.3.1.parquet and metadata.v2025.1.0.1.parquet.
# pseudobulk atlas version: hca_2024/0.4.1 and hca_2025/0.1.1

library(dplyr)
library(cellNexus)
library(zellkonverter)
library(tidySingleCellExperiment)
cache <- "/vast/scratch/users/shen.m/cellNexus"

metadata <- get_metadata(cloud_metadata = get_metadata_url(c("metadata.v2024.2.3.1.parquet",
                                                           "metadata.v2025.1.0.1.parquet")),
                         cache_directory = cache) |>
  keep_quality_cells()

metadata <- metadata |>
  # This threshold return samples sharing at least 15000 genes
  dplyr::filter(feature_count >= 30000)

cols_dropped <- c(
  "run_from_cell_id", "cell_count", "default_embedding", "feature_count", 
  "filesize", "mean_genes_per_cell", "primary_cell_count", "suspension_type",
  "url", "x_approximate_distribution", "tissue_groups", "suspension_type", "tissue_type"
)

metadata <- metadata |> select(-any_of(cols_dropped))

pb <- metadata |> get_pseudobulk(cache_directory = cache)
# pb |> dim()
# 12956 333218

priority_cols <- c(".aggregated_cells", "sample_id", "dataset_id", 
                   "cell_type_unified_ensemble", "disease", 
                   "tissue", "age_days", "assay")

colData(pb) <- pb |> colData() |> as.data.frame() %>%
  select(all_of(priority_cols), sort(setdiff(names(.), priority_cols))) |> 
  mutate(.aggregated_cells = as.integer(.aggregated_cells), 
         across(where(is.factor), as.character)) |> DataFrame()

colData(pb)$published_at <- as.character(colData(pb)$published_at)
colData(pb)$revised_at   <- as.character(colData(pb)$revised_at)

pb |> writeH5AD("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/hca_2024_2025_pseudobulk_se.h5ad", compression = "gzip")
