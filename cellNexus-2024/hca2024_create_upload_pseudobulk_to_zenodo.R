# This script generates pseudobulk SummarizedExperiment data for cellNexus samples 
#   sharing at least 15,000 genes. 
# The output is saved as AnnData format and uploaded to Zenodo, with the dataset DOI 
#   linked to the associated preprint article.
# cellNexus metadata version: hca2024_v2.4.0 pseudobulk atlas version: hca_2024/0.5.0

library(dplyr)
library(cellNexus)
library(zellkonverter)
library(tidySingleCellExperiment)
cache <- "/vast/scratch/users/shen.m/cellNexus"

metadata <- get_metadata(cloud_metadata = get_metadata_url("hca_2024"),
                         cache_directory = cache) |>
  keep_quality_cells(min_features = 30000)   # This threshold return samples sharing at least 15000 genes
census_metadata <- cellNexus:::get_census_metadata("2024-07-01")
con <- dbplyr::remote_con(metadata)
duckdb::duckdb_register_arrow(con, "census_metadata", census_metadata)

metadata <- metadata |>
  dplyr::left_join(tbl(con, "census_metadata") |> 
                     dplyr::select(observation_joinid, dataset_id, tissue, 
                                   self_reported_ethnicity, assay, disease))

se <- metadata |> get_pseudobulk(cache_directory = cache)

priority_cols <- c(".aggregated_cells", "sample_id", "dataset_id", 
                   "cell_type_unified_ensemble", "disease", 
                   "tissue", "age_days", "assay")

colData(se) <- se |> colData() |> as.data.frame() %>%
  select(all_of(priority_cols), sort(setdiff(names(.), priority_cols))) |> 
  mutate(.aggregated_cells = as.integer(.aggregated_cells), 
         across(where(is.factor), as.character)) |> DataFrame()

cols_dropped <- c(
  "run_from_cell_id", "cell_count", "default_embedding", "feature_count", 
  "filesize", "mean_genes_per_cell", "primary_cell_count", "suspension_type",
  "url", "x_approximate_distribution", "max_lt_10.min_lt_0.rounding_error"
)

se <- se |> select(-any_of(cols_dropped))

job::job({
  se |> writeH5AD("/vast/scratch/users/shen.m/cellNexus/hca_2024/pseudobulk_se.h5ad", compression = "gzip")
})

# Validate saved file
x = readH5AD("/vast/scratch/users/shen.m/cellNexus/hca_2024/pseudobulk_se.h5ad", reader = "R", use_hdf5 = T)
colData(x) |> dim()


file.copy("/vast/scratch/users/shen.m/cellNexus/hca_2024/pseudobulk_se.h5ad",
          "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/hca2024_pseudobulk_se.h5ad",
          overwrite = T)

# Git clone https://github.com/jhpoelen/zenodo-upload.git
Sys.setenv(ZENODO_TOKEN = Sys.getenv("ZENODO_TOKEN"))
system("echo $ZENODO_TOKEN")

job::job({system("/home/users/allstaff/shen.m/git_control/zenodo-upload/zenodo_upload.sh 22668580 /vast/scratch/users/shen.m/cellNexus/hca_2024/pseudobulk_se.h5ad -v")})
#job::job({system("/home/users/allstaff/shen.m/git_control/zenodo-upload/zenodo_upload.sh 22668580 /home/users/allstaff/shen.m/git_control/cellNexus_development/hca2024_create_upload_pseudobulk_to_zenodo.R -v")})

