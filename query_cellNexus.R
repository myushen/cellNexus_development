# Build a combined cell-metadata table from HCA, the prostate atlas, and HTA.
# query_quality_hca_primary_data() joins HCA 2024 and 2025 metadata to
# CELLxGENE Census obs and dataset fields (2025 only for datasets absent
# from 2024), keeps quality cells, and drops two samples with conflicting sex.
# query_prostate(), query_hta_lung(), and query_hta_full_TEMPORARY() load the
# other atlases. The script unions those tables onto HCA, keeps the HCA
# columns, then checks that each sample_id has one sex and one ethnicity.

library(cellNexus)
library(dplyr)

query_quality_hca_primary_data <- function() {
  metadata_tbl <- get_metadata(
    cloud_metadata = get_metadata_url(c("hca_2024", "hca_2025"))
  )

  con <- dbplyr::remote_con(metadata_tbl)

  read_census <- function(version) {
    census <- cellxgene.census::open_soma(census_version = version)
    on.exit(census$close(), add = TRUE)

    list(
      obs = census$
        get("census_data")$
        get("homo_sapiens")$
        get("obs")$
        read(
          column_names = c(
            "observation_joinid",
            "dataset_id",
            "sex",
            "tissue",
            "assay",
            "disease",
            "self_reported_ethnicity",
            "cell_type",
            "is_primary_data"
          ),
          value_filter = "is_primary_data == 'TRUE'"
        )$
        concat(),

      datasets = census$
        get("census_info")$
        get("datasets")$
        read(
          column_names = c(
            "dataset_id",
            "dataset_title",
            "collection_id"
          )
        )$
        concat()
    )
  }

  census_2024 <- read_census("2024-07-01")
  census_2025 <- read_census("2025-11-08")

  duckdb::duckdb_register_arrow(con, "census_obs_2024", census_2024$obs)
  duckdb::duckdb_register_arrow(con, "census_datasets_2024", census_2024$datasets)
  duckdb::duckdb_register_arrow(con, "census_obs_2025", census_2025$obs)
  duckdb::duckdb_register_arrow(con, "census_datasets_2025", census_2025$datasets)

  # All datasets already included in hca_2024.
  old_dataset_ids <- dplyr::tbl(con, "census_datasets_2024") |>
    dplyr::select(dataset_id) |>
    dplyr::distinct()

  # Only datasets genuinely added in hca_2025.
  new_dataset_ids <- dplyr::tbl(con, "census_datasets_2025") |>
    dplyr::select(dataset_id) |>
    dplyr::distinct() |>
    dplyr::anti_join(old_dataset_ids, by = "dataset_id")

  # The two sides of this union are now disjoint by dataset_id.
  census_obs_lookup <- dplyr::union_all(
    dplyr::tbl(con, "census_obs_2024") |>
      dplyr::select(-is_primary_data) |>
      dplyr::distinct(
        dataset_id,
        observation_joinid,
        .keep_all = TRUE
      ),

    dplyr::tbl(con, "census_obs_2025") |>
      dplyr::semi_join(new_dataset_ids, by = "dataset_id") |>
      dplyr::select(-is_primary_data) |>
      dplyr::distinct(
        dataset_id,
        observation_joinid,
        .keep_all = TRUE
      )
  )

  census_dataset_lookup <- dplyr::union_all(
    dplyr::tbl(con, "census_datasets_2024") |>
      dplyr::distinct(dataset_id, .keep_all = TRUE),

    dplyr::tbl(con, "census_datasets_2025") |>
      dplyr::semi_join(new_dataset_ids, by = "dataset_id") |>
      dplyr::distinct(dataset_id, .keep_all = TRUE)
  ) |>
    dplyr::rename(title = dataset_title)

  metadata_tbl |>
    dplyr::left_join(
      census_obs_lookup,
      by = c("dataset_id", "observation_joinid")
    ) |>
    dplyr::left_join(
      census_dataset_lookup,
      by = "dataset_id"
    ) |>
    keep_quality_cells() |>
    # TEMPORARY REMOVE: these two sample match different sex for some reasons
    dplyr::filter(!sample_id %in% c(
      "e6b21291ca8ed20cb62541c2a1c7b2ac",
      "20b3bf1fa3f5f28eb62130efbe464f94"
    ))
}

query_wehi_prostate_primary_data <- function() {
  get_metadata(
    local_metadata = "/vast/projects/cellxgene_curated/prostate_atlas/prostate_metadata.v0.3.3.parquet",
    cloud_metadata = NULL
  ) |>
    dplyr::filter(is_primary_data == TRUE)
}

query_hta_lung <- function() {
  get_metadata(cloud_metadata = get_metadata_url("metadata_hta_lung.v0.1.0.parquet"))
}

query_hta_full_TEMPORARY <- function() {
  get_metadata(cloud_metadata = get_metadata_url("hta_2026.v0.2.0.parquet"))
}

my_tbl <- query_quality_hca_primary_data()
cols_to_keep <- colnames(my_tbl)

my_tbl <- my_tbl |>
  dplyr::union_all(query_wehi_prostate_primary_data(), copy = TRUE) |>
  dplyr::union_all(query_hta_lung(), copy = TRUE) |>
  dplyr::select(dplyr::all_of(cols_to_keep))

# Sanity check
my_tbl |> distinct(sample_id, self_reported_ethnicity) |> dplyr::count(sample_id) |> filter(n > 1) |> dplyr::count()
my_tbl |> distinct(sample_id, sex) |> dplyr::count(sample_id) |> filter(n > 1)
