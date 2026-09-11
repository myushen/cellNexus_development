library(dplyr)
library(tibble)
library(glue)
library(purrr)
library(stringr)
library(HPCell)
library(arrow)
library(targets)
library(crew)
library(crew.cluster)
library(duckdb)

# ── Paths (MODIFY HERE) ────────────────────────────────────────────────────────
directory <- "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2025-11-08/"
metadata_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/2025-11-08_census_samples_to_download.parquet"
sample_summary_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/sample_distribution_summary.parquet"
sample_tbl_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/updated_transform_sample_tbl_2025_Nov.parquet"

prep_store <- "/vast/scratch/users/shen.m/cellNexus/2025-11-08/step6_sample_tbl_prep_store"

# ── Targets pipeline 1: assemble and classify sample_tbl ──────────────────────
tar_script(
  {
    library(dplyr)
    library(glue)
    library(arrow)
    library(HPCell)
    library(targets)

    directory <- "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2025-11-08/"
    metadata_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/2025-11-08_census_samples_to_download.parquet"
    sample_summary_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/sample_distribution_summary.parquet"
    sample_tbl_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/updated_transform_sample_tbl_2025_Nov.parquet"
    prep_store <- "/vast/scratch/users/shen.m/cellNexus/2025-11-08/step6_sample_tbl_prep_store"
    # Samples should be excluded in STEP_5 (no observation_joinid in sample colData, unable to slice from dataset)
    failed_samples <- c(
      "1d507aba40ec4e7307749428d4bc5c49", "221407d6715567ce06f8ec2efcf271d2", "24a0ab92f770ce82fc4bdbcb13ae6981", "29e98fa31c0ba4ed1b1e80a3c5aadccd", "2cdc34efd4a9c2ed7f80d1467d86f668", "2d70a0f016494707fae85a71de749dc7", "2ff166a68fa0b1ef5d0d8f9a91ac5923", "4e1b4b0af1998a6be416ab683a079764", "6081de2a72b2ea4bf7be8546aaf1d28c", "82e03a6cf864facd80211f76fbb1c85d", "8922d56becfeed3dc5919ce97b16c779", "9209879be3ca213fd2e73ebab31af71b", "9af2aafd1092aa19cbe10b2a63c8002c", "a42337b23e661eba31c7cf33e6e27e26", "a9f96f268e56e43adfa07f530bcbd268", "aca3b66fa8f730e4ba51af989b849be5", "ad739b85daca78093534a977699d8949", "b7251dabd2b20cc8f72b7c6e9344129f", "b806903916d5397c6f04335f5d4fa129", "c6fd0bb7191baf46f7d768dc86c00755", "c8d5fe854a5773a43ffdfa883b0a9cb0", "c95bb81d7e8d6d53251247233a01ce0e", "c96ca36a6fb86de531f23dabc9f015ca", "d0ba2ed2cb64609f200bd14dbf3cc36e", "d4cd0d2b0f603474faafe976de9a2492", "d58b0e1bccdb7ff916358b130fdf68fa", "d87d8f55a6a82b33ac23337fe74fa1e8", "d914e292ae1ee5cd01e5deca69193643", "e6195dc957219843cc470ee0b0960d93", "e6cd4914d5a9aa0bea4971b9b344911c", "ebac20fd555a13fa3d1e6ed638cb84c9", "ebdaf704fb553d7bd11becbe7ea6a7d3", "ef1780945cc1674f6eba884046f4fe6e", "f7ae602b0489415c1e4ce21981de6959"
    )
    list(
      tar_target(
        sample_summary_df,
        arrow::read_parquet(sample_summary_parquet),
        deployment = "main"
      ),
      tar_target(
        sample_tbl,
        {
          downloaded <- arrow::read_parquet(metadata_parquet) |>
            dplyr::rename(cell_number = list_length)

          tbl <- downloaded |>
            dplyr::filter(!sample_id %in% failed_samples) |>
            dplyr::left_join(
              cellxgenedp::datasets() |>
                dplyr::select(dataset_id, x_approximate_distribution) |>
                dplyr::distinct(),
              by = "dataset_id", copy = TRUE
            ) |>
            dplyr::mutate(
              cell_number    = as.integer(cell_number),
              file_name      = glue("{directory}{sample_id}.h5ad") |> as.character(),
              feature_thresh = ifelse(assay == "BD Rhapsody Targeted mRNA", 11, 200)
            )

          sample_summary_classified <- sample_summary_df |>
            HPCell::impute_x_approximate_distribution(
              counts_gap_threshold = 0.25,
              pos_mode_threshold   = 1
            ) |>
            dplyr::mutate(
              count_upper_bound = 10,
              method_to_apply = dplyr::case_when(
                inferred_distribution == "double_log1p" ~ "safe_expm1",
                inferred_distribution == "log1p" ~ "expm1",
                inferred_distribution == "log_negative_max_10" ~ "exp",
                inferred_distribution %in% c("raw", "raw_scaled", "raw_negative_scaled") ~ "identity"
              )
            )

          tbl |>
            dplyr::left_join(
              sample_summary_classified |>
                dplyr::mutate(sample_id = stringr::str_remove(sample_id, ".h5ad")) |>
                dplyr::select(sample_id, method_to_apply, dataset_id, count_upper_bound, inferred_distribution),
              by = c("sample_id", "dataset_id")
            ) |>
            dplyr::filter(!sample_id %in% failed_samples) |>
            dplyr::select(
              file_name, cell_number, dataset_id, sample_id, inferred_distribution,
              method_to_apply, assay, count_upper_bound, feature_thresh
            )
        },
        deployment = "main"
      ),
      tar_target(
        sample_tbl_parquet_file,
        {
          arrow::write_parquet(sample_tbl, sample_tbl_parquet)
          sample_tbl_parquet
        },
        format = "file",
        deployment = "main"
      )
    )
  },
  ask = FALSE,
  script = glue("{prep_store}/_targets.R")
)

job::job({
  tar_make(
    reporter = "summary",
    script   = glue("{prep_store}/_targets.R"),
    store    = glue("{prep_store}/_targets")
  )
})

# ── Read assembled sample_tbl for HPCell ──────────────────────────────────────
sample_tbl <- read_parquet(sample_tbl_parquet)
sample_names <- sample_tbl |>
  pull(file_name) |>
  set_names(sample_tbl |> pull(sample_id))
functions <- sample_tbl |> pull(method_to_apply)
feature_thresh <- sample_tbl |> pull(feature_thresh)
count_upper_bound <- sample_tbl |> pull(count_upper_bound)


my_store <- "/vast/scratch/users/shen.m/cellNexus_target_store_2025-11-08" # MODIFY HERE: HPCell targets store (used throughout this script)

new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 1, backup = NULL) {
  crew_controller_slurm(
    name = name,
    workers = workers,
    crashes_max = crashes_max,
    seconds_idle = 30,
    options_cluster = crew_options_slurm(
      memory_gigabytes_required = mem_gb,
      cpus_per_task = cpus_per_task,
      time_minutes = time_min
    ),
    backup = backup
  )
}
elastic_300 <- new_elastic("elastic_300", 300, 60 * 24, workers = 8, crashes_max = 2)
elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 10, crashes_max = 2)
elastic_120 <- new_elastic("elastic_120", 120, 60 * 4, workers = 24, crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
elastic_80 <- new_elastic("elastic_80", 80, 60 * 4, workers = 35, crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
elastic_40 <- new_elastic("elastic_40", 40, 60 * 4, workers = 70, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
elastic_20 <- new_elastic("elastic_20", 20, 60 * 4, workers = 140, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
elastic_10 <- new_elastic("elastic_10", 10, 60 * 4, workers = 290, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
elastic_5_minimal <- new_elastic("elastic_5_minimal", 5, 60 * 4, workers = 440, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)

controllers <- crew_controller_group(
  elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_300, elastic_5_minimal
)

job::job({
  library(HPCell)

  sample_names |>
    initialise_hpc(
      store = my_store,
      gene_nomenclature = "ensembl",
      data_container_type = "anndata",
      computing_resources = list(
        elastic_5_minimal, elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_300
      ),
      default_controller = "elastic_20",
      verbosity = "summary",
      update = "never",
      # update = "thorough",
      error = "continue",
      garbage_collection = 100,
      workspace_on_error = TRUE
    ) |>
    transform_assay(fx = functions, target_output = "sce_transformed", scale_max = count_upper_bound) |>
    # Sanity-check non-sensical expression samples (flags max_lt_10, min_lt_0, rounding_error)
    sanity_check_transform_samples(target_input = "sce_transformed", target_output = "sanity_checked_sample_stats") |>
    # Remove empty outliers based on RNA count threshold per cell
    remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = feature_thresh) |>
    # Annotation
    annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |>
    # Cell type harmonisation
    celltype_consensus_constructor(
      target_input = "sce_transformed",
      target_output = "cell_type_concensus_tbl"
    ) |>
    # Alive identification
    remove_dead_scuttle(
      target_input = "sce_transformed", target_annotation = "cell_type_concensus_tbl",
      group_by = "cell_type_unified_ensemble"
    ) |>
    # Doublets identification
    remove_doublets_scDblFinder(target_input = "sce_transformed") |>
    # SCT
    normalise_abundance_seurat_SCT(target_input = "sce_transformed", factors_to_regress = c(
      "subsets_Mito_percent",
      "subsets_Ribo_percent"
    )) |>
    # Pseudobulk
    calculate_pseudobulk(
      target_input = "sce_transformed",
      group_by = "cell_type_unified_ensemble"
    ) |>
    # # metacell
    # cluster_metacell(target_input = "sce_transformed",  group_by = "cell_type_unified_ensemble") |>
    #
    # # Cell Chat
    # ligand_receptor_cellchat(target_input = "sce_transformed",
    #                          group_by = "cell_type_unified_ensemble") |>

    print()
})


# Sample metadata
# ── Targets pipeline 2: assemble cell-level annotation ────────────────────────
metadata_assembly_store <- "/vast/scratch/users/shen.m/cellNexus/2025-11-08/step6_cell_metadata_assembly_store"

tar_script(
  {
    library(dplyr)
    library(duckdb)
    library(targets)
    library(stringr)
    library(crew)
    library(crew.cluster)

    my_store <- "/vast/scratch/users/shen.m/cellNexus_target_store_2025-11-08" # MODIFY HERE: HPCell targets store (must match my_store above)
    cell_metadata_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet"
    cell_annotation_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_annotation.parquet"

    elastic_500 <- crew_controller_slurm(
      name = "elastic_500",
      workers = 2,
      crashes_max = 1,
      seconds_idle = 30,
      options_cluster = crew_options_slurm(
        memory_gigabytes_required = 500,
        cpus_per_task             = 1,
        time_minutes              = 60 * 24
      )
    )

    tar_option_set(
      memory             = "transient",
      garbage_collection = 100,
      error              = "continue",
      format             = "qs",
      controller         = crew_controller_group(elastic_500)
    )

    list(
      # Joins cell metadata with all HPCell outputs and writes the annotation parquet.
      # Keeps a single DuckDB connection alive across all copy=TRUE left_joins.
      tar_target(
        cell_annotation_parquet_file,
        {
          con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
          on.exit(DBI::dbDisconnect(con), add = TRUE)

          cell_metadata <- tbl(
            con,
            dplyr::sql(paste0("SELECT * FROM read_parquet('", cell_metadata_parquet, "')"))
          ) |>
            dplyr::mutate(cell_ = paste0(cell_, "___", dataset_id)) |>
            dplyr::select(
              cell_, observation_joinid, dplyr::contains("cell_type"), dataset_id,
              self_reported_ethnicity, tissue, donor_id, sample_id, is_primary_data, assay
            )

          empty_droplet <- tar_read(empty_tbl, store = my_store) |>
            dplyr::bind_rows() |>
            dplyr::rename(cell_ = .cell)

          alive_cells <- tar_read(alive_tbl, store = my_store) |>
            dplyr::bind_rows() |>
            dplyr::select(-dplyr::any_of(c("cell_type_unified_ensemble", "observation_originalid"))) |>
            dplyr::rename(cell_ = .cell)

          doublet_cells <- tar_read(doublet_tbl, store = my_store) |>
            dplyr::bind_rows() |>
            dplyr::rename(cell_ = .cell)

          cell_type_concensus_tbl <- tar_read(cell_type_concensus_tbl, store = my_store) |>
            dplyr::bind_rows() |>
            dplyr::rename(cell_ = .cell) |>
            dplyr::mutate(cell_type_unified_ensemble = ifelse(
              is.na(cell_type_unified_ensemble), "Unknown", cell_type_unified_ensemble
            ))

          cell_metadata_joined <- cell_metadata |>
            dplyr::left_join(empty_droplet, copy = TRUE) |>
            dplyr::left_join(cell_type_concensus_tbl, copy = TRUE) |>
            dplyr::left_join(alive_cells, copy = TRUE) |>
            dplyr::left_join(doublet_cells, copy = TRUE)
          # |>
          #   dplyr::left_join(metacell, copy = TRUE)

          cell_metadata_labels_filled <- cell_metadata_joined |>
            dplyr::mutate(
              cell_type_unified_ensemble    = dplyr::coalesce(cell_type_unified_ensemble, "Unknown"),
              data_driven_ensemble          = dplyr::coalesce(data_driven_ensemble, "Unknown"),
              blueprint_first_labels_fine   = dplyr::coalesce(blueprint_first_labels_fine, "Other"),
              monaco_first_labels_fine      = dplyr::coalesce(monaco_first_labels_fine, "Other"),
              azimuth_predicted_celltype_l2 = dplyr::coalesce(azimuth_predicted_celltype_l2, "Other"),
              azimuth                       = dplyr::coalesce(azimuth, "Other"),
              blueprint                     = dplyr::coalesce(blueprint, "Other"),
              monaco                        = dplyr::coalesce(monaco, "Other")
            )

          final_sql <- dbplyr::sql_render(cell_metadata_labels_filled)
          DBI::dbExecute(con, sprintf(
            "COPY (%s) TO '%s' (FORMAT PARQUET, COMPRESSION 'zstd')",
            final_sql, cell_annotation_parquet
          ))
          cell_annotation_parquet
        },
        format = "file",
        resources = tar_resources(
          crew = tar_resources_crew(controller = "elastic_500")
        )
      )
    )
  },
  ask = FALSE,
  script = glue("{metadata_assembly_store}/_targets.R")
)

job::job({
  tar_make(
    reporter = "summary",
    script   = glue("{metadata_assembly_store}/_targets.R"),
    store    = glue("{metadata_assembly_store}/_targets")
  )
})
