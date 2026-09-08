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
directory              <- "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01/"
metadata_parquet       <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/census_samples_to_download_groups_MODIFIED.parquet"
sample_summary_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/sample_distribution_summary.parquet"
sample_tbl_parquet     <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/updated_transform_sample_tbl_2024_Jul.parquet"
cellnexus_cache        <- "/vast/scratch/users/shen.m/cellNexus"
prep_store             <- "/vast/scratch/users/shen.m/cellNexus/2024-07-01/step6_sample_tbl_prep_store"

# ── Targets pipeline: assemble and classify sample_tbl ────────────────────────
tar_script({
  library(dplyr)
  library(glue)
  library(arrow)
  library(HPCell)
  library(targets)
  
  directory              <- "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01/"
  metadata_parquet       <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/census_samples_to_download_groups_MODIFIED.parquet"
  sample_summary_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/sample_distribution_summary.parquet"
  sample_tbl_parquet     <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/updated_transform_sample_tbl_2024_Jul.parquet"
  cellnexus_cache        <- "/vast/scratch/users/shen.m/cellNexus"
  
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
          dplyr::rename(cell_number = list_length) |>
          dplyr::mutate(
            cell_number = as.integer(cell_number),
            file_name   = glue("{directory}{sample_2}.h5ad") |> as.character()
          )
        
        tbl <- downloaded |>
          dplyr::filter(!dataset_id %in% c(
            "99950e99-2758-41d2-b2c9-643edcdf6d82",
            "9fcb0b73-c734-40a5-be9c-ace7eea401c9"
          )) |>
          dplyr::left_join(
            cellxgenedp::datasets() |>
              dplyr::select(dataset_id, x_approximate_distribution) |>
              dplyr::distinct(),
            by = "dataset_id", copy = TRUE
          ) |>
          dplyr::mutate(
            cell_number = as.integer(cell_number),
            file_name   = glue("{directory}{sample_2}.h5ad") |> as.character()
          ) |>
          dplyr::left_join(
            cellNexus::get_metadata(
              cellNexus::get_metadata_url("metadata.v2024.2.3.1.parquet"),
              cache_directory = cellnexus_cache
            ) |> dplyr::distinct(sample_id, assay),
            by = c("sample_2" = "sample_id"), copy = TRUE
          ) |>
          # Propositional set up expressed genes threshold for panel technologies
          dplyr::mutate(feature_thresh = ifelse(assay == "BD Rhapsody Targeted mRNA", 11, 200))
        
        sample_summary_classified <- sample_summary_df |>
          HPCell::impute_x_approximate_distribution(
            counts_gap_threshold = 0.25,
            pos_mode_threshold   = 1
          ) |>
          dplyr::mutate(
            count_upper_bound = 10,
            method_to_apply   = dplyr::case_when(
              inferred_distribution == "double_log1p" ~ "safe_expm1",
              inferred_distribution == "log1p" ~ "expm1",
              inferred_distribution == "log_negative_max_10" ~ "exp",
              inferred_distribution %in% c("raw", "raw_scaled", "raw_negative_scaled") ~ "identity"
            )
          )
        
        tbl |>
          dplyr::left_join(
            sample_summary_classified |>
              dplyr::select(sample_2, method_to_apply, dataset_id, count_upper_bound, inferred_distribution),
            by = c("sample_2", "dataset_id")
          ) |>
          dplyr::select(file_name, cell_number, dataset_id, sample_2, inferred_distribution,
                        method_to_apply, assay, count_upper_bound, feature_thresh)
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
  
}, ask = FALSE, script = glue("{prep_store}/_targets.R"))

job::job({
  tar_make(
    reporter = "summary",
    script   = glue("{prep_store}/_targets.R"),
    store    = glue("{prep_store}/_targets")
  )
})

# ── Read assembled sample_tbl for HPCell ──────────────────────────────────────
sample_tbl        <- read_parquet(sample_tbl_parquet)
sample_names      <- sample_tbl |> pull(file_name) |> set_names(sample_tbl |> pull(sample_2))
functions         <- sample_tbl |> pull(method_to_apply)
feature_thresh    <- sample_tbl |> pull(feature_thresh)
count_upper_bound <- sample_tbl |> pull(count_upper_bound)

my_store = "/vast/scratch/users/shen.m/cellNexus/2024-07-01/process_updated_samples_transform_hpcell_target_store_v1" # MODIFY HERE: HPCell targets store (used throughout this script)

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
elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 10,  crashes_max = 2)
elastic_120  <- new_elastic("elastic_120",  120,  60 * 4,  workers = 24, crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
elastic_80  <- new_elastic("elastic_80",   80,  60 * 4,  workers = 35, crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 70, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 140, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 290, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)

elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 440, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)

# Group for targets (small → large)
controllers <- crew_controller_group(
  elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_5_minimal
)

job::job({
  
  library(HPCell)
  
  sample_names |>
    initialise_hpc(
      store = my_store,
      gene_nomenclature = "ensembl",
      data_container_type = "anndata",
      computing_resources = list(
        elastic_5_minimal, elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160
      ),
      default_controller = "elastic_5_minimal", 
      verbosity = "summary",
      update = "never", 
      #update = "thorough", 
      error = "continue",
      garbage_collection = 100, 
      workspace_on_error = TRUE
      
    ) |> 
    transform_assay(fx = functions, target_output = "sce_transformed", scale_max = count_upper_bound) |>
    
    # sanity_check_ non sensical expression samples
    sanity_check_transform_samples(target_input = "sce_transformed", target_output = "sanity_checked_sample_stats") |>
    
    # # Remove empty outliers based on RNA count threshold per cell
    remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = feature_thresh) |>

    # Annotation
    annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |>

    # Cell type harmonisation
    celltype_consensus_constructor(target_input = "sce_transformed",
                                   target_output = "cell_type_concensus_tbl") |>

    # Alive identification
    remove_dead_scuttle(target_input = "sce_transformed", target_annotation = "cell_type_concensus_tbl",
                        group_by = "cell_type_unified_ensemble") |>

    # Doublets identification
    remove_doublets_scDblFinder(target_input = "sce_transformed") |>

    # SCT
    normalise_abundance_seurat_SCT(target_input = "sce_transformed", factors_to_regress = c(
      "subsets_Mito_percent",
      "subsets_Ribo_percent")) |>

    # Pseudobulk
    calculate_pseudobulk(target_input = "sce_transformed",
                         group_by = "cell_type_unified_ensemble") |>
    
    # # metacell
    # cluster_metacell(target_input = "sce_transformed",  group_by = "cell_type_unified_ensemble") |>
    
    # # Cell Chat
    # ligand_receptor_cellchat(target_input = "sce_transformed",
    #                          group_by = "cell_type_unified_ensemble") |>
    
    print()
  
  
})

#' # # View target metadata if needed
#' # tar_meta(store = my_store) |> filter(!is.na(error)) |>  arrange(desc(time)) |> View()
#' # tar_meta(store = my_store) |> filter(!is.na(error)) |> distinct(name, error)
#' # tar_meta(starts_with("annotation_tbl_"), store = "/vast/scratch/users/shen.m/cellNexus_target_store") |> 
#' #   filter(!data |> is.na()) |> arrange(desc(time)) |> select(error, name)
#' # 
#' # # Debug cellchat
#' # tar_workspace(ligand_receptor_tbl_f1dcd76261dd9c86, store = my_store, 
#' #               script = paste0(my_store,".R"))
#' # debugonce(cell_communication)
#' # cell_communication(sce_transformed, empty_droplets_tbl = empty_tbl, alive_identification_tbl = alive_tbl,
#' #                   doublet_identification_tbl = doublet_tbl, cell_type_tbl = cell_type_concensus_tbl,
#' #                   cell_type_column = "cell_type_unified_ensemble",
#' #                   feature_nomenclature = gene_nomenclature)
#' 
#' 
#' #' Pipeline for Lightening Annotations in High-Performance Computing Environment
#' #' 
#' #' This pipeline is designed to read, process, and "lighten" large annotation tables in an HPC environment.
#' #' It uses the `targets` package for reproducibility and `crew` for efficient job scheduling on a Slurm cluster.
#' #' The `lighten_annotation` function selects and processes specific columns from large tables to reduce memory usage.
#' #' 
#' #' The pipeline consists of:
#' #' - **Crew Controllers**: Four tiers of Slurm controllers with varying memory allocations to optimize resource usage.
#' #' - **Targets**:
#' #'   - `my_store`: Defines the path to the target storage directory, ensuring all targets use the correct storage location.
#' #'   - `target_name`: Retrieves metadata to identify branch targets for annotation.
#' #'   - `annotation_tbl_light`: Applies `lighten_annotation` to process each target name, optimally running with `tier_1` resources.
#' #' 
#' #' @libraries:
#' #'   - `dplyr`, `magrittr`, `tibble`, `targets`, `tarchetypes` for data manipulation and pipeline structure.
#' #'   - `crew`, `crew.cluster` for parallel computation and cluster scheduling in an HPC environment.
#' #' 
#' #' @options:
#' #'   - Memory settings, garbage collection frequency, and error handling are set to handle large data efficiently.
#' #'   - The `cue` option is set to `never` for forced target updates if needed.
#' #'   - `controller` is a group of Slurm controllers to manage computation across memory tiers.
#' #' 
#' #' @function `lighten_annotation`: Processes each annotation table target, unnesting and selecting specific columns to reduce data size.
#' #'
#' #' @example Usage:
#' #'   The pipeline script is saved as `/vast/scratch/users/shen.m/lighten_annotation_tbl_target.R` by tar_script and can be run using `tar_make()`.
#' tar_script({
#'   library(dplyr)
#'   library(magrittr)
#'   library(tibble)
#'   library(targets)
#'   library(tarchetypes)
#'   library(crew)
#'   library(crew.cluster)
#'   
#'   
#'   # Helper (optional) to avoid repetition
#'   new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 2, backup = NULL) {
#'     crew_controller_slurm(
#'       name = name,
#'       workers = workers,
#'       crashes_max = crashes_max,
#'       seconds_idle = 30,
#'       options_cluster = crew_options_slurm(
#'         memory_gigabytes_required = mem_gb,
#'         cpus_per_task = cpus_per_task,
#'         time_minutes = time_min
#'       ),
#'       backup = backup
#'     )
#'   }
#'   
#'   # Small → large, with fallbacks to the next size up
#'   elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2)
#'   elastic_120  <- new_elastic("elastic_120",  120,  60 * 4,  workers = 16, crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
#'   elastic_80  <- new_elastic("elastic_80",   80,  60 * 4,  workers = 24, crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
#'   elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 32, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
#'   elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 48, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
#'   elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
#'   
#'   elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)
#'   
#'   # Group for targets (small → large)
#'   controllers <- crew_controller_group(
#'     elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_5_minimal
#'   )
#'   
#'   tar_option_set(
#'     memory = "transient", 
#'     garbage_collection = 100, 
#'     storage = "worker", 
#'     retrieval = "worker", 
#'     error = "continue", 
#'     cue = tar_cue(mode = "never"), 
#'     resources = tar_resources(
#'       crew = tar_resources_crew(controller = "elastic_5_minimal")
#'     ),
#'     controller = controllers, 
#'     trust_object_timestamps = TRUE
#'   )
#'   
#'   lighten_annotation = function(target_name, my_store ){
#'     annotation_tbl = tar_read_raw( target_name,  store = my_store )
#'     if(annotation_tbl |> is.null()) { 
#'       warning("this annotation is null -> ", target_name)
#'       return(NULL) 
#'     }
#'     
#'     annotation_tbl |> 
#'       unnest(blueprint_scores_fine) |> 
#'       select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), monaco_scores_fine, contains("macro"), contains("CD4") ) |> 
#'       unnest(monaco_scores_fine) |> 
#'       select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), contains("macro") , contains("CD4"), contains("helper"), contains("Th")) |> 
#'       rename(cell_ = .cell)
#'   }
#'   
#'   list(
#'     
#'     # The input DO NOT DELETE
#'     tar_target(my_store, "/vast/scratch/users/shen.m/cellNexus/2024-07-01/process_updated_samples_transform_hpcell_target_store_v1", deployment = "main"), # MODIFY HERE: HPCell targets store (must match my_store above)
#'     
#'     tar_target(
#'       target_name,
#'       tar_meta(
#'         starts_with("annotation_tbl_"), 
#'         store = my_store) |> 
#'         filter(type=="branch") |> 
#'         pull(name),
#'       deployment = "main"
#'     )    ,
#'     
#'     tar_target(
#'       annotation_tbl_light,
#'       lighten_annotation(target_name, my_store),
#'       packages = c("dplyr", "tidyr"),
#'       pattern = map(target_name),
#'       resources = tar_resources(
#'         crew = tar_resources_crew(controller = "elastic_5_minimal")
#'       )
#'     )
#'   )
#'   
#'   
#' }, script = "/vast/scratch/users/shen.m/lighten_annotation_tbl_target_2024_Jul.R", ask = FALSE) # MODIFY HERE: output path for the tar_script file
#' 
#' job::job({
#'   
#'   tar_make(
#'     script = "/vast/scratch/users/shen.m/lighten_annotation_tbl_target_2024_Jul.R", # MODIFY HERE: must match the script path above
#'     store = "/vast/scratch/users/shen.m/lighten_annotation_tbl_target_2024_Jul", # MODIFY HERE: targets store for the lighten-annotation pipeline
#'     reporter = "summary"
#'   )
#'   
#' })

# Sample metadata
# ── Paths (MODIFY HERE) ────────────────────────────────────────────────────────
cell_metadata_parquet   <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet"
cell_annotation_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_annotation_2024_Jul.parquet"
lr_duckdb               <- "/vast/scratch/users/shen.m/cellNexus/cellNexus_lr_signaling_pathway_strength.duckdb"
metadata_assembly_store <- "/vast/scratch/users/shen.m/cellNexus/2024-07-01/step6_cell_metadata_assembly_store"

# ── Targets pipeline: assemble cell-level annotation ─────────────────────────
tar_script({
  library(dplyr)
  library(duckdb)
  library(targets)
  library(stringr)
  library(crew)
  library(crew.cluster)
  
  my_store                <- "/vast/scratch/users/shen.m/cellNexus/2024-07-01/process_updated_samples_transform_hpcell_target_store_v1" # MODIFY HERE: HPCell targets store (must match my_store above)
  cell_metadata_parquet   <- "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet"
  cell_annotation_parquet <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_annotation_2024_Jul.parquet"
  #lr_duckdb               <- "/vast/scratch/users/shen.m/cellNexus/cellNexus_lr_signaling_pathway_strength.duckdb"
  
  elastic_500 <- crew_controller_slurm(
    name         = "elastic_500",
    workers      = 2,
    crashes_max  = 1,
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
          dplyr::select(cell_, observation_joinid, dplyr::contains("cell_type"), dataset_id,
                        self_reported_ethnicity, tissue, donor_id, sample_id, is_primary_data, assay)
        
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
          dplyr::left_join(empty_droplet,           copy = TRUE) |>
          dplyr::left_join(cell_type_concensus_tbl, copy = TRUE) |>
          dplyr::left_join(alive_cells,             copy = TRUE) |>
          dplyr::left_join(doublet_cells,           copy = TRUE)
        # |>
        #   dplyr::left_join(metacell, copy = TRUE)
        
        cell_metadata_joined2 <- cell_metadata_joined |>
          dplyr::mutate(
            cell_type_unified_ensemble    = dplyr::coalesce(cell_type_unified_ensemble,    "Unknown"),
            data_driven_ensemble          = dplyr::coalesce(data_driven_ensemble,          "Unknown"),
            blueprint_first_labels_fine   = dplyr::coalesce(blueprint_first_labels_fine,   "Other"),
            monaco_first_labels_fine      = dplyr::coalesce(monaco_first_labels_fine,      "Other"),
            azimuth_predicted_celltype_l2 = dplyr::coalesce(azimuth_predicted_celltype_l2, "Other"),
            azimuth                       = dplyr::coalesce(azimuth,                       "Other"),
            blueprint                     = dplyr::coalesce(blueprint,                     "Other"),
            monaco                        = dplyr::coalesce(monaco,                        "Other")
          )
        
        final_sql <- dbplyr::sql_render(cell_metadata_joined2)
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
    ),
    
    # Cellchat output
    tar_target(
      lr_duckdb_file,
      {
        ligand_receptor_tbl <- tar_read(ligand_receptor_tbl, store = my_store) |>
          dplyr::bind_rows()
        con <- DBI::dbConnect(duckdb::duckdb(), dbdir = lr_duckdb)
        on.exit(DBI::dbDisconnect(con), add = TRUE)
        duckdb::dbWriteTable(con, "lr_pathway_table", ligand_receptor_tbl, overwrite = TRUE)
        lr_duckdb
      },
      format = "file",
      deployment = "main"
    )
  )
  
}, ask = FALSE, script = glue("{metadata_assembly_store}/_targets.R"))

job::job({
  tar_make(
    reporter = "summary",
    script   = glue("{metadata_assembly_store}/_targets.R"),
    store    = glue("{metadata_assembly_store}/_targets")
  )
})

