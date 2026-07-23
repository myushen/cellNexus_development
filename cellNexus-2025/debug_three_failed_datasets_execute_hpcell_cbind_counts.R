library(dplyr)
library(duckdb)
ds <- c("21722308-5091-4b63-9c07-4f116ab9a7b3", "145dcf6a-2461-4fa3-a0af-7fc56db0bd33", "7f7faf6b-f11d-4f07-bc1c-188a4472748d")
tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/2025-11-08_census_samples_to_download.parquet')")
) |> filter(dataset_id %in% ds) |>
  distinct(sample_id)

# Get sample counts summary -----------------------------------------------
library(targets)
library(dplyr)
library(stringr)
library(glue)
library(arrow)
library(tidySingleCellExperiment)
set.seed(12345)
# Identify raw counts range from new census data
summary_store = "/vast/scratch/users/shen.m/debug_three_failed_dataset_target_store"
tar_script({
  library(dplyr)
  library(SummarizedExperiment)
  library(zellkonverter)
  library(crew)
  library(crew.cluster)
  library(duckdb)
  
  # Helper to avoid repetition
  new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 2, backup = NULL) {
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
  elastic_500      <- new_elastic("elastic_500",      500, 60 * 24, workers = 3,   crashes_max = 2)
  elastic_300      <- new_elastic("elastic_300",      300, 60 * 24, workers = 8,   crashes_max = 2, backup = elastic_500)
  elastic_160      <- new_elastic("elastic_160",      160, 60 * 24, workers = 8,   crashes_max = 2, backup = elastic_300)
  elastic_120      <- new_elastic("elastic_120",      120, 60 * 4,  workers = 16,  crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
  elastic_80       <- new_elastic("elastic_80",        80, 60 * 4,  workers = 24,  crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
  elastic_40       <- new_elastic("elastic_40",        40, 60 * 4,  workers = 32,  crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
  elastic_20       <- new_elastic("elastic_20",        20, 60 * 4,  workers = 48,  crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
  elastic_10       <- new_elastic("elastic_10",        10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
  elastic_5_minimal <- new_elastic("elastic_5_minimal", 5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)
  
  controllers <- crew_controller_group(
    elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_300, elastic_500, elastic_5_minimal
  )
  
  tar_option_set(
    memory             = "transient",
    garbage_collection = 100,
    storage            = "worker",
    retrieval          = "worker",
    error              = "continue",
    cue                = tar_cue(mode = "thorough"),
    format             = "qs",
    workspace_on_error = TRUE,
    controller         = controllers,
    trust_object_timestamps = TRUE,
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    )
  )
  
  # ── helpers ────────────────────────────────────────────────────────────────
  
  pos_min_med_ratio <- function(x) {
    pos <- x[x > 0]
    min(pos) / median(pos)
  }
  
  pos_min_mean_ratio <- function(x) {
    pos <- x[x > 0]
    min(pos) / mean(pos)
  }
  
  get_positive_mode <- function(x) {
    sort(table(x[x > 0]), decreasing = TRUE)[1] |> names() |> as.numeric()
  }
  
  # Stage 1 – read SCE and return a minimal list with only what downstream needs.
  # Keeping the raw matrix avoids re-reading the file in the metrics stage, while
  # returning a plain list (not a full SCE) keeps the serialised object small.
  read_sce_counts <- function(file) {
    sce        <- readH5AD(file, reader = "R", use_hdf5 = TRUE)
    if (ncol(sce) == 0) return(NULL)
    
    assay_name <- names(sce@assays)[1]
    counts_mat <- as.matrix(assay(sce, assay_name))
    
    list(
      sample_id  = basename(file),
      counts_mat = counts_mat,
      n_cells    = ncol(counts_mat),
      n_genes    = nrow(counts_mat)
    )
  }
  
  # Stage 2 – pure computation on the pre-loaded matrix; no I/O.
  calc_counts_metrics <- function(sce_data) {
    if (is.null(sce_data)) return(NULL)
    
    counts_vec <- as.numeric(sce_data$counts_mat)
    tol        <- 1e-4
    all_int    <- all(counts_vec == floor(counts_vec), na.rm = TRUE)
    
    tibble::tibble(
      sample_id          = sce_data$sample_id,
      min_val            = min(counts_vec,    na.rm = TRUE),
      median_val         = median(counts_vec, na.rm = TRUE),
      max_val            = max(counts_vec,    na.rm = TRUE),
      counts_gap_min_med = pos_min_med_ratio(counts_vec),
      counts_gap_min_mean= pos_min_mean_ratio(counts_vec),
      positive_mode      = get_positive_mode(counts_vec),
      has_negative       = min(counts_vec, na.rm = TRUE) < 0,
      max_gt_10          = max(counts_vec, na.rm = TRUE) > 10,
      all_integer        = all_int,
      has_floating       = !all_int && all(abs(counts_vec - round(counts_vec)) < tol, na.rm = TRUE),
      n_cells            = sce_data$n_cells,
      n_genes            = sce_data$n_genes
    ) |>
      left_join(
        tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/2025-11-08_census_samples_to_download.parquet')")
        ) |>
          select(-observation_joinid, list_length) |>
          distinct() |>
          mutate(sample_id = paste0(sample_id, ".h5ad")) |>
          collect(),
        by   = "sample_id",
        copy = TRUE
      )
  }
  
  # ── pipeline ───────────────────────────────────────────────────────────────
  
  list(
    tar_target(
      files,
      c("/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2025-11-08/f365e988ae4449b991b093d6e458466d.h5ad", 
        "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2025-11-08/6513fc618af9ab8a79fd05e0321f9033.h5ad", 
        "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2025-11-08/4152c6c633d769bc4e14d6b7022e568c.h5ad"),
      deployment = "main"
    ),
    
    # Stage 1 – I/O-bound; needs memory for the full matrix
    tar_target(
      sce_counts,
      read_sce_counts(files),
      pattern   = map(files),
      iteration = "list",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_20")
      )
    ),
    
    # Stage 2 – CPU-bound, no disk I/O; can run on smaller workers
    tar_target(
      sample_summary_df,
      calc_counts_metrics(sce_counts),
      pattern   = map(sce_counts),
      iteration = "list",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_20")
      )
    )
  )
  
}, ask = FALSE, script = glue("{summary_store}/_targets.R"))

job::job({
  
  tar_make(
    # callr_function = NULL,
    reporter = "summary",
    script = glue("{summary_store}/_targets.R"),
    store = glue("{summary_store}/_targets")
  )
  
})

sample_summary_df = tar_read(sample_summary_df,  store = glue("{summary_store}/_targets")) |> bind_rows() |> 
  mutate(max_gt_20 = ifelse(max_val > 20, TRUE, FALSE))
sample_summary_df

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

version <- "2025-11-08"
directory = glue::glue("/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/{version}/")

downloaded_samples_tbl <- 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/2025-11-08_census_samples_to_download.parquet')")
  ) |>
  dplyr::rename(cell_number = list_length) |>
  collect() |>
  mutate(file_name = file.path(directory, paste0(sample_id, ".h5ad")) |> as.character())

sample_meta <- tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/dataset.parquet')")
) |> collect()

sample_tbl <- downloaded_samples_tbl |> left_join(
  cellxgenedp::datasets() |>
    dplyr::select(dataset_id, x_approximate_distribution) |>
    distinct(), by = "dataset_id", copy = TRUE) |>
  mutate(cell_number = cell_number |> as.integer(),
         file_name = glue("{directory}{sample_id}.h5ad") |> as.character()) |>
  mutate(feature_thresh = ifelse(assay == "BD Rhapsody Targeted mRNA", 11, 200))


sample_summary_df = tar_read(sample_summary_df,  store = glue("{summary_store}/_targets")) |> bind_rows() |> 
  mutate(max_gt_20 = ifelse(max_val > 20, TRUE, FALSE))

# Impute distribution decision tree
impute_x_approximate_distribution <- function(df,
                                              counts_gap_threshold,
                                              pos_mode_threshold) {
  df |>
    dplyr::mutate(
      inferred_distribution = dplyr::case_when(
        
        # 0) When counts gap between 0 and next min value >= threshold
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          (counts_gap_min_mean >= counts_gap_threshold) & (positive_mode > pos_mode_threshold) ~ "double_log1p",
        
        # 1) Small counts gap
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          !(
            (counts_gap_min_mean >= counts_gap_threshold) &
              (positive_mode > pos_mode_threshold)
            
          ) ~ "log1p",
        
        # 2) No negatives, has large values
        !has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 3) Large values, integer counts
        !has_negative & max_gt_20 & all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 4) Has negatives, compressed range
        has_negative & !max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 5) Has negatives and large values
        has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # fallback
        TRUE ~ NA_character_
      )
    )
}

sample_summary_df = sample_summary_df |> impute_x_approximate_distribution(0.25, 1) |> 
  mutate(count_upper_bound = case_when(
    # 0) When counts gap between 0 and next min value >= 0.25, double log. Max value before exp is 10.
    inferred_distribution == "double_log1p" ~ 10,
    
    # 1) make 10 as max before exp
    inferred_distribution == "log1p" ~ 10,
    
    # 2,3,5) transform algo picks up negative value. should always scale max to 10
    # 4) Has negatives, no large values, no integer, no floating. Counts peak at 10
    inferred_distribution == "raw_limit_max_to_10" ~ 10
    
  )) |>
  # Inverse distribution
  mutate(method_to_apply = case_when(inferred_distribution == "double_log1p" ~ "safe_expm1",
                                     inferred_distribution == "log1p" ~ "expm1",
                                     inferred_distribution == "raw_limit_max_to_10" ~ "identity_with_max_limit"))

sample_tbl = sample_tbl |> left_join(sample_summary_df |>
                                       mutate(sample_id = sample_id |> str_remove(".h5ad"))|>
                                       select(sample_id,
                                              method_to_apply,
                                              dataset_id,
                                              count_upper_bound),
                                     by = c("sample_id", "dataset_id")) |>
  
  filter(dataset_id %in% x)

sample_tbl = sample_tbl |>
  
  select(file_name, cell_number, dataset_id, sample_id, method_to_apply, assay, count_upper_bound, feature_thresh)

sample_tbl |> write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/updated_transform_sample_tbl_2025_Nov_rest_three_dataset.parquet")
sample_tbl <- read_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/updated_transform_sample_tbl_2025_Nov_rest_three_dataset.parquet")  # MODIFY HERE: output path for sliced_sample_tbl RDS


# # Datasets missing colData are excluded. ??????
# dataset_to_exclude <- sample_tbl |> dplyr::filter(inferred_distribution |> is.na(), 
#                                                   is.na(x_approximate_distribution),  
#                                                   is.na(has_negative),  is.na(max_gt_20),  
#                                                   is.na(all_integer), is.na(has_floating)) |>
#   pull(dataset_id) |> unique()


sample_names <-
  sample_tbl |> 
  pull(file_name) |> 
  str_replace("/home/users/allstaff/shen.m/scratch", "/vast/scratch/users/shen.m") |> 
  set_names(sample_tbl |> pull(sample_id))
functions = sample_tbl |> pull(method_to_apply)
feature_thresh = sample_tbl |> pull(feature_thresh)
count_upper_bound = sample_tbl |> pull(count_upper_bound)


my_store = "/vast/scratch/users/shen.m/cellNexus_target_store_2025-11-08_rest_three_dataset"

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
elastic_300 <- new_elastic("elastic_300", 300, 60 * 24, workers = 4,  crashes_max = 2)
elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2, backup = elastic_300)
elastic_120  <- new_elastic("elastic_120",  120,  60 * 8,  workers = 16, crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
elastic_80  <- new_elastic("elastic_80",   80,  60 * 8,  workers = 24, crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 32, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 48, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)

elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)

# Group for targets (small → large)
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
      default_controller = "elastic_5_minimal", #"elastic_40", 
      verbosity = "summary",
      update = "never", 
      #update = "thorough", 
      error = "continue",
      garbage_collection = 100, 
      workspace_on_error = TRUE
      
    ) |> 
    transform_assay(fx = functions, target_output = "sce_transformed", scale_max = count_upper_bound) |>
    
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
    # 
    # # Cell Chat
    # ligand_receptor_cellchat(target_input = "sce_transformed",
    #                          group_by = "cell_type_unified_ensemble") |>
    
    print()
  
  
})

tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  # Helper (optional) to avoid repetition
  new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 2, backup = NULL) {
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
  
  # Small → large, with fallbacks to the next size up
  elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2)
  elastic_120  <- new_elastic("elastic_120",  120,  60 * 4,  workers = 16, crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
  elastic_80  <- new_elastic("elastic_80",   80,  60 * 4,  workers = 24, crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
  elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 32, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
  elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 48, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
  elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
  
  elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)
  
  # Group for targets (small → large)
  controllers <- crew_controller_group(
    elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_5_minimal
  )
  
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    cue = tar_cue(mode = "never"), 
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    ),
    controller = controllers, 
    trust_object_timestamps = TRUE
  )
  
  lighten_annotation = function(target_name, my_store ){
    annotation_tbl = tar_read_raw( target_name,  store = my_store )
    if(annotation_tbl |> is.null()) { 
      warning("this annotation is null -> ", target_name)
      return(NULL) 
    }
    
    annotation_tbl |> 
      unnest(blueprint_scores_fine) |> 
      select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), monaco_scores_fine, contains("macro"), contains("CD4") ) |> 
      unnest(monaco_scores_fine) |> 
      select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), contains("macro") , contains("CD4"), contains("helper"), contains("Th")) |> 
      rename(cell_ = .cell)
  }
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/scratch/users/shen.m/cellNexus_target_store_2025-11-08_rest_three_dataset", deployment = "main"),
    
    tar_target(
      target_name,
      tar_meta(
        starts_with("annotation_tbl_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    )    ,
    
    tar_target(
      annotation_tbl_light,
      lighten_annotation(target_name, my_store),
      packages = c("dplyr", "tidyr"),
      pattern = map(target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      )
    )
  )
  
  
}, script = glue::glue("/vast/scratch/users/shen.m/lighten_annotation_tbl_target_{version}_rest_three_dataset.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = glue::glue("/vast/scratch/users/shen.m/lighten_annotation_tbl_target_{version}_rest_three_dataset.R"),
    store = glue::glue("/vast/scratch/users/shen.m/lighten_annotation_tbl_target_{version}_rest_three_dataset"), 
    reporter = "summary"
  )
  
})

# Sample metadata
library(arrow)
library(dplyr)
library(duckdb)
library(targets)

# Write annotation light
cell_metadata <- 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet')")
  ) |>
  mutate(cell_ = paste0(cell_, "___", dataset_id)) |> 
  filter(dataset_id %in% ds) |>
  select(cell_, observation_joinid, contains("cell_type"), dataset_id,  self_reported_ethnicity, tissue, donor_id,  sample_id, is_primary_data, assay)


cell_annotation = 
  tar_read(annotation_tbl_light, store = glue::glue("/vast/scratch/users/shen.m/lighten_annotation_tbl_target_{version}_rest_three_dataset")) |> 
  dplyr::rename(
    blueprint_first_labels_fine = blueprint_first.labels.fine, 
    monaco_first_labels_fine = monaco_first.labels.fine, 
    azimuth_predicted_celltype_l2 = azimuth_predicted.celltype.l2
  ) 

cell_annotation = cell_annotation |> mutate(
  blueprint_first_labels_fine = ifelse(is.na(blueprint_first_labels_fine), "Other", blueprint_first_labels_fine),
  monaco_first_labels_fine = ifelse(is.na(monaco_first_labels_fine), "Other", monaco_first_labels_fine),
  azimuth_predicted_celltype_l2=ifelse(is.na(azimuth_predicted_celltype_l2), "Other", azimuth_predicted_celltype_l2))

# # cell_annotation |> arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/annotation_tbl_light.parquet",
# #                                         compression = "zstd")
# cell_annotation |> arrow::write_parquet("~/scratch/cache_temp/annotation_tbl_light.parquet",
#                                         compression = "zstd")

empty_droplet = 
  tar_read(empty_tbl, store = "/vast/scratch/users/shen.m/cellNexus_target_store_2025-11-08_rest_three_dataset") |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

alive_cells = 
  tar_read(alive_tbl, store = "/vast/scratch/users/shen.m/cellNexus_target_store_2025-11-08_rest_three_dataset") |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell) |>
  select(-any_of("sample_id"))

doublet_cells =
  tar_read(doublet_tbl, store = "/vast/scratch/users/shen.m/cellNexus_target_store_2025-11-08_rest_three_dataset") |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell) |>
  select(-any_of("sample_id"))

# metacell = 
#   tar_read(metacell_tbl, store = glue::glue("/vast/scratch/users/shen.m/lighten_annotation_tbl_target_{version}_rest_three_dataset")) |> 
#   bind_rows() |> 
#   dplyr::rename(cell_ = cell) |> 
#   dplyr::rename_with(
#     ~ stringr::str_replace(.x, "^gamma", "metacell_"),
#     starts_with("gamma")
#   )

# Save cell type concensus tbl from HPCell output to disk
cell_type_concensus_tbl = tar_read(cell_type_concensus_tbl, store = "/vast/scratch/users/shen.m/cellNexus_target_store_2025-11-08_rest_three_dataset") |>  
  bind_rows() |> 
  dplyr::rename(cell_ = .cell)|>
  select(-any_of("sample_id"))

cell_type_concensus_tbl = cell_type_concensus_tbl |> mutate(cell_type_unified_ensemble = 
                                                              ifelse(is.na(cell_type_unified_ensemble),
                                                                     "Unknown",
                                                                     cell_type_unified_ensemble)) 

# cell_type_concensus_tbl |> arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_type_concensus_tbl_from_hpcell.parquet",
#                                                 compression = "zstd")


# cell_type_concensus_tbl |> arrow::write_parquet("~/scratch/cache_temp/cell_type_concensus_tbl_from_hpcell.parquet",
#                                                 compression = "zstd")
# This command needs a big memory machine
cell_metadata_joined = cell_metadata |> 
  left_join(empty_droplet, copy=TRUE) |>  
  left_join(cell_type_concensus_tbl, copy=TRUE) |>
  #left_join(cell_annotation, copy=TRUE) |>  
  left_join(alive_cells, copy=TRUE) |> 
  left_join(doublet_cells, copy=TRUE)
# |>
#   left_join(metacell, copy=TRUE)

cell_metadata_joined |> filter(is.na(blueprint_first_labels_fine))

cell_metadata_joined2 = cell_metadata_joined |> as_tibble() |> 
  # Match to how pseudobulk annotations get parsed in HPCell/R/functions preprocessing_output()
  mutate(cell_type_unified_ensemble = ifelse(cell_type_unified_ensemble |> is.na(), "Unknown", cell_type_unified_ensemble)) |>
  mutate(data_driven_ensemble = ifelse(data_driven_ensemble |> is.na(), "Unknown", data_driven_ensemble))   |>
  mutate(blueprint_first_labels_fine = ifelse(blueprint_first_labels_fine |> is.na(), "Other", blueprint_first_labels_fine)) |> 
  mutate(monaco_first_labels_fine = ifelse(monaco_first_labels_fine |> is.na(), "Other", monaco_first_labels_fine)) |> 
  mutate(azimuth_predicted_celltype_l2 = ifelse(azimuth_predicted_celltype_l2 |> is.na(), "Other", azimuth_predicted_celltype_l2)) |> 
  mutate(azimuth = ifelse(azimuth |> is.na(), "Other", azimuth)) |> 
  mutate(blueprint = ifelse(blueprint |> is.na(), "Other", blueprint)) |> 
  mutate(monaco = ifelse(monaco |> is.na(), "Other", monaco))

cell_metadata_joined2 |>
  arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_annotation_rest_three_dataset.parquet",
                       compression = "zstd")

# Downstream
library(arrow)
library(dplyr)
library(duckdb)

job::job({
  
  get_file_ids = function(cell_annotation ){
    sample_chunk_df = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue::glue("SELECT * FROM read_parquet('{cell_annotation}')"))
      ) |> 
      # Define chunks
      dplyr::count(dataset_id, sample_id, name = "cell_count") |>  # Ensure unique dataset_id and sample_id combinations
      distinct(dataset_id, sample_id, cell_count) |>  # Ensure unique dataset_id and sample_id combinations
      group_by(dataset_id) |> 
      dbplyr::window_order(dataset_id, cell_count, sample_id) |>  # Ensure order. Note: order cell_count only is not enough because it needs a secondary tie-breaker
      mutate(sample_index = row_number()) |>  # Create sequential index within each dataset
      mutate(sample_chunk = (sample_index - 1) %/% 1000 + 1) |>  # Assign chunks (up to 1000 samples per chunk)
      mutate(sample_pseudobulk_chunk = (sample_index - 1) %/% 250 + 1) |> # Max combination of dataset_id, sample_pseudobulk_chunk and file_id_pseudobulk up to 10000
      mutate(cell_chunk = cumsum(cell_count) %/% 100000 + 1) |> # max 20K cells per sample
      ungroup() 
    
    # Test whether cell_chunk and sample_chunk are unique for this sample
    run_chunk_once <- function(column_name, id) {
      sample_chunk_df |> filter(sample_id == id) |> pull(!!column_name)
    }
    
    sample_chunk_results <- replicate(20, run_chunk_once("sample_chunk", "d6e942a09a140ee8bb6f0c3da8defea4___exp7-human-150well."), simplify = FALSE)
    sample_chunk_identical <- all(sapply(sample_chunk_results[-1], function(x) identical(x, sample_chunk_results[[1]])))
    if (!sample_chunk_identical) {
      stop("Inconsistent sample chunk value was generated in multiple runs, this will lead to file id changes")
    }
    
    cell_chunk_results <- replicate(20, run_chunk_once("cell_chunk",  "d6e942a09a140ee8bb6f0c3da8defea4___exp7-human-150well."), simplify = FALSE)
    cell_chunk_identical <- all(sapply(cell_chunk_results[-1], function(x) identical(x, cell_chunk_results[[1]])))
    if (!cell_chunk_identical) {
      stop("Inconsistent cell chunk value was generated in multiple runs, this will lead to file id changes")
    }
    
    
    tbl(
      dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
      sql(glue::glue("SELECT * FROM read_parquet('{cell_annotation}')"))
    ) |> 
      # Cells in cell_annotation could be more than cells in cell_consensus. In order to avoid NA happens in cell_consensus cell_type column
      mutate(cell_type_unified_ensemble = ifelse(cell_type_unified_ensemble |> is.na(),
                                                 "Unknown",
                                                 cell_type_unified_ensemble)) |>
      
      left_join(sample_chunk_df |> select(dataset_id, sample_chunk, sample_pseudobulk_chunk, cell_chunk, sample_id), copy=TRUE) |> 
      
      # Define chunks
      group_by(dataset_id, sample_chunk, cell_chunk, sample_pseudobulk_chunk, cell_type, sample_id) |>
      summarise(cell_count = n(), .groups = "drop") |>
      group_by(dataset_id, sample_chunk, cell_chunk, cell_type) |>
      dbplyr::window_order(desc(cell_count), sample_id) |> # Important!
      mutate(chunk = cumsum(cell_count) %/% 20000 + 1) |> # max 20K cells per sample
      ungroup() |> 
      as_tibble() |> 
      
      # Single cell file ID
      mutate(file_id_cellNexus_single_cell = 
               glue::glue("{dataset_id}___{sample_chunk}___{cell_chunk}___{cell_type}") |> 
               sapply(digest::digest) |> 
               paste0("___", chunk, ".h5ad") 
      ) |> 
      
      # Pseudobulk file id
      mutate(file_id_cellNexus_pseudobulk = 
               glue::glue("{dataset_id}___{sample_pseudobulk_chunk}") |> 
               sapply(digest::digest) |>
               paste0("___", chunk, ".h5ad"))
    
  }
  
  
  # FOR MENGYUAN CELL_METADATA COULD BE BIGGER THAN CELL_ANNOTATION
  
  get_file_ids(
    "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_annotation_rest_three_dataset.parquet"
  )  |> 
    write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/file_id_cellNexus_single_cell_rest_three_dataset.parquet")
  
  gc()
  
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Create a view for cell_annotation in DuckDB
  dbExecute(con, "
  CREATE VIEW cell_metadata AS
  SELECT 
    CONCAT(cell_, '___', dataset_id) AS cell_,
    * EXCLUDE (cell_, dataset_id_1, cell_type)  -- drop original cell_ and dataset_id_1  
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet')
  WHERE dataset_id IN ('21722308-5091-4b63-9c07-4f116ab9a7b3', '145dcf6a-2461-4fa3-a0af-7fc56db0bd33', '7f7faf6b-f11d-4f07-bc1c-188a4472748d')
")
  
  dbExecute(con, "
  CREATE VIEW hpcell_output_metadata AS
   SELECT * EXCLUDE (
    observation_joinid,
    cell_type_ontology_term_id,
    assay,
    donor_id,
    is_primary_data,
    self_reported_ethnicity,
    tissue,
    azimuth,
    blueprint,
    monaco,
    subsets_Mito_sum,
    subsets_Mito_detected,
    ensemble_joinid,
    cell_type_unified,
    data_driven_ensemble,
    observation_originalid
)

  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_annotation_rest_three_dataset.parquet')
")

dbExecute(con, "
  CREATE VIEW file_id_cellNexus_single_cell AS
  SELECT 
    dataset_id,
    sample_chunk,
    cell_chunk,
    sample_pseudobulk_chunk,
    cell_type,
    sample_id,
    file_id_cellNexus_single_cell,
    file_id_cellNexus_pseudobulk
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/file_id_cellNexus_single_cell_rest_three_dataset.parquet')
")

# MODIFY HERE: transformation data frame
dbExecute(con, "
  CREATE VIEW sample_distribution_method_tbl AS
  SELECT 
    sample_id,
    count_upper_bound,
    feature_thresh AS nfeature_expressed_thresh,
    method_to_apply AS inverse_transform
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/updated_transform_sample_tbl_2025_Nov_rest_three_dataset.parquet')
")

# Perform the left join and save to Parquet
copy_query <- "
  COPY (
     SELECT 
        cell_metadata.cell_ AS cell_id, -- Rename cell_ to cell_id
        COALESCE(hpcell_output_metadata.alive, FALSE) AS alive, -- Set alive column NULL to FALSE
        cell_metadata.* EXCLUDE (cell_),          -- drop cell_ since it's already aliased as cell_id
        hpcell_output_metadata.* EXCLUDE (cell_, dataset_id, sample_id, alive), -- Deduplicate join keys, and aliased column
        file_id_cellNexus_single_cell.* EXCLUDE (sample_id, dataset_id, cell_type), -- Deduplicate join keys
        sample_distribution_method_tbl.* EXCLUDE (sample_id) -- Deduplicate join keys
      FROM cell_metadata

      LEFT JOIN hpcell_output_metadata
        ON hpcell_output_metadata.cell_ = cell_metadata.cell_
        AND hpcell_output_metadata.dataset_id = cell_metadata.dataset_id
    
      LEFT JOIN file_id_cellNexus_single_cell
        ON file_id_cellNexus_single_cell.sample_id = hpcell_output_metadata.sample_id
        AND file_id_cellNexus_single_cell.dataset_id = hpcell_output_metadata.dataset_id
        AND file_id_cellNexus_single_cell.cell_type = hpcell_output_metadata.cell_type
      
      LEFT JOIN sample_distribution_method_tbl
        ON sample_distribution_method_tbl.sample_id = cell_metadata.sample_id
        
      -- (THESE DATASETS DOESNT contain meaningful data - no observation_joinid etc), thus was excluded in the final metadata.
      WHERE cell_metadata.dataset_id NOT IN ('99950e99-2758-41d2-b2c9-643edcdf6d82', '9fcb0b73-c734-40a5-be9c-ace7eea401c9', '60a29d0b-1a37-4447-ac32-00d701580b47', '09b518f9-da64-44cc-aec8-70a89d55611f', 'cb252df6-6e49-4553-abd1-495a00006fb1') 
         
  ) TO  '/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata_cell_type_consensus_v1_0_0_rest_three_dataset_mengyuan.parquet'
  (FORMAT PARQUET, COMPRESSION 'gzip');
"

# Execute the final query to write the result to a Parquet file
dbExecute(con, copy_query)

# Disconnect from the database
dbDisconnect(con, shutdown = TRUE)

print("Done.")
})

# We decided to make cell_id lighter without re-run everything in HPCell pipeline. Here to swap cell_id in the metadata
# cell_map is processed in a separate target script step6_supp_dataset_cell_map.R
job::job({
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Create a view for cell_annotation in DuckDB
  # MODIFY HERE: v1_2_2 merged metadata parquet path inside the SQL string below (should match the COPY TO output above)
  dbExecute(con, "
  CREATE VIEW cell_metadata AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata_cell_type_consensus_v1_0_0_rest_three_dataset_mengyuan.parquet')
")
  
  # MODIFY HERE: cell_id dictionary parquet path inside the SQL string below
  dbExecute(con, "
  CREATE VIEW cell_map AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/file_id_cell_id_dict_v1_0_0.parquet')
")
  
  # Perform the left join and save to Parquet
  copy_query <- "
  COPY (
     SELECT cell_metadata.*,
            cell_map.* EXCLUDE (cell_id, file_id_cellNexus_single_cell)
      FROM cell_metadata
    
      LEFT JOIN cell_map
        ON cell_metadata.cell_id = cell_map.cell_id
        AND cell_metadata.file_id_cellNexus_single_cell = cell_map.file_id_cellNexus_single_cell

  ) TO  '/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata_cell_type_consensus_v1_0_1_mengyuan.parquet' -- MODIFY HERE: output final metadata parquet with new cell IDs (v1_3_2)
  (FORMAT PARQUET, COMPRESSION 'gzip');
"
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
  print("Done.")
  
  
})





