library(targets)

target_store <- "/vast/projects/cellxgene_curated/prostate_atlas/"

tar_script({
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(glue)
  library(arrow)
  library(Seurat)
  library(tidyseurat)
  library(SingleCellExperiment)
  library(tidySingleCellExperiment)
  library(cellNexus)
  library(HPCell)
  
  # ---- Config ------------------------------------------------------------------
  atlas_dir       <- "/vast/projects/cellxgene_curated/prostate_atlas/"
  atlas_id        <- "prostate_2026/0.3.3" # v0.3.3 to match rds version
  collection_id   <- "hutchison_prostate_2026"
  seurat_rds_path <- "/vast/scratch/users/hutchison.w/prostate_atlas_transfer/prostate_atlas_sc_v3.3.rds"
  counts_dir      <- file.path(atlas_dir, atlas_id,  "counts")
  cpm_dir         <- file.path(atlas_dir, atlas_id, "cpm")
  metadata_path   <- file.path(atlas_dir, "prostate_metadata.v0.3.3.parquet")
  smoke_test_cache <- "/home/users/allstaff/shen.m/cellxgene_curated/prostate_atlas"
  
  # Columns placed first in the final metadata; the remaining columns are sorted.
  first_cols <- c(
    "cell_id", "cell_type", "sample_id", "donor_id", "dataset_id", "disease",
    "tissue", "collection_id", "age_days", "assay", "sex",
    "self_reported_ethnicity", "feature_count", "is_primary_data",
    "cell_annotation_blueprint_singler", "cell_annotation_monaco_singler",
    "cell_annotation_azimuth_l2", "cell_type_unified_ensemble", "is_immune",
    "file_id_cellNexus_single_cell", "atlas_id"
  )
  
  # ---- Slurm controllers -------------------------------------------------------
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
  
  # Each controller falls back to the next larger one when a worker crashes.
  elastic_80         <- new_elastic("elastic_80",         80, 60 * 4,  35, crashes_max = 1, cpus_per_task = 1)
  elastic_40         <- new_elastic("elastic_40",         40, 60 * 4,  70, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
  elastic_20         <- new_elastic("elastic_20",         20, 60 * 4, 140, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
  elastic_10         <- new_elastic("elastic_10",         10, 60 * 4, 290, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
  elastic_5_minimal  <- new_elastic("elastic_5_minimal",   5, 60 * 4, 440, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)
  
  # Group for targets (small → large)
  controllers <- crew_controller_group(
    elastic_10, elastic_20, elastic_40, elastic_80, elastic_5_minimal
  )
  
  res_10 <- tar_resources(crew = tar_resources_crew(controller = "elastic_10"))
  res_20 <- tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
  
  tar_option_set(
    memory = "transient",
    garbage_collection = 100,
    storage = "worker",
    retrieval = "worker",
    error = "continue",
    cue = tar_cue(mode = "never"),
    format = "qs",
    workspace_on_error = TRUE,
    controller = controllers,
    trust_object_timestamps = TRUE,
    resources = tar_resources(crew = tar_resources_crew(controller = "elastic_5_minimal"))
  )
  
  list(
    
    # ---- 1. Load atlas and harmonise metadata (main process: large objects) ------
    
    tar_target(
      seurat_atlas,
      {
        se <- readRDS(seurat_rds_path)
        se[["RNA_decontx"]] <- NULL
        se
      },
      deployment = "main"
    ),
    
    tar_target(
      gene_counts_by_study,
      seurat_atlas |>
        group_by(study_id) |>
        group_split() |>
        map_dfr(\(x) {
          tibble(study_id = unique(x$study_id),
                 feature_count = length(rownames(x)))
        }, .progress = TRUE),
      deployment = "main"
    ),
    
    # Retrieve metadata
    tar_target(
      harmonised_meta,
      seurat_atlas[[]] |>
        rownames_to_column("cell_id") |>
        as_tibble() |>
        mutate(
          annotation_manual_collapsed = case_when(
            str_starts(annotation_manual_fine, "unidentified") ~ "unidentified",
            str_starts(annotation_manual_fine, "macrophage") ~ "macrophage",
            .default = annotation_manual_fine
          ),
          cell_type = annotation_manual_collapsed,
          age_days = sample_age_years * 365,
          sex = "male",
          tissue = if_else(is.na(sample_location), "prostate", sample_location),
          collection_id = collection_id,
          # The study GSE172301 was accessed using cellNexus and all others were accessed by other methods.
          # is_primary_data FALSE is about whether those same cells are already represented elsewhere
          is_primary_data = ifelse(str_detect(study_data, "cellxgene.cziscience.com"), 
                                   FALSE, TRUE),
          dataset_id = ifelse(
            str_detect(study_data, "cellxgene.cziscience.com"),
            str_extract(study_data, "(?<=/e/)[0-9a-f-]+(?=\\.cxg)"),
            study_id
          )
        ) |>
        dplyr::rename(
          self_reported_ethnicity = sample_ethnicity,
          assay = sample_technology_sequencer,
          disease = sample_disease
        ) |>
        left_join(
          gene_counts_by_study,
          by = "study_id"
        ),
      deployment = "main"
    ),
    
    # Add harmonised metadata back to the Seurat object, then split by sample.
    # split by sample for now. Split by study_id and cell type is more ideal in the future.
    # Done in one target so the full Seurat object is not stored a second time.
    tar_target(
      seurat_sample_split_list,
      {
        se <- seurat_atlas
        meta <- harmonised_meta |> column_to_rownames("cell_id")
        stopifnot(setequal(rownames(meta), colnames(se)))
        se@meta.data <- meta[colnames(se), , drop = FALSE]
        SplitObject(se, split.by = "sample_id")
      },
      deployment = "main"
    ),
    
    # ---- 2. Per-sample processing ------------------------------------------------
    
    tar_target(
      sample_ids,
      unique(harmonised_meta$sample_id),
      deployment = "main"
    ),
    
    # Process each sample independently
    tar_target(
      sce_by_sample,
      {
        seu_obj <- seurat_sample_split_list[[sample_ids]]
        SingleCellExperiment(
          assays = list(counts = GetAssayData(seu_obj, layer = "counts", assay = "RNA")),
          colData = DataFrame(seu_obj[[]][, "sample_id", drop = FALSE])
        )
      },
      pattern = map(sample_ids),
      resources = res_20
    ),
    
    tar_target(
      save_counts,
      {
        dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
        HPCell::save_experiment_data(sce_by_sample, dir = file.path(counts_dir, sample_ids))
      },
      pattern = map(sce_by_sample, sample_ids)
    ),
    
    tar_target(
      save_cpm,
      {
        dir.create(cpm_dir, recursive = TRUE, showWarnings = FALSE)
        cellNexus::get_counts_per_million(
          sce_by_sample,
          file.path(cpm_dir, paste0(sample_ids, ".h5ad")))
      },
      pattern = map(sce_by_sample, sample_ids)
    ),
    
    tar_target(
      cell_type_annotation_raw,
      HPCell::annotation_label_transfer(
        sce_by_sample,
        feature_nomenclature = "symbol",
        reference_azimuth = "pbmcref"
      ),
      pattern = map(sce_by_sample),
      resources = res_20
    ),
    
    # ---- 3. Cell type harmonisation (main process: all samples combined) ---------
    
    # Natural joins onto celltype_unification_maps: if a join key changes
    # upstream, check the "Joining with `by = ...`" message for this target.
    tar_target(
      cell_type_annotation,
      cell_type_annotation_raw |>
        bind_rows() |>
        dplyr::rename(
          blueprint_first_labels_fine = blueprint_first.labels.fine,
          blueprint_first_labels_coarse = blueprint_first.labels.coarse,
          monaco_first_labels_fine = monaco_first.labels.fine,
          monaco_first_labels_coarse = monaco_first.labels.coarse,
          azimuth_predicted_celltype_l2 = azimuth_predicted.celltype.l2
        ) |>
        left_join(celltype_unification_maps$azimuth, copy = TRUE) |>
        left_join(celltype_unification_maps$blueprint, copy = TRUE) |>
        left_join(celltype_unification_maps$monaco, copy = TRUE) |>
        # harmonised_meta$cell_type is already the collapsed manual annotation
        left_join(
          harmonised_meta |> select(cell_id, sample_id, cell_type_unified = cell_type),
          by = c(".cell" = "cell_id"),
          copy = TRUE
        ) |>
        mutate(ensemble_joinid = paste(azimuth, blueprint, monaco, cell_type_unified, sep = "_")),
      deployment = "main"
    ),
    
    # One row per unique label combination, so the ensemble is computed once per
    # combination rather than once per cell.
    tar_target(
      ensemble_map,
      cell_type_annotation |>
        dplyr::count(ensemble_joinid, azimuth, blueprint, monaco, cell_type_unified, name = "NCells") |>
        as_tibble() |>
        mutate(
          cellxgene = if_else(cell_type_unified %in% nonimmune_cellxgene, "non immune", cell_type_unified),
          data_driven_ensemble = ensemble_annotation(
            cbind(azimuth, blueprint, monaco),
            override_celltype = c("non immune", "nkt", "mast")
          ),
          cell_type_unified_ensemble = ensemble_annotation(
            cbind(azimuth, blueprint, monaco, cellxgene),
            method_weights = c(1, 1, 1, 2),
            override_celltype = c("non immune", "nkt", "mast")
          ),
          cell_type_unified_ensemble = case_when(
            cell_type_unified_ensemble == "non immune" & cellxgene == "non immune" ~ cell_type_unified,
            cell_type_unified_ensemble == "non immune" & cellxgene != "non immune" ~ "other",
            .default = cell_type_unified_ensemble
          ),
          is_immune = !cell_type_unified_ensemble %in% nonimmune_cellxgene
        ) |>
        select(ensemble_joinid, data_driven_ensemble, cell_type_unified_ensemble, is_immune),
      deployment = "main"
    ),
    
    tar_target(
      cell_type_annotation_final,
      cell_type_annotation |>
        left_join(ensemble_map, by = join_by(ensemble_joinid), copy = TRUE) |>
        mutate(cell_type_unified_ensemble = ifelse(is.na(cell_type_unified_ensemble), "Unknown", cell_type_unified_ensemble)) |>
        select(
          .cell,
          cell_annotation_blueprint_singler = blueprint_first_labels_fine,
          cell_annotation_monaco_singler = monaco_first_labels_fine,
          cell_annotation_azimuth_l2 = azimuth_predicted_celltype_l2,
          sample_id,
          cell_type_unified_ensemble,
          is_immune
        ),
      deployment = "main"
    ),
    
    # ---- 4. Final metadata -------------------------------------------------------
    
    tar_target(
      final_metadata,
      {
        meta <- harmonised_meta |>
          left_join(cell_type_annotation_final, by = c("cell_id" = ".cell", "sample_id"), copy = TRUE) |>
          mutate(
            file_id_cellNexus_single_cell = paste0(sample_id, ".h5ad"),
            atlas_id = atlas_id
          )
        select(meta, any_of(first_cols), all_of(sort(setdiff(names(meta), first_cols))))
      },
      deployment = "main"
    ),
    
    tar_target(
      metadata_parquet,
      {
        write_parquet(final_metadata, metadata_path)
        metadata_path
      },
      format = "file",
      deployment = "main"
    )
  )
  
}, script = paste0(target_store, "_target_script.R"), ask = FALSE)

# Run the pipeline as a background RStudio job
job::job({
  targets::tar_make(
    script = paste0(target_store, "_target_script.R"),
    store = target_store,
    reporter = "summary"
  )
})

sce = get_metadata(local_metadata = "/vast/projects/cellxgene_curated/prostate_atlas/prostate_metadata.v0.3.3.parquet",
                   cloud_metadata = NULL) |>
  get_single_cell_experiment(cache_directory = "/home/users/allstaff/shen.m/cellxgene_curated/prostate_atlas",
                             assays = c("counts","cpm"))

