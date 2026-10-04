library(cellNexus)
library(tidyseurat)
library(arrow)
library(purrr)
library(SingleCellExperiment)
library(glue)
library(stringr)
library(targets)
library(dplyr)
prostate_se = readRDS("/vast/scratch/users/hutchison.w/prostate_atlas_transfer/prostate_atlas_sc_v3.3.rds")
prostate_se[["RNA_decontx"]] <- NULL

metadata_df <- read.csv("/vast/scratch/users/hutchison.w/prostate_atlas_transfer/prostate_atlas_metadata_v3.3.csv",
                        header = TRUE, na.strings = c("NA",""))

# Retrieve metadata
prostate_se_meta <- prostate_se[[]] |> tibble::rownames_to_column("cell_id") |> as_tibble()

harmonised_meta <- prostate_se_meta |>
  mutate(
    annotation_manual_collapsed = case_when(
      str_starts(annotation_manual_fine, "unidentified") ~ "unidentified",
      str_starts(annotation_manual_fine, "macrophage") ~ "macrophage",
      .default = annotation_manual_fine
    )
  ) |>
  mutate(cell_type = annotation_manual_collapsed,
         # file_id_cellNexus_single_cell = paste(study_id, cell_type, sep = "_") |> 
         #   sapply(digest::digest) |> 
         #       paste0(".h5ad") ,
         age_days = sample_age_years * 365,
         sex = "male",
         # organism = "Homo sapiens",
         tissue = if_else(is.na(sample_location), "prostate", sample_location),
         collection_id = "hutchison_prostate_2026"
         # is_primary_data_x = ifelse(dataset_id == "human_cell_atlas_prostate", "FALSE", "TRUE")
  ) |> 
  dplyr::rename(self_reported_ethnicity = sample_ethnicity,
                assay = sample_technology_sequencer,
                disease = sample_disease,
                cell_annotation_blueprint_singler = annotation_blueprint_fine
  )

first_cols <- c(
  "cell_id",
  "cell_type",
  "sample_id",
  "donor_id",
  "disease",
  "tissue",
  "collection_id",
  "age_days",
  "assay",
  "sex",
  "cell_annotation_blueprint_singler"
)

harmonised_meta = harmonised_meta|>
  select(
    all_of(first_cols),
    all_of(sort(setdiff(names(harmonised_meta), first_cols)))
  )


# Add harmonised prostate_se_meta back to seurat object
harmonised_meta2 <- harmonised_meta |>
  tibble::column_to_rownames("cell_id")

stopifnot(setequal(rownames(harmonised_meta2), colnames(prostate_se)))
harmonised_meta2 <- harmonised_meta2[colnames(prostate_se), , drop = FALSE]
prostate_se@meta.data <- harmonised_meta2

# split by sample for now. Split by study_id and cell type is more ideal in the future.
seurat_sample_split_list <-  Seurat::SplitObject(prostate_se, split.by = "sample_id")

saved <- seurat_sample_split_list |> saveRDS("/vast/projects/cellxgene_curated/prostate_atlas/prostate_split_by_sample_id.rds")

# Processing 
target_store = "/vast/projects/cellxgene_curated/prostate_atlas/"
tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  library(Seurat)
  library(tidyseurat)
  library(SingleCellExperiment)
  library(glue)
  library(stringr)
  library(tidySingleCellExperiment)
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

elastic_80  <- new_elastic("elastic_80",   80,  60 * 4,  workers = 35, crashes_max = 1, cpus_per_task = 1)
elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 70, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 140, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 290, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 440, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)

# Group for targets (small → large)
controllers <- crew_controller_group(
  elastic_10, elastic_20, elastic_40, elastic_80, elastic_5_minimal
)
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
  resources = tar_resources(
    crew = tar_resources_crew(controller = "elastic_5_minimal")
  ) 
)

list(
  tar_target(
    counts_path,
    "/vast/projects/cellxgene_curated/prostate_atlas/counts/",
    deployment = "main"
  ),
  
  tar_target(
    cpm_path,
    "/vast/projects/cellxgene_curated/prostate_atlas/cpm/",
    deployment = "main"
  ),
  
  tar_target(
    seurat_split_list,
    readRDS("/vast/projects/cellxgene_curated/prostate_atlas/prostate_split_by_sample_id.rds")
  ),
  
  tar_target(
    sample_id,
    names(seurat_split_list),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_10")
    )
  ),
  
  # Process each file_id independently
  tar_target(
    processed_single_cell_sce,
    {
      name <- sample_id
      
      seu_obj <- seurat_split_list[[name]]
      
      mat <- GetAssayData(
        seu_obj,
        layer = "counts",
        assay = "RNA"
      )
      
      sce <- SingleCellExperiment(
        assays = list(counts = mat),
        colData = DataFrame(
          seu_obj[[]][, "sample_id", drop = FALSE]
        )
      )
      sce
      },
    pattern = map(sample_id),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_20")
    )
    ),
    
    tar_target(
      save_sce,
      {
        
        .name <- processed_single_cell_sce |> pull(sample_id) |> unique()
        
        counts_dir <- file.path(counts_path, .name)

        
        HPCell::save_experiment_data( processed_single_cell_sce,  dir = counts_dir)
        
        },
      pattern = map(processed_single_cell_sce)
      ),
      
      tar_target(
        save_cpm,
        {
          .name = processed_single_cell_sce |> pull(sample_id) |> unique() |> paste0(".h5ad")
          cpm_file <- file.path(cpm_path, .name) 
          cellNexus::get_counts_per_million(processed_single_cell_sce,cpm_file)
        },
        pattern = map(processed_single_cell_sce)
      ),
  
   tar_target(
     cell_type_annotation_df,
     HPCell::annotation_label_transfer(processed_single_cell_sce, feature_nomenclature = "symbol", reference_azimuth = "pbmcref"),
     pattern = map(processed_single_cell_sce),
     packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "HPCell", "digest", "scater", "arrow", "dplyr"),
     resources = tar_resources(
       crew = tar_resources_crew(controller = "elastic_20")
     )
   )
)}, script = paste0(target_store, "_target_script.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = paste0(target_store, "_target_script.R"), 
    store = target_store, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

cell_type_annotation_df = tar_read(cell_type_annotation_df,store=target_store) |> bind_rows()

cell_type_annotation_df = cell_type_annotation_df |>
  dplyr::rename(
    blueprint_first_labels_fine = blueprint_first.labels.fine, 
    blueprint_first_labels_coarse = blueprint_first.labels.coarse,
    monaco_first_labels_fine = monaco_first.labels.fine,
    monaco_first_labels_coarse = monaco_first.labels.coarse,
    azimuth_predicted_celltype_l2 = azimuth_predicted.celltype.l2
  ) |>

cell_type_annotation_df = cell_type_annotation_df |>
  left_join(celltype_unification_maps$azimuth, copy = TRUE) |>
  left_join(celltype_unification_maps$blueprint, copy = TRUE) |>
  left_join(celltype_unification_maps$monaco, copy = TRUE) |>
  left_join(prostate_se[[]] |> tibble::rownames_to_column("cell_id") |>
              mutate(
                annotation_manual_collapsed = case_when(
                  str_starts(cell_type, "unidentified") ~ "unidentified",
                  str_starts(cell_type, "macrophage") ~ "macrophage",
                  .default = cell_type
                )
              ) |>
              select(cell_id, sample_id,
                     cell_type_unified = annotation_manual_collapsed), by = c(".cell"="cell_id"), copy = TRUE) |>
  mutate(ensemble_joinid = paste(azimuth, blueprint, monaco, cell_type_unified, sep = "_"))

df_map = cell_type_annotation_df |>
  dplyr::count(ensemble_joinid, azimuth, blueprint, monaco, cell_type_unified, name = "NCells") |>
  as_tibble() |>
  mutate(
    cellxgene = if_else(cell_type_unified %in% nonimmune_cellxgene, "non immune", cell_type_unified),
    data_driven_ensemble = ensemble_annotation(cbind(azimuth, blueprint, monaco), override_celltype = c("non immune", "nkt", "mast")),
    cell_type_unified_ensemble = ensemble_annotation(cbind(azimuth, blueprint, monaco, cellxgene), method_weights = c(1, 1, 1, 2), override_celltype = c("non immune", "nkt", "mast")),
    cell_type_unified_ensemble = case_when(
      cell_type_unified_ensemble == "non immune" & cellxgene == "non immune" ~ cell_type_unified,
      cell_type_unified_ensemble == "non immune" & cellxgene != "non immune" ~ "other",
      .default = cell_type_unified_ensemble
    ),
    is_immune = !cell_type_unified_ensemble %in% nonimmune_cellxgene
  ) |>
  select(
    ensemble_joinid,
    data_driven_ensemble,
    cell_type_unified_ensemble,
    is_immune
  )

cell_type_annotation_df2 = cell_type_annotation_df |>
  left_join(df_map, by = join_by(ensemble_joinid), copy = TRUE) |> 
  mutate(cell_type_unified_ensemble = ifelse(cell_type_unified_ensemble |> is.na(), "Unknown", cell_type_unified_ensemble))

cell_type_annotation_df2 = cell_type_annotation_df2 |> select(.cell, monaco_first_labels_fine, 
                                                              azimuth_predicted_celltype_l2, sample_id, 
                                                              cell_type_unified_ensemble, is_immune)

final_metadata = harmonised_meta |> left_join(cell_type_annotation_df2, by = c("cell_id" = ".cell", "sample_id"), copy=TRUE) |>
  mutate(file_id_cellNexus_single_cell = paste0(sample_id, ".h5ad"),
         # v0.3.3 to match rds version
         atlas_id = "prostate_2026/0.3.3")

final_metadata |>
  arrow::write_parquet("/vast/projects/cellxgene_curated/prostate_atlas/prostate_metadata.v0.3.3.parquet")

x = get_metadata(local_metadata = "/vast/projects/cellxgene_curated/prostate_atlas/prostate_metadata.v0.3.3.parquet",
             cloud_metadata = NULL) |>
  get_single_cell_experiment(cache_directory = "/home/users/allstaff/shen.m/cellxgene_curated/prostate_atlas",
                             assays = c("counts","cpm"))






