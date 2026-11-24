library(dplyr)
library(ggplot2)
library(patchwork)

read_rds_from_box <- function(remote_path) {
  tmp <- tempfile(fileext = ".rds")
  
  system2(
    command = "/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone",
    args = c("copyto", remote_path, tmp)
  )
  
  readRDS(tmp)
}

cell_metadata_1_0_15_sccomp_input_counts <- read_rds_from_box(
  "box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/cell_metadata_1_0_15_sccomp_input_counts.rds"
)

cell_metadata_1_0_14_sccomp_input_counts <- read_rds_from_box(
  "box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/cell_metadata_1_0_14_sccomp_input_counts.rds"
)


# Healthy donors only
healthy_15 <- cell_metadata_1_0_15_sccomp_input_counts |>
  filter(disease_groups == "Normal") |>
  mutate(age_years = age_days / 365) |>
  filter(!is.na(age_years))

healthy_14 <- cell_metadata_1_0_14_sccomp_input_counts |>
  filter(disease_groups == "Normal") |>
  mutate(age_years = age_days / 365) |>
  filter(!is.na(age_years))


calculate_age_range <- function(data, group_col) {
  data |>
    group_by(across(all_of(group_col))) |>
    summarise(
      n = n(),
      min_age = min(age_years),
      median_age = median(age_years),
      max_age = max(age_years),
      .groups = "drop"
    )
}

plot_age_range <- function(data, group_col, title) {
  
  age_summary <- calculate_age_range(data, group_col)
  
  ggplot(
    age_summary,
    aes(
      y = reorder(.data[[group_col]], median_age),
      x = median_age
    )
  ) +
    geom_errorbar(
      aes(
        xmin = min_age,
        xmax = max_age
      ),
      orientation = "y",
      width = 0.3,
      linewidth = 0.7
    ) +
    geom_point(
      size = 2.5,
      colour = "steelblue"
    ) +
    labs(
      x = "Age (years)",
      y = NULL,
      title = title
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank()
    )
}

# # Plot 1: v15 by tissue
# p1 <- plot_age_range(
#   healthy_15,
#   "tissue",
#   "X: Healthy donor age ranges by tissue"
# )
# 
# # Plot 2: v14 by tissue
# p2 <- plot_age_range(
#   healthy_14,
#   "tissue",
#   "Y: Healthy donor age ranges by tissue"
# )

# by tissue_groups
p3 <- plot_age_range(
  healthy_15,
  "tissue_groups",
  "Healthy donor age ranges by tissue group cell_metadata_1_0_15_sccomp_input_counts"
)

# by tissue_groups
p4 <- plot_age_range(
  healthy_14,
  "tissue_groups",
  "Healthy donor age ranges by tissue group cell_metadata_1_0_14_sccomp_input_counts"
)

p3_p4 <- p3/p4

ggsave(
  filename = "/home/users/allstaff/shen.m/git_control/cellNexus_development/sccomp_immune_analysis/healthy_donor_age_ranges_by_tissue_group_v1.0.14.png",
  plot = p4,
  width = 12,
  height = 8,
  dpi = 300
)


# Unique donor-tissue combinations in x
healthy_15_donors <- healthy_15 |>
  distinct(donor_id, tissue)

# Additional donor-tissue combinations in y
extra <- healthy_14 |>
  anti_join(
    healthy_15_donors,
    by = c("donor_id", "tissue")
  ) |>
  distinct(
    donor_id,
    tissue,
    tissue_groups,
    age_years
  )

# Plot 6: Additional donors by tissue_groups
p6 <- plot_age_range(
  extra,
  "tissue_groups",
  "Additional healthy donors in cell metadata 1.0.15: Age ranges by tissue group"
)
