# Load shared paths; run analysis scripts from the project root.
# Input: config/config.yaml. Output: path objects and output directories.
project_config <- yaml::read_yaml("config/config.yaml")
project_paths <- project_config$paths
cohort_files <- unlist(project_config$cohort_results, use.names = TRUE)
for (directory in c(
  project_paths$ukb_processed, project_paths$ukb_results,
  project_paths$cross_cohort, project_paths$plots,
  file.path(project_paths$ukb_processed, c("bootstrap_rec", "bootstrap_2i"))
)) {
  if (!dir.exists(directory) &&
      !dir.create(directory, recursive = TRUE)) {
    stop("Cannot create output directory: ", directory)
  }
}
