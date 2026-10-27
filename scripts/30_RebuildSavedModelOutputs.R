# Reuse full fitted models; run their existing reporting code and reduced fits.
# This avoids overwriting expensive full models solely to regenerate outputs.
project_dir <- getwd()
output_scripts <- c("05_GAMRec.R", "12_GAMImg2i.R", "13_GAMlogLF.R")
for (script in output_scripts) {
  environment <- new.env(parent = globalenv())
  expressions <- parse(file.path(project_dir, "scripts", script))
  for (expression in expressions) {
    code <- paste(deparse(expression), collapse = " ")
    assignment <- is.call(expression) &&
      identical(expression[[1]], as.name("<-"))
    target <- if (assignment) as.character(expression[[2]]) else ""
    if (assignment && target %in% c(
      "slopebmregdf", "slopebmregdf_2i", "slopebmregdf_2ilf"
    ) &&
      grepl("mod1 = map", code, fixed = TRUE)) {
      filename <- paste0(target, ".rds")
      assign(target, readRDS(file.path(
        project_dir, "data", "UKB", "processed_data", filename
      )), envir = environment)
      next
    }
    if (assignment && target == "slopebmregdf_2ilf" &&
          grepl("bind_rows", code, fixed = TRUE)) next
    eval(expression, envir = environment)
  }
  # Save inputs to each composite, not only the assembled image.
  plot_names <- c("bivardistrib_plot", "varex_plot", "varexp_plot",
                  "fai_plot", "relfai_plot", "pdiff_plot")
  plot_names <- plot_names[vapply(
    plot_names, exists, logical(1), envir = environment, inherits = FALSE
  )]
  plot_inputs <- mget(plot_names, envir = environment)
  saveRDS(plot_inputs, file.path(
    project_dir, "data", "UKB", "processed_data",
    paste0("figure_panels_", sub(".R$", "", script), ".rds")
  ))
}
