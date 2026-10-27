# Shared style for the main-figure scripts 43-47: one theme, display labels
# for the traits and a function that saves a figure at print size.
library(ggplot2)

paper_theme <- theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom", panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    plot.tag = element_text(face = "bold", size = 11),
    legend.key.width = grid::unit(1.2, "lines")
  )

trait_labels <- c(
  ALT = "ALT", AST = "AST", GGT = "GGT", HDL = "HDL-C", LDL = "LDL-C",
  TG = "Triglycerides", SBP = "SBP", DBP = "DBP", GLU = "Glucose",
  CRP = "CRP", SCR = "Creatinine", BF_perc = "Body fat %",
  AndF_perc = "Android fat %", GynF_perc = "Gynoid fat %",
  ArmsF_perc = "Arms fat %", LegsF_perc = "Legs fat %",
  TrunkF_perc = "Trunk fat %", logLiverF_perc = "Liver fat (log)"
)

# Figures are drawn at 7 inches wide (print width) and 300 dpi.
save_figure <- function(figure, name, width, height) {
  ggsave(
    file.path("plots", name), figure, width = width, height = height,
    units = "in", dpi = 300, bg = "white"
  )
  cat("Saved plots/", name, "\n", sep = "")
}
