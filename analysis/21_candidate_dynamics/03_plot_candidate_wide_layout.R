Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")

.libPaths(c(file.path(getwd(), "Rlib"), .libPaths()))

library(ggplot2)

library(dplyr)

library(patchwork)

input_dir <- "work/pathogen_dual_axis_trend/input"

cached <- "work/pathogen_dual_axis_trend/output"

read_cached <- function(name) {
  d <- read.csv(file.path(cached, name))
  d$Date <- as.Date(d$Date)
  d
}

sample_metrics <- read_cached("Pathogen_richness_relative_abundance_sample_data.csv")

richness_predictions <- read_cached("Pathogen_richness_GAMM_predictions.csv")

relative_predictions <- read_cached("Pathogen_relative_abundance_GAMM_predictions.csv")

axis_multiplier <- relative_predictions$scale_for_plot[1]

richness_axis_max <- 250

format_p <- function(p) ifelse(p < 0.001, "<0.001", paste0("=", formatC(p, digits = 3, format = "f")))

richness_fit <- relative_fit <- list(statistics = list(p_smooth = 0.001))

src <- readLines(file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/21_candidate_dynamics/02_plot_candidate_dual_axis.R"), encoding = "UTF-8")

start <- grep("^water_observed <-", src)[1]

end <- grep("^combined_plot <-", src)[1] - 1

eval(parse(text = src[start:end]))

p_main <- p_main + labs(title = NULL, subtitle = NULL, x = NULL) + scale_x_date(limits = date_limits_plot,
  breaks = as.Date(c("2024-12-01", "2025-02-01", "2025-04-01", "2025-06-01", "2025-08-01", "2025-10-01",
    "2026-01-01", "2026-04-01")), date_labels = "%Y-%m", expand = expansion(mult = c(0.025, 0.025))) +
  theme(legend.position = "none", axis.text = element_text(size = 11), axis.text.x = element_text(size = 13,
    angle = 45, hjust = 1), axis.title = element_text(size = 12, face = "bold"), plot.margin = margin(5,
    8, 5, 8))

p_water_inset <- p_water_inset + theme(axis.text = element_text(size = 7.5), axis.text.x = element_text(size = 7,
  angle = 45, hjust = 1), axis.title = element_text(size = 7, face = "bold"), plot.title = element_text(size = 8,
  face = "bold"), plot.subtitle = element_text(size = 6.5))

combined_plot <- p_main + inset_element(p_water_inset, left = 0.075, bottom = 0.53, right = 0.47, top = 0.965,
  align_to = "panel", on_top = TRUE, clip = TRUE)

out <- "Pathogen_trend_wider_20260921"

dir.create(out, showWarnings = FALSE)

prefix <- file.path(out, "Pathogen_temporal_trend_wider")

ggsave(paste0(prefix, ".png"), combined_plot, width = 10.5, height = 5.5, dpi = 600, device = ragg::agg_png,
  bg = "white")

ggsave(paste0(prefix, ".pdf"), combined_plot, width = 10.5, height = 5.5, device = cairo_pdf, bg = "white")

ggsave(paste0(prefix, ".svg"), combined_plot, width = 10.5, height = 5.5, device = svglite::svglite,
  bg = "white")

ggsave(paste0(prefix, "_preview.png"), combined_plot, width = 10.5, height = 5.5, dpi = 120, device = ragg::agg_png,
  bg = "white")
