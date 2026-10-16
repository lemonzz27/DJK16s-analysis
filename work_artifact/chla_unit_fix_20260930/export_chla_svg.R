try(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"), silent = TRUE)

.libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "Rlib"), .libPaths()))

source("01_分析项目/环境变化曲线/plot_monthly_environmental_curves_CNS4panel.R", encoding = "UTF-8")

stopifnot("Chla" %in% names(plot_list))

ggplot2::ggsave(filename = "work_artifact/chla_unit_fix_20260930/replot/individual/Chla_monthly_curve.svg",
  plot = plot_list[["Chla"]], width = 9.6, height = 4.8, bg = "white", device = svglite::svglite)
