.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")

source(file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/22_indicator_thresholds/03_plot_indicator_ridgelines.R"),
  encoding = "UTF-8")
