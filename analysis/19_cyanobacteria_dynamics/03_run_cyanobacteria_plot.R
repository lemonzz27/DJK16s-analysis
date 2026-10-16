.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")

source(file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/19_cyanobacteria_dynamics/02_plot_monthly_cyanobacteria.R"),
  encoding = "UTF-8")
