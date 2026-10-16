.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")

source(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/36_TITAN2致病菌ASV水平_水位/分析代码/运行代码/TITAN2_WL_top20_不同种代表稳健ASV山脊图_最新注释.R"),
  encoding = "UTF-8")
