.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")

source(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/34_蓝藻门结构方程模型/分析代码/03_蓝藻门月均相对丰度纵向曲线.R"),
  encoding = "UTF-8")
