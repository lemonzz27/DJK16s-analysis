ROOT <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

PROJECT <- file.path(ROOT, "01_分析项目", "34_蓝藻门结构方程模型")

SAMPLE_DATA_FILE <- file.path(PROJECT, "原始数据", "蓝藻门SEM_每样本原始数据_120样本.csv")

RESULT_TABLE_DIR <- file.path(PROJECT, "结果表格", "三组多框架SEM")

RESULT_FIGURE_DIR <- file.path(PROJECT, "结果图片", "三组多框架SEM")

ANALYSIS_GROUPS <- c("all", "large", "small")

GROUP_LABELS <- c(all = "All samples (n = 120)", large = "Large fraction (n = 60)", small = "Small fraction (n = 60)")

RUN_LAVAAN <- TRUE

RUN_PLS <- TRUE

RUN_BAYES <- FALSE

LAVAAN_ESTIMATOR <- "MLR"

PLS_BOOTSTRAP_REPS <- 1000L

BAYES_CHAINS <- 2L

BAYES_CORES <- 2L

BAYES_ITER <- 2000L

BAYES_WARMUP <- 1000L

BAYES_SEED <- 20260820L

VARIABLE_SETTINGS <- data.frame(model_name = c("WL", "Q", "WT", "pH", "DO", "CODmn", "BOD5", "NH3N",
  "TN", "TP", "TNTP", "Cu", "Zn", "F", "Se", "As", "Pb", "Chla", "Cyano"), source_column = c("WL",
  "Q", "WT", "pH", "DO", "CODmn", "BOD5", "NH3N", "TN", "TP", "TN_TP_mass_ratio", "Cu", "Zn", "F",
  "Se", "As", "Pb", "Chla", "cyanobacteria_relative"), transform = c("log10", "log10_half_min", "log10",
  "none", "log10", "log10", "log10_half_min", "log10_half_min", "log10", "log10", "log10", "log10_half_min",
  "log10_half_min", "log10_half_min", "log10_half_min", "log10_half_min", "log10_half_min", "log10",
  "log10_cyano"), standardize = TRUE, stringsAsFactors = FALSE)

CYANO_PSEUDOCOUNT <- 0.5/21000

MODEL_PRESET <- "latent_nutrients"

if (MODEL_PRESET == "observed_path") {
  MEASUREMENT_BLOCKS <- list(WL = "WL", WT = "WT", TN = "TN", TP = "TP", Chla = "Chla", Cyano = "Cyano")
  STRUCTURAL_PATHS <- list(TN = c("WL", "WT"), TP = c("WL", "WT"), Chla = c("WT", "TN", "TP"), Cyano = c("WL",
    "WT", "TN", "TP", "Chla"))
  CORRELATED_RESIDUALS <- list(c("TN", "TP"))
}

if (MODEL_PRESET == "latent_nutrients") {
  MEASUREMENT_BLOCKS <- list(WL = "WL", WT = "WT", Nutrients = c("TN", "TP", "NH3N"), Chla = "Chla",
    Cyano = "Cyano")
  STRUCTURAL_PATHS <- list(Nutrients = c("WL", "WT"), Chla = c("WT", "Nutrients"), Cyano = c("WL",
    "WT", "Nutrients", "Chla"))
  CORRELATED_RESIDUALS <- list()
}

if (!MODEL_PRESET %in% c("observed_path", "latent_nutrients")) {
  stop("未知 MODEL_PRESET：", MODEL_PRESET)
}

CONTROL_SITE <- FALSE

CONTROL_SEASONAL_CYCLE <- FALSE

CONTROL_SIZE_IN_ALL <- FALSE

SHOW_CONTROL_PATHS_IN_DIAGRAM <- FALSE

GLOBAL_SEED <- 20260820L

SAVE_PREPARED_GROUP_DATA <- TRUE

PLOT_WIDTH <- 11

PLOT_HEIGHT <- 7.5

PLOT_DPI <- 220
