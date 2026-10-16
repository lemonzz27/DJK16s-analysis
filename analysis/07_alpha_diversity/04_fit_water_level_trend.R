required_packages <- c("ggplot2", "dplyr", "mgcv", "nlme")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
}

library(ggplot2)

library(dplyr)

library(mgcv)

library(nlme)

env_data <- read.csv("env.csv", stringsAsFactors = FALSE, check.names = FALSE)

group_month <- read.csv("groupmonth.csv", stringsAsFactors = FALSE)

colnames(group_month)[1:2] <- c("sample", "month_group")

matched <- group_month %>%
  inner_join(env_data %>%
    select(sample, WL), by = "sample") %>%
  mutate(WL = as.numeric(WL), Date = as.Date(paste0("20", month_group, "01"), format = "%Y%m%d")) %>%
  filter(!is.na(Date), is.finite(WL))

if (nrow(matched) != nrow(group_month)) {
  stop("水位与月份匹配不完整：匹配", nrow(matched), "行，应为", nrow(group_month),
    "行。")
}

monthly_check <- matched %>%
  group_by(Date, month_group) %>%
  summarise(n_sample_rows = n(), n_unique_WL = n_distinct(WL), WL = first(WL), .groups = "drop") %>%
  arrange(Date)

if (any(monthly_check$n_unique_WL != 1L)) {
  stop("同一月份出现多个不同水位值，不能直接建立单月水位序列。")
}

monthly_data <- monthly_check %>%
  mutate(time_num = as.numeric(Date - min(Date)), series = factor("Water level"))

if (nrow(monthly_data) != 15L) {
  stop("水位月份数不是预期的15个月，实际为", nrow(monthly_data), "个月。")
}

model_warnings <- character()

water_fit <- withCallingHandlers(mgcv::gamm(WL ~ s(time_num, bs = "cr", k = 6), data = monthly_data,
  method = "REML", correlation = nlme::corCAR1(form = ~time_num | series), control = nlme::lmeControl(opt = "optim",
    maxIter = 200, msMaxIter = 200, niterEM = 0, returnObject = TRUE)), warning = function(w) {
  model_warnings <<- c(model_warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
})

smooth_table <- summary(water_fit$gam)$s.table

smooth_edf <- unname(smooth_table[1, "edf"])

smooth_statistic <- if ("F" %in% colnames(smooth_table)) {
  unname(smooth_table[1, "F"])
} else {
  unname(smooth_table[1, "Chi.sq"])
}

smooth_p <- unname(smooth_table[1, "p-value"])

car1_phi <- tryCatch(unname(coef(water_fit$lme$modelStruct$corStruct, unconstrained = FALSE)[1]), error = function(e) NA_real_)

observed_fit <- as.numeric(predict(water_fit$gam, type = "response"))

pseudo_r2 <- cor(monthly_data$WL, observed_fit, use = "complete.obs")^2

prediction_data <- data.frame(Date = seq(min(monthly_data$Date), max(monthly_data$Date), length.out = 300))

prediction_data$time_num <- as.numeric(prediction_data$Date - min(monthly_data$Date))

prediction_data$series <- factor("Water level", levels = levels(monthly_data$series))

predicted <- predict(water_fit$gam, newdata = prediction_data, se.fit = TRUE)

prediction_data$fit <- as.numeric(predicted$fit)

prediction_data$se <- as.numeric(predicted$se.fit)

prediction_data$lower <- prediction_data$fit - 1.96 * prediction_data$se

prediction_data$upper <- prediction_data$fit + 1.96 * prediction_data$se

statistics <- data.frame(metric = "WL", model = "Gaussian GAMM", smoother = "cubic regression spline",
  method = "REML", k = 6L, n_months = nrow(monthly_data), first_month = min(monthly_data$Date), last_month = max(monthly_data$Date),
  smooth_edf = smooth_edf, smooth_statistic = smooth_statistic, smooth_p = smooth_p, car1_phi = car1_phi,
  pseudo_r2 = pseudo_r2, AIC_lme = AIC(water_fit$lme), model_warnings = paste(unique(model_warnings),
    collapse = " | "), stringsAsFactors = FALSE)

format_p <- function(p_value) {
  if (!is.finite(p_value))
    return("P=NA")
  if (p_value < 0.001)
    return("P<0.001")
  sprintf("P=%.3f", p_value)
}

subtitle_text <- paste0("GAMM smooth: ", format_p(smooth_p), "; k=6; REML; continuous-time CAR(1); n=15 months")

water_plot <- ggplot(monthly_data, aes(x = Date, y = WL)) + geom_ribbon(data = prediction_data, aes(x = Date,
  ymin = lower, ymax = upper), inherit.aes = FALSE, fill = "#BDBDBD", alpha = 0.45) + geom_line(data = prediction_data,
  aes(x = Date, y = fit), inherit.aes = FALSE, color = "#666666", linewidth = 1) + geom_point(color = "#666666",
  size = 2.2, alpha = 0.9) + scale_x_date(date_labels = "%Y-%m", date_breaks = "1 month") + theme_classic() +
  labs(title = "Trend of Water Level over Time", subtitle = subtitle_text, x = "Month", y = "Water level (m)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.subtitle = element_text(size = 8.5,
    lineheight = 1.05))

write.csv(monthly_data %>%
  select(Date, month_group, WL, n_sample_rows), "WaterLevel_Monthly_Data.csv", row.names = FALSE)

write.csv(statistics, "GAMM_WaterLevel_Model_Statistics.csv", row.names = FALSE)

write.csv(prediction_data %>%
  select(Date, time_num, fit, se, lower, upper), "GAMM_WaterLevel_Predictions.csv", row.names = FALSE)

ggsave("WaterLevel_Trend.png", plot = water_plot, width = 8, height = 6, dpi = 300)

ggsave("WaterLevel_Trend.pdf", plot = water_plot, width = 8, height = 6)

if (interactive()) print(water_plot)
