options(stringsAsFactors = FALSE, digits = 16)

set.seed(20260921)

args <- commandArgs(trailingOnly = TRUE)

out <- if (length(args)) args[1] else getwd()

d <- read.csv(file.path(out, "monthly_inputs.csv"), fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))

grid <- data.frame(month = format(seq(as.Date("1958-01-01"), as.Date("2026-04-01"), by = "month"), "%Y-%m"))

d <- merge(grid, d, by = "month", all.x = TRUE, sort = TRUE)

d$source[is.na(d$source)] <- "not_supplied"

d$status <- ifelse(is.na(d$Q), "missing", "observed")

historical <- d$month <= "2018-12" & is.finite(d$Q)

stopifnot(sum(historical) == 492L, all(d$Q[historical] > 0))

fit_gamma <- function(x) {
  x <- x/mean(x)
  a <- -mean(log(x))
  shape <- uniroot(function(k) log(k) - digamma(k) - a, c(1e-05, 1e+06), tol = 1e-11)$root
  shape
}

ks_stat <- function(x, shape, scale) {
  u <- sort(pgamma(x, shape = shape, scale = scale))
  n <- length(x)
  max(seq_len(n)/n - u, u - (seq_len(n) - 1)/n)
}

d$SRI <- d$SRI_empirical <- d$probability <- NA_real_

d$outside_reference_range <- NA

fits <- list()

for (m in 1:12) {
  ix <- as.integer(substr(d$month, 6, 7)) == m
  x <- d$Q[ix & historical]
  n <- length(x)
  stopifnot(n == 41)
  shape <- fit_gamma(x)
  scale <- mean(x)/shape
  q <- d$Q[ix]
  p <- pgamma(q, shape = shape, scale = scale)
  p <- pmin(1 - 1e-08, pmax(1e-08, p))
  d$probability[ix] <- p
  d$SRI[ix] <- qnorm(p)
  ep <- vapply(q, function(v) if (!is.finite(v))
    NA_real_
  else (sum(x < v) + 0.5 * sum(x == v) + 0.5)/(n + 1), numeric(1))
  d$SRI_empirical[ix] <- qnorm(ep)
  d$outside_reference_range[ix] <- q < min(x) | q > max(x)
  D <- ks_stat(x, shape, scale)
  boot <- replicate(499, {
    z <- rgamma(n, shape = shape, scale = scale)
    k <- fit_gamma(z)
    ks_stat(z, k, mean(z)/k)
  })
  bp <- (1 + sum(boot >= D))/500
  fits[[m]] <- data.frame(month = m, n = n, min_Q = min(x), max_Q = max(x), mean_Q = mean(x), gamma_shape = shape,
    gamma_scale = scale, KS_D = D, KS_bootstrap_p = bp, diagnostic = ifelse(bp < 0.05, "review_gamma_fit",
      "no_flag_not_proof"))
}

rdfai <- function(a, b) (b - a) * (abs(a) + abs(b)) * 3.2^(-abs(a + b))

stopifnot(rdfai(-1, 1) == 4, rdfai(1, -1) == -4, rdfai(1, 1) == 0, is.na(rdfai(NA, 1)))

d$previous_SRI <- c(NA, head(d$SRI, -1))

d$delta_SRI <- d$SRI - d$previous_SRI

d$RDFAI <- rdfai(d$previous_SRI, d$SRI)

d$RDFAI_empirical <- rdfai(c(NA, head(d$SRI_empirical, -1)), d$SRI_empirical)

d$state <- ifelse(is.na(d$SRI), NA, ifelse(d$SRI <= -0.5, "dry", ifelse(d$SRI >= 0.5, "wet", "normal")))

d$state_reversal <- ifelse(is.na(d$RDFAI), NA, ifelse(d$previous_SRI <= -0.5 & d$SRI >= 0.5, "dry_to_wet",
  ifelse(d$previous_SRI >= 0.5 & d$SRI <= -0.5, "wet_to_dry", "no_reversal")))

d$index_candidate <- ifelse(is.na(d$RDFAI), NA, ifelse(d$RDFAI > 1, "positive_gt1", ifelse(d$RDFAI <
  -1, "negative_lt_minus1", "within_threshold")))

d$event <- ifelse(is.na(d$RDFAI), NA, ifelse(d$state_reversal == "dry_to_wet" & d$RDFAI > 1, "dry_to_wet",
  ifelse(d$state_reversal == "wet_to_dry" & d$RDFAI < -1, "wet_to_dry", "not_event")))

d$metadata_status <- "historical_station_code_and_unit_unverified"

write.csv(d, file.path(out, "all_months.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

write.csv(d[is.finite(d$SRI), ], file.path(out, "available_months.csv"), row.names = FALSE, na = "",
  fileEncoding = "UTF-8")

write.csv(d[d$month >= "2024-11", ], file.path(out, "recent_months.csv"), row.names = FALSE, na = "",
  fileEncoding = "UTF-8")

write.csv(d[d$event %in% c("dry_to_wet", "wet_to_dry"), ], file.path(out, "state_reversal_events.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8")

write.csv(do.call(rbind, fits), file.path(out, "monthly_fit_diagnostics.csv"), row.names = FALSE, na = "",
  fileEncoding = "UTF-8")

capture.output(sessionInfo(), file = file.path(out, "R_session.txt"))

cat("SRI months:", sum(is.finite(d$SRI)), "RDFAI pairs:", sum(is.finite(d$RDFAI)), "\n")

print(d[d$month >= "2025-07" & d$month <= "2025-11", c("month", "Q", "SRI", "RDFAI", "event")])
