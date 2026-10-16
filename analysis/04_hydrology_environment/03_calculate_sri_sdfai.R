options(stringsAsFactors = FALSE, digits = 17)

Sys.setlocale("LC_ALL", "Chinese_China.utf8")

.libPaths(c(normalizePath("outputs/baihe_p3_all_data_20260922/Rlib", winslash = "/"), normalizePath("Rlib",
  winslash = "/"), .libPaths()))

library(lmom)

library(jsonlite)

out <- "outputs/baihe_sri_matrix_20260922"

d <- read.csv("outputs/baihe_monthly_indices_20260921/all_months.csv")

d <- d[, c("month", "Q", "source", "source_cell", "days_valid", "days_expected")]

d$calendar_month <- as.integer(substr(d$month, 6, 7))

stopifnot(nrow(d) == 820, sum(is.finite(d$Q)) == 510, all(d$Q[is.finite(d$Q)] > 0), identical(d$month,
  format(seq(as.Date("1958-01-01"), as.Date("2026-04-01"), by = "month"), "%Y-%m")))

meta <- data.frame(id = c("Gamma_MLE", "Gamma_Thom", "Gamma_Lmom", "P3_Lmom", "P3_MOM", "Lognormal_MLE",
  "Weibull_MLE", "Loglogistic_MLE", "GEV_Lmom", "Emp_Gringorten", "Emp_Weibull"), name = c("Gamma 数值MLE",
  "Gamma Thom近似", "Gamma L矩", "P-III L矩", "P-III 普通矩", "对数正态 MLE", "Weibull MLE",
  "对数Logistic MLE", "GEV L矩", "经验 Gringorten", "经验 Weibull秩"), family = c(rep("Gamma",
  3), rep("Pearson-III", 2), "Lognormal", "Weibull", "Loglogistic", "GEV", "Empirical", "Empirical"),
  estimation = c("数值极大似然", "Thom近似极大似然", "L矩", "L矩", "普通中心矩 分母n",
    "极大似然", "极大似然", "极大似然", "L矩", "(平均秩-0.44)/(n+0.12)", "平均秩/(n+1)"),
  params = c("shape;scale", "shape;scale", "shape;scale", "mean;sd;skew", "mean;sd;skew", "meanlog;sdlog",
    "shape;scale", "log_location;log_scale", "location;scale;shape(lmom)", "none", "none"))

fitdist <- function(x, id) {
  if (grepl("Gamma", id)) {
    A <- log(mean(x)) - mean(log(x))
    if (id == "Gamma_Lmom")
      return(unname(pelgam(samlmu(x, nmom = 2))))
    k <- if (id == "Gamma_Thom")
      (1 + sqrt(1 + 4 * A/3))/(4 * A)
    else uniroot(function(k) log(k) - digamma(k) - A, c(1e-06, 1e+07), tol = 1e-10)$root
    return(c(k, mean(x)/k))
  }
  if (id == "P3_Lmom")
    return(unname(pelpe3(samlmu(x, nmom = 3))))
  if (id == "P3_MOM") {
    mu <- mean(x)
    v <- mean((x - mu)^2)
    return(c(mu, sqrt(v), mean((x - mu)^3)/v^1.5))
  }
  if (id == "Lognormal_MLE") {
    lx <- log(x)
    return(c(mean(lx), sqrt(mean((lx - mean(lx))^2))))
  }
  if (id == "Loglogistic_MLE") {
    lx <- log(x)
    st <- c(median(lx), log(sd(lx) * sqrt(3)/pi))
    f <- optim(st, function(p) -sum(dlogis(lx, p[1], exp(p[2]), log = TRUE)), method = "BFGS", control = list(reltol = 1e-12,
      maxit = 1000))
    if (f$convergence != 0)
      stop("Loglogistic convergence")
    return(c(f$par[1], exp(f$par[2])))
  }
  if (id == "Weibull_MLE") {
    lx <- log(x)
    lc <- lx - max(lx)
    eq <- function(k) 1/k + mean(lx) - sum(lx * exp(k * lc))/sum(exp(k * lc))
    k <- uniroot(eq, c(1e-05, 10000), tol = 1e-10)$root
    lambda <- exp(max(lx) + log(mean(exp(k * lc)))/k)
    return(c(k, lambda))
  }
  if (id == "GEV_Lmom")
    return(unname(pelgev(samlmu(x, nmom = 3))))
  stop(id)
}

tails <- function(x, id, p) {
  if (grepl("Gamma", id))
    return(list(lp = pgamma(x, p[1], scale = p[2], log.p = TRUE), lq = pgamma(x, p[1], scale = p[2],
      log.p = TRUE, lower.tail = FALSE)))
  if (grepl("P3", id)) {
    if (abs(p[3]) < 1e-10)
      return(list(lp = pnorm(x, p[1], p[2], log.p = TRUE), lq = pnorm(x, p[1], p[2], log.p = TRUE,
        lower.tail = FALSE)))
    sh <- 4/p[3]^2
    sc <- p[2] * abs(p[3])/2
    loc <- p[1] - 2 * p[2]/p[3]
    y <- sign(p[3]) * (x - loc)
    return(list(lp = pgamma(y, sh, scale = sc, lower.tail = p[3] > 0, log.p = TRUE), lq = pgamma(y,
      sh, scale = sc, lower.tail = p[3] < 0, log.p = TRUE)))
  }
  if (id == "Lognormal_MLE")
    return(list(lp = plnorm(x, p[1], p[2], log.p = TRUE), lq = plnorm(x, p[1], p[2], log.p = TRUE,
      lower.tail = FALSE)))
  if (id == "Weibull_MLE")
    return(list(lp = pweibull(x, p[1], p[2], log.p = TRUE), lq = pweibull(x, p[1], p[2], log.p = TRUE,
      lower.tail = FALSE)))
  if (id == "Loglogistic_MLE")
    return(list(lp = plogis(log(x), p[1], p[2], log.p = TRUE), lq = plogis(log(x), p[1], p[2], log.p = TRUE,
      lower.tail = FALSE)))
  if (id == "GEV_Lmom") {
    z <- (x - p[1])/p[2]
    k <- p[3]
    if (abs(k) < 1e-10)
      y <- z
    else {
      t <- 1 - k * z
      y <- rep(NA_real_, length(z))
      ok <- t > 0
      y[ok] <- -log(t[ok])/k
      y[!ok] <- if (k > 0)
        Inf
      else -Inf
    }
    lp <- -exp(-y)
    lq <- log(-expm1(lp))
    return(list(lp = lp, lq = lq))
  }
  stop(id)
}

cdf <- function(x, id, p) exp(tails(x, id, p)$lp)

rmodel <- function(n, id, p) {
  if (grepl("Gamma", id))
    return(rgamma(n, p[1], scale = p[2]))
  if (grepl("P3", id))
    return(quape3(runif(n), p))
  switch(id, Lognormal_MLE = rlnorm(n, p[1], p[2]), Weibull_MLE = rweibull(n, p[1], p[2]), Loglogistic_MLE = exp(rlogis(n,
    p[1], p[2])), GEV_Lmom = quagev(runif(n), p))
}

logdensity <- function(x, id, p) switch(id, Gamma_MLE = dgamma(x, p[1], scale = p[2], log = TRUE), Lognormal_MLE = dlnorm(x,
  p[1], p[2], log = TRUE), Weibull_MLE = dweibull(x, p[1], p[2], log = TRUE), Loglogistic_MLE = dlogis(log(x),
  p[1], p[2], log = TRUE) - log(x), NULL)

KS <- function(x, id, p) {
  v <- cdf(sort(x), id, p)
  n <- length(x)
  max(seq_len(n)/n - v, v - (seq_len(n) - 1)/n)
}

zfrom <- function(t) ifelse(t$lp <= log(0.5), qnorm(t$lp, log.p = TRUE), qnorm(t$lq, lower.tail = FALSE,
  log.p = TRUE))

prev <- function(x) c(NA, head(x, -1))

direction <- function(a, b) {
  v <- rep(NA_integer_, length(a))
  ok <- is.finite(a) & is.finite(b)
  v[ok] <- 0L
  v[which(ok & a <= -0.5 & b >= 0.5)] <- 1L
  v[which(ok & a >= 0.5 & b <= -0.5)] <- -1L
  v
}

label <- function(x) ifelse(is.na(x), "不可计算", ifelse(x > 0, "枯转丰", ifelse(x < 0, "丰转枯",
  "未识别")))

event_threshold <- function(x, inclusive = FALSE) label(ifelse(is.finite(x), ifelse(if (inclusive) x >=
  1 else x > 1, 1, ifelse(if (inclusive) x <= -1 else x < -1, -1, 0)), NA_real_))

stopifnot(identical(event_threshold(c(-1, 1)), c("未识别", "未识别")), identical(event_threshold(c(-1,
  1), TRUE), c("丰转枯", "枯转丰")))

rsdfai <- function(a, b) {
  v <- rep(NA_real_, length(a))
  ok <- which(is.finite(a) & is.finite(b))
  v[ok] <- 0
  k <- ok[a[ok] != b[ok] & a[ok] != 0 & b[ok] != 0]
  aa <- abs(a[k])
  bb <- abs(b[k])
  f1 <- b[k] - a[k]
  f2 <- aa + bb
  I <- (abs(f1)/f2)^abs(a[k] + b[k]) * pmin(aa, bb)
  M <- pmax(aa, bb)
  m <- pmin(aa, bb)
  D <- abs(f1) + f2
  ip <- (I/0.5)^(M^2/D) * (I^(M/D) + I^(m/D))/2
  ip0 <- (0.5^(M/D) + 0.5^(m/D))/2
  v[k] <- sign(f1) * ((ip/ip0) * (I/0.5))^((M/D)^(1 - M/D))
  v
}

stopifnot(abs(rsdfai(-0.5, 0.5) - 1) < 1e-12, abs(rsdfai(-1, 1) - 1.444) < 0.001)

B <- 499L

set.seed(20260922)

fits <- list()

SRI <- list()

indices <- list()

summaries <- list()

events <- list()

failure <- list()

nr <- 0

for (j in seq_len(nrow(meta))) {
  id <- meta$id[j]
  parametric <- j <= 9
  z <- prob <- rep(NA_real_, nrow(d))
  status <- rep("原始流量缺测", nrow(d))
  for (m in 1:12) {
    ix <- which(d$calendar_month == m & is.finite(d$Q))
    x <- d$Q[ix]
    n <- length(x)
    if (!parametric) {
      pp <- if (id == "Emp_Gringorten")
        (rank(x, ties.method = "average") - 0.44)/(n + 0.12)
      else rank(x, ties.method = "average")/(n + 1)
      prob[ix] <- pp
      z[ix] <- qnorm(pp)
      status[ix] <- "可计算"
      next
    }
    p <- fitdist(x, id)
    tt <- tails(x, id, p)
    zz <- zfrom(tt)
    prob[ix] <- exp(tt$lp)
    z[ix] <- zz
    status[ix] <- ifelse(is.finite(zz), "可计算", ifelse(is.infinite(zz), "分布支持域边界 SRI无穷",
      "拟合或数值失败"))
    if (grepl("P3", id))
      stopifnot(max(abs(cdf(x, id, p) - cdfpe3(x, p))) < 1e-07)
    if (id == "GEV_Lmom")
      stopifnot(max(abs(cdf(x, id, p) - cdfgev(x, p))) < 1e-07)
    obs <- KS(x, id, p)
    bs <- replicate(B, tryCatch({
      sim <- rmodel(n, id, p)
      ps <- fitdist(sim, id)
      KS(sim, id, ps)
    }, error = function(e) NA_real_))
    bp <- if (all(is.finite(bs)))
      (1 + sum(bs >= obs))/(B + 1)
    else NA_real_
    v <- sort(prob[ix])
    cvm <- 1/(12 * n) + sum((v - (seq_len(n) - 0.5)/n)^2)
    ld <- logdensity(x, id, p)
    ll <- if (is.null(ld))
      NA_real_
    else sum(ld)
    aicc <- if (is.finite(ll))
      -2 * ll + 4 + 12/(n - 3)
    else NA_real_
    nr <- nr + 1
    fits[[nr]] <- data.frame(method = id, name = meta$name[j], calendar_month = m, n = n, param1 = p[1],
      param2 = p[2], param3 = if (length(p) >= 3)
        p[3]
      else NA_real_, parameter_definition = meta$params[j], KS = obs, CvM = cvm, bootstrap_p = bp,
      bootstrap_failures = sum(!is.finite(bs)), negative_Q_probability = cdf(0, id, p), nonfinite_SRI = sum(!is.finite(zz)),
      loglik = ll, AICc = aicc)
  }
  finite <- is.finite(z)
  nonfinite_obs <- which(is.finite(d$Q) & !finite)
  if (length(nonfinite_obs))
    failure[[length(failure) + 1]] <- data.frame(method = id, name = meta$name[j], month = d$month[nonfinite_obs],
      Q = d$Q[nonfinite_obs], CDF = prob[nonfinite_obs], raw_SRI = as.character(z[nonfinite_obs]),
      reason = status[nonfinite_obs])
  sa <- if (length(nonfinite_obs) == 0)
    diff(range(z, na.rm = TRUE))
  else NA_real_
  clean <- z
  clean[!finite] <- NA_real_
  a <- prev(clean)
  b <- clean
  di <- direction(a, b)
  ok <- is.finite(a) & is.finite(b)
  sw <- rep(NA_real_, length(z))
  if (is.finite(sa))
    sw[ok] <- ifelse(di[ok] == 0, 0, (b[ok] - a[ok])/2 * (1 + abs(a[ok] + b[ok])/sa))
  sf <- (b - a) * (abs(a) + abs(b)) * 3.2^(-abs(a + b))
  rs <- rsdfai(a, b)
  swlab <- if (is.finite(sa))
    label(di)
  else rep("不可计算", length(z))
  SRI[[j]] <- data.frame(method = id, name = meta$name[j], month = d$month, Q = d$Q, CDF = prob, SRI = clean,
    state = ifelse(!finite, "不可计算", ifelse(z <= -0.5, "枯", ifelse(z >= 0.5, "丰", "正常"))),
    status = status)
  ii <- data.frame(method = id, name = meta$name[j], month = d$month, previous_SRI = a, SRI = b, delta_SRI = b -
    a, Sa = sa, SDWAI = sw, SDFAI = sf, RSDFAI = rs, SDWAI_event = swlab, SDFAI_event = event_threshold(sf),
    RSDFAI_event = event_threshold(rs, TRUE), endpoints = label(di), status = ifelse(!is.finite(d$Q),
      "原始流量缺测", ifelse(!finite, "本月SRI不可计算", ifelse(!is.finite(a), "上月SRI不可计算",
        "相邻SRI可计算"))), SDWAI_status = if (is.finite(sa))
      "Sa有效"
    else "全序列有非有限SRI，Sa不可定义")
  indices[[j]] <- ii
  for (iname in c("SDWAI", "SDFAI", "RSDFAI")) {
    lab <- ii[[paste0(iname, "_event")]]
    e <- which(lab %in% c("枯转丰", "丰转枯"))
    recent <- d$month >= "2024-11"
    if (length(e))
      events[[length(events) + 1]] <- data.frame(method = id, name = meta$name[j], index = iname,
        start_month = d$month[e - 1], end_month = d$month[e], previous_SRI = a[e], SRI = b[e],
        value = ii[[iname]][e], direction = lab[e], endpoints = ii$endpoints[e], method_warning = if (length(nonfinite_obs))
          "本版部分SRI支持域失败，仅局部值可计算"
        else "")
    summaries[[length(summaries) + 1]] <- data.frame(method = id, name = meta$name[j], index = iname,
      finite_SRI = sum(finite), failed_SRI = length(nonfinite_obs), valid_pairs = sum(is.finite(ii[[iname]])),
      dry_to_wet = sum(lab == "枯转丰"), wet_to_dry = sum(lab == "丰转枯"), recent_valid_pairs = sum(is.finite(ii[[iname]]) &
        recent), recent_dry_to_wet = sum(lab == "枯转丰" & recent), recent_wet_to_dry = sum(lab ==
        "丰转枯" & recent), Sa = sa)
  }
  cat("Completed", id, "finite SRI", sum(finite), "invalid", length(nonfinite_obs), "\n")
  flush.console()
}

fits <- do.call(rbind, fits)

fits$Holm_p <- p.adjust(fits$bootstrap_p, "holm")

SRI <- do.call(rbind, SRI)

indices <- do.call(rbind, indices)

summaries <- do.call(rbind, summaries)

events <- do.call(rbind, events)

failure <- do.call(rbind, failure)

summaries[summaries$valid_pairs == 0, c("dry_to_wet", "wet_to_dry")] <- NA_integer_

summaries[summaries$recent_valid_pairs == 0, c("recent_dry_to_wet", "recent_wet_to_dry")] <- NA_integer_

events <- events[order(events$end_month, decreasing = TRUE), ]

quality <- do.call(rbind, lapply(seq_len(nrow(meta)), function(j) {
  id <- meta$id[j]
  s <- SRI[SRI$method == id, ]
  f <- fits[fits$method == id, ]
  data.frame(method = id, name = meta$name[j], family = meta$family[j], estimation = meta$estimation[j],
    finite_SRI = sum(is.finite(s$SRI)), failed_SRI = sum(is.finite(s$Q) & !is.finite(s$SRI)), min_SRI = min(s$SRI,
      na.rm = TRUE), max_SRI = max(s$SRI, na.rm = TRUE), sum_AICc = if (nrow(f) == 12 && all(is.finite(f$AICc)))
      sum(f$AICc)
    else NA_real_, mean_KS = if (nrow(f))
      mean(f$KS)
    else NA_real_, KS_Holm_lt_05 = if (nrow(f))
      sum(f$Holm_p < 0.05, na.rm = TRUE)
    else NA_integer_, max_negative_Q_probability = if (nrow(f))
      max(f$negative_Q_probability)
    else NA_real_)
}))

wide <- d

for (id in meta$id) {
  s <- SRI[SRI$method == id, ]
  wide[[id]] <- s$SRI
  wide[[paste0(id, "_status")]] <- s$status
}

write_csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")

write_csv(wide, "SRI_all_months.csv")

write_csv(wide[wide$month >= "2024-11", ], "SRI_recent.csv")

write_csv(SRI, "SRI_long.csv")

write_csv(indices, "all_indices.csv")

write_csv(events, "events.csv")

write_csv(fits, "monthly_fits.csv")

write_csv(failure, "support_failures.csv")

write_csv(summaries, "comparison.csv")

write_csv(quality, "SRI_quality.csv")

write_json(list(meta = meta, data = d, sri = SRI, wide = wide, indices = indices, events = events, fits = fits,
  failures = failure, comparison = summaries, quality = quality, bootstrap = B, seed = 20260922, R = R.version.string,
  lmom = as.character(packageVersion("lmom"))), file.path(out, "analysis.json"), dataframe = "rows",
  auto_unbox = TRUE, na = "null", digits = 16, pretty = TRUE)

capture.output(sessionInfo(), file = file.path(out, "R_session.txt"))

print(quality, row.names = FALSE)

print(wide[wide$month %in% c("2025-07", "2025-08", "2025-09"), c("month", meta$id)], row.names = FALSE)

print(summaries[, c("name", "index", "valid_pairs", "recent_dry_to_wet", "recent_wet_to_dry")], row.names = FALSE)
