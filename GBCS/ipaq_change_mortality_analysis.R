#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(survival)
})

options(stringsAsFactors = FALSE, scipen = 999)

main_path <- Sys.getenv("GBCS_MAIN_DTA", "/Users/linxu/Documents/GBCS/gbcs_main.dta")
mort_path <- Sys.getenv("GBCS_MORTALITY_DTA", "/Users/linxu/Documents/GBCS/newvital_2023.dta")
out_dir <- Sys.getenv("IPAQ_CHANGE_MORTALITY_OUTDIR", file.path("output", "ipaq_change_mortality"))

if (!file.exists(main_path)) stop("Main dataset not found: ", main_path)
if (!file.exists(mort_path)) stop("Mortality dataset not found: ", mort_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(paste0(...), "\n")
  flush.console()
}

as_num <- function(x) {
  if (inherits(x, "haven_labelled")) x <- unclass(x)
  suppressWarnings(as.numeric(x))
}

normalize_id <- function(x) {
  if (is.numeric(x)) {
    y <- ifelse(is.na(x), NA_character_, sprintf("%.0f", x))
  } else {
    y <- as.character(x)
    y <- sub("\\.0+$", "", y)
  }
  y <- trimws(y)
  y[y == ""] <- NA_character_
  y
}

recode_ipaq_zero_based <- function(x) {
  x <- as_num(x)
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x) & x %in% c(0, 1, 2)
  out[ok] <- x[ok] + 1
  out
}

recode_ipaq_one_based <- function(x) {
  x <- as_num(x)
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x) & x %in% c(1, 2, 3)
  out[ok] <- x[ok]
  out
}

compute_ipaq_category_from_v7 <- function(df, suffix = "") {
  getv <- function(stem) as_num(df[[paste0(stem, suffix)]])

  v7_1a <- getv("v7_1a")
  v7_1b <- getv("v7_1b")
  v7_2a <- getv("v7_2a")
  v7_2b <- getv("v7_2b")
  v7_3a <- getv("v7_3a")
  v7_3b <- getv("v7_3b")
  v7_4a <- getv("v7_4a")
  v7_4b <- getv("v7_4b")
  v7_5a <- getv("v7_5a")
  v7_5b <- getv("v7_5b")
  v7_6a <- getv("v7_6a")
  v7_6b <- getv("v7_6b")

  any_info <- rowSums(!is.na(cbind(
    v7_1a, v7_1b, v7_2a, v7_2b,
    v7_3a, v7_3b, v7_4a, v7_4b,
    v7_5a, v7_5b, v7_6a, v7_6b
  ))) > 0

  use_vig <- !is.na(v7_1a) & v7_1a == 1
  use_mod <- !is.na(v7_3a) & v7_3a == 1
  use_walk <- !is.na(v7_5a) & v7_5a == 1

  calc_minutes_week <- function(days, min_day, min_week, eligible) {
    out <- rep(0, length(days))
    out[eligible & !is.na(min_week)] <- min_week[eligible & !is.na(min_week)]
    derive <- eligible & is.na(min_week) & !is.na(days) & !is.na(min_day)
    out[derive] <- days[derive] * min_day[derive]
    out[eligible & out < 0] <- NA_real_
    out
  }

  calc_minutes_day <- function(days, min_day, min_week, eligible) {
    out <- rep(0, length(days))
    out[eligible & !is.na(min_day)] <- min_day[eligible & !is.na(min_day)]
    derive <- eligible & is.na(min_day) & !is.na(min_week) & !is.na(days) & days > 0
    out[derive] <- min_week[derive] / days[derive]
    out[eligible & out < 0] <- NA_real_
    out
  }

  d_vig <- ifelse(use_vig & !is.na(v7_1b), v7_1b, 0)
  d_mod <- ifelse(use_mod & !is.na(v7_3b), v7_3b, 0)
  d_walk <- ifelse(use_walk & !is.na(v7_5b), v7_5b, 0)

  mw_vig <- calc_minutes_week(v7_1b, v7_2a, v7_2b, use_vig)
  mw_mod <- calc_minutes_week(v7_3b, v7_4a, v7_4b, use_mod)
  mw_walk <- calc_minutes_week(v7_5b, v7_6a, v7_6b, use_walk)

  md_vig <- calc_minutes_day(v7_1b, v7_2a, v7_2b, use_vig)
  md_mod <- calc_minutes_day(v7_3b, v7_4a, v7_4b, use_mod)
  md_walk <- calc_minutes_day(v7_5b, v7_6a, v7_6b, use_walk)

  total_days <- d_vig + d_mod + d_walk
  total_met <- 8 * mw_vig + 4 * mw_mod + 3.3 * mw_walk
  vig_met <- 8 * mw_vig

  is_high <- (d_vig >= 3 & vig_met >= 1500) | (total_days >= 7 & total_met >= 3000)
  is_moderate <- (d_vig >= 3 & md_vig >= 20) |
    (((d_mod >= 5 & md_mod >= 30) | (d_walk >= 5 & md_walk >= 30))) |
    (total_days >= 5 & total_met >= 600)

  out <- rep(NA_real_, nrow(df))
  out[any_info] <- 1
  out[is_moderate] <- 2
  out[is_high] <- 3
  out
}

compute_met_wave <- function(df, suffix = "") {
  getv <- function(stem) as_num(df[[paste0(stem, suffix)]])

  v7_1a <- getv("v7_1a")
  v7_1b <- getv("v7_1b")
  v7_2a <- getv("v7_2a")
  v7_2b <- getv("v7_2b")
  v7_3a <- getv("v7_3a")
  v7_3b <- getv("v7_3b")
  v7_4a <- getv("v7_4a")
  v7_4b <- getv("v7_4b")
  v7_5a <- getv("v7_5a")
  v7_5b <- getv("v7_5b")
  v7_6a <- getv("v7_6a")
  v7_6b <- getv("v7_6b")

  any_info <- rowSums(!is.na(cbind(
    v7_1a, v7_1b, v7_2a, v7_2b,
    v7_3a, v7_3b, v7_4a, v7_4b,
    v7_5a, v7_5b, v7_6a, v7_6b
  ))) > 0

  use_vig <- !is.na(v7_1a) & v7_1a == 1
  use_mod <- !is.na(v7_3a) & v7_3a == 1
  use_walk <- !is.na(v7_5a) & v7_5a == 1

  calc_minutes_week <- function(days, min_day, min_week, eligible) {
    out <- rep(0, length(days))
    out[eligible & !is.na(min_week)] <- min_week[eligible & !is.na(min_week)]
    derive <- eligible & is.na(min_week) & !is.na(days) & !is.na(min_day)
    out[derive] <- days[derive] * min_day[derive]
    out[eligible & out < 0] <- NA_real_
    out
  }

  mw_vig <- calc_minutes_week(v7_1b, v7_2a, v7_2b, use_vig)
  mw_mod <- calc_minutes_week(v7_3b, v7_4a, v7_4b, use_mod)
  mw_walk <- calc_minutes_week(v7_5b, v7_6a, v7_6b, use_walk)

  total_met <- 8 * mw_vig + 4 * mw_mod + 3.3 * mw_walk
  total_met[!any_info] <- NA_real_
  total_met
}

extract_cox <- function(fit, model_id, interval_id, exposure_type) {
  s <- summary(fit)
  coef_mat <- s$coefficients
  ci_mat <- s$conf.int
  data.frame(
    interval = interval_id,
    exposure_type = exposure_type,
    model_id = model_id,
    term = rownames(coef_mat),
    coef = coef_mat[, "coef"],
    se = coef_mat[, "se(coef)"],
    z = coef_mat[, "z"],
    p_value = coef_mat[, "Pr(>|z|)"],
    hr = ci_mat[, "exp(coef)"],
    ci_lower = ci_mat[, "lower .95"],
    ci_upper = ci_mat[, "upper .95"],
    n = fit$n,
    nevent = fit$nevent,
    row.names = NULL
  )
}

fit_stats <- function(fit, model_id, interval_id, exposure_type) {
  s <- summary(fit)
  data.frame(
    interval = interval_id,
    exposure_type = exposure_type,
    model_id = model_id,
    n = fit$n,
    nevent = fit$nevent,
    loglik_null = fit$loglik[1],
    loglik_fit = fit$loglik[2],
    aic = AIC(fit),
    concordance = if (!is.null(s$concordance)) s$concordance[1] else NA_real_,
    concordance_se = if (!is.null(s$concordance)) s$concordance[2] else NA_real_,
    row.names = NULL
  )
}

ph_table <- function(fit, model_id, interval_id, exposure_type) {
  z <- cox.zph(fit)
  tab <- as.data.frame(z$table, stringsAsFactors = FALSE)
  tab$term <- rownames(tab)
  names(tab) <- sub("^p$", "p_value", names(tab))
  tab$interval <- interval_id
  tab$exposure_type <- exposure_type
  tab$model_id <- model_id
  tab
}

fit_cox_safe <- function(formula, data, model_id, interval_id, exposure_type) {
  out <- list(coef = NULL, stats = NULL, ph = NULL, error = NA_character_)
  fit <- tryCatch(
    coxph(formula, data = data, x = TRUE, y = TRUE, model = TRUE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    out$error <- fit$message
    return(out)
  }
  out$coef <- extract_cox(fit, model_id, interval_id, exposure_type)
  out$stats <- fit_stats(fit, model_id, interval_id, exposure_type)
  if (model_id == "M3") {
    ph <- tryCatch(ph_table(fit, model_id, interval_id, exposure_type), error = function(e) NULL)
    out$ph <- ph
  }
  out
}

fixed_effect_pool <- function(df, effect_terms) {
  rows <- list()
  i <- 1L
  for (term_i in effect_terms) {
    sub <- df %>% filter(term == term_i, !is.na(coef), !is.na(se), se > 0)
    if (nrow(sub) == 0) next
    w <- 1 / (sub$se^2)
    beta <- sum(w * sub$coef) / sum(w)
    se <- sqrt(1 / sum(w))
    rows[[i]] <- data.frame(
      term = term_i,
      n_intervals = nrow(sub),
      coef = beta,
      se = se,
      hr = exp(beta),
      ci_lower = exp(beta - 1.96 * se),
      ci_upper = exp(beta + 1.96 * se),
      z = beta / se,
      p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
      row.names = NULL
    )
    i <- i + 1L
  }
  bind_rows(rows)
}

log_msg("Step 1: loading selected columns")
main_meta <- read_dta(main_path, n_max = 0)
main_names <- names(main_meta)

v7_base <- c("v7_1a", "v7_1b", "v7_2a", "v7_2b", "v7_3a", "v7_3b", "v7_4a", "v7_4b", "v7_5a", "v7_5b", "v7_6a", "v7_6b")
need_main <- unique(c(
  "obje_id",
  "reg_date", "reg_date_f", "reg_date_f2", "reg_date_f3",
  "ipaq", "ipaq_f2", "ipaq_g3", "ipaq_f3",
  "agec", "sex", "phase", "edu", "incomeh", "marital", "smk", "drk1",
  "hyp", "diab", "hxchd", "hxstroke", "v12_1a",
  v7_base,
  paste0(v7_base, "_f"),
  paste0(v7_base, "_f2"),
  paste0(v7_base, "_f3")
))
need_main <- need_main[need_main %in% main_names]

main_df <- read_dta(main_path, col_select = all_of(need_main))
mort_df <- read_dta(mort_path, col_select = all_of(c("obje_id", "dead", "futime", "ddate", "reg_date")))

main_df$obje_id <- normalize_id(main_df$obje_id)
mort_df$obje_id <- normalize_id(mort_df$obje_id)

if (sum(duplicated(main_df$obje_id)) > 0) stop("Duplicate obje_id in main file.")
if (sum(duplicated(mort_df$obje_id)) > 0) stop("Duplicate obje_id in mortality file.")

log_msg("Step 2: merging files and deriving exposure variables")
dat <- main_df %>%
  left_join(mort_df, by = "obje_id", suffix = c("", "_mort")) %>%
  mutate(
    reg_date = as.Date(reg_date),
    reg_date_f = as.Date(reg_date_f),
    reg_date_f2 = as.Date(reg_date_f2),
    reg_date_f3 = as.Date(reg_date_f3),
    ddate = as.Date(ddate),
    dead_num = as_num(dead),
    event_allcause = case_when(dead_num == 1 ~ 1L, dead_num == 0 ~ 0L, TRUE ~ NA_integer_),
    futime = as_num(futime),
    IPAQ_cat_bl = recode_ipaq_zero_based(ipaq)
  ) %>%
  mutate(
    agec = as_num(agec),
    sex = as_num(sex),
    phase = as_num(phase),
    edu = as_num(edu),
    incomeh = as_num(incomeh),
    marital = as_num(marital),
    smk = as_num(smk),
    drk1 = as_num(drk1),
    hyp = as_num(hyp),
    diab = as_num(diab),
    hxchd = as_num(hxchd),
    hxstroke = as_num(hxstroke),
    v12_1a = as_num(v12_1a),
    agec_5 = agec / 5,
    sex_fct = factor(ifelse(sex %in% c(0, 1), sex, NA), levels = c(0, 1), labels = c("female", "male")),
    phase_fct = factor(ifelse(phase %in% c(1, 2, 3), phase, NA), levels = c(1, 2, 3)),
    edu_fct = factor(ifelse(edu %in% 1:6, edu, NA), levels = c(6, 5, 4, 3, 2, 1),
                     labels = c("college", "junior_college", "senior_middle", "junior_middle", "primary", "lt_primary")),
    incomeh_fct = factor(ifelse(incomeh %in% 1:7, incomeh, NA), levels = c(6, 5, 4, 3, 2, 1, 7),
                         labels = c("ge50000", "30000_band", "20000_band", "10000_band", "5000_band", "lt5000", "notknow")),
    marital_fct = factor(ifelse(marital %in% c(0, 1, 2, 3), marital, NA), levels = c(0, 1, 2, 3),
                         labels = c("married", "separate", "widow", "never_married")),
    smk_fct = factor(ifelse(smk %in% c(0, 1, 2, 3), smk, NA), levels = c(0, 1, 2, 3),
                     labels = c("never", "ex_occasional", "ex_daily", "current")),
    drk1_fct = factor(ifelse(drk1 %in% 0:5, drk1, NA), levels = c(0, 1, 2, 3, 4, 5),
                      labels = c("never", "lt1_month", "lt1_week", "1to4_week", "ge5_week", "ex_drinker")),
    hyp_fct = factor(ifelse(hyp %in% c(0, 1), hyp, NA), levels = c(0, 1), labels = c("no", "yes")),
    diab_fct = factor(ifelse(diab %in% c(0, 1), diab, NA), levels = c(0, 1), labels = c("no", "yes")),
    hxchd_fct = factor(ifelse(hxchd %in% c(0, 1), hxchd, NA), levels = c(0, 1), labels = c("no", "yes")),
    hxstroke_fct = factor(ifelse(hxstroke %in% c(0, 1), hxstroke, NA), levels = c(0, 1), labels = c("no", "yes")),
    cancerhx_bin = case_when(v12_1a == 1 ~ 0L, v12_1a == 2 ~ 1L, TRUE ~ NA_integer_),
    cancerhx_fct = factor(cancerhx_bin, levels = c(0, 1), labels = c("no", "yes"))
  )

dat$IPAQ_cat_f1 <- compute_ipaq_category_from_v7(dat, suffix = "_f")
ipaq_f2_one <- recode_ipaq_one_based(dat$ipaq_f2)
ipaq_f2_zero <- recode_ipaq_zero_based(dat$ipaq_g3)
dat$IPAQ_cat_f2 <- ifelse(!is.na(ipaq_f2_one), ipaq_f2_one, ipaq_f2_zero)
ipaq_f3_zero <- recode_ipaq_zero_based(dat$ipaq_f3)
ipaq_f3_from_v7 <- compute_ipaq_category_from_v7(dat, suffix = "_f3")
dat$IPAQ_cat_f3 <- ifelse(!is.na(ipaq_f3_zero), ipaq_f3_zero, ipaq_f3_from_v7)
dat$MET_bl <- compute_met_wave(dat, suffix = "")
dat$MET_f1 <- compute_met_wave(dat, suffix = "_f")
dat$MET_f2 <- compute_met_wave(dat, suffix = "_f2")
dat$MET_f3 <- compute_met_wave(dat, suffix = "_f3")

flow <- data.frame(
  step = c(
    "main_rows",
    "matched_to_mortality",
    "valid_event_allcause",
    "valid_futime_positive",
    "valid_baseline_ipaq",
    "valid_baseline_met"
  ),
  n = c(
    nrow(main_df),
    sum(!is.na(dat$dead_num)),
    sum(!is.na(dat$event_allcause)),
    sum(!is.na(dat$futime) & dat$futime > 0),
    sum(!is.na(dat$IPAQ_cat_bl)),
    sum(!is.na(dat$MET_bl))
  )
)
write.csv(flow, file.path(out_dir, "cohort_flow.csv"), row.names = FALSE)

wave_summary <- bind_rows(lapply(c("bl", "f1", "f2", "f3"), function(w) {
  ipaq_v <- dat[[paste0("IPAQ_cat_", w)]]
  met_v <- dat[[paste0("MET_", w)]]
  cat_tab <- table(factor(ipaq_v, levels = c(1, 2, 3), labels = c("low", "moderate", "high")), useNA = "no")
  data.frame(
    wave = w,
    n_ipaq_nonmissing = sum(!is.na(ipaq_v)),
    n_ipaq_low = if ("low" %in% names(cat_tab)) as.integer(cat_tab[["low"]]) else 0L,
    n_ipaq_moderate = if ("moderate" %in% names(cat_tab)) as.integer(cat_tab[["moderate"]]) else 0L,
    n_ipaq_high = if ("high" %in% names(cat_tab)) as.integer(cat_tab[["high"]]) else 0L,
    n_met_nonmissing = sum(!is.na(met_v)),
    met_mean = ifelse(sum(!is.na(met_v)) > 0, mean(met_v, na.rm = TRUE), NA_real_),
    met_sd = ifelse(sum(!is.na(met_v)) > 1, sd(met_v, na.rm = TRUE), NA_real_),
    met_median = ifelse(sum(!is.na(met_v)) > 0, median(met_v, na.rm = TRUE), NA_real_),
    met_q1 = ifelse(sum(!is.na(met_v)) > 0, as.numeric(quantile(met_v, 0.25, na.rm = TRUE, names = FALSE)), NA_real_),
    met_q3 = ifelse(sum(!is.na(met_v)) > 0, as.numeric(quantile(met_v, 0.75, na.rm = TRUE, names = FALSE)), NA_real_),
    stringsAsFactors = FALSE
  )
}))
write.csv(wave_summary, file.path(out_dir, "wave_activity_summary.csv"), row.names = FALSE)

build_transition <- function(df, bl_var, fu_var, label) {
  tmp <- tibble(
    bl = df[[bl_var]],
    fu = df[[fu_var]]
  ) %>%
    mutate(
      bl = ifelse(bl %in% 1:3, bl, NA_real_),
      fu = ifelse(fu %in% 1:3, fu, NA_real_)
    ) %>%
    filter(!is.na(bl), !is.na(fu))
  if (nrow(tmp) == 0) return(tibble())
  tab <- table(
    factor(tmp$bl, levels = 1:3, labels = c("low", "moderate", "high")),
    factor(tmp$fu, levels = 1:3, labels = c("low", "moderate", "high"))
  )
  as.data.frame(tab, stringsAsFactors = FALSE) %>%
    rename(from = Var1, to = Var2, n = Freq) %>%
    group_by(from) %>%
    mutate(row_pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    mutate(interval = label, total_pct = 100 * n / sum(n))
}

transition_tbl <- bind_rows(
  build_transition(dat, "IPAQ_cat_bl", "IPAQ_cat_f1", "BL_to_F1"),
  build_transition(dat, "IPAQ_cat_f1", "IPAQ_cat_f2", "F1_to_F2"),
  build_transition(dat, "IPAQ_cat_f2", "IPAQ_cat_f3", "F2_to_F3")
)
write.csv(transition_tbl, file.path(out_dir, "ipaq_transition_tables.csv"), row.names = FALSE)

missingness_vars <- c(
  "event_allcause", "futime",
  "IPAQ_cat_bl", "IPAQ_cat_f1", "IPAQ_cat_f2", "IPAQ_cat_f3",
  "MET_bl", "MET_f1", "MET_f2", "MET_f3",
  "agec", "sex", "phase", "edu", "incomeh", "marital", "smk", "drk1", "hyp", "diab", "hxchd", "hxstroke", "v12_1a"
)
missingness <- bind_rows(lapply(missingness_vars, function(v) {
  x <- dat[[v]]
  data.frame(
    variable = v,
    n = length(x),
    n_missing = sum(is.na(x)),
    pct_missing = 100 * mean(is.na(x)),
    n_non_missing = sum(!is.na(x))
  )
}))
write.csv(missingness, file.path(out_dir, "variable_missingness.csv"), row.names = FALSE)

build_interval_data <- function(df, prior_cat, follow_cat, prior_met, follow_met, landmark_date_var, interval_id) {
  out <- df %>%
    mutate(
      landmark_date = as.Date(.data[[landmark_date_var]]),
      prior_ipaq = .data[[prior_cat]],
      follow_ipaq = .data[[follow_cat]],
      prior_met = .data[[prior_met]],
      follow_met = .data[[follow_met]],
      years_to_landmark = as.numeric(landmark_date - reg_date) / 365.25,
      age_landmark = agec + years_to_landmark,
      age_landmark_5 = age_landmark / 5,
      time_landmark = as.numeric(ddate - landmark_date) / 365.25
    ) %>%
    filter(!is.na(event_allcause), !is.na(ddate), !is.na(landmark_date), !is.na(reg_date)) %>%
    filter(landmark_date >= reg_date) %>%
    filter(ddate > landmark_date) %>%
    filter(!is.na(time_landmark), time_landmark > 0) %>%
    mutate(
      event_landmark = ifelse(event_allcause == 1, 1L, 0L),
      ipaq_delta = ifelse(!is.na(prior_ipaq) & !is.na(follow_ipaq), follow_ipaq - prior_ipaq, NA_real_),
      ipaq_direction = case_when(
        is.na(ipaq_delta) ~ NA_character_,
        ipaq_delta > 0 ~ "improved",
        ipaq_delta == 0 ~ "stable",
        ipaq_delta < 0 ~ "worsened"
      ),
      ipaq_direction = factor(ipaq_direction, levels = c("stable", "improved", "worsened")),
      prior_label = factor(prior_ipaq, levels = 1:3, labels = c("low", "moderate", "high")),
      follow_label = factor(follow_ipaq, levels = 1:3, labels = c("low", "moderate", "high")),
      ipaq_transition = ifelse(!is.na(prior_label) & !is.na(follow_label), paste0(prior_label, "->", follow_label), NA_character_),
      met_change = ifelse(!is.na(prior_met) & !is.na(follow_met), follow_met - prior_met, NA_real_),
      interval = interval_id
    )
  sd_met <- sd(out$met_change, na.rm = TRUE)
  if (is.finite(sd_met) && sd_met > 0) {
    out$met_change_sd <- out$met_change / sd_met
  } else {
    out$met_change_sd <- NA_real_
  }
  out
}

log_msg("Step 3: constructing landmark interval datasets")
int_a <- build_interval_data(dat, "IPAQ_cat_bl", "IPAQ_cat_f1", "MET_bl", "MET_f1", "reg_date_f", "BL_to_F1")
int_b <- build_interval_data(dat, "IPAQ_cat_f1", "IPAQ_cat_f2", "MET_f1", "MET_f2", "reg_date_f2", "F1_to_F2")
int_c <- build_interval_data(dat, "IPAQ_cat_f2", "IPAQ_cat_f3", "MET_f2", "MET_f3", "reg_date_f3", "F2_to_F3")
interval_list <- list(BL_to_F1 = int_a, F1_to_F2 = int_b, F2_to_F3 = int_c)

interval_summary <- bind_rows(lapply(names(interval_list), function(nm) {
  x <- interval_list[[nm]]
  data.frame(
    interval = nm,
    n_landmark_valid = nrow(x),
    n_events = sum(x$event_landmark == 1, na.rm = TRUE),
    n_ipaq_change_nonmissing = sum(!is.na(x$ipaq_direction)),
    n_met_change_nonmissing = sum(!is.na(x$met_change)),
    followup_mean_years = ifelse(nrow(x) > 0, mean(x$time_landmark, na.rm = TRUE), NA_real_),
    followup_median_years = ifelse(nrow(x) > 0, median(x$time_landmark, na.rm = TRUE), NA_real_),
    stringsAsFactors = FALSE
  )
}))
write.csv(interval_summary, file.path(out_dir, "interval_eligibility_summary.csv"), row.names = FALSE)

interval_characteristics <- bind_rows(lapply(names(interval_list), function(nm) {
  x <- interval_list[[nm]]
  if (nrow(x) == 0) return(data.frame())
  x %>%
    summarise(
      interval = nm,
      n = n(),
      events = sum(event_landmark == 1, na.rm = TRUE),
      age_landmark_mean = mean(age_landmark, na.rm = TRUE),
      age_landmark_sd = sd(age_landmark, na.rm = TRUE),
      prop_male = mean(sex == 1, na.rm = TRUE),
      met_change_mean = mean(met_change, na.rm = TRUE),
      met_change_sd = sd(met_change, na.rm = TRUE),
      met_change_median = median(met_change, na.rm = TRUE)
    )
}))
write.csv(interval_characteristics, file.path(out_dir, "interval_characteristics.csv"), row.names = FALSE)

ipaq_dir_dist <- bind_rows(lapply(names(interval_list), function(nm) {
  x <- interval_list[[nm]]
  tb <- table(x$ipaq_direction, useNA = "no")
  if (length(tb) == 0) return(data.frame())
  data.frame(
    interval = nm,
    ipaq_direction = names(tb),
    n = as.integer(tb),
    pct = 100 * as.integer(tb) / sum(tb),
    stringsAsFactors = FALSE
  )
}))
write.csv(ipaq_dir_dist, file.path(out_dir, "ipaq_direction_distribution_by_interval.csv"), row.names = FALSE)

log_msg("Step 4: fitting landmark Cox models")
coef_rows <- list()
stats_rows <- list()
ph_rows <- list()
fail_rows <- list()
idx <- 1L
st_idx <- 1L
ph_idx <- 1L
fl_idx <- 1L

for (nm in names(interval_list)) {
  x <- interval_list[[nm]]
  if (nrow(x) == 0) next

  m1_ipaq <- as.formula(Surv(time_landmark, event_landmark) ~ ipaq_direction + age_landmark_5 + sex_fct + phase_fct)
  m3_ipaq <- as.formula(Surv(time_landmark, event_landmark) ~ ipaq_direction + age_landmark_5 + sex_fct + phase_fct +
                          edu_fct + incomeh_fct + marital_fct + smk_fct + drk1_fct +
                          hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct)
  m1_met <- as.formula(Surv(time_landmark, event_landmark) ~ met_change_sd + age_landmark_5 + sex_fct + phase_fct)
  m3_met <- as.formula(Surv(time_landmark, event_landmark) ~ met_change_sd + age_landmark_5 + sex_fct + phase_fct +
                         edu_fct + incomeh_fct + marital_fct + smk_fct + drk1_fct +
                         hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct)

  model_specs <- list(
    list(formula = m1_ipaq, model_id = "M1", exposure_type = "ipaq_direction"),
    list(formula = m3_ipaq, model_id = "M3", exposure_type = "ipaq_direction"),
    list(formula = m1_met, model_id = "M1", exposure_type = "met_change_sd"),
    list(formula = m3_met, model_id = "M3", exposure_type = "met_change_sd")
  )

  for (sp in model_specs) {
    res <- fit_cox_safe(sp$formula, x, sp$model_id, nm, sp$exposure_type)
    if (!is.na(res$error)) {
      fail_rows[[fl_idx]] <- data.frame(
        interval = nm,
        exposure_type = sp$exposure_type,
        model_id = sp$model_id,
        error = res$error,
        stringsAsFactors = FALSE
      )
      fl_idx <- fl_idx + 1L
      next
    }
    coef_rows[[idx]] <- res$coef
    stats_rows[[st_idx]] <- res$stats
    idx <- idx + 1L
    st_idx <- st_idx + 1L
    if (!is.null(res$ph)) {
      ph_rows[[ph_idx]] <- res$ph
      ph_idx <- ph_idx + 1L
    }
  }
}

coef_df <- bind_rows(coef_rows)
stats_df <- bind_rows(stats_rows)
ph_df <- bind_rows(ph_rows)
fail_df <- bind_rows(fail_rows)

if (nrow(fail_df) == 0) {
  fail_df <- data.frame(
    interval = character(),
    exposure_type = character(),
    model_id = character(),
    error = character(),
    stringsAsFactors = FALSE
  )
}

write.csv(coef_df, file.path(out_dir, "cox_landmark_coefficients.csv"), row.names = FALSE)
write.csv(stats_df, file.path(out_dir, "cox_landmark_model_stats.csv"), row.names = FALSE)
write.csv(ph_df, file.path(out_dir, "cox_landmark_ph_tests.csv"), row.names = FALSE)
write.csv(fail_df, file.path(out_dir, "cox_landmark_failed_models.csv"), row.names = FALSE)

main_effects <- coef_df %>%
  filter(
    (exposure_type == "ipaq_direction" & term %in% c("ipaq_directionimproved", "ipaq_directionworsened")) |
      (exposure_type == "met_change_sd" & term == "met_change_sd")
  )
write.csv(main_effects, file.path(out_dir, "cox_landmark_main_effects.csv"), row.names = FALSE)

pooled <- bind_rows(
  fixed_effect_pool(
    main_effects %>% filter(exposure_type == "ipaq_direction", model_id == "M3"),
    c("ipaq_directionimproved", "ipaq_directionworsened")
  ) %>% mutate(exposure_type = "ipaq_direction", model_id = "M3"),
  fixed_effect_pool(
    main_effects %>% filter(exposure_type == "met_change_sd", model_id == "M3"),
    c("met_change_sd")
  ) %>% mutate(exposure_type = "met_change_sd", model_id = "M3")
)
write.csv(pooled, file.path(out_dir, "cox_landmark_pooled_main_effects.csv"), row.names = FALSE)

log_msg("Step 5: writing markdown summary")
fmt_num <- function(x, d = 3) ifelse(is.na(x), "NA", format(round(x, d), nsmall = d))
fmt_pct <- function(x, d = 1) ifelse(is.na(x), "NA", paste0(format(round(x, d), nsmall = d), "%"))

summary_lines <- c(
  "# Physical Activity Change and All-Cause Mortality: Results Summary",
  "",
  "## 1) Cohort flow",
  paste0("- Main rows: ", flow$n[flow$step == "main_rows"]),
  paste0("- Matched to mortality: ", flow$n[flow$step == "matched_to_mortality"]),
  paste0("- Valid all-cause event coding: ", flow$n[flow$step == "valid_event_allcause"]),
  paste0("- Valid positive futime: ", flow$n[flow$step == "valid_futime_positive"]),
  "",
  "## 2) Wave exposure completeness",
  paste0("- Baseline IPAQ non-missing: ", wave_summary$n_ipaq_nonmissing[wave_summary$wave == "bl"]),
  paste0("- F1 IPAQ non-missing: ", wave_summary$n_ipaq_nonmissing[wave_summary$wave == "f1"]),
  paste0("- F2 IPAQ non-missing: ", wave_summary$n_ipaq_nonmissing[wave_summary$wave == "f2"]),
  paste0("- F3 IPAQ non-missing: ", wave_summary$n_ipaq_nonmissing[wave_summary$wave == "f3"]),
  paste0("- Baseline MET non-missing: ", wave_summary$n_met_nonmissing[wave_summary$wave == "bl"]),
  paste0("- F1 MET non-missing: ", wave_summary$n_met_nonmissing[wave_summary$wave == "f1"]),
  paste0("- F2 MET non-missing: ", wave_summary$n_met_nonmissing[wave_summary$wave == "f2"]),
  paste0("- F3 MET non-missing: ", wave_summary$n_met_nonmissing[wave_summary$wave == "f3"]),
  "",
  "## 3) Landmark interval samples",
  paste0("- BL->F1: N=", interval_summary$n_landmark_valid[interval_summary$interval == "BL_to_F1"],
         ", events=", interval_summary$n_events[interval_summary$interval == "BL_to_F1"]),
  paste0("- F1->F2: N=", interval_summary$n_landmark_valid[interval_summary$interval == "F1_to_F2"],
         ", events=", interval_summary$n_events[interval_summary$interval == "F1_to_F2"]),
  paste0("- F2->F3: N=", interval_summary$n_landmark_valid[interval_summary$interval == "F2_to_F3"],
         ", events=", interval_summary$n_events[interval_summary$interval == "F2_to_F3"]),
  "",
  "## 4) Primary model effects (M3)",
  "Per-interval and pooled effects are in `cox_landmark_main_effects.csv` and `cox_landmark_pooled_main_effects.csv`."
)

if (nrow(pooled) > 0) {
  summary_lines <- c(summary_lines, "", "### Pooled fixed-effect estimates (M3)")
  for (i in seq_len(nrow(pooled))) {
    rr <- pooled[i, ]
    summary_lines <- c(
      summary_lines,
      paste0(
        "- ", rr$term, ": HR=", fmt_num(rr$hr, 3),
        " (95% CI ", fmt_num(rr$ci_lower, 3), ", ", fmt_num(rr$ci_upper, 3), "), p=", fmt_num(rr$p_value, 4),
        ", intervals=", rr$n_intervals
      )
    )
  }
}

if (nrow(fail_df) > 0) {
  summary_lines <- c(summary_lines, "", "## 5) Model failures", paste0("- Failed models: ", nrow(fail_df), " (see `cox_landmark_failed_models.csv`)"))
} else {
  summary_lines <- c(summary_lines, "", "## 5) Model failures", "- No model fitting failures.")
}

writeLines(summary_lines, con = file.path(out_dir, "ipaq_change_mortality_results_summary.md"), useBytes = TRUE)

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, con = file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)

log_msg("Analysis completed. Outputs written to: ", normalizePath(out_dir))
