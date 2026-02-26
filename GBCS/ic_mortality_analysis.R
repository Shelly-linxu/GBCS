#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(survival)
  library(splines)
})

main_path <- Sys.getenv("GBCS_MAIN_DTA", "/Users/linxu/Documents/GBCS/gbcs_main.dta")
mort_path <- Sys.getenv("GBCS_MORTALITY_DTA", "/Users/linxu/Documents/GBCS/newvital_2023.dta")
out_dir <- Sys.getenv("IC_MORTALITY_OUTDIR", file.path("output", "ic_mortality"))

if (!file.exists(main_path)) stop("Main dataset not found: ", main_path)
if (!file.exists(mort_path)) stop("Mortality dataset not found: ", mort_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  flush.console()
}

write_text <- function(path, lines) {
  writeLines(lines, con = path, useBytes = TRUE)
}

normalize_id <- function(x) {
  y <- trimws(as.character(x))
  y[y == ""] <- NA_character_
  y
}

recode_binary <- function(x, one_val, two_val) {
  out <- rep(NA_real_, length(x))
  out[x == 1] <- one_val
  out[x == 2] <- two_val
  out
}

row_sum_strict <- function(df) {
  m <- as.data.frame(df)
  vals <- suppressWarnings(as.matrix(m))
  mode(vals) <- "numeric"
  rowSums(vals, na.rm = FALSE)
}

pair_mean_strict <- function(a, b) {
  out <- rep(NA_real_, length(a))
  ok <- !is.na(a) & !is.na(b)
  out[ok] <- (a[ok] + b[ok]) / 2
  out
}

calc_bmi <- function(weight_kg, height_cm) {
  out <- rep(NA_real_, length(weight_kg))
  ok <- !is.na(weight_kg) & !is.na(height_cm) & height_cm > 0
  out[ok] <- weight_kg[ok] / (height_cm[ok] / 100)^2
  out
}

score_locomotion <- function(gugt) {
  out <- rep(NA_real_, length(gugt))
  out[!is.na(gugt) & gugt < 4.7] <- 2
  out[!is.na(gugt) & gugt >= 4.7 & gugt <= 5.3] <- 1
  out[!is.na(gugt) & gugt > 5.3] <- 0
  out
}

score_vitality <- function(rgripmax) {
  out <- rep(NA_real_, length(rgripmax))
  out[!is.na(rgripmax) & rgripmax < 0.87] <- 0
  out[!is.na(rgripmax) & rgripmax >= 0.87 & rgripmax <= 1.15] <- 1
  out[!is.na(rgripmax) & rgripmax > 1.15] <- 2
  out
}

score_psy <- function(gds) {
  out <- rep(NA_real_, length(gds))
  out[!is.na(gds) & gds < 5] <- 2
  out[!is.na(gds) & gds >= 5 & gds <= 10] <- 1
  out[!is.na(gds) & gds > 10] <- 0
  out
}

score_cognition <- function(mmse) {
  out <- rep(NA_real_, length(mmse))
  out[!is.na(mmse) & mmse >= 27] <- 2
  out[!is.na(mmse) & mmse >= 21 & mmse <= 26] <- 1
  out[!is.na(mmse) & mmse <= 20] <- 0
  out
}

score_sensory_part <- function(x) {
  out <- rep(NA_real_, length(x))
  out[x %in% c(4, 5)] <- 0
  out[x == 3] <- 0.5
  out[x %in% c(1, 2)] <- 1
  out
}

mmse_suffixes <- c("a", "b", "c", "e", "f", "g", "h", "i", "j", "k", "l")
gds_keys <- letters[1:15]
gds_map <- list(
  a = c(0, 1), b = c(1, 0), c = c(1, 0), d = c(1, 0), e = c(0, 1),
  f = c(1, 0), g = c(0, 1), h = c(1, 0), i = c(1, 0), j = c(1, 0),
  k = c(0, 1), l = c(1, 0), m = c(0, 1), n = c(1, 0), o = c(1, 0)
)

derive_baseline_ic <- function(df) {
  tug1 <- df$v13_7d
  tug2 <- df$v13_7e
  tug1[tug1 > 50] <- NA
  tug2[tug2 > 50] <- NA
  gugt <- pair_mean_strict(tug1, tug2)
  locomotion <- score_locomotion(gugt)

  grip_left_1 <- df$v13_7f
  grip_left_2 <- df$v13_7g
  grip_right_1 <- df$v13_7h
  grip_right_2 <- df$v13_7i
  for (obj_name in c("grip_left_1", "grip_left_2", "grip_right_1", "grip_right_2")) {
    obj <- get(obj_name)
    obj[obj > 60] <- NA
    assign(obj_name, obj)
  }
  gripl <- pair_mean_strict(grip_left_1, grip_left_2)
  gripr <- pair_mean_strict(grip_right_1, grip_right_2)
  gripmax <- pmax(gripl, gripr, na.rm = TRUE)
  gripmax[is.infinite(gripmax)] <- NA

  bmi_wave <- df$bmi
  bmi_hw <- calc_bmi(df$weight, df$height)
  bmi_wave[is.na(bmi_wave)] <- bmi_hw[is.na(bmi_wave)]
  rgripmax <- gripmax / bmi_wave
  vitality <- score_vitality(rgripmax)

  mmse_vars <- paste0("v20_1", mmse_suffixes)
  mmse_inputs <- as.data.frame(df[, mmse_vars]) - 1
  mmse <- row_sum_strict(mmse_inputs)
  cognition <- score_cognition(mmse)

  gds_inputs <- vector("list", length(gds_keys))
  names(gds_inputs) <- gds_keys
  for (k in gds_keys) {
    rec <- gds_map[[k]]
    gds_inputs[[k]] <- recode_binary(df[[paste0("v20_3", k)]], one_val = rec[1], two_val = rec[2])
  }
  gds <- row_sum_strict(as.data.frame(gds_inputs))
  psy <- score_psy(gds)

  vision <- score_sensory_part(df$v10_19c)
  hearing <- score_sensory_part(df$v10_20c)
  sensory <- vision + hearing

  ic5 <- locomotion + vitality + cognition + psy + sensory
  ic3_classic <- rep(NA_real_, nrow(df))
  ic3_classic[!is.na(ic5) & ic5 >= 9] <- 0
  ic3_classic[!is.na(ic5) & ic5 < 9 & ic5 >= 6] <- 1
  ic3_classic[!is.na(ic5) & ic5 < 6] <- 2

  data.frame(
    obje_id = df$obje_id,
    id = if ("id" %in% names(df)) df$id else NA_character_,
    phase = if ("phase" %in% names(df)) df$phase else NA_real_,
    locomotion = locomotion,
    vitality = vitality,
    cognition = cognition,
    psy = psy,
    sensory = sensory,
    ic5 = ic5,
    ic3_classic = ic3_classic,
    gugt = gugt,
    gripmax = gripmax,
    bmi_ic = bmi_wave,
    mmse = mmse,
    gds = gds,
    stringsAsFactors = FALSE
  )
}

summary_num_by_group <- function(data, var, group) {
  x <- data[[var]]
  g <- data[[group]]
  levs <- levels(g)
  rows <- list()
  i <- 1L
  for (lv in levs) {
    xv <- x[g == lv]
    ok <- !is.na(xv)
    if (!any(ok)) {
      rows[[i]] <- data.frame(variable = var, group = lv, n = 0, missing = length(xv), mean = NA_real_, sd = NA_real_,
                              median = NA_real_, p25 = NA_real_, p75 = NA_real_, stringsAsFactors = FALSE)
    } else {
      q <- quantile(xv[ok], probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
      rows[[i]] <- data.frame(variable = var, group = lv, n = sum(ok), missing = sum(!ok), mean = mean(xv[ok]),
                              sd = sd(xv[ok]), median = q[2], p25 = q[1], p75 = q[3], stringsAsFactors = FALSE)
    }
    i <- i + 1L
  }
  bind_rows(rows)
}

summary_cat_by_group <- function(data, var, group) {
  g <- data[[group]]
  x <- data[[var]]
  levs <- levels(g)
  rows <- list()
  i <- 1L
  for (lv in levs) {
    idx <- which(g == lv)
    xv <- x[idx]
    denom <- sum(!is.na(xv))
    tab <- table(xv, useNA = "ifany")
    tab_names <- names(tab)
    for (j in seq_along(tab)) {
      nm <- tab_names[j]
      is_na_level <- identical(nm, "<NA>")
      n_val <- as.integer(tab[j])
      pct <- if (is_na_level || denom == 0) NA_real_ else 100 * n_val / denom
      rows[[i]] <- data.frame(
        variable = var,
        group = lv,
        level = nm,
        n = n_val,
        denom_non_missing = denom,
        pct = pct,
        stringsAsFactors = FALSE
      )
      i <- i + 1L
    }
  }
  bind_rows(rows)
}

extract_cox <- function(fit, model_id, outcome = "all_cause") {
  s <- summary(fit)
  coef_mat <- s$coefficients
  ci_mat <- s$conf.int
  out <- data.frame(
    outcome = outcome,
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
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out
}

fit_stats <- function(fit, model_id, outcome = "all_cause") {
  data.frame(
    outcome = outcome,
    model_id = model_id,
    n = fit$n,
    nevent = fit$nevent,
    loglik_null = fit$loglik[1],
    loglik_fit = fit$loglik[2],
    aic = AIC(fit),
    concordance = if (!is.null(summary(fit)$concordance)) summary(fit)$concordance[1] else NA_real_,
    concordance_se = if (!is.null(summary(fit)$concordance)) summary(fit)$concordance[2] else NA_real_,
    stringsAsFactors = FALSE
  )
}

ph_table <- function(fit, model_id, outcome = "all_cause") {
  z <- cox.zph(fit)
  tab <- as.data.frame(z$table, stringsAsFactors = FALSE)
  tab$term <- rownames(tab)
  names(tab) <- sub("^p$", "p_value", names(tab))
  tab$model_id <- model_id
  tab$outcome <- outcome
  tab
}

first_icd_letter <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y == ""] <- NA_character_
  y <- gsub("\\s+", "", y)
  substr(y, 1, 1)
}

icd_group_broad <- function(und_icd_clean, event_allcause) {
  letter <- first_icd_letter(und_icd_clean)
  out <- rep(NA_character_, length(und_icd_clean))
  out[event_allcause == 1 & !is.na(letter) & letter == "C"] <- "cancer"
  out[event_allcause == 1 & !is.na(letter) & letter == "I"] <- "circulatory"
  out[event_allcause == 1 & !is.na(letter) & letter == "J"] <- "respiratory"
  out[event_allcause == 1 & !is.na(letter) & !(letter %in% c("C", "I", "J"))] <- "other"
  out[event_allcause == 1 & is.na(letter)] <- "missing_icd"
  out
}

log_msg("Reading selected columns from main and mortality datasets")

main_vars <- unique(c(
  "obje_id", "id", "phase",
  "agec", "sex", "edu", "incomeh", "marital",
  "smk", "syn", "sm_stat", "sm_stat2", "drk1", "v7_5a", "bmi",
  "hyp", "hxht", "diab", "hxdm", "hxchd", "hxstroke", "v12_1a", "cxcopd",
  "weight", "height", "v10_19c", "v10_20c",
  "v13_7d", "v13_7e", "v13_7f", "v13_7g", "v13_7h", "v13_7i",
  paste0("v20_1", c("a","b","c","e","f","g","h","i","j","k","l")),
  paste0("v20_3", letters[1:15])
))

main_meta <- read_dta(main_path, n_max = 0)
main_vars <- main_vars[main_vars %in% names(main_meta)]
main_df <- read_dta(main_path, col_select = all_of(main_vars))

mort_vars <- c("obje_id", "dead", "vstatus", "futime", "und_icd", "ddate", "dyear", "reg_date", "fu_date", "phase")
mort_df <- read_dta(mort_path, col_select = all_of(mort_vars))

main_df$obje_id <- normalize_id(main_df$obje_id)
mort_df$obje_id <- normalize_id(mort_df$obje_id)

dup_main <- sum(duplicated(main_df$obje_id))
dup_mort <- sum(duplicated(mort_df$obje_id))
if (dup_main > 0 || dup_mort > 0) {
  stop("Duplicate obje_id detected (main=", dup_main, ", mortality=", dup_mort, ").")
}

log_msg("Deriving baseline IC scores using existing scoring rules")
ic_base <- derive_baseline_ic(main_df)

main_keep <- main_df %>%
  select(obje_id, agec, sex, phase, edu, incomeh, marital,
         smk, syn, sm_stat, sm_stat2, drk1, v7_5a, bmi,
         hyp, hxht, diab, hxdm, hxchd, hxstroke, v12_1a, cxcopd)

merged <- ic_base %>%
  left_join(main_keep, by = "obje_id", suffix = c("", "_main")) %>%
  left_join(mort_df, by = "obje_id", suffix = c("", "_mort"))

# Prefer phase from main if available; fill from mortality otherwise.
if ("phase_main" %in% names(merged)) {
  merged$phase <- ifelse(is.na(merged$phase_main), merged$phase, merged$phase_main)
  merged$phase_main <- NULL
}

log_msg("Creating analysis variables")
merged <- merged %>%
  mutate(
    event_allcause = case_when(dead == 1 ~ 1L, dead == 0 ~ 0L, TRUE ~ NA_integer_),
    time_allcause = futime,
    und_icd_clean = toupper(trimws(as.character(und_icd))),
    und_icd_clean = gsub("\\s+", "", und_icd_clean),
    und_icd_clean = ifelse(und_icd_clean == "", NA_character_, und_icd_clean),
    icd3 = ifelse(is.na(und_icd_clean), NA_character_, sub("\\..*$", "", und_icd_clean)),
    icd_group = icd_group_broad(und_icd_clean, event_allcause),
    ic5_per1lower = -1 * ic5,
    ic3_classic_fct = factor(ic3_classic, levels = c(0, 1, 2), labels = c("high", "middle", "low")),
    agec_5 = agec / 5,
    sex_fct = factor(ifelse(sex %in% c(0,1), sex, NA), levels = c(0,1), labels = c("female", "male")),
    phase_fct = factor(ifelse(phase %in% c(1,2,3), phase, NA), levels = c(1,2,3)),
    edu_fct = factor(ifelse(edu %in% 1:6, edu, NA), levels = c(6,5,4,3,2,1),
                     labels = c("college","junior_college","senior_middle","junior_middle","primary","lt_primary")),
    incomeh_fct = factor(ifelse(incomeh %in% 1:7, incomeh, NA), levels = c(6,5,4,3,2,1,7),
                         labels = c("ge50000","30000_band","20000_band","10000_band","5000_band","lt5000","notknow")),
    marital_clean = ifelse(marital %in% c(0,1,2,3), marital, NA),
    marital_fct = factor(marital_clean, levels = c(0,1,2,3),
                         labels = c("married","separate","widow","never_married")),
    marital5_fct = addNA(marital_fct),
    smk_fct = factor(ifelse(smk %in% c(0,1,2,3), smk, NA), levels = c(0,1,2,3),
                     labels = c("never","ex_occasional","ex_daily","current")),
    syn_fct = factor(ifelse(syn %in% c(0,1), syn, NA), levels = c(0,1), labels = c("never","ever")),
    drk1_fct = factor(ifelse(drk1 %in% 0:5, drk1, NA), levels = c(0,1,2,3,4,5),
                      labels = c("never","lt1_month","lt1_week","1to4_week","ge5_week","ex_drinker")),
    walk_any = case_when(v7_5a %in% c(1,2) ~ 1L, v7_5a == 3 ~ 0L, TRUE ~ NA_integer_),
    walk_any_fct = factor(walk_any, levels = c(0,1), labels = c("no","yes")),
    walk3_fct = factor(ifelse(v7_5a %in% c(1,2,3), v7_5a, NA), levels = c(3,2,1),
                       labels = c("no","yes_lt10min","yes_ge10min")),
    hyp_fct = factor(ifelse(hyp %in% c(0,1), hyp, NA), levels = c(0,1), labels = c("no","yes")),
    hxht_fct = factor(ifelse(hxht %in% c(0,1), hxht, NA), levels = c(0,1), labels = c("no","yes")),
    diab_fct = factor(ifelse(diab %in% c(0,1), diab, NA), levels = c(0,1), labels = c("no","yes")),
    hxdm_fct = factor(ifelse(hxdm %in% c(0,1), hxdm, NA), levels = c(0,1), labels = c("no","yes")),
    hxchd_fct = factor(ifelse(hxchd %in% c(0,1), hxchd, NA), levels = c(0,1), labels = c("no","yes")),
    hxstroke_fct = factor(ifelse(hxstroke %in% c(0,1), hxstroke, NA), levels = c(0,1), labels = c("no","yes")),
    cancerhx_bin = case_when(v12_1a == 1 ~ 0L, v12_1a == 2 ~ 1L, TRUE ~ NA_integer_),
    cancerhx_fct = factor(cancerhx_bin, levels = c(0,1), labels = c("no","yes")),
    cxcopd_any = case_when(cxcopd == 30 ~ 0L, cxcopd %in% c(31,32,33,34) ~ 1L, TRUE ~ NA_integer_),
    cxcopd_any_fct = factor(cxcopd_any, levels = c(0,1), labels = c("none","any_abnormal"))
  )
levels(merged$marital5_fct)[is.na(levels(merged$marital5_fct))] <- "missing_or_other"

flow <- data.frame(
  step = c(
    "main_rows",
    "matched_to_mortality",
    "valid_ic5",
    "valid_ic3_classic",
    "valid_event_allcause",
    "valid_time_allcause",
    "time_allcause_positive",
    "analysis_base_allcause"
  ),
  n = c(
    nrow(main_df),
    sum(!is.na(merged$dead)),
    sum(!is.na(merged$ic5)),
    sum(!is.na(merged$ic3_classic)),
    sum(!is.na(merged$event_allcause)),
    sum(!is.na(merged$time_allcause)),
    sum(!is.na(merged$time_allcause) & merged$time_allcause > 0),
    sum(!is.na(merged$ic5) & !is.na(merged$event_allcause) & !is.na(merged$time_allcause) & merged$time_allcause > 0)
  ),
  stringsAsFactors = FALSE
)
write.csv(flow, file.path(out_dir, "cohort_flow.csv"), row.names = FALSE)

analysis_base <- merged %>%
  filter(!is.na(ic5), !is.na(ic3_classic_fct), !is.na(event_allcause), !is.na(time_allcause), time_allcause > 0)

key_vars <- c(
  "ic5","ic3_classic","event_allcause","time_allcause","und_icd_clean",
  "agec","sex","phase","edu","incomeh","marital_clean","smk","syn","sm_stat","sm_stat2","drk1","v7_5a","walk_any",
  "bmi","hyp","hxht","diab","hxdm","hxchd","hxstroke","v12_1a","cxcopd"
)
missingness <- bind_rows(lapply(key_vars, function(v) {
  x <- merged[[v]]
  data.frame(
    variable = v,
    n = length(x),
    n_missing = sum(is.na(x)),
    pct_missing = 100 * mean(is.na(x)),
    n_non_missing = sum(!is.na(x)),
    stringsAsFactors = FALSE
  )
}))
write.csv(missingness, file.path(out_dir, "variable_missingness.csv"), row.names = FALSE)

lifestyle_cmp <- data.frame(
  variable = c("smk","syn","sm_stat","sm_stat2","drk1","v7_5a","walk_any","hxht","hyp","hxdm","diab"),
  n_non_missing = sapply(merged[c("smk","syn","sm_stat","sm_stat2","drk1","v7_5a","walk_any","hxht","hyp","hxdm","diab")], function(x) sum(!is.na(x))),
  n_missing = sapply(merged[c("smk","syn","sm_stat","sm_stat2","drk1","v7_5a","walk_any","hxht","hyp","hxdm","diab")], function(x) sum(is.na(x))),
  stringsAsFactors = FALSE
)
write.csv(lifestyle_cmp, file.path(out_dir, "candidate_completeness_comparison.csv"), row.names = FALSE)

simple_num_summary <- function(x, metric) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  data.frame(
    metric = metric,
    n = sum(!is.na(x)),
    missing = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    p25 = q[1],
    median = q[2],
    p75 = q[3],
    max = max(x, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
ic_dist <- bind_rows(
  simple_num_summary(analysis_base$ic5, "ic5"),
  simple_num_summary(analysis_base$ic5_per1lower, "ic5_per1lower")
)
write.csv(ic_dist, file.path(out_dir, "ic_distribution_summary.csv"), row.names = FALSE)

cont_tbl <- bind_rows(
  summary_num_by_group(analysis_base, "agec", "ic3_classic_fct"),
  summary_num_by_group(analysis_base, "bmi", "ic3_classic_fct"),
  summary_num_by_group(analysis_base, "time_allcause", "ic3_classic_fct"),
  summary_num_by_group(analysis_base, "ic5", "ic3_classic_fct")
)
write.csv(cont_tbl, file.path(out_dir, "baseline_continuous_by_ic3.csv"), row.names = FALSE)

cat_vars_for_table <- c("sex_fct","phase_fct","edu_fct","incomeh_fct","marital_fct","smk_fct","drk1_fct","walk_any_fct",
                        "hyp_fct","diab_fct","hxchd_fct","hxstroke_fct","cancerhx_fct")
cat_tbl <- bind_rows(lapply(cat_vars_for_table, function(v) summary_cat_by_group(analysis_base, v, "ic3_classic_fct")))
write.csv(cat_tbl, file.path(out_dir, "baseline_categorical_by_ic3.csv"), row.names = FALSE)

deaths_pt <- analysis_base %>%
  group_by(ic3_classic_fct) %>%
  summarise(
    n = n(),
    deaths = sum(event_allcause == 1, na.rm = TRUE),
    person_time = sum(time_allcause, na.rm = TRUE),
    death_rate_per_1000_py = 1000 * deaths / person_time,
    .groups = "drop"
  )
write.csv(deaths_pt, file.path(out_dir, "deaths_person_time_by_ic3.csv"), row.names = FALSE)

log_msg("Saving Kaplan-Meier plot")
km_fit <- survfit(Surv(time_allcause, event_allcause) ~ ic3_classic_fct, data = analysis_base)
png(file.path(out_dir, "km_allcause_by_ic3.png"), width = 1200, height = 900, res = 140)
plot(km_fit, col = c("#1b9e77", "#d95f02", "#7570b3"), lwd = 2,
     xlab = "Follow-up time (futime)", ylab = "Survival probability",
     main = "All-cause mortality by baseline IC category")
legend("bottomleft", legend = levels(analysis_base$ic3_classic_fct),
       col = c("#1b9e77", "#d95f02", "#7570b3"), lwd = 2, bty = "n", title = "IC category")
grid()
dev.off()

fit_and_collect <- function(formula, data, model_id, outcome = "all_cause") {
  fit <- coxph(formula, data = data, x = TRUE, y = TRUE, model = TRUE)
  list(
    fit = fit,
    coef = extract_cox(fit, model_id = model_id, outcome = outcome),
    stats = fit_stats(fit, model_id = model_id, outcome = outcome)
  )
}

log_msg("Fitting all-cause Cox models (stepwise)")
f_m1_ic5 <- as.formula(Surv(time_allcause, event_allcause) ~ ic5_per1lower + agec_5 + sex_fct + phase_fct)
f_m2_ic5 <- as.formula(Surv(time_allcause, event_allcause) ~ ic5_per1lower + agec_5 + sex_fct + phase_fct +
                         edu_fct + incomeh_fct + marital_fct)
f_m3_ic5 <- as.formula(Surv(time_allcause, event_allcause) ~ ic5_per1lower + agec_5 + sex_fct + phase_fct +
                         edu_fct + incomeh_fct + marital_fct +
                         smk_fct + drk1_fct + walk_any_fct + bmi +
                         hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct)
f_m3_ic3 <- as.formula(Surv(time_allcause, event_allcause) ~ ic3_classic_fct + agec_5 + sex_fct + phase_fct +
                         edu_fct + incomeh_fct + marital_fct +
                         smk_fct + drk1_fct + walk_any_fct + bmi +
                         hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct)

m1_ic5 <- fit_and_collect(f_m1_ic5, analysis_base, "M1_ic5")
m2_ic5 <- fit_and_collect(f_m2_ic5, analysis_base, "M2_ic5")
m3_ic5 <- fit_and_collect(f_m3_ic5, analysis_base, "M3_ic5")
m3_ic3 <- fit_and_collect(f_m3_ic3, analysis_base, "M3_ic3")

all_cox_coef <- bind_rows(m1_ic5$coef, m2_ic5$coef, m3_ic5$coef, m3_ic3$coef)
all_cox_stats <- bind_rows(m1_ic5$stats, m2_ic5$stats, m3_ic5$stats, m3_ic3$stats)
write.csv(all_cox_coef, file.path(out_dir, "cox_allcause_coefficients.csv"), row.names = FALSE)
write.csv(all_cox_stats, file.path(out_dir, "cox_allcause_model_stats.csv"), row.names = FALSE)

main_effect_rows <- c("ic5_per1lower", "ic3_classic_fctmiddle", "ic3_classic_fctlow")
main_effects <- all_cox_coef %>% filter(term %in% main_effect_rows)
write.csv(main_effects, file.path(out_dir, "cox_allcause_main_effects.csv"), row.names = FALSE)

log_msg("Running proportional hazards diagnostics")
ph_out <- bind_rows(
  ph_table(m3_ic5$fit, "M3_ic5"),
  ph_table(m3_ic3$fit, "M3_ic3")
)
write.csv(ph_out, file.path(out_dir, "cox_ph_tests.csv"), row.names = FALSE)

log_msg("Fitting sensitivity models")
sens_results <- list()
sens_stats <- list()
si <- 1L

analysis_lag1y <- analysis_base %>% filter(time_allcause > 1)
if (nrow(analysis_lag1y) > 0 && sum(analysis_lag1y$event_allcause) > 0) {
  fit <- fit_and_collect(f_m3_ic5, analysis_lag1y, "S1_M3_ic5_lag1y")
  sens_results[[si]] <- fit$coef; sens_stats[[si]] <- fit$stats; si <- si + 1L
}

f_m3_ic5_syn <- as.formula(Surv(time_allcause, event_allcause) ~ ic5_per1lower + agec_5 + sex_fct + phase_fct +
                             edu_fct + incomeh_fct + marital_fct +
                             syn_fct + drk1_fct + walk_any_fct + bmi +
                             hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct)
fit <- fit_and_collect(f_m3_ic5_syn, analysis_base, "S2_M3_ic5_syn")
sens_results[[si]] <- fit$coef; sens_stats[[si]] <- fit$stats; si <- si + 1L

f_m3_ic5_hx <- as.formula(Surv(time_allcause, event_allcause) ~ ic5_per1lower + agec_5 + sex_fct + phase_fct +
                            edu_fct + incomeh_fct + marital_fct +
                            smk_fct + drk1_fct + walk_any_fct + bmi +
                            hxht_fct + hxdm_fct + hxchd_fct + hxstroke_fct + cancerhx_fct)
fit <- fit_and_collect(f_m3_ic5_hx, analysis_base, "S3_M3_ic5_hxht_hxdm")
sens_results[[si]] <- fit$coef; sens_stats[[si]] <- fit$stats; si <- si + 1L

f_m3_ic5_maritalmiss <- as.formula(Surv(time_allcause, event_allcause) ~ ic5_per1lower + agec_5 + sex_fct + phase_fct +
                                     edu_fct + incomeh_fct + marital5_fct +
                                     smk_fct + drk1_fct + walk_any_fct + bmi +
                                     hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct)
fit <- fit_and_collect(f_m3_ic5_maritalmiss, analysis_base, "S4_M3_ic5_marital_missingcat")
sens_results[[si]] <- fit$coef; sens_stats[[si]] <- fit$stats; si <- si + 1L

f_m3_ic5_copd <- as.formula(Surv(time_allcause, event_allcause) ~ ic5_per1lower + agec_5 + sex_fct + phase_fct +
                              edu_fct + incomeh_fct + marital_fct +
                              smk_fct + drk1_fct + walk_any_fct + bmi +
                              hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct +
                              cxcopd_any_fct)
fit <- fit_and_collect(f_m3_ic5_copd, analysis_base, "S5_M3_ic5_plus_cxcopd")
sens_results[[si]] <- fit$coef; sens_stats[[si]] <- fit$stats; si <- si + 1L

sens_coef_df <- bind_rows(sens_results)
sens_stats_df <- bind_rows(sens_stats)
write.csv(sens_coef_df, file.path(out_dir, "cox_sensitivity_coefficients.csv"), row.names = FALSE)
write.csv(sens_stats_df, file.path(out_dir, "cox_sensitivity_model_stats.csv"), row.names = FALSE)
write.csv(sens_coef_df %>% filter(term == "ic5_per1lower"),
          file.path(out_dir, "cox_sensitivity_main_ic_effects.csv"), row.names = FALSE)

log_msg("Preparing cause-specific mortality analyses (broad ICD groups)")
icd_counts <- analysis_base %>%
  mutate(icd_group = ifelse(event_allcause == 0, "alive_or_censored", ifelse(is.na(icd_group), "missing_icd", icd_group))) %>%
  count(icd_group, sort = TRUE, name = "n")
write.csv(icd_counts, file.path(out_dir, "icd_group_counts_broad.csv"), row.names = FALSE)

cause_groups <- c("circulatory", "cancer", "respiratory", "other")
cs_coef <- list()
cs_stats <- list()
ci <- 1L
for (grp in cause_groups) {
  dat_cs <- analysis_base
  dat_cs$event_cs <- ifelse(dat_cs$event_allcause == 1 & is.na(dat_cs$icd_group), NA_integer_,
                            ifelse(dat_cs$event_allcause == 1 & dat_cs$icd_group == grp, 1L, 0L))
  n_events_grp <- sum(dat_cs$event_cs == 1, na.rm = TRUE)
  n_valid_grp <- sum(!is.na(dat_cs$event_cs))
  if (n_events_grp < 50) {
    log_msg("Skipping cause-specific model for ", grp, " due to low events (", n_events_grp, ")")
    next
  }
  f_cs <- as.formula(Surv(time_allcause, event_cs) ~ ic5_per1lower + agec_5 + sex_fct + phase_fct +
                       edu_fct + incomeh_fct + marital_fct +
                       smk_fct + drk1_fct + walk_any_fct + bmi +
                       hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct)
  fit <- fit_and_collect(f_cs, dat_cs, paste0("CS_", grp, "_M3_ic5"), outcome = paste0("cause_specific_", grp))
  fit$coef$cause_group <- grp
  fit$coef$n_valid_eventcoding <- n_valid_grp
  fit$coef$n_events_target <- n_events_grp
  fit$stats$cause_group <- grp
  fit$stats$n_valid_eventcoding <- n_valid_grp
  fit$stats$n_events_target <- n_events_grp
  cs_coef[[ci]] <- fit$coef
  cs_stats[[ci]] <- fit$stats
  ci <- ci + 1L
}
if (length(cs_coef) > 0) {
  cs_coef_df <- bind_rows(cs_coef)
  cs_stats_df <- bind_rows(cs_stats)
} else {
  cs_coef_df <- data.frame()
  cs_stats_df <- data.frame()
}
write.csv(cs_coef_df, file.path(out_dir, "cox_cause_specific_coefficients.csv"), row.names = FALSE)
write.csv(cs_stats_df, file.path(out_dir, "cox_cause_specific_model_stats.csv"), row.names = FALSE)
if (nrow(cs_coef_df) > 0) {
  write.csv(cs_coef_df %>% filter(term == "ic5_per1lower"),
            file.path(out_dir, "cox_cause_specific_main_ic_effects.csv"), row.names = FALSE)
}

log_msg("Creating results summary markdown")

get_term <- function(df, model_id, term) {
  x <- df %>% filter(model_id == model_id, term == term)
  x
}

fmt_hr <- function(df, model_id, term) {
  x <- df[df$model_id == model_id & df$term == term, , drop = FALSE]
  if (nrow(x) == 0) return("NA")
  sprintf("%.3f (%.3f, %.3f), p=%.3g", x$hr[1], x$ci_lower[1], x$ci_upper[1], x$p_value[1])
}

fmt_model_n <- function(stats_df, model_id) {
  x <- stats_df[stats_df$model_id == model_id, , drop = FALSE]
  if (nrow(x) == 0) return("NA")
  paste0("n=", x$n[1], ", events=", x$nevent[1])
}

ph_global_m3_ic5 <- ph_out %>% filter(model_id == "M3_ic5", term == "GLOBAL")
ph_global_txt <- if (nrow(ph_global_m3_ic5) == 1) sprintf("GLOBAL PH test p=%.3g", ph_global_m3_ic5$p_value[1]) else "GLOBAL PH test unavailable"

summary_lines <- c(
  "# IC-Mortality Analysis Results (Stepwise)",
  "",
  "## Data and cohort",
  paste0("- Main rows: ", nrow(main_df)),
  paste0("- Mortality rows: ", nrow(mort_df)),
  paste0("- Analytic base cohort (valid baseline IC + all-cause follow-up): ", nrow(analysis_base)),
  paste0("- All-cause deaths in analytic base cohort: ", sum(analysis_base$event_allcause == 1, na.rm = TRUE)),
  "",
  "## IC category counts",
  paste(capture.output(print(table(analysis_base$ic3_classic_fct, useNA = "ifany"))), collapse = "\n"),
  "",
  "## Primary Cox models (IC effect)",
  paste0("- M1 (continuous ic5_per1lower): ", fmt_hr(all_cox_coef, "M1_ic5", "ic5_per1lower"), " [", fmt_model_n(all_cox_stats, "M1_ic5"), "]"),
  paste0("- M2 (continuous ic5_per1lower): ", fmt_hr(all_cox_coef, "M2_ic5", "ic5_per1lower"), " [", fmt_model_n(all_cox_stats, "M2_ic5"), "]"),
  paste0("- M3 (continuous ic5_per1lower): ", fmt_hr(all_cox_coef, "M3_ic5", "ic5_per1lower"), " [", fmt_model_n(all_cox_stats, "M3_ic5"), "]"),
  paste0("- M3 (IC category middle vs high): ", fmt_hr(all_cox_coef, "M3_ic3", "ic3_classic_fctmiddle")),
  paste0("- M3 (IC category low vs high): ", fmt_hr(all_cox_coef, "M3_ic3", "ic3_classic_fctlow")),
  "",
  "## Diagnostics",
  paste0("- ", ph_global_txt),
  "",
  "## Sensitivity models (continuous ic5_per1lower)",
  if (nrow(sens_coef_df) > 0) paste(capture.output(print(sens_coef_df %>% filter(term == "ic5_per1lower") %>% select(model_id, hr, ci_lower, ci_upper, p_value, n, nevent))), collapse = "\n") else "- None",
  "",
  "## Cause-specific Cox models (broad ICD groups; continuous ic5_per1lower)",
  if (nrow(cs_coef_df) > 0) paste(capture.output(print(cs_coef_df %>% filter(term == "ic5_per1lower") %>% select(cause_group, model_id, hr, ci_lower, ci_upper, p_value, n, nevent))), collapse = "\n") else "- No cause-specific models fit (event threshold not met)",
  "",
  "## Notes",
  "- Smoking primary confounder used in fitted models: `smk` (lower missingness than `sm_stat`/`sm_stat2`).",
  "- Physical activity primary confounder used in fitted models: binary `walk_any` from `v7_5a`.",
  "- `marital` has high missingness; both primary complete-case and missing-category sensitivity models were run.",
  "- `cxcopd` is sensitivity-only due high missingness."
)
write_text(file.path(out_dir, "ic_mortality_results_summary.md"), summary_lines)

session_lines <- capture.output(sessionInfo())
write_text(file.path(out_dir, "sessionInfo.txt"), session_lines)

log_msg("Analysis complete. Outputs written to ", normalizePath(out_dir))
