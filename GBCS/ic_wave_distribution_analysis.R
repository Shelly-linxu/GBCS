#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

data_path <- Sys.getenv("GBCS_MAIN_DTA", "/Users/linxu/Documents/GBCS/gbcs_main.dta")
out_dir <- Sys.getenv("IC_WAVE_OUTDIR", file.path("output", "ic_wave"))
phase_filter_env <- Sys.getenv("IC_PHASE_FILTER", "")

if (!file.exists(data_path)) {
  stop("Input dataset not found: ", data_path)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_label <- function(x) {
  lab <- attr(x, "label", exact = TRUE)
  if (is.null(lab)) return("")
  out <- as.character(lab)
  iconv(out, from = "", to = "UTF-8", sub = "byte")
}

get_label <- function(meta, var_name) {
  idx <- match(var_name, meta$name)
  if (is.na(idx)) "" else meta$label[idx]
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
  out <- rowSums(vals, na.rm = FALSE)
  out
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

summary_numeric <- function(x) {
  ok <- !is.na(x)
  if (!any(ok)) {
    return(data.frame(
      n = 0, missing = length(x), mean = NA_real_, sd = NA_real_,
      median = NA_real_, p25 = NA_real_, p75 = NA_real_,
      min = NA_real_, max = NA_real_
    ))
  }
  q <- quantile(x[ok], probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  data.frame(
    n = sum(ok),
    missing = sum(!ok),
    mean = mean(x[ok]),
    sd = sd(x[ok]),
    median = q[2],
    p25 = q[1],
    p75 = q[3],
    min = min(x[ok]),
    max = max(x[ok])
  )
}

mmse_suffixes <- c("a", "b", "c", "e", "f", "g", "h", "i", "j", "k", "l")
gds_keys <- letters[1:15]

gds_map <- list(
  a = c(0, 1),
  b = c(1, 0),
  c = c(1, 0),
  d = c(1, 0),
  e = c(0, 1),
  f = c(1, 0),
  g = c(0, 1),
  h = c(1, 0),
  i = c(1, 0),
  j = c(1, 0),
  k = c(0, 1),
  l = c(1, 0),
  m = c(0, 1),
  n = c(1, 0),
  o = c(1, 0)
)

wave_cfg <- list(
  BL = list(
    suffix = "",
    date_var = "examdate",
    bmi_mode = "baseline_bmi",
    grip_pairs = list(left = c("v13_7f", "v13_7g"), right = c("v13_7h", "v13_7i")),
    sensory_available = TRUE
  ),
  F1 = list(
    suffix = "_f",
    date_var = "reg_date_f",
    bmi_mode = "weight_f_plus_baseline_height",
    grip_pairs = list(left = c("v13_7f_f", "v13_7h_f"), right = c("v13_7g_f", "v13_7i_f")),
    sensory_available = FALSE
  ),
  F2 = list(
    suffix = "_f2",
    date_var = "reg_date_f2",
    bmi_mode = "f2_height_weight",
    grip_pairs = list(left = c("v13_7f_f2", "v13_7g_f2"), right = c("v13_7h_f2", "v13_7i_f2")),
    sensory_available = FALSE
  ),
  F3 = list(
    suffix = "_f3",
    date_var = "reg_date_f3",
    bmi_mode = "f3_height_weight",
    grip_pairs = list(left = c("v13_7f_f3", "v13_7h_f3"), right = c("v13_7g_f3", "v13_7i_f3")),
    sensory_available = FALSE
  )
)

meta_src <- read_dta(data_path, n_max = 0)
meta <- data.frame(
  name = names(meta_src),
  label = vapply(meta_src, safe_label, character(1)),
  stringsAsFactors = FALSE
)

build_var_map <- function(meta, wave_cfg) {
  rows <- list()
  i <- 1L
  for (wave in names(wave_cfg)) {
    cfg <- wave_cfg[[wave]]
    s <- cfg$suffix
    component_vars <- c(
      setNames(paste0(c("v13_7d", "v13_7e"), s), c("tug1", "tug2")),
      setNames(unlist(cfg$grip_pairs), c("grip_left_1", "grip_left_2", "grip_right_1", "grip_right_2")),
      setNames(paste0("v20_1", mmse_suffixes, s), paste0("mmse_", mmse_suffixes)),
      setNames(paste0("v20_3", gds_keys, s), paste0("gds_", gds_keys))
    )
    if (wave == "BL") {
      component_vars <- c(
        component_vars,
        vision = "v10_19c",
        hearing = "v10_20c",
        bmi = "bmi",
        height = "height",
        weight = "weight"
      )
    }
    if (wave == "F1") {
      component_vars <- c(component_vars, weight = "weight_f", height_proxy = "height")
    }
    if (wave == "F2") {
      component_vars <- c(component_vars, weight = "v13_9n_f2", height = "v13_3_f2")
    }
    if (wave == "F3") {
      component_vars <- c(component_vars, weight = "weight_f3", height = "height_f3")
    }
    component_vars <- c(component_vars, exam_date = cfg$date_var)
    for (role in names(component_vars)) {
      var_name <- component_vars[[role]]
      rows[[i]] <- data.frame(
        wave = wave,
        role = role,
        variable = var_name,
        label = get_label(meta, var_name),
        present = var_name %in% meta$name,
        stringsAsFactors = FALSE
      )
      i <- i + 1L
    }
  }
  bind_rows(rows)
}

var_map <- build_var_map(meta, wave_cfg)
write.csv(var_map, file.path(out_dir, "ic_variable_map.csv"), row.names = FALSE)

vars_to_read <- unique(c(
  "id", "phase", "examdate", "reg_date", "reg_date_f", "reg_date_f2", "reg_date_f3",
  "height", "weight", "bmi", "weight_f", "v13_3_f2", "v13_9n_f2", "height_f3", "weight_f3",
  "v10_19c", "v10_20c",
  paste0(rep(c("v13_7d", "v13_7e", "v13_7f", "v13_7g", "v13_7h", "v13_7i"), each = 4), c("", "_f", "_f2", "_f3")),
  as.vector(outer(paste0("v20_1", c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l")), c("", "_f", "_f2", "_f3"), paste0)),
  as.vector(outer(paste0("v20_3", gds_keys), c("", "_f", "_f2", "_f3"), paste0))
))
vars_to_read <- vars_to_read[vars_to_read %in% meta$name]

message("Reading selected columns from gbcs_main.dta: ", length(vars_to_read))
raw_df <- read_dta(data_path, col_select = all_of(vars_to_read))

if (nzchar(phase_filter_env)) {
  phase_keep <- suppressWarnings(as.numeric(phase_filter_env))
  if (!is.na(phase_keep) && "phase" %in% names(raw_df)) {
    raw_df <- raw_df %>% filter(.data$phase == phase_keep)
    message("Applied phase filter: phase == ", phase_keep)
  }
}

baseline_vision <- score_sensory_part(raw_df$v10_19c)
baseline_hearing <- score_sensory_part(raw_df$v10_20c)
baseline_sensory <- baseline_vision + baseline_hearing

derive_bmi_for_wave <- function(df, wave) {
  if (wave == "BL") {
    bmi_from_hw <- calc_bmi(df$weight, df$height)
    out <- df$bmi
    out[is.na(out)] <- bmi_from_hw[is.na(out)]
    return(out)
  }
  if (wave == "F1") return(calc_bmi(df$weight_f, df$height))
  if (wave == "F2") return(calc_bmi(df$v13_9n_f2, df$v13_3_f2))
  if (wave == "F3") return(calc_bmi(df$weight_f3, df$height_f3))
  stop("Unknown wave: ", wave)
}

derive_wave <- function(df, wave, cfg, baseline_sensory) {
  s <- cfg$suffix

  tug1 <- df[[paste0("v13_7d", s)]]
  tug2 <- df[[paste0("v13_7e", s)]]
  tug1[tug1 > 50] <- NA
  tug2[tug2 > 50] <- NA
  gugt <- pair_mean_strict(tug1, tug2)
  locomotion <- score_locomotion(gugt)

  grip_left_1 <- df[[cfg$grip_pairs$left[1]]]
  grip_left_2 <- df[[cfg$grip_pairs$left[2]]]
  grip_right_1 <- df[[cfg$grip_pairs$right[1]]]
  grip_right_2 <- df[[cfg$grip_pairs$right[2]]]
  for (obj_name in c("grip_left_1", "grip_left_2", "grip_right_1", "grip_right_2")) {
    obj <- get(obj_name)
    obj[obj > 60] <- NA
    assign(obj_name, obj)
  }
  gripl <- pair_mean_strict(grip_left_1, grip_left_2)
  gripr <- pair_mean_strict(grip_right_1, grip_right_2)
  gripmax <- pmax(gripl, gripr, na.rm = TRUE)
  gripmax[is.infinite(gripmax)] <- NA

  bmi_wave <- derive_bmi_for_wave(df, wave)
  rgripmax <- gripmax / bmi_wave
  vitality <- score_vitality(rgripmax)

  mmse_vars <- paste0("v20_1", mmse_suffixes, s)
  mmse_inputs <- as.data.frame(df[, mmse_vars]) - 1
  mmse <- row_sum_strict(mmse_inputs)
  cognition <- score_cognition(mmse)

  gds_inputs <- vector("list", length(gds_keys))
  names(gds_inputs) <- gds_keys
  for (k in gds_keys) {
    raw_name <- paste0("v20_3", k, s)
    rec <- gds_map[[k]]
    gds_inputs[[k]] <- recode_binary(df[[raw_name]], one_val = rec[1], two_val = rec[2])
  }
  gds_df <- as.data.frame(gds_inputs)
  gds <- row_sum_strict(gds_df)
  psy <- score_psy(gds)

  if (cfg$sensory_available) {
    vision <- score_sensory_part(df$v10_19c)
    hearing <- score_sensory_part(df$v10_20c)
    sensory <- vision + hearing
  } else {
    vision <- rep(NA_real_, nrow(df))
    hearing <- rep(NA_real_, nrow(df))
    sensory <- rep(NA_real_, nrow(df))
  }

  ic3_repeat <- vitality + psy + cognition
  ic4 <- locomotion + vitality + psy + cognition
  ic5 <- ic4 + sensory
  ic5_cf_sensory <- ic4 + baseline_sensory

  ic3_repeat_pct <- (ic3_repeat / 6) * 100
  ic4_pct <- (ic4 / 8) * 100
  ic5_cf_pct <- (ic5_cf_sensory / 10) * 100

  ic3_classic <- rep(NA_real_, nrow(df))
  ic3_classic[!is.na(ic5) & ic5 >= 9] <- 0
  ic3_classic[!is.na(ic5) & ic5 < 9 & ic5 >= 6] <- 1
  ic3_classic[!is.na(ic5) & ic5 < 6] <- 2

  data.frame(
    id = df$id,
    phase = if ("phase" %in% names(df)) df$phase else NA_real_,
    wave = wave,
    wave_date = as.character(df[[cfg$date_var]]),
    has_wave_date = !is.na(df[[cfg$date_var]]),
    gugt = gugt,
    locomotion = locomotion,
    gripl = gripl,
    gripr = gripr,
    gripmax = gripmax,
    bmi_wave = bmi_wave,
    rgripmax = rgripmax,
    vitality = vitality,
    mmse = mmse,
    cognition = cognition,
    gds = gds,
    psy = psy,
    vision = vision,
    hearing = hearing,
    sensory = sensory,
    baseline_sensory = baseline_sensory,
    ic3_repeat = ic3_repeat,
    ic3_repeat_pct = ic3_repeat_pct,
    ic4 = ic4,
    ic4_pct = ic4_pct,
    ic5 = ic5,
    ic5_cf_sensory = ic5_cf_sensory,
    ic5_cf_pct = ic5_cf_pct,
    ic3_classic = ic3_classic,
    stringsAsFactors = FALSE
  )
}

wave_rows <- lapply(names(wave_cfg), function(w) derive_wave(raw_df, w, wave_cfg[[w]], baseline_sensory))
ic_long <- bind_rows(wave_rows)

avail_summary <- ic_long %>%
  group_by(wave) %>%
  summarise(
    n_rows = n(),
    n_with_wave_date = sum(has_wave_date, na.rm = TRUE),
    n_locomotion = sum(!is.na(locomotion)),
    n_vitality = sum(!is.na(vitality)),
    n_cognition = sum(!is.na(cognition)),
    n_psy = sum(!is.na(psy)),
    n_sensory = sum(!is.na(sensory)),
    n_ic3_repeat = sum(!is.na(ic3_repeat)),
    n_ic4 = sum(!is.na(ic4)),
    n_ic5 = sum(!is.na(ic5)),
    n_ic5_cf_sensory = sum(!is.na(ic5_cf_sensory)),
    .groups = "drop"
  )
write.csv(avail_summary, file.path(out_dir, "ic_wave_availability_summary.csv"), row.names = FALSE)

score_vars <- c("ic3_repeat", "ic3_repeat_pct", "ic4", "ic4_pct", "ic5", "ic5_cf_sensory", "ic5_cf_pct")
dist_rows <- list()
idx <- 1L
for (w in unique(ic_long$wave)) {
  sub <- ic_long[ic_long$wave == w, , drop = FALSE]
  for (sv in score_vars) {
    stats <- summary_numeric(sub[[sv]])
    dist_rows[[idx]] <- cbind(data.frame(wave = w, score = sv, stringsAsFactors = FALSE), stats)
    idx <- idx + 1L
  }
}
dist_summary <- bind_rows(dist_rows)
write.csv(dist_summary, file.path(out_dir, "ic_wave_distribution_summary.csv"), row.names = FALSE)

component_vars <- c("locomotion", "vitality", "cognition", "psy", "sensory")
component_dist <- bind_rows(lapply(unique(ic_long$wave), function(w) {
  sub <- ic_long[ic_long$wave == w, , drop = FALSE]
  bind_rows(lapply(component_vars, function(v) {
    tbl <- as.data.frame(table(sub[[v]], useNA = "ifany"), stringsAsFactors = FALSE)
    names(tbl) <- c("value", "n")
    tbl$wave <- w
    tbl$component <- v
    tbl
  }))
}))
write.csv(component_dist, file.path(out_dir, "ic_wave_component_value_counts.csv"), row.names = FALSE)

baseline_ic3_counts <- as.data.frame(table(ic_long$ic3_classic[ic_long$wave == "BL"], useNA = "ifany"), stringsAsFactors = FALSE)
names(baseline_ic3_counts) <- c("ic3_classic", "n")
write.csv(baseline_ic3_counts, file.path(out_dir, "ic_baseline_ic3_classic_counts.csv"), row.names = FALSE)

write.csv(
  ic_long %>% select(id, phase, wave, wave_date, locomotion, vitality, cognition, psy, sensory, ic3_repeat, ic4, ic5, ic5_cf_sensory),
  file.path(out_dir, "ic_wave_derived_scores.csv"),
  row.names = FALSE
)

message("Wrote outputs to: ", normalizePath(out_dir))
message("Availability summary:")
print(avail_summary)
message("Key score summaries (n by wave):")
print(dist_summary %>% select(wave, score, n, missing) %>% arrange(score, wave))
