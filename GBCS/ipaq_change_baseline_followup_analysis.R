#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit[length(hit)] == length(args)) return(default)
  args[hit[length(hit)] + 1L]
}

infile <- get_arg("--in", "/Users/linxu/Documents/GBCS/gbcs_main.dta")
outdir <- get_arg("--outdir", "output/ipaq_change")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(infile)) {
  stop("Input file not found: ", infile)
}

cat("Scanning column names from: ", infile, "\n", sep = "")
hdr <- read_dta(infile, n_max = 0)
nms <- names(hdr)

find_col <- function(pattern, required = TRUE) {
  hit <- grep(pattern, nms, value = TRUE)
  if (!length(hit)) {
    if (required) stop("Missing required column pattern: ", pattern)
    return(NA_character_)
  }
  hit[1]
}

nm_obje <- find_col("^obje_id$")
nm_ipaq_bl <- find_col("^ipaq$")
nm_ipaq_f2 <- find_col("^ipaq_f2$", required = FALSE)
nm_ipaq_g3 <- find_col("^ipaq_g3$", required = FALSE)
nm_ipaq_f3 <- find_col("^ipaq_f3\\s*$", required = FALSE) # legacy file may contain trailing space

v7_f_cols <- c(
  "v7_1a_f", "v7_1b_f", "v7_2a_f", "v7_2b_f",
  "v7_3a_f", "v7_3b_f", "v7_4a_f", "v7_4b_f",
  "v7_5a_f", "v7_5b_f", "v7_6a_f", "v7_6b_f"
)
missing_v7_f <- setdiff(v7_f_cols, nms)
if (length(missing_v7_f)) {
  stop("Missing F1 v7 columns: ", paste(missing_v7_f, collapse = ", "))
}

sel <- unique(na.omit(c(nm_obje, nm_ipaq_bl, v7_f_cols, nm_ipaq_f2, nm_ipaq_g3, nm_ipaq_f3)))
cat("Reading ", length(sel), " required columns...\n", sep = "")
df <- read_dta(infile, col_select = all_of(sel))

as_numeric_vec <- function(x) {
  if (inherits(x, "haven_labelled")) x <- unclass(x)
  suppressWarnings(as.numeric(x))
}

recode_ipaq_zero_based_num <- function(x) {
  x_num <- as_numeric_vec(x)
  out <- rep(NA_real_, length(x_num))
  ok <- !is.na(x_num) & x_num %in% c(0, 1, 2)
  out[ok] <- x_num[ok] + 1
  out
}

recode_ipaq_one_based_num <- function(x) {
  x_num <- as_numeric_vec(x)
  out <- rep(NA_real_, length(x_num))
  ok <- !is.na(x_num) & x_num %in% c(1, 2, 3)
  out[ok] <- x_num[ok]
  out
}

compute_ipaq_category_from_v7_num <- function(dat, suffix = "") {
  getv <- function(stem) {
    nm <- paste0(stem, suffix)
    if (!(nm %in% names(dat))) return(rep(NA_real_, nrow(dat)))
    as_numeric_vec(dat[[nm]])
  }

  v7_1a <- getv("v7_1a"); v7_1b <- getv("v7_1b"); v7_2a <- getv("v7_2a"); v7_2b <- getv("v7_2b")
  v7_3a <- getv("v7_3a"); v7_3b <- getv("v7_3b"); v7_4a <- getv("v7_4a"); v7_4b <- getv("v7_4b")
  v7_5a <- getv("v7_5a"); v7_5b <- getv("v7_5b"); v7_6a <- getv("v7_6a"); v7_6b <- getv("v7_6b")

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
    need_derive <- eligible & is.na(min_week) & !is.na(days) & !is.na(min_day)
    out[need_derive] <- days[need_derive] * min_day[need_derive]
    out[eligible & out < 0] <- NA_real_
    out
  }

  calc_minutes_day <- function(days, min_day, min_week, eligible) {
    out <- rep(0, length(days))
    out[eligible & !is.na(min_day)] <- min_day[eligible & !is.na(min_day)]
    need_derive <- eligible & is.na(min_day) & !is.na(min_week) & !is.na(days) & days > 0
    out[need_derive] <- min_week[need_derive] / days[need_derive]
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

  out <- rep(NA_real_, nrow(dat))
  out[any_info] <- 1
  out[is_moderate] <- 2
  out[is_high] <- 3
  out
}

df$IPAQ_cat_bl <- recode_ipaq_zero_based_num(df[[nm_ipaq_bl]])
df$IPAQ_cat_f1 <- compute_ipaq_category_from_v7_num(df, suffix = "_f")
df$IPAQ_cat_f2 <- if (!is.na(nm_ipaq_f2) && nm_ipaq_f2 %in% names(df)) {
  recode_ipaq_one_based_num(df[[nm_ipaq_f2]])
} else if (!is.na(nm_ipaq_g3) && nm_ipaq_g3 %in% names(df)) {
  recode_ipaq_zero_based_num(df[[nm_ipaq_g3]])
} else {
  rep(NA_real_, nrow(df))
}
df$IPAQ_cat_f3 <- if (!is.na(nm_ipaq_f3) && nm_ipaq_f3 %in% names(df)) {
  recode_ipaq_zero_based_num(df[[nm_ipaq_f3]])
} else {
  rep(NA_real_, nrow(df))
}

labs <- c("Low", "Moderate", "High")

summarize_wave <- function(dat, fu_var, wave_label) {
  tmp <- tibble(
    bl = dat$IPAQ_cat_bl,
    fu = dat[[fu_var]]
  ) %>%
    mutate(
      bl = ifelse(bl %in% 1:3, bl, NA_real_),
      fu = ifelse(fu %in% 1:3, fu, NA_real_)
    )

  paired <- tmp %>% filter(!is.na(bl) & !is.na(fu))
  b_nonmiss <- sum(!is.na(tmp$bl))
  f_nonmiss <- sum(!is.na(tmp$fu))

  trans_long <- tibble()
  paired_dist <- tibble()

  if (nrow(paired) > 0) {
    trans_tbl <- table(
      factor(paired$bl, levels = 1:3, labels = labs),
      factor(paired$fu, levels = 1:3, labels = labs)
    )

    trans_long <- as.data.frame(trans_tbl, stringsAsFactors = FALSE) %>%
      rename(baseline = Var1, followup = Var2, n = Freq) %>%
      group_by(baseline) %>%
      mutate(row_pct = 100 * n / sum(n)) %>%
      ungroup() %>%
      mutate(total_pct = 100 * n / sum(n), wave = wave_label)

    bl_tab <- table(factor(paired$bl, levels = 1:3, labels = labs))
    fu_tab <- table(factor(paired$fu, levels = 1:3, labels = labs))
    paired_dist <- bind_rows(
      tibble(wave = wave_label, timepoint = "baseline", category = names(bl_tab), n = as.integer(bl_tab)),
      tibble(wave = wave_label, timepoint = "followup", category = names(fu_tab), n = as.integer(fu_tab))
    ) %>%
      group_by(wave, timepoint) %>%
      mutate(pct = 100 * n / sum(n)) %>%
      ungroup()
  }

  delta <- paired$fu - paired$bl
  improved_n <- sum(delta > 0)
  no_change_n <- sum(delta == 0)
  worsened_n <- sum(delta < 0)

  summary_row <- tibble(
    wave = wave_label,
    fu_var = fu_var,
    n_total = nrow(tmp),
    n_baseline_nonmissing = b_nonmiss,
    n_followup_nonmissing = f_nonmiss,
    n_paired = nrow(paired),
    paired_pct_of_baseline_nonmissing = ifelse(b_nonmiss > 0, 100 * nrow(paired) / b_nonmiss, NA_real_),
    improved_n = improved_n,
    improved_pct = ifelse(nrow(paired) > 0, 100 * improved_n / nrow(paired), NA_real_),
    no_change_n = no_change_n,
    no_change_pct = ifelse(nrow(paired) > 0, 100 * no_change_n / nrow(paired), NA_real_),
    worsened_n = worsened_n,
    worsened_pct = ifelse(nrow(paired) > 0, 100 * worsened_n / nrow(paired), NA_real_),
    mean_delta = ifelse(nrow(paired) > 0, mean(delta), NA_real_),
    median_delta = ifelse(nrow(paired) > 0, median(delta), NA_real_)
  )

  list(summary = summary_row, transitions = trans_long, paired_dist = paired_dist)
}

waves <- list(
  list(var = "IPAQ_cat_f1", label = "F1 (first follow-up; derived from v7_*_f)"),
  list(var = "IPAQ_cat_f2", label = "F2"),
  list(var = "IPAQ_cat_f3", label = "F3")
)

res <- lapply(waves, function(w) summarize_wave(df, w$var, w$label))
summary_tbl <- bind_rows(lapply(res, `[[`, "summary"))
trans_tbl <- bind_rows(lapply(res, `[[`, "transitions"))
paired_dist_tbl <- bind_rows(lapply(res, `[[`, "paired_dist"))

write.csv(summary_tbl, file.path(outdir, "ipaq_change_summary_baseline_to_followups.csv"), row.names = FALSE)
write.csv(trans_tbl, file.path(outdir, "ipaq_transition_tables_baseline_to_followups.csv"), row.names = FALSE)
write.csv(paired_dist_tbl, file.path(outdir, "ipaq_paired_distributions_baseline_vs_followups.csv"), row.names = FALSE)

fmt_pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f%%", x))
md_lines <- c("# IPAQ Changes from Baseline to Follow-up", "")
for (i in seq_len(nrow(summary_tbl))) {
  s <- summary_tbl[i, ]
  md_lines <- c(
    md_lines,
    paste0("## ", s$wave),
    paste0("- Paired N: ", s$n_paired, " (", fmt_pct(s$paired_pct_of_baseline_nonmissing), " of baseline non-missing)"),
    paste0("- Improved: ", s$improved_n, " (", fmt_pct(s$improved_pct), ")"),
    paste0("- No change: ", s$no_change_n, " (", fmt_pct(s$no_change_pct), ")"),
    paste0("- Worsened: ", s$worsened_n, " (", fmt_pct(s$worsened_pct), ")"),
    paste0("- Mean delta (follow-up - baseline): ", sprintf("%.3f", s$mean_delta)),
    ""
  )
}
writeLines(md_lines, file.path(outdir, "ipaq_change_summary_baseline_to_followups.md"))

cat("\nSummary table:\n")
print(summary_tbl, width = Inf)
cat("\nWrote outputs to: ", normalizePath(outdir), "\n", sep = "")
