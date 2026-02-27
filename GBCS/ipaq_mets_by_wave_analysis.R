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
outdir <- get_arg("--outdir", "output/ipaq_mets")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(infile)) {
  stop("Input file not found: ", infile)
}

cat("Scanning columns in: ", infile, "\n", sep = "")
hdr <- read_dta(infile, n_max = 0)
nms <- names(hdr)

need_base <- c("v7_1a", "v7_1b", "v7_2a", "v7_2b", "v7_3a", "v7_3b", "v7_4a", "v7_4b", "v7_5a", "v7_5b", "v7_6a", "v7_6b")
need_all <- unique(c(
  "obje_id",
  need_base,
  paste0(need_base, "_f"),
  paste0(need_base, "_f2"),
  paste0(need_base, "_f3")
))

missing_cols <- setdiff(need_all, nms)
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("Reading ", length(need_all), " columns...\n", sep = "")
dat <- read_dta(infile, col_select = all_of(need_all))

as_numeric_vec <- function(x) {
  if (inherits(x, "haven_labelled")) x <- unclass(x)
  suppressWarnings(as.numeric(x))
}

compute_met_wave <- function(df, suffix = "") {
  getv <- function(stem) as_numeric_vec(df[[paste0(stem, suffix)]])

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

out <- tibble(
  obje_id = as.character(dat$obje_id),
  MET_bl = compute_met_wave(dat, ""),
  MET_f1 = compute_met_wave(dat, "_f"),
  MET_f2 = compute_met_wave(dat, "_f2"),
  MET_f3 = compute_met_wave(dat, "_f3")
)

wave_names <- c("MET_bl", "MET_f1", "MET_f2", "MET_f3")
summary_tbl <- bind_rows(lapply(wave_names, function(vn) {
  x <- out[[vn]]
  x_valid <- x[!is.na(x)]
  tibble(
    wave = vn,
    n_nonmissing = length(x_valid),
    n_missing = sum(is.na(x)),
    mean = ifelse(length(x_valid) > 0, mean(x_valid), NA_real_),
    sd = ifelse(length(x_valid) > 1, sd(x_valid), NA_real_),
    median = ifelse(length(x_valid) > 0, median(x_valid), NA_real_),
    q1 = ifelse(length(x_valid) > 0, as.numeric(quantile(x_valid, 0.25, na.rm = TRUE, names = FALSE)), NA_real_),
    q3 = ifelse(length(x_valid) > 0, as.numeric(quantile(x_valid, 0.75, na.rm = TRUE, names = FALSE)), NA_real_),
    p10 = ifelse(length(x_valid) > 0, as.numeric(quantile(x_valid, 0.10, na.rm = TRUE, names = FALSE)), NA_real_),
    p90 = ifelse(length(x_valid) > 0, as.numeric(quantile(x_valid, 0.90, na.rm = TRUE, names = FALSE)), NA_real_),
    min = ifelse(length(x_valid) > 0, min(x_valid), NA_real_),
    max = ifelse(length(x_valid) > 0, max(x_valid), NA_real_)
  )
}))

write.csv(out, file.path(outdir, "ipaq_mets_per_participant_by_wave.csv"), row.names = FALSE)
write.csv(summary_tbl, file.path(outdir, "ipaq_mets_wave_summary.csv"), row.names = FALSE)

md <- c(
  "# IPAQ MET-min/week by Examination Wave",
  "",
  "Derived using IPAQ scoring formula: `8 * vigorous + 4 * moderate + 3.3 * walking` (minutes/week).",
  ""
)
for (i in seq_len(nrow(summary_tbl))) {
  s <- summary_tbl[i, ]
  md <- c(
    md,
    paste0("## ", s$wave),
    paste0("- Non-missing N: ", s$n_nonmissing),
    paste0("- Mean (SD): ", sprintf("%.1f", s$mean), " (", sprintf("%.1f", s$sd), ")"),
    paste0("- Median (Q1, Q3): ", sprintf("%.1f", s$median), " (", sprintf("%.1f", s$q1), ", ", sprintf("%.1f", s$q3), ")"),
    paste0("- P10, P90: ", sprintf("%.1f", s$p10), ", ", sprintf("%.1f", s$p90)),
    paste0("- Min, Max: ", sprintf("%.1f", s$min), ", ", sprintf("%.1f", s$max)),
    ""
  )
}
writeLines(md, file.path(outdir, "ipaq_mets_wave_summary.md"))

cat("\nWave summary:\n")
print(summary_tbl, width = Inf)
cat("\nWrote outputs to: ", normalizePath(outdir), "\n", sep = "")
