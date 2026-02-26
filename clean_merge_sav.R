#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(readr)
})

base_dir <- "/Users/linxu/Library/CloudStorage/OneDrive-Personal/广州乳腺癌队列/乳腺癌队列/任泽舫课题组数据库/病例队列基线数据库by李岳霖2019年8月更新/病例队列库（2008年10月-2018年1月）/病例队列总库（2008年10月-2018年1月）"
output_dir <- "/Users/linxu/Documents/breast cancer/output"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

files <- list.files(base_dir, pattern = "\\.sav$", full.names = TRUE)
if (length(files) != 4) {
  stop("Expected 4 .sav files but found: ", length(files))
}

id_candidates <- c("ADNO", "QNO", "QNO2", "SNO", "SNO2", "A1")

normalize_id <- function(x) {
  y <- as.character(x)
  y <- trimws(y)
  y[y %in% c("", "NA", "NULL")] <- NA_character_
  y
}

source_priority <- function(fname) {
  if (grepl("补充临床病理特征", fname, fixed = TRUE)) return(4L)
  if (grepl("病因总库", fname, fixed = TRUE)) return(3L)
  if (grepl("恶性病例病因库", fname, fixed = TRUE)) return(2L)
  if (grepl("良性病例病因库", fname, fixed = TRUE)) return(1L)
  0L
}

make_record_id <- function(df, source_name) {
  n <- nrow(df)
  record_id <- rep(NA_character_, n)
  record_id_type <- rep(NA_character_, n)

  for (id_name in id_candidates) {
    if (!id_name %in% names(df)) next
    v <- normalize_id(df[[id_name]])
    idx <- is.na(record_id) & !is.na(v)
    record_id[idx] <- paste0(id_name, ":", v[idx])
    record_id_type[idx] <- id_name
  }

  noid <- is.na(record_id)
  if (any(noid)) {
    record_id[noid] <- paste0("NOID:", source_name, ":", which(noid))
    record_id_type[noid] <- "NOID"
  }

  list(record_id = record_id, record_id_type = record_id_type)
}

count_nonmissing <- function(df) {
  sum_cols <- vapply(df, function(col) {
    if (inherits(col, "POSIXct") || inherits(col, "Date")) {
      !is.na(col)
    } else if (is.character(col)) {
      !is.na(col) & trimws(col) != ""
    } else {
      !is.na(col)
    }
  }, logical(nrow(df)))
  rowSums(sum_cols)
}

message("Reading 4 SAV files...")
all_data <- lapply(files, function(f) {
  df <- read_sav(f)
  src <- basename(f)
  ids <- make_record_id(df, src)
  df$source_file <- src
  df$source_priority <- source_priority(src)
  df$record_id <- ids$record_id
  df$record_id_type <- ids$record_id_type
  df$non_missing_count <- count_nonmissing(df)
  df
}) %>% bind_rows()

message("Computing duplicate summary...")
dup_stats <- all_data %>%
  count(record_id, name = "n_rows") %>%
  filter(n_rows > 1)

duplicate_ids <- nrow(dup_stats)
duplicate_rows <- if (duplicate_ids == 0) 0 else sum(dup_stats$n_rows)

message("Selecting best row per ID by non-missing values...")
deduped <- all_data %>%
  group_by(record_id) %>%
  arrange(desc(non_missing_count), desc(source_priority), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

dup_detail <- all_data %>%
  semi_join(dup_stats, by = "record_id") %>%
  group_by(record_id) %>%
  arrange(desc(non_missing_count), desc(source_priority), .by_group = TRUE) %>%
  mutate(rank_within_id = row_number()) %>%
  ungroup() %>%
  select(record_id, rank_within_id, source_file, record_id_type, non_missing_count)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
sav_out <- file.path(output_dir, paste0("case-cohort_merged_cleaned_dedup_", ts, ".sav"))
csv_out <- file.path(output_dir, paste0("case-cohort_merged_cleaned_dedup_", ts, ".csv"))
dup_out <- file.path(output_dir, paste0("case-cohort_duplicate_detail_", ts, ".csv"))
summary_out <- file.path(output_dir, paste0("case-cohort_merge_summary_", ts, ".txt"))

message("Writing outputs...")
write_sav(deduped, sav_out)
write_csv(deduped, csv_out, na = "")
write_csv(dup_detail, dup_out, na = "")

summary_lines <- c(
  paste0("Input files: ", length(files)),
  paste0("Rows before merge/dedup: ", nrow(all_data)),
  paste0("Rows after dedup: ", nrow(deduped)),
  paste0("Duplicate IDs found: ", duplicate_ids),
  paste0("Rows involved in duplicate IDs: ", duplicate_rows),
  paste0("ID priority used: ", paste(id_candidates, collapse = " -> ")),
  paste0("Output SAV: ", sav_out),
  paste0("Output CSV: ", csv_out),
  paste0("Duplicate detail CSV: ", dup_out)
)
writeLines(summary_lines, summary_out)

message("Done.")
message(paste(summary_lines, collapse = "\n"))
