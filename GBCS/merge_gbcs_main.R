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
  args[hit[length(hit)] + 1]
}

data_dir <- get_arg("--data-dir", "/Users/linxu/Documents/GBCS")
out_file <- get_arg("--out", file.path(data_dir, "gbcs_main.dta"))
qc_dir <- get_arg("--qc-dir", file.path(data_dir, "gbcs_merge_qc"))
seed_val <- as.integer(get_arg("--seed", "20260222"))
if (is.na(seed_val)) seed_val <- 20260222L

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

normalize_id <- function(x) {
  if (inherits(x, "haven_labelled")) {
    x <- unclass(x)
  }
  if (is.numeric(x)) {
    out <- ifelse(is.na(x), NA_character_, sprintf("%.0f", x))
  } else {
    out <- as.character(x)
    out <- trimws(out)
    out[out == ""] <- NA_character_
    out <- sub("\\.0+$", "", out)
  }
  out <- trimws(out)
  out[out %in% c("", "0")] <- NA_character_
  out
}

normalize_visit_id <- function(x) {
  out <- as.character(x)
  out <- trimws(out)
  out[out == ""] <- NA_character_
  out
}

row_nonmissing_count <- function(df, exclude = character()) {
  cols <- setdiff(names(df), exclude)
  if (!length(cols)) {
    return(rep.int(0L, nrow(df)))
  }
  mats <- lapply(df[cols], function(v) {
    if (is.character(v)) {
      !(is.na(v) | trimws(v) == "")
    } else {
      !is.na(v)
    }
  })
  Reduce(`+`, lapply(mats, as.integer))
}

append_unique_suffix <- function(base, used, fallback_suffix) {
  cand <- base
  if (!(cand %in% used)) return(cand)
  cand <- paste0(base, fallback_suffix)
  if (!(cand %in% used)) return(cand)
  i <- 1L
  repeat {
    cand_i <- paste0(base, fallback_suffix, "_", i)
    if (!(cand_i %in% used)) return(cand_i)
    i <- i + 1L
  }
}

rename_followup_cols <- function(df, from_suffix, to_suffix, dataset_tag, key = "obje_id") {
  old_names <- names(df)
  stage1 <- sub(paste0(from_suffix, "$"), to_suffix, old_names)
  names(df) <- stage1

  current <- names(df)
  new_names <- character(length(current))
  used <- character()
  rename_map <- tibble(
    old_name = old_names,
    after_suffix_swap = current,
    final_name = NA_character_
  )

  for (i in seq_along(current)) {
    nm <- current[i]
    final_nm <- nm
    if (nm != key && !grepl("_(f2|f3)$", nm)) {
      desired <- paste0(nm, to_suffix)
      final_nm <- if (!(desired %in% used)) {
        desired
      } else {
        append_unique_suffix(nm, used, paste0("_", dataset_tag))
      }
    } else if (nm %in% used) {
      final_nm <- append_unique_suffix(nm, used, paste0("_", dataset_tag))
    }
    new_names[i] <- final_nm
    used <- c(used, final_nm)
    rename_map$final_name[i] <- final_nm
  }

  names(df) <- new_names
  list(data = df, rename_map = rename_map)
}

avoid_name_collisions <- function(df, existing_names, dataset_tag, key = "obje_id") {
  old_names <- names(df)
  new_names <- old_names
  used <- existing_names
  rename_map <- tibble(old_name = old_names, final_name = old_names, reason = NA_character_)

  for (i in seq_along(old_names)) {
    nm <- old_names[i]
    if (nm == key) {
      used <- c(used, nm)
      next
    }
    if (nm %in% existing_names || nm %in% used[seq_along(used) < length(used)]) {
      new_nm <- append_unique_suffix(nm, used, paste0("_", dataset_tag))
      new_names[i] <- new_nm
      rename_map$final_name[i] <- new_nm
      rename_map$reason[i] <- "collision_with_existing"
      used <- c(used, new_nm)
    } else {
      used <- c(used, nm)
    }
  }
  names(df) <- new_names
  list(data = df, rename_map = rename_map)
}

resolve_duplicates_by_nonmissing <- function(df, dataset_name, key = "obje_id", seed = 1L) {
  if (!(key %in% names(df))) stop("Missing key column: ", key)
  df$.rowid_input <- seq_len(nrow(df))
  key_vals <- df[[key]]
  dup_keys <- unique(key_vals[!is.na(key_vals) & duplicated(key_vals)])

  if (!length(dup_keys)) {
    return(list(
      data = df %>% select(-.rowid_input),
      log = tibble(),
      summary = tibble(
        dataset = dataset_name,
        duplicate_key_groups = 0L,
        duplicate_rows = 0L,
        rows_dropped = 0L,
        tie_groups = 0L
      )
    ))
  }

  set.seed(seed)
  drop_rowids <- integer()
  log_rows <- list()
  tie_groups <- 0L
  group_idx <- 0L

  for (k in dup_keys) {
    group_idx <- group_idx + 1L
    idx <- which(!is.na(df[[key]]) & df[[key]] == k)
    subdf <- df[idx, , drop = FALSE]
    nonmiss <- row_nonmissing_count(subdf, exclude = c(key))
    max_nonmiss <- max(nonmiss)
    tied <- which(nonmiss == max_nonmiss)
    keep_local <- if (length(tied) == 1L) tied else sample(tied, size = 1L)
    if (length(tied) > 1L) tie_groups <- tie_groups + 1L
    keep_global <- idx[keep_local]
    drop_global <- setdiff(idx, keep_global)
    drop_rowids <- c(drop_rowids, df$.rowid_input[drop_global])

    log_rows[[length(log_rows) + 1L]] <- tibble(
      dataset = dataset_name,
      group_index = group_idx,
      obje_id = as.character(k),
      rowid_input = subdf$.rowid_input,
      nonmissing_count = as.integer(nonmiss),
      selected_to_keep = subdf$.rowid_input == df$.rowid_input[keep_global],
      tie_nonmissing = length(tied) > 1L,
      resolution_rule = ifelse(length(tied) > 1L, "random_keep_tie", "keep_max_nonmissing")
    )
  }

  log_df <- bind_rows(log_rows)
  out_df <- df %>% filter(!(.rowid_input %in% drop_rowids)) %>% select(-.rowid_input)

  list(
    data = out_df,
    log = log_df,
    summary = tibble(
      dataset = dataset_name,
      duplicate_key_groups = length(dup_keys),
      duplicate_rows = nrow(log_df),
      rows_dropped = length(drop_rowids),
      tie_groups = tie_groups
    )
  )
}

write_csv_safe <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "")
}

cat("Reading source datasets...\n")
fu_path <- file.path(data_dir, "fu_may2014.dta")
g3_path <- file.path(data_dir, "gbcs3+端粒.dta")
g4_path <- file.path(data_dir, "gbcs4+端粒.dta")

fu <- read_dta(fu_path)
g3 <- read_dta(g3_path)
g4 <- read_dta(g4_path)

cat("Normalizing canonical merge keys (obje_id)...\n")
fu$obje_id <- normalize_id(fu$obje_id)
g3$obje_id <- normalize_id(g3$obje_id)
g4$obje_id <- normalize_id(g4$obje_id)

for (nm in intersect(c("obje_id_2f", "obje_id_3f"), names(g3))) {
  g3[[nm]] <- normalize_visit_id(g3[[nm]])
}
for (nm in intersect(c("obje_id_2f", "obje_id_3f"), names(g4))) {
  g4[[nm]] <- normalize_visit_id(g4[[nm]])
}

cat("Applying suffix renaming rules (_2f -> _f2, _3f -> _f3) and wave-tagging follow-up variables...\n")
g3_ren <- rename_followup_cols(g3, from_suffix = "_2f", to_suffix = "_f2", dataset_tag = "g3")
g3 <- g3_ren$data
g3_rename_map <- g3_ren$rename_map

g4_ren <- rename_followup_cols(g4, from_suffix = "_3f", to_suffix = "_f3", dataset_tag = "g4")
g4 <- g4_ren$data
g4_rename_map <- g4_ren$rename_map

cat("Resolving duplicates using non-missing-count rule...\n")
fu_dups <- resolve_duplicates_by_nonmissing(fu, "fu_may2014", seed = seed_val)
fu <- fu_dups$data
g3_dups <- resolve_duplicates_by_nonmissing(g3, "gbcs3", seed = seed_val)
g3 <- g3_dups$data
g4_dups <- resolve_duplicates_by_nonmissing(g4, "gbcs4", seed = seed_val)
g4 <- g4_dups$data

dup_summary <- bind_rows(fu_dups$summary, g3_dups$summary, g4_dups$summary)
dup_log <- bind_rows(fu_dups$log, g3_dups$log, g4_dups$log)

cat("Preparing follow-up flags and avoiding column collisions with anchor...\n")
g3$matched_g3 <- TRUE
g3$id_source_g3 <- ifelse(!is.na(g3$obje_id), "obje_id", NA_character_)

g4$matched_g4 <- TRUE
g4$id_source_g4 <- ifelse(!is.na(g4$obje_id), "obje_id", NA_character_)

g3_colfix <- avoid_name_collisions(g3, existing_names = names(fu), dataset_tag = "g3")
g3 <- g3_colfix$data
g3_collision_map <- g3_colfix$rename_map %>% filter(!is.na(reason))

tmp_join1 <- fu %>% left_join(g3, by = "obje_id")
g4_colfix <- avoid_name_collisions(g4, existing_names = names(tmp_join1), dataset_tag = "g4")
g4 <- g4_colfix$data
g4_collision_map <- g4_colfix$rename_map %>% filter(!is.na(reason))

cat("Generating unmatched ID QC tables...\n")
unmatched_g3 <- anti_join(g3 %>% distinct(obje_id), fu %>% distinct(obje_id), by = "obje_id")
unmatched_g4 <- anti_join(g4 %>% distinct(obje_id), fu %>% distinct(obje_id), by = "obje_id")

cat("Merging datasets (left joins on fu anchor)...\n")
merged <- fu %>%
  left_join(g3, by = "obje_id") %>%
  left_join(g4, by = "obje_id")

if ("matched_g3" %in% names(merged)) merged$matched_g3[is.na(merged$matched_g3)] <- FALSE
if ("matched_g4" %in% names(merged)) merged$matched_g4[is.na(merged$matched_g4)] <- FALSE

qc_summary <- tibble(
  metric = c(
    "fu_rows",
    "fu_unique_obje_id",
    "g3_rows_after_dedup",
    "g3_unique_obje_id_after_dedup",
    "g4_rows_after_dedup",
    "g4_unique_obje_id_after_dedup",
    "merged_rows",
    "merged_unique_obje_id",
    "matched_g3_rows_in_merged",
    "matched_g4_rows_in_merged",
    "unmatched_g3_ids_vs_fu",
    "unmatched_g4_ids_vs_fu"
  ),
  value = c(
    nrow(fu),
    dplyr::n_distinct(fu$obje_id, na.rm = TRUE),
    nrow(g3),
    dplyr::n_distinct(g3$obje_id, na.rm = TRUE),
    nrow(g4),
    dplyr::n_distinct(g4$obje_id, na.rm = TRUE),
    nrow(merged),
    dplyr::n_distinct(merged$obje_id, na.rm = TRUE),
    sum(isTRUE(merged$matched_g3) | (!is.na(merged$matched_g3) & merged$matched_g3)),
    sum(isTRUE(merged$matched_g4) | (!is.na(merged$matched_g4) & merged$matched_g4)),
    nrow(unmatched_g3),
    nrow(unmatched_g4)
  )
)

if (nrow(merged) != nrow(fu)) {
  stop("Merged row count changed unexpectedly: ", nrow(merged), " vs fu ", nrow(fu))
}
if (anyDuplicated(merged$obje_id[!is.na(merged$obje_id)])) {
  stop("Merged dataset has duplicate non-missing obje_id values.")
}

cat("Writing QC outputs...\n")
write_csv_safe(dup_summary, file.path(qc_dir, "duplicate_summary.csv"))
write_csv_safe(dup_log, file.path(qc_dir, "duplicate_resolution_log.csv"))
write_csv_safe(unmatched_g3, file.path(qc_dir, "unmatched_g3_vs_fu.csv"))
write_csv_safe(unmatched_g4, file.path(qc_dir, "unmatched_g4_vs_fu.csv"))
write_csv_safe(g3_rename_map, file.path(qc_dir, "g3_rename_map.csv"))
write_csv_safe(g4_rename_map, file.path(qc_dir, "g4_rename_map.csv"))
write_csv_safe(g3_collision_map, file.path(qc_dir, "g3_collision_renames.csv"))
write_csv_safe(g4_collision_map, file.path(qc_dir, "g4_collision_renames.csv"))
write_csv_safe(qc_summary, file.path(qc_dir, "merge_qc_summary.csv"))

cat("Writing merged dataset to:", out_file, "\n")
write_dta(merged, out_file, version = 14)

cat("Done.\n")
print(qc_summary)
