#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

case_control_file <- "/Users/linxu/Library/CloudStorage/OneDrive-Personal/广州乳腺癌队列/乳腺癌队列/任泽舫课题组数据库/病例队列基线数据库by李岳霖2019年8月更新/病例对照库（2008年10月-2012年3月）/总库（1551对1605）/2008年10月-2012年3月-病例对照总库（1551病例+1605对照）.sav"
questionnaire_index_candidates <- c(
  "/Users/linxu/Documents/breast cancer/output/questionnaire_code_index.csv",
  "/Users/linxu/Documents/breast cancer/repo/data/metadata/questionnaire_code_index.csv"
)
questionnaire_text_dir <- "/Users/linxu/Documents/breast cancer/output/questionnaire_text"
output_dir <- "/Users/linxu/Documents/breast cancer/output"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!file.exists(case_control_file)) {
  stop("Case-control SAV file not found: ", case_control_file)
}

questionnaire_index_file <- questionnaire_index_candidates[file.exists(questionnaire_index_candidates)][1]
if (is.na(questionnaire_index_file)) {
  questionnaire_index_file <- ""
}

message("Using case-control file: ", basename(case_control_file))
if (questionnaire_index_file != "") {
  message("Using questionnaire index: ", questionnaire_index_file)
} else {
  message("Questionnaire index not found; questionnaire-assisted mapping disabled.")
}

df <- read_sav(case_control_file)

extract_var_label <- function(x) {
  lb <- attr(x, "label")
  if (is.null(lb) || is.na(lb)) "" else as.character(lb)
}

parse_code_pairs_from_label <- function(label_text) {
  if (is.na(label_text) || label_text == "") return("")

  out <- character(0)

  # Pattern 1: "1=病例" / "0:对照"
  p1 <- str_match_all(
    label_text,
    "(-?\\d+(?:\\.\\d+)?)\\s*[=：:]\\s*([^；;，,。()（）\\n]+)"
  )[[1]]
  if (nrow(p1) > 0) {
    out <- c(out, paste0(trimws(p1[, 2]), "=", trimws(p1[, 3])))
  }

  # Pattern 2: "恶性肿瘤.....1" / "是...1"
  p2 <- str_match_all(
    label_text,
    "([^；;，,。()（）\\n\\d]{1,20})\\s*[.…·]{1,8}\\s*(-?\\d+(?:\\.\\d+)?)"
  )[[1]]
  if (nrow(p2) > 0) {
    out <- c(out, paste0(trimws(p2[, 3]), "=", trimws(p2[, 2])))
  }

  out <- out[out != ""]
  out <- unique(out)
  if (length(out) == 0) "" else paste(out, collapse = " | ")
}

as_char_value <- function(x) {
  if (inherits(x, "POSIXct") || inherits(x, "Date")) {
    return(format(x, "%Y-%m-%d"))
  }
  as.character(x)
}

value_summary <- function(x, max_show = 12L) {
  xc <- as_char_value(x)
  xc <- trimws(xc)
  miss <- is.na(x) | xc == "" | xc == "NA"
  xnm <- xc[!miss]
  if (length(xnm) == 0) return("all_missing")

  tab <- sort(table(xnm), decreasing = TRUE)
  n_dist <- length(tab)
  if (n_dist <= max_show) {
    vals <- names(tab)
    counts <- as.integer(tab)
    return(paste0(vals, " (n=", counts, ")", collapse = " | "))
  }

  top_n <- min(10L, n_dist)
  vals <- names(tab)[seq_len(top_n)]
  counts <- as.integer(tab)[seq_len(top_n)]
  paste0(vals, " (n=", counts, ")", collapse = " | ")
}

numeric_range <- function(x) {
  if (!(is.numeric(x) || inherits(x, "Date") || inherits(x, "POSIXct"))) return("")
  if (all(is.na(x))) return("")
  if (inherits(x, "Date") || inherits(x, "POSIXct")) {
    mn <- min(x, na.rm = TRUE)
    mx <- max(x, na.rm = TRUE)
    return(paste0(format(mn, "%Y-%m-%d"), " ~ ", format(mx, "%Y-%m-%d")))
  }
  mn <- suppressWarnings(min(x, na.rm = TRUE))
  mx <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(mn) || !is.finite(mx)) return("")
  paste0(mn, " ~ ", mx)
}

questionnaire_idx <- tibble(
  code_upper = character(),
  questionnaire_file = character(),
  questionnaire_page = character(),
  questionnaire_snippet = character()
)

if (questionnaire_index_file != "") {
  questionnaire_idx <- read_csv(
    questionnaire_index_file,
    show_col_types = FALSE
  ) %>%
    transmute(
      code_upper = toupper(trimws(code_upper)),
      questionnaire_file = as.character(questionnaire_file),
      questionnaire_page = as.character(page),
      questionnaire_snippet = as.character(snippet)
    ) %>%
    group_by(code_upper) %>%
    summarize(
      questionnaire_file = first(questionnaire_file),
      questionnaire_page = first(questionnaire_page),
      questionnaire_snippet = first(questionnaire_snippet),
      .groups = "drop"
  )
}

questionnaire_text_files <- character(0)
if (dir.exists(questionnaire_text_dir)) {
  questionnaire_text_files <- list.files(
    questionnaire_text_dir,
    pattern = "\\.pdf\\.txt$",
    full.names = TRUE
  )
}

q_texts <- map(questionnaire_text_files, function(fp) {
  txt <- read_file(fp)
  list(
    file = basename(fp),
    text = txt,
    text_upper = toupper(txt)
  )
})

match_questionnaire_text <- function(var_name) {
  var_u <- toupper(var_name)
  if (length(q_texts) == 0) {
    return(list(file = "", page = "", snippet = "", mode = ""))
  }

  ord <- order(!str_detect(map_chr(q_texts, "file"), "女性健康调查表"))
  q_ord <- q_texts[ord]

  for (qi in q_ord) {
    txt_u <- qi$text_upper
    txt <- qi$text

    m <- str_locate(txt_u, fixed(paste0("(", var_u, ")")))[1, ]
    if (!is.na(m[1])) {
      s <- max(1, m[1] - 80)
      e <- min(nchar(txt), m[2] + 120)
      return(list(
        file = qi$file,
        page = "",
        snippet = str_sub(txt, s, e),
        mode = "text_parentheses"
      ))
    }

    m2 <- str_locate(txt_u, regex(paste0("(^|[^A-Z0-9])", var_u, "([^A-Z0-9]|$)")))[1, ]
    if (!is.na(m2[1])) {
      s <- max(1, m2[1] - 80)
      e <- min(nchar(txt), m2[2] + 120)
      return(list(
        file = qi$file,
        page = "",
        snippet = str_sub(txt, s, e),
        mode = "text_token"
      ))
    }
  }

  list(file = "", page = "", snippet = "", mode = "")
}

vars <- names(df)
var_labels <- map_chr(df, extract_var_label)

message("Building dictionary for ", length(vars), " variables...")

dict <- map_dfr(seq_along(vars), function(i) {
  v <- vars[i]
  x <- df[[v]]
  lb <- var_labels[i]

  xc <- as_char_value(x)
  xc <- trimws(xc)
  miss <- is.na(x) | xc == "" | xc == "NA"
  nonmiss_n <- sum(!miss)
  miss_n <- sum(miss)
  distinct_n <- length(unique(xc[!miss]))

  base_var <- str_remove(v, "_(new|old)$")
  base_upper <- toupper(base_var)
  q_hit <- questionnaire_idx %>% filter(code_upper == base_upper)
  q_text_hit <- list(file = "", page = "", snippet = "", mode = "")
  if (nrow(q_hit) == 0) {
    q_text_hit <- match_questionnaire_text(base_var)
  }

  coded <- parse_code_pairs_from_label(lb)
  unclear <- str_detect(lb, "代码|→|不祥|不详|其它|其他") ||
    str_detect(lb, "^[岁年月天分钟/\\.\\-\\sA-Za-z0-9]+$")

  source_group <- case_when(
    str_detect(v, "_old$") ~ "old_only_variable",
    str_detect(v, "_new$") ~ "new_only_variable",
    TRUE ~ "merged_variable"
  )

  tibble(
    variable_name = v,
    base_variable = base_var,
    source_group = source_group,
    variable_label = lb,
    data_type = class(x)[1],
    non_missing_n = nonmiss_n,
    missing_n = miss_n,
    distinct_non_missing_n = distinct_n,
    value_range = numeric_range(x),
    coded_values_from_label = coded,
    observed_values_summary = value_summary(x),
    questionnaire_file = ifelse(
      nrow(q_hit) > 0,
      q_hit$questionnaire_file[1],
      q_text_hit$file
    ),
    questionnaire_page = ifelse(
      nrow(q_hit) > 0,
      q_hit$questionnaire_page[1],
      q_text_hit$page
    ),
    questionnaire_match_mode = ifelse(
      nrow(q_hit) > 0,
      "index_exact",
      q_text_hit$mode
    ),
    questionnaire_snippet = ifelse(
      nrow(q_hit) > 0,
      q_hit$questionnaire_snippet[1],
      q_text_hit$snippet
    ),
    interpretation_flag = ifelse(
      unclear & ifelse(nrow(q_hit) > 0, FALSE, q_text_hit$file == ""),
      "manual_review_recommended",
      ""
    )
  )
})

prefix <- str_extract(dict$variable_name, "^[A-Za-z]+")
prefix[is.na(prefix)] <- "OTHER"
prefix_summary <- dict %>%
  mutate(section_prefix = prefix) %>%
  count(section_prefix, sort = TRUE, name = "n_variables")

manual_review <- dict %>%
  filter(interpretation_flag != "") %>%
  arrange(variable_name)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dict_csv <- file.path(output_dir, paste0("case-control_variable_dictionary_", ts, ".csv"))
manual_csv <- file.path(output_dir, paste0("case-control_variable_dictionary_manual-review_", ts, ".csv"))
prefix_csv <- file.path(output_dir, paste0("case-control_variable_prefix_summary_", ts, ".csv"))
summary_md <- file.path(output_dir, paste0("case-control_variable_dictionary_summary_", ts, ".md"))

write_csv(dict, dict_csv, na = "")
write_csv(manual_review, manual_csv, na = "")
write_csv(prefix_summary, prefix_csv, na = "")

summary_lines <- c(
  "# Case-Control Variable Dictionary Summary",
  "",
  paste0("- Source SAV: ", case_control_file),
  paste0("- Number of variables: ", nrow(dict)),
  paste0("- Number of records: ", nrow(df)),
  paste0("- Variables with questionnaire code match: ", sum(dict$questionnaire_file != "")),
  paste0("-   - by index exact match: ", sum(dict$questionnaire_match_mode == "index_exact")),
  paste0("-   - by questionnaire text match: ", sum(dict$questionnaire_match_mode %in% c("text_parentheses", "text_token"))),
  paste0("- Variables with coded values parsed from labels: ", sum(dict$coded_values_from_label != "")),
  paste0("- Variables flagged for manual review: ", nrow(manual_review)),
  "",
  "## Variable Source Group",
  paste0("- Merged variables (no suffix): ", sum(dict$source_group == "merged_variable")),
  paste0("- Old-only variables (`_old`): ", sum(dict$source_group == "old_only_variable")),
  paste0("- New-only variables (`_new`): ", sum(dict$source_group == "new_only_variable")),
  "",
  "## Output Files",
  paste0("- Dictionary CSV: ", dict_csv),
  paste0("- Manual-review CSV: ", manual_csv),
  paste0("- Prefix summary CSV: ", prefix_csv),
  "",
  "## Top Prefixes",
  paste(apply(head(prefix_summary, 15), 1, function(r) paste0("- ", r[[1]], ": ", r[[2]])), collapse = "\n"),
  "",
  "## Notes",
  "- Variable suffix rules from source note:",
  "  - `_old`: variable only in old subset",
  "  - `_new`: variable only in new subset",
  "  - no suffix: merged variable"
)

writeLines(summary_lines, summary_md)

message("Dictionary CSV: ", dict_csv)
message("Manual-review CSV: ", manual_csv)
message("Prefix summary CSV: ", prefix_csv)
message("Summary MD: ", summary_md)
