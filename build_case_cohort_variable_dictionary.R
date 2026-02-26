#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

cohort_dir <- "/Users/linxu/Library/CloudStorage/OneDrive-Personal/广州乳腺癌队列/乳腺癌队列/任泽舫课题组数据库/病例队列基线数据库by李岳霖2019年8月更新/病例队列库（2008年10月-2018年1月）/病例队列总库（2008年10月-2018年1月）"
questionnaire_text_dir <- "/Users/linxu/Documents/breast cancer/output/questionnaire_text"
output_dir <- "/Users/linxu/Documents/breast cancer/output"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

sav_files <- list.files(cohort_dir, pattern = "\\.sav$", full.names = TRUE)
if (length(sav_files) == 0) stop("No .sav files found in cohort directory.")

preferred_file <- sav_files[grepl("病因总库-共9029例", sav_files, fixed = TRUE)][1]
if (is.na(preferred_file)) preferred_file <- sav_files[1]

message("Using main cohort file: ", basename(preferred_file))
df <- read_sav(preferred_file)

extract_var_label <- function(x) {
  lb <- attr(x, "label")
  if (is.null(lb) || is.na(lb)) "" else as.character(lb)
}

parse_code_pairs_from_label <- function(label_text) {
  if (is.na(label_text) || label_text == "") return("")

  out <- character(0)
  # Pattern: 1=病例 / 0:对照
  p1 <- str_match_all(label_text, "(-?\\d+(?:\\.\\d+)?)\\s*[=：:]\\s*([^；;，,。()（）\\n]+)")[[1]]
  if (nrow(p1) > 0) {
    out <- c(out, paste0(trimws(p1[, 2]), "=", trimws(p1[, 3])))
  }

  # Pattern: 是…1 否…2 / 恶性肿瘤.....1
  p2 <- str_match_all(label_text, "([^；;，,。()（）\\n\\d]{1,20})\\s*[.…·]{1,8}\\s*(-?\\d+(?:\\.\\d+)?)")[[1]]
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
    paste0(vals, " (n=", counts, ")", collapse = " | ")
  } else {
    top_n <- min(10L, n_dist)
    vals <- names(tab)[seq_len(top_n)]
    counts <- as.integer(tab)[seq_len(top_n)]
    paste0(vals, " (n=", counts, ")", collapse = " | ")
  }
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

questionnaire_files <- list.files(questionnaire_text_dir, pattern = "\\.pdf\\.txt$", full.names = TRUE)
q_texts <- map(questionnaire_files, function(fp) {
  list(
    file = basename(fp),
    text = read_file(fp),
    text_upper = toupper(read_file(fp))
  )
})

match_questionnaire <- function(var_name) {
  var_u <- toupper(var_name)
  if (length(q_texts) == 0) return(list(file = "", snippet = "", match_mode = ""))

  # Priority: baseline questionnaire first.
  ord <- order(!str_detect(map_chr(q_texts, "file"), "女性健康调查表"))
  q_ord <- q_texts[ord]

  for (qi in q_ord) {
    txt_u <- qi$text_upper
    txt <- qi$text

    # 1) Exact in parentheses, e.g., (AdNo)
    m <- str_locate(txt_u, fixed(paste0("(", var_u, ")")))[1, ]
    if (!is.na(m[1])) {
      s <- max(1, m[1] - 80)
      e <- min(nchar(txt), m[2] + 120)
      return(list(file = qi$file, snippet = str_sub(txt, s, e), match_mode = "parentheses"))
    }

    # 2) Exact token-like occurrence
    m2 <- str_locate(txt_u, regex(paste0("(^|[^A-Z0-9])", var_u, "([^A-Z0-9]|$)")))[1, ]
    if (!is.na(m2[1])) {
      s <- max(1, m2[1] - 80)
      e <- min(nchar(txt), m2[2] + 120)
      return(list(file = qi$file, snippet = str_sub(txt, s, e), match_mode = "token"))
    }
  }
  list(file = "", snippet = "", match_mode = "")
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

  q <- match_questionnaire(v)
  coded <- parse_code_pairs_from_label(lb)
  unclear <- str_detect(lb, "代码|→|不祥|不详|其它|其他") ||
    str_detect(lb, "^[岁年月天分钟/\\.\\-\\sA-Za-z0-9]+$")

  tibble(
    variable_name = v,
    variable_label = lb,
    data_type = class(x)[1],
    non_missing_n = nonmiss_n,
    missing_n = miss_n,
    distinct_non_missing_n = distinct_n,
    value_range = numeric_range(x),
    coded_values_from_label = coded,
    observed_values_summary = value_summary(x),
    questionnaire_file = q$file,
    questionnaire_match_mode = q$match_mode,
    questionnaire_snippet = q$snippet,
    interpretation_flag = ifelse(unclear & q$file == "", "manual_review_recommended", "")
  )
})

prefix <- str_extract(dict$variable_name, "^[A-Za-z]+")
prefix[is.na(prefix)] <- "OTHER"
prefix_summary <- dict %>%
  mutate(section_prefix = prefix) %>%
  count(section_prefix, sort = TRUE, name = "n_variables")

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dict_csv <- file.path(output_dir, paste0("case-cohort_variable_dictionary_", ts, ".csv"))
prefix_csv <- file.path(output_dir, paste0("case-cohort_variable_prefix_summary_", ts, ".csv"))
summary_md <- file.path(output_dir, paste0("case-cohort_variable_dictionary_summary_", ts, ".md"))

write_csv(dict, dict_csv, na = "")
write_csv(prefix_summary, prefix_csv, na = "")

summary_lines <- c(
  "# Case-Cohort Variable Dictionary Summary",
  "",
  paste0("- Source directory: ", cohort_dir),
  paste0("- Main file profiled for values: ", basename(preferred_file)),
  paste0("- Number of SAV files detected: ", length(sav_files)),
  paste0("- Number of variables: ", nrow(dict)),
  paste0("- Variables with questionnaire match: ", sum(dict$questionnaire_file != "")),
  paste0("- Variables flagged for manual review: ", sum(dict$interpretation_flag != "")),
  "",
  "## Output Files",
  paste0("- Dictionary CSV: ", dict_csv),
  paste0("- Prefix summary CSV: ", prefix_csv),
  "",
  "## Top Prefixes",
  paste(apply(head(prefix_summary, 15), 1, function(r) paste0("- ", r[[1]], ": ", r[[2]])), collapse = "\n")
)

writeLines(summary_lines, summary_md)

message("Dictionary CSV: ", dict_csv)
message("Prefix summary CSV: ", prefix_csv)
message("Summary MD: ", summary_md)
