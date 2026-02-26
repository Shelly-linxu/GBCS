#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
min_events <- if (length(args) >= 1) as.integer(args[[1]]) else 10L
if (is.na(min_events) || min_events < 1L) min_events <- 10L

out_dir <- file.path("output", "ic_mortality")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Sourcing ic_mortality_analysis.R to build analysis_base and model helpers...")
source("ic_mortality_analysis.R", local = TRUE)

if (!exists("analysis_base")) stop("`analysis_base` not found after sourcing ic_mortality_analysis.R")
if (!exists("fit_and_collect")) stop("`fit_and_collect` not found after sourcing ic_mortality_analysis.R")

normalize_icd3 <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y == ""] <- NA_character_
  y <- gsub("\\s+", "", y)
  y <- sub("\\..*$", "", y)
  y <- ifelse(is.na(y), NA_character_, substr(y, 1, 3))
  y[!grepl("^[A-Z][0-9][0-9]$", y)] <- NA_character_
  y
}

icd_letter <- function(code3) substr(code3, 1, 1)
icd_num2 <- function(code3) suppressWarnings(as.integer(substr(code3, 2, 3)))

in_range <- function(code3, letter, lo, hi) {
  !is.na(code3) & icd_letter(code3) == letter & !is.na(icd_num2(code3)) & icd_num2(code3) >= lo & icd_num2(code3) <= hi
}

# Compact labels for plot annotations. Unknown codes fall back to ICD-10 code.
icd3_name_lookup <- c(
  C11 = "Nasopharyngeal cancer",
  C15 = "Esophageal cancer",
  C16 = "Gastric cancer",
  C18 = "Colon cancer",
  C20 = "Rectal cancer",
  C22 = "Liver cancer",
  C24 = "Biliary tract cancer",
  C25 = "Pancreatic cancer",
  C32 = "Laryngeal cancer",
  C34 = "Lung cancer",
  C50 = "Breast cancer",
  C56 = "Ovarian cancer",
  C61 = "Prostate cancer",
  C67 = "Bladder cancer",
  C90 = "Multiple myeloma",
  C92 = "Myeloid leukemia",
  E11 = "Type 2 diabetes mellitus",
  E46 = "Malnutrition (unspecified)",
  I10 = "Essential hypertension",
  I11 = "Hypertensive heart disease",
  I21 = "Acute myocardial infarction",
  I24 = "Other acute ischemic heart disease",
  I25 = "Chronic ischemic heart disease",
  I46 = "Cardiac arrest",
  I50 = "Heart failure",
  I60 = "Subarachnoid hemorrhage",
  I61 = "Intracerebral hemorrhage",
  I63 = "Cerebral infarction",
  I64 = "Stroke, unspecified type",
  I67 = "Other cerebrovascular disease",
  I69 = "Sequelae of cerebrovascular disease",
  I71 = "Aortic aneurysm/dissection",
  J15 = "Bacterial pneumonia",
  J18 = "Pneumonia, unspecified organism",
  J44 = "COPD",
  J47 = "Bronchiectasis",
  J98 = "Other respiratory disorders",
  K74 = "Fibrosis/cirrhosis of liver",
  K92 = "Other digestive system diseases",
  W01 = "Fall from slipping/tripping (same level)",
  W19 = "Unspecified fall",
  X80 = "Intentional self-harm by jumping",
  Y86 = "Sequelae of external causes"
)

icd3_to_name <- function(code3) {
  out <- unname(icd3_name_lookup[code3])
  out[is.na(out)] <- code3[is.na(out)]
  out
}

message("Preparing ICD-10 3-character cause codes...")
dat <- analysis_base
dat$icd3_std <- normalize_icd3(dat$und_icd_clean)

cause_counts <- dat %>%
  filter(event_allcause == 1, !is.na(icd3_std)) %>%
  count(icd3_std, sort = TRUE, name = "n_events_target")

eligible_codes <- cause_counts %>%
  filter(n_events_target >= min_events) %>%
  arrange(icd3_std) %>%
  pull(icd3_std)

message("Eligible ICD-10 3-char causes (>= ", min_events, " deaths): ", length(eligible_codes))
write.csv(cause_counts, file.path(out_dir, paste0("icd3_counts_analysis_base_min", min_events, ".csv")), row.names = FALSE)

f_cs <- as.formula(Surv(time_allcause, event_cs) ~ ic5_per1lower + agec_5 + sex_fct + phase_fct +
                     edu_fct + incomeh_fct + marital_fct +
                     smk_fct + drk1_fct + walk_any_fct + bmi +
                     hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct)

safe_fit_icd3 <- function(code3) {
  dat_cs <- dat
  dat_cs$event_cs <- ifelse(dat_cs$event_allcause == 1 & is.na(dat_cs$icd3_std), NA_integer_,
                            ifelse(dat_cs$event_allcause == 1 & dat_cs$icd3_std == code3, 1L, 0L))
  n_events_target <- sum(dat_cs$event_cs == 1, na.rm = TRUE)
  n_valid_eventcoding <- sum(!is.na(dat_cs$event_cs))
  if (n_events_target < min_events) return(NULL)

  fit <- tryCatch(
    fit_and_collect(f_cs, dat_cs, paste0("CS_ICD3_", code3, "_M3_ic5"), outcome = paste0("cause_specific_icd3_", code3)),
    error = function(e) {
      message("Skipping ", code3, " due to model error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  fit$coef$cause_code <- code3
  fit$coef$cause_level <- "icd3"
  fit$coef$n_valid_eventcoding <- n_valid_eventcoding
  fit$coef$n_events_target <- n_events_target

  fit$stats$cause_code <- code3
  fit$stats$cause_level <- "icd3"
  fit$stats$n_valid_eventcoding <- n_valid_eventcoding
  fit$stats$n_events_target <- n_events_target
  fit
}

message("Fitting cause-specific Cox models across ICD-10 3-char causes...")
fits <- vector("list", length(eligible_codes))
for (i in seq_along(eligible_codes)) {
  code3 <- eligible_codes[[i]]
  if (i %% 10 == 1 || i == length(eligible_codes)) {
    message("  Progress: ", i, "/", length(eligible_codes), " (", code3, ")")
  }
  fits[[i]] <- safe_fit_icd3(code3)
}
fits <- Filter(Negate(is.null), fits)

if (length(fits) == 0) stop("No ICD-10 3-char cause-specific models were fit.")

coef_df <- bind_rows(lapply(fits, `[[`, "coef"))
stats_df <- bind_rows(lapply(fits, `[[`, "stats"))

main_ic_df <- coef_df %>%
  filter(term == "ic5_per1lower", cause_level == "icd3") %>%
  mutate(
    p_value = as.numeric(p_value),
    hr = as.numeric(hr),
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
    significant = p_value < 0.05,
    cause_name = icd3_to_name(cause_code)
  ) %>%
  left_join(cause_counts, by = c("cause_code" = "icd3_std")) %>%
  arrange(cause_code)

coef_out <- file.path(out_dir, paste0("cox_cause_specific_icd3_coefficients_min", min_events, ".csv"))
stats_out <- file.path(out_dir, paste0("cox_cause_specific_icd3_model_stats_min", min_events, ".csv"))
main_out <- file.path(out_dir, paste0("cox_cause_specific_icd3_main_ic_effects_min", min_events, ".csv"))
sig_out <- file.path(out_dir, paste0("ic_continuous_cause_specific_manhattan_icd3_min", min_events, "_significant_causes.csv"))
png_out <- file.path(out_dir, paste0("ic_continuous_cause_specific_manhattan_icd3_min", min_events, "_p005.png"))

write.csv(coef_df, coef_out, row.names = FALSE)
write.csv(stats_df, stats_out, row.names = FALSE)
write.csv(main_ic_df, main_out, row.names = FALSE)
write.csv(main_ic_df %>% filter(significant), sig_out, row.names = FALSE)

## Selected grouped examples requested by user (IHD/COPD/DM etc.)
group_defs <- list(
  IHD = function(code3) in_range(code3, "I", 20, 25),
  Stroke = function(code3) in_range(code3, "I", 60, 69),
  Hypertensive = function(code3) in_range(code3, "I", 10, 15),
  Heart_Failure = function(code3) code3 == "I50",
  COPD = function(code3) (in_range(code3, "J", 40, 44) | code3 == "J47"),
  Pneumonia = function(code3) in_range(code3, "J", 12, 18),
  Diabetes_Mellitus = function(code3) in_range(code3, "E", 10, 14),
  CKD = function(code3) code3 == "N18",
  Dementia = function(code3) in_range(code3, "F", 1, 3) | code3 == "G30"
)

fit_group <- function(group_name, predicate) {
  dat_cs <- dat
  flag <- predicate(dat_cs$icd3_std)
  dat_cs$event_cs <- ifelse(dat_cs$event_allcause == 1 & is.na(dat_cs$icd3_std), NA_integer_,
                            ifelse(dat_cs$event_allcause == 1 & flag, 1L, 0L))
  n_events_target <- sum(dat_cs$event_cs == 1, na.rm = TRUE)
  n_valid_eventcoding <- sum(!is.na(dat_cs$event_cs))
  if (n_events_target < min_events) return(NULL)
  fit <- tryCatch(
    fit_and_collect(f_cs, dat_cs, paste0("CS_GROUP_", group_name, "_M3_ic5"), outcome = paste0("cause_specific_group_", group_name)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  out <- fit$coef %>%
    filter(term == "ic5_per1lower") %>%
    mutate(
      cause_group_name = group_name,
      n_valid_eventcoding = n_valid_eventcoding,
      n_events_target = n_events_target
    )
  out
}

selected_group_df <- bind_rows(lapply(names(group_defs), function(nm) fit_group(nm, group_defs[[nm]])))
if (nrow(selected_group_df) > 0) {
  selected_group_df <- selected_group_df %>%
    mutate(
      p_value = as.numeric(p_value),
      hr = as.numeric(hr),
      neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
      significant = p_value < 0.05
    ) %>%
    arrange(p_value)
  write.csv(selected_group_df, file.path(out_dir, paste0("cox_cause_specific_selected_groups_main_ic_effects_min", min_events, ".csv")), row.names = FALSE)
}

plot_df <- main_ic_df %>% filter(is.finite(p_value))
plot_df$x <- seq_len(nrow(plot_df))
plot_df$chapter_letter <- substr(plot_df$cause_code, 1, 1)
chapter_levels <- unique(plot_df$chapter_letter)
chapter_palette <- c("#4E79A7", "#F28E2B")
plot_df$point_col <- ifelse(plot_df$significant, "#C62828", chapter_palette[(match(plot_df$chapter_letter, chapter_levels) %% 2) + 1])
plot_df$point_bg <- ifelse(plot_df$significant, "#EF9A9A", "#DCEAF7")

sig_threshold <- -log10(0.05)
ylim_top <- max(c(plot_df$neg_log10_p, sig_threshold), na.rm = TRUE) * 1.12

png(png_out, width = max(2200, 35 * nrow(plot_df)), height = 1400, res = 170)
par(mar = c(12, 6, 5, 1) + 0.1)
plot(
  plot_df$x, plot_df$neg_log10_p,
  type = "n",
  xaxt = "n",
  xlab = "",
  ylab = expression(-log[10](italic(P))),
  ylim = c(0, ylim_top),
  main = paste0("Continuous IC vs Cause-Specific Mortality (ICD-10 3-char causes, n≥", min_events, " deaths)"),
  sub = "Cause-specific Cox M3 models for each ICD-10 3-character cause; red points are P < 0.05"
)

usr <- par("usr")
letter_blocks <- split(plot_df$x, plot_df$chapter_letter)
for (j in seq_along(letter_blocks)) {
  idx <- letter_blocks[[j]]
  if (j %% 2 == 1) {
    rect(min(idx) - 0.5, usr[3], max(idx) + 0.5, usr[4], col = rgb(0, 0, 0, 0.03), border = NA)
  }
}
grid(nx = NA, ny = NULL, col = "#EEEEEE", lty = 1)
abline(h = sig_threshold, col = "#D32F2F", lty = 2, lwd = 2)
text(usr[1] + 0.01 * diff(usr[1:2]), sig_threshold, labels = "P = 0.05", pos = 3, cex = 0.9, col = "#D32F2F")

points(plot_df$x, plot_df$neg_log10_p, pch = 21, cex = 1.25, lwd = 1.1, col = plot_df$point_col, bg = plot_df$point_bg)
axis(1, at = plot_df$x, labels = plot_df$cause_code, las = 2, cex.axis = ifelse(nrow(plot_df) > 80, 0.55, 0.7))

sig_idx <- which(plot_df$significant)
if (length(sig_idx) > 0) {
  sig_labels <- sprintf("%s (%s, p=%.2g)", plot_df$cause_name[sig_idx], plot_df$cause_code[sig_idx], plot_df$p_value[sig_idx])
  text(plot_df$x[sig_idx], plot_df$neg_log10_p[sig_idx], labels = sig_labels, pos = 3, cex = 0.65, col = "#7F0000", xpd = NA)
}

legend(
  "topright",
  legend = c("P < 0.05", "P >= 0.05"),
  pt.bg = c("#EF9A9A", "#DCEAF7"),
  col = c("#C62828", "#4E79A7"),
  pch = 21,
  pt.cex = 1.2,
  bty = "n"
)

dev.off()

message("Saved ICD-10 3-char cause-specific outputs:")
message("  ", main_out)
message("  ", sig_out)
message("  ", png_out)
if (nrow(selected_group_df) > 0) {
  message("Saved selected grouped causes (IHD/COPD/DM etc.) table:")
  message("  ", file.path(out_dir, paste0("cox_cause_specific_selected_groups_main_ic_effects_min", min_events, ".csv")))
}
