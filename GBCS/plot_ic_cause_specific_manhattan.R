#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

infile <- if (length(args) >= 1) args[[1]] else file.path("output", "ic_mortality", "cox_cause_specific_main_ic_effects.csv")
outfile <- if (length(args) >= 2) args[[2]] else file.path("output", "ic_mortality", "ic_continuous_cause_specific_manhattan_p005.png")

if (!file.exists(infile)) {
  stop("Input file not found: ", infile)
}

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

df <- read.csv(infile, stringsAsFactors = FALSE)

required_cols <- c("p_value", "hr")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Keep the continuous IC effect rows if the file includes additional terms.
if ("term" %in% names(df)) {
  df <- df[df$term %in% c("ic5_per1lower", "ic_continuous"), , drop = FALSE]
}

if (nrow(df) == 0) {
  stop("No continuous IC rows found in input file.")
}

if (!("cause_group" %in% names(df))) {
  if ("outcome" %in% names(df)) {
    df$cause_group <- sub("^cause_specific_", "", df$outcome)
  } else {
    stop("Could not determine cause labels (need `cause_group` or `outcome`).")
  }
}

df$p_value <- suppressWarnings(as.numeric(df$p_value))
df$hr <- suppressWarnings(as.numeric(df$hr))
df <- df[is.finite(df$p_value) & !is.na(df$cause_group), , drop = FALSE]

if (nrow(df) == 0) {
  stop("No valid rows with p-values available for plotting.")
}

df$cause_label <- gsub("_", " ", df$cause_group)
df$cause_label <- tools::toTitleCase(df$cause_label)
df$significant <- df$p_value < 0.05
df$neg_log10_p <- -log10(pmax(df$p_value, .Machine$double.xmin))

# Order by p-value (smallest first) for a Manhattan-like ranking view.
df <- df[order(df$p_value, df$cause_label), , drop = FALSE]
df$x <- seq_len(nrow(df))

sig_threshold <- -log10(0.05)
point_cols <- ifelse(df$significant, "#C62828", "#607D8B")
point_bg <- ifelse(df$significant, "#EF9A9A", "#CFD8DC")

png(outfile, width = 1400, height = 850, res = 150)
par(mar = c(8.5, 5, 4.5, 1) + 0.1)

ylim_top <- max(c(df$neg_log10_p, sig_threshold), na.rm = TRUE) * 1.18
plot(
  df$x, df$neg_log10_p,
  type = "n",
  xaxt = "n",
  xlab = "",
  ylab = expression(-log[10](italic(P))),
  ylim = c(0, ylim_top),
  main = "Continuous IC and Cause-Specific Mortality",
  sub = "Cause-specific Cox models (IC term: ic5_per1lower); red points are P < 0.05"
)

abline(h = sig_threshold, col = "#D32F2F", lty = 2, lwd = 2)
text(
  x = par("usr")[1] + 0.02 * diff(par("usr")[1:2]),
  y = sig_threshold,
  labels = "P = 0.05",
  pos = 3,
  cex = 0.85,
  col = "#D32F2F"
)

grid(nx = NA, ny = NULL, col = "#EEEEEE", lty = 1)

points(df$x, df$neg_log10_p, pch = 21, cex = 1.9, lwd = 1.3, col = point_cols, bg = point_bg)
axis(1, at = df$x, labels = df$cause_label, las = 2, cex.axis = 0.9)

sig_idx <- which(df$significant)
if (length(sig_idx) > 0) {
  sig_labels <- sprintf("%s (p=%.3g)", df$cause_label[sig_idx], df$p_value[sig_idx])
  text(df$x[sig_idx], df$neg_log10_p[sig_idx], labels = sig_labels, pos = 3, cex = 0.78, col = "#7F0000", xpd = NA)
}

legend(
  "topright",
  legend = c("P < 0.05", "P >= 0.05"),
  pt.bg = c("#EF9A9A", "#CFD8DC"),
  col = c("#C62828", "#607D8B"),
  pch = 21,
  pt.cex = 1.4,
  bty = "n"
)

dev.off()

sig_df <- df[df$significant, c("cause_group", "hr", "p_value", "neg_log10_p"), drop = FALSE]
sig_outfile <- sub("\\.png$", "_significant_causes.csv", outfile, ignore.case = TRUE)
write.csv(sig_df, sig_outfile, row.names = FALSE)

message("Saved plot: ", outfile)
message("Saved significant causes table: ", sig_outfile)
