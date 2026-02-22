# GBCS Merge Planning

This folder stores merge planning notes and scripts for building a single cohort-wide main dataset from:

- `/Users/linxu/Documents/GBCS/fu_may2014.dta` (anchor cohort; baseline + first follow-up content)
- `/Users/linxu/Documents/GBCS/gbcs3+端粒.dta` (second repeated measurement / wave 2)
- `/Users/linxu/Documents/GBCS/gbcs4+端粒.dta` (third repeated measurement / wave 3)

Current scope:

- Merge script implemented and executed (see `GBCS/merge_gbcs_main.R`)
- Duplicate-resolution rule implemented in the merge script
- Git commits/pushes include only code and markdown documentation

Data handling guardrails:

- Do not commit `.dta` source files
- Do not commit merged outputs (`.csv`, `.rds`, `.RData`, `.sav`)
- Use `obje_id` as the canonical person-level merge key
- Treat `obje_id_2f` / `obje_id_3f` as visit-specific IDs for audit/QC

## Merge Script

- Script: `/Users/linxu/Documents/breast cancer/GBCS/merge_gbcs_main.R`
- Default output dataset: `/Users/linxu/Documents/GBCS/gbcs_main.dta`
- Default QC folder: `/Users/linxu/Documents/GBCS/gbcs_merge_qc/`

Core behavior:

- Renames `_2f -> _f2` in `gbcs3+端粒.dta`
- Renames `_3f -> _f3` in `gbcs4+端粒.dta`
- Harmonizes `obje_id` to character for matching
- Left-joins onto `fu_may2014.dta` by canonical `obje_id`
- Resolves duplicates by keeping the row with more non-missing values; tied rows are randomly resolved with a fixed seed

## Latest Run (2026-02-22)

- Generated merged dataset: `/Users/linxu/Documents/GBCS/gbcs_main.dta` (about `1.3G`)
- QC outputs written to: `/Users/linxu/Documents/GBCS/gbcs_merge_qc/`
- Summary:
  - `fu_rows = 30518`
  - `g3_rows_after_dedup = 10580`
  - `g4_rows_after_dedup = 5387` (1 duplicate row dropped after tie randomization)
  - `merged_rows = 30518`
  - `matched_g3_rows_in_merged = 10566`
  - `matched_g4_rows_in_merged = 5379`
