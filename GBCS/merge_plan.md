# GBCS Cohort Merge Plan (Planning Only)

## Objective

Create one person-level main dataset for the GBCS cohort by merging:

- `/Users/linxu/Documents/GBCS/fu_may2014.dta` (anchor)
- `/Users/linxu/Documents/GBCS/gbcs3+端粒.dta` (wave 2 / repeated measurement)
- `/Users/linxu/Documents/GBCS/gbcs4+端粒.dta` (wave 3 / repeated measurement)

No merge is executed in this plan step.

## Confirmed Rules From User

1. Canonical merge key is `obje_id`.
2. `obje_id_2f` and `obje_id_3f` are visit-specific IDs (string IDs ending in `03` / `04`) and should be retained for audit/QC, not used as the primary merge key by default.
3. In `gbcs3+端粒.dta`, rename variables ending in `_2f` to `_f2`.
4. In `gbcs4+端粒.dta`, rename variables ending in `_3f` to `_f3`.
5. Exception: exact canonical `obje_id` stays `obje_id`.
6. Use `fu_may2014.dta` as the anchor; unmatched rows from later-wave files should not create new rows in the main dataset (left-join behavior).
7. Duplicate handling rule:
   - Check duplicates by canonical `obje_id` before merge.
   - Keep the row with more non-missing values.
   - If duplicate rows have the same non-missing count, randomly drop one and keep one.

## Dataset Profiling Summary (already checked)

- `fu_may2014.dta`
  - `30,518` rows
  - `obje_id` unique
  - anchor dataset
- `gbcs3+端粒.dta`
  - `10,580` rows
  - `obje_id` unique
- `gbcs4+端粒.dta`
  - `5,388` rows
  - `obje_id` has `1` duplicate (`1000002315601`)
  - duplicate rows differ in telomere / beta-globin assay values

## Step-by-Step Merge Plan

### Step 1. Create a merge script and QC outputs (local only)

- Implement in `R` (using `haven`, `dplyr`, `stringr`, optional `readr`).
- Script outputs are local only (not committed), for example:
  - duplicate logs
  - unmatched ID lists
  - merge summary counts

### Step 2. Read all 3 datasets and preserve raw IDs

- Read:
  - `/Users/linxu/Documents/GBCS/fu_may2014.dta`
  - `/Users/linxu/Documents/GBCS/gbcs3+端粒.dta`
  - `/Users/linxu/Documents/GBCS/gbcs4+端粒.dta`
- Preserve original ID fields for audit before any transformation:
  - `obje_id`
  - `obje_id_2f` (if present)
  - `obje_id_3f` (if present)

### Step 3. Rename wave suffixes in follow-up datasets

- `gbcs3+端粒.dta`: rename all column names ending in `_2f` to `_f2`
- `gbcs4+端粒.dta`: rename all column names ending in `_3f` to `_f3`
- Keep exact `obje_id` unchanged

Implementation note:

- Apply renaming after reading, before joining.
- Confirm whether `obje_id_2f` / `obje_id_3f` should also be renamed to `obje_id_f2` / `obje_id_f3`.
  - Current plan assumes they will follow the suffix rule unless you later choose to keep original names.

### Step 4. Harmonize merge-key data type (`obje_id`) across all datasets

This is the key requirement for successful matching.

- Convert canonical `obje_id` to a common type in all three datasets (recommended: character).
- Normalize formatting:
  - trim whitespace
  - remove trailing `.0` artifacts if created during conversion
  - avoid scientific notation
  - convert `""`, `"0"`, and numeric `0` to missing if treated as invalid IDs
- Keep raw/original ID columns for audit (`*_raw` copies optional).

Important:

- Non-key variables do not need global type harmonization to perform a merge by `obje_id`.
- String/numeric differences matter only when:
  - the same column name exists in multiple datasets and would collide, or
  - you intentionally coalesce/combine values into one column.

### Step 5. (Optional) Canonical-ID recovery using visit-specific IDs

Default matching uses canonical `obje_id` only.

If a row in `gbcs3`/`gbcs4` has missing or zero canonical `obje_id`, optionally recover a candidate canonical ID using visit-specific IDs:

- `obje_id_2f` / `obje_id_f2` ends with `03`
- `obje_id_3f` / `obje_id_f3` ends with `04`
- Derive canonical candidate ID ending with `01`

QC requirements for this rescue step:

- Record `id_source` (`obje_id` vs `derived_from_visit_id`)
- Flag rows where derived canonical ID conflicts with provided `obje_id`
- Keep unresolved rows in unmatched-ID QC tables

### Step 6. Check duplicates by canonical `obje_id` and resolve before merge

Perform duplicate checks in each dataset after key harmonization (and after optional ID recovery):

- `fu` (expected no duplicates)
- `gbcs3` (expected no duplicates)
- `gbcs4` (known duplicate exists)

Duplicate resolution rule (user-specified):

1. Group by canonical `obje_id`.
2. For each duplicate row, count non-missing values (recommended on all non-key columns; log whether key columns are excluded).
3. Keep the row with the largest non-missing count.
4. If there is a tie in non-missing counts:
   - randomly keep one row
   - drop the others
   - log the dropped rows and tie condition

Reproducibility recommendation:

- Set a fixed random seed (for example `set.seed(...)`) before tie-breaking so the random drop is reproducible and auditable.

### Step 7. Rename or isolate non-key overlapping columns to avoid collisions

Before joining, inspect overlapping names between datasets and prevent overwrites.

- Keep wave-specific variables wave-specific.
- If unsuffixed variables in `gbcs3` or `gbcs4` are wave-specific (for example telomere/lab fields), rename them to include wave identity (`_f2`, `_f3`, or clear prefixes).
- Do not overwrite columns already present in `fu`.

This step also avoids string/numeric type conflicts caused by accidental same-name collisions.

### Step 8. Perform left joins using canonical `obje_id`

Merge sequence:

1. `fu_may2014` LEFT JOIN cleaned `gbcs3` by canonical `obje_id`
2. result LEFT JOIN cleaned `gbcs4` by canonical `obje_id`

Expected behavior:

- Row count remains the same as `fu`
- `obje_id` remains unique
- Rows in `fu` without later-wave matches keep missing values for appended variables
- Later-wave rows not found in `fu` are not added to the main dataset

### Step 9. Produce QC tables and merge diagnostics (local only)

Generate and review:

- duplicate summary table (before/after resolution)
- tie-break log for random duplicate drops
- unmatched `gbcs3` IDs vs `fu`
- unmatched `gbcs4` IDs vs `fu`
- row counts and unique-ID counts before/after each join
- match-rate summary by wave

### Step 10. Validate final merged dataset structure

Validation checklist:

- row count equals `nrow(fu_may2014)`
- canonical `obje_id` unique
- `_f2` and `_f3` renamed columns present as expected
- no unexpected overwrites of anchor variables
- provenance flags added (recommended):
  - `matched_g3`
  - `matched_g4`
  - `id_source_g3`
  - `id_source_g4`

### Step 11. Save outputs locally and keep Git clean

- Save merged dataset locally only (not committed)
- Save QC outputs locally only (not committed)
- Commit and push only:
  - scripts (`.R`, `.py`) if created
  - markdown documentation (`README`, merge plan)

## Known Issues Already Identified

- `gbcs4+端粒.dta` has one duplicated canonical `obje_id` and will require the duplicate rule above.
- Some `gbcs3`/`gbcs4` IDs do not match `fu` by canonical `obje_id`; these should remain in QC unmatched lists and not be appended as new rows to the anchor dataset.

## Implementation Notes for Data Types

- Join key (`obje_id`) type must be harmonized across datasets before merge.
- Other variables can remain as native types unless combining/coalescing into the same output column.
- `haven_labelled` columns should be handled intentionally:
  - preserve labelled class if not transforming, or
  - convert explicitly to numeric/character with documentation if needed later.
