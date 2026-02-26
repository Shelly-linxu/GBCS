# Intrinsic Capacity (IC) Across Waves: Variable Review and Analysis Plan

## Objective

Use the merged dataset `gbcs_main.dta` and the legacy IC scoring logic in `IC计算.R` to prepare a reproducible analysis workflow for comparing IC-related distributions across:

- Baseline (`BL`, no suffix)
- Follow-up 1 (`F1`, `_f`)
- Follow-up 2 (`F2`, `_f2`)
- Follow-up 3 (`F3`, `_f3`)

## Variable Review (from names + labels in `gbcs_main.dta`)

The IC components in `IC计算.R` can be mapped to repeated variables in `gbcs_main.dta` as follows:

- Locomotion (TUG): `v13_7d`, `v13_7e` with wave suffixes through `F3`
- Vitality (grip / BMI):
  - Grip trials repeat through `F3` (`v13_7f` to `v13_7i`, with suffixes)
  - Baseline/F2 grip trial order is consistent with the legacy code
  - F1/F3 labels show different left/right ordering, so pairing must follow labels (not only variable letters)
  - BMI is directly available only at baseline (`bmi`)
  - Follow-up BMI must be derived from wave-specific height/weight variables:
    - F1: `weight_f` + baseline `height` (no standing height field found for F1)
    - F2: `v13_9n_f2` + `v13_3_f2`
    - F3: `weight_f3` + `height_f3`
- Cognition (MMSE): `v20_1a-l` with suffixes through `F3` (legacy scoring excludes `v20_1d`)
- Psychological (GDS-15): `v20_3a-o` with suffixes through `F3`
- Sensory (vision/hearing): baseline only (`v10_19c`, `v10_20c`); no matching `*_f`, `*_f2`, `*_f3` variables found

## Key Feasibility Findings (derived with `ic_wave_distribution_analysis.R`)

Rows in merged dataset: `30,518`.

Wave-level derived availability (`/Users/linxu/Documents/breast cancer/GBCS/output/ic_wave/ic_wave_availability_summary.csv`):

| Wave | With wave date | Locomotion | Vitality | Cognition | Psy | Sensory | IC (3 repeated domains) | IC4 (incl. locomotion) | Full IC5 | IC5 + baseline sensory carry-forward |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BL | 30,518 | 10,007 | 6,981 | 9,615 | 9,991 | 30,068 | 6,672 | 6,657 | 6,632 | 6,632 |
| F1 | 18,129 | 44 | 9,187 | 17,893 | 18,007 | 0 | 9,020 | 38 | 0 | 37 |
| F2 | 10,566 | 13 | 9,875 | 10,452 | 10,480 | 0 | 9,711 | 13 | 0 | 13 |
| F3 | 5,379 | 6 | 5,321 | 5,311 | 5,353 | 0 | 5,248 | 6 | 0 | 6 |

Implications:

- Full 5-domain IC (legacy definition) is only feasible at baseline because sensory variables are not repeated.
- A strict 4-domain IC including locomotion is technically derivable at follow-ups, but sample size is too small for meaningful wave comparison (F1/F2/F3: `38/13/6`).
- The most feasible repeated composite uses the 3 domains with adequate follow-up coverage: vitality + cognition + psychological.

## Recommended Analysis Strategy

### 1. Primary longitudinal comparison (feasible)

Analyze the distribution across waves of a repeated-domain IC composite:

- `IC3_repeat = vitality + cognition + psy` (range `0-6`)
- Also report standardized percentage scale: `IC3_repeat_pct = IC3_repeat / 6 * 100`

Outputs to report by wave:

- `n`, missing count
- mean, SD
- median, IQR
- min/max
- density/violin/box plot distributions

Rationale:

- Preserves repeated domains available in all waves
- Avoids unsupported sensory carry-forward assumptions
- Avoids the severe locomotion bottleneck

### 2. Baseline-only full IC (legacy-compatible)

At baseline, reproduce the original IC scoring logic from `IC计算.R`:

- 5-domain IC total (`ic5`, range `0-10`)
- Legacy 3-category IC (`ic3_classic`, thresholds: `>=9`, `6-<9`, `<6`)

Report:

- Distribution summary of `ic5`
- Frequency table of `ic3_classic`
- Baseline component distributions (locomotion, vitality, cognition, psy, sensory)

### 3. Sensitivity analyses (explicitly labeled as non-primary)

#### A. `IC4` including locomotion (very small N at follow-ups)

- `IC4 = locomotion + vitality + cognition + psy` (range `0-8`)
- Report as feasibility/sensitivity only due follow-up Ns (`38/13/6`)

#### B. `IC5_cf_sensory` (carry-forward baseline sensory)

- `IC5_cf_sensory = IC4 + baseline_sensory`
- Use only as a sensitivity analysis
- State assumption clearly: sensory status is fixed from baseline (likely unrealistic over long follow-up)

#### C. Vitality denominator sensitivity (F1 only)

Because F1 lacks a standing height variable, BMI uses `weight_f` + baseline `height`.

Sensitivity options:

- Main: use baseline height (implemented in script)
- Alternative: use baseline BMI directly for F1 grip normalization
- Compare resulting vitality distribution and IC3_repeat distribution

### 4. Missingness and selection diagnostics

Before interpreting wave differences, quantify selection:

- Wave participation (`wave_date` present)
- Component availability by wave
- Compare baseline characteristics of:
  - all participants with wave date
  - participants contributing to `IC3_repeat`
  - participants contributing to `IC4`

This is critical because locomotion is sparse at follow-ups and may reflect a non-random substudy.

### 5. Optional inferential extension (after descriptive review)

If the descriptive distributions are acceptable and the outcome definition is finalized:

- Mixed-effects model or GEE for `IC3_repeat` / `IC3_repeat_pct` over wave
- Robust SEs and sensitivity to complete-case vs available-case definitions

Do not start modeling until the score definition (IC3_repeat vs sensitivity variants) is agreed.

## Implementation Notes

The companion script:

- `/Users/linxu/Documents/breast cancer/GBCS/ic_wave_distribution_analysis.R`

What it does:

- inspects IC-related variable names and labels in `gbcs_main.dta`
- harmonizes wave-specific variables (including F1/F3 grip left/right pairing by labels)
- derives baseline full IC and repeated-domain follow-up-compatible scores
- writes summary outputs to:
  - `/Users/linxu/Documents/breast cancer/GBCS/output/ic_wave/`

Run (default path points to `/Users/linxu/Documents/GBCS/gbcs_main.dta`):

```bash
Rscript /Users/linxu/Documents/breast\ cancer/GBCS/ic_wave_distribution_analysis.R
```

Optional environment variables:

- `GBCS_MAIN_DTA` to override dataset path
- `IC_WAVE_OUTDIR` to override output folder
- `IC_PHASE_FILTER` (e.g., `3`) to match the legacy phase restriction if needed
