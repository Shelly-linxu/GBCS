# Physical Activity Change and All-Cause Mortality Analysis Plan (Planning Only)

## Scope

This document defines a step-by-step analysis plan to evaluate whether changes in physical activity across GBCS examination waves are associated with all-cause mortality.

Data sources:

- Mortality follow-up: `/Users/linxu/Documents/GBCS/newvital_2023.dta`
- GBCS main examination data: `/Users/linxu/Documents/GBCS/gbcs_main.dta`

No inferential analysis is performed in this step.

## Research Question

Among GBCS participants with repeated physical activity assessments, is worsening or improvement in physical activity over time (IPAQ category and MET-min/week) associated with subsequent all-cause mortality risk?

## Outcome Definition (All-Cause Mortality)

From `newvital_2023.dta`:

- Follow-up time: `futime`
- Event indicator: `dead` (`1=dead`, `0=alive`)
- Exclude unknown status: `dead==2`

Planned event coding:

- event = 1 if `dead==1`
- event = 0 if `dead==0`
- exclude missing/invalid `futime` or missing `dead`

## Exposure Definitions (Physical Activity Change)

### A. IPAQ categorical activity level across waves

Use 3-level IPAQ categories: Low, Moderate, High (coded internally as 1/2/3).

Planned wave construction aligned with existing scripts:

- Baseline: `ipaq` (0/1/2 recoded to 1/2/3)
- Follow-up 1 (F1): derive from `v7_*_f` items using IPAQ rule-based scoring
- Follow-up 2 (F2): `ipaq_f2` (1/2/3) if available, otherwise `ipaq_g3` (0/1/2 recoded)
- Follow-up 3 (F3): `ipaq_f3` (0/1/2 recoded; handle legacy trailing-space variable name if present)

Change parameterizations:

- Transition groups by adjacent waves (e.g., Low->Low, Low->Moderate/High, High->Low/Moderate, etc.)
- Directional change: improved / unchanged / worsened
- Net category change score per interval (`follow-up - prior wave`)
- Trajectory patterns across multiple waves (stable low, stable moderate/high, improving, worsening, fluctuating)

### B. MET-min/week across waves

Derive wave-specific MET-min/week from IPAQ items using:

- `8 * vigorous_minutes_week + 4 * moderate_minutes_week + 3.3 * walking_minutes_week`

Planned variables (aligned with current derivation):

- `MET_bl`, `MET_f1`, `MET_f2`, `MET_f3`

Change parameterizations:

- Absolute change per interval (`MET_t - MET_t-1`)
- Relative change per interval (`(MET_t - MET_t-1) / MET_t-1`)
- Directional classes (increase / stable / decrease using prespecified thresholds)
- Cumulative average MET up to each wave (secondary)
- Time-varying MET (continuous, with spline form if needed)

## Cohort Construction

1. Load `gbcs_main.dta` and `newvital_2023.dta`.
2. Harmonize linkage key `obje_id` to character in both files.
3. Merge 1:1 on `obje_id`.
4. Keep participants with valid mortality information (`dead` in {0,1} and non-missing `futime`).
5. Build wave-specific exposure datasets for:
   - Baseline reference analyses
   - Baseline->F1 change
   - F1->F2 change
   - F2->F3 change
6. For each change analysis, require non-missing exposure at both interval endpoints and required covariates.
7. Produce cohort-flow counts for each interval and final analytic samples.

## Primary Analysis Strategy

### 1) Landmark Cox models (primary)

Use landmark analyses to avoid immortal-time bias when defining change:

- Interval A: define change from Baseline to F1; start risk time at F1 exam date
- Interval B: define change from F1 to F2; start risk time at F2 exam date
- Interval C: define change from F2 to F3; start risk time at F3 exam date

For each interval:

- Outcome: time from landmark to death/censoring
- Model: Cox proportional hazards
- Exposures:
  - IPAQ change metrics (categorical transition and improved/unchanged/worsened)
  - MET change metrics (continuous and categorized)

Pool interval-specific estimates by meta-analytic combination (secondary) if effect sizes are directionally consistent and heterogeneity is acceptable.

### 2) Extended Cox with time-varying exposure (secondary)

Build start-stop person-period data with wave-updated physical activity and covariates:

- Interval rows per participant bounded by exam dates
- Time-varying IPAQ category and MET values
- Estimate hazard ratios for current activity level and/or recent change

This model provides a complementary estimate using all repeated measurements jointly.

## Confounder Adjustment

Core confounders (minimum set):

- `agec`
- `sex`
- `phase`

Socioeconomic/social:

- `edu`
- `incomeh`
- `marital`

Lifestyle:

- smoking (prefer `sm_stat`, fallback `sm_stat2` based on completeness)
- alcohol (`drk1` if coding/completeness acceptable)

Comorbidity history (baseline):

- hypertension (`hxht` or `hyp`)
- diabetes (`hxdm` or `diab`)
- CHD (`hxchd`)
- stroke (`hxstroke`)
- cancer history (`v12_1a`)
- COPD/chronic respiratory (`cxcopd`)

Model sequence:

1. Model 1: age/sex/phase
2. Model 2: Model 1 + SES/social
3. Model 3 (primary adjusted): Model 2 + lifestyle + comorbidities

## Functional Form and Contrast Definitions

IPAQ category change:

- Reference: stable High/Moderate group (primary) or stable High only (sensitivity)
- Report HRs for worsening and improving patterns separately

MET change:

- Continuous per SD decrease/increase
- Categorized into meaningful cutoffs (for example tertiles or clinically motivated thresholds)
- Nonlinearity check with restricted cubic splines

## Diagnostics

1. Proportional hazards testing (Schoenfeld residuals).
2. Nonlinearity checks for continuous MET and age.
3. Collinearity checks for covariate sets.
4. Influence and outlier review for extreme follow-up times and MET values.

## Missing Data Plan

Primary:

- Complete-case analysis within each model and interval.

Sensitivity:

- Multiple imputation for covariates and exposure components where appropriate.
- Compare complete-case vs imputed estimates.

## Sensitivity Analyses

1. Exclude deaths within first 1-2 years after each landmark (reverse causation check).
2. Restrict to participants with at least 2 (or 3) valid waves of activity data.
3. Alternative IPAQ change grouping (3-level directional vs detailed transitions).
4. Alternative MET scaling (log(MET+1), winsorized MET).
5. Stratified analyses by sex and baseline age group if events are sufficient.
6. Joint-model style sensitivity using baseline health status interactions (effect modification).

## Planned Outputs

- Cohort-flow tables for each interval and final analysis sets.
- Wave-specific descriptive summaries of activity category and MET distributions.
- Transition matrices and directional change summaries.
- Cox model result tables:
  - main effects (HR, 95% CI, p-value)
  - model diagnostics
  - sensitivity models
- Figures:
  - Kaplan-Meier curves by selected change groups (landmark-specific)
  - forest plot of interval-specific change effects
  - spline plot for MET change vs hazard (if nonlinear terms used)
- Reproducible analysis summary markdown in `output/` with session info.

## Implementation Notes

Suggested script workflow:

1. Data preparation script:
   - Merge mortality and GBCS data
   - Construct wave variables and interval datasets
2. Modeling script:
   - Landmark Cox models by interval
   - Time-varying Cox secondary model
3. Reporting script:
   - Tables, figures, and markdown summary

Existing scripts that can be reused:

- `ipaq_change_baseline_followup_analysis.R`
- `ipaq_mets_by_wave_analysis.R`
- `ic_mortality_analysis.R` (modeling/output structure reference)
