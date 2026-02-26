# Intrinsic Capacity (IC) and Mortality Analysis Plan (Planning Only)

## Scope

This document defines a step-by-step analysis plan to examine the association between intrinsic capacity (IC) and mortality using:

- Mortality follow-up file: `/Users/linxu/Documents/GBCS/newvital_2023.dta`
- IC/covariate source file: `/Users/linxu/Documents/GBCS/gbcs_main.dta`

No data analysis is performed in this step.

## Confirmed Mortality Variables (from `newvital_2023.dta`)

- `futime`: length of follow-up (planned survival time variable)
- `dead`: vital status as on 21Nov2023 (`0=Alive`, `1=Dead`, `2=Unknown`)
- `und_icd`: ICD-10 underlying cause code (available for deaths; some missing among deaths)
- `ddate`: date of death/censoring
- `obje_id`: registry identifier for linkage

## Confirmed Linkage Feasibility

- Linkage key: `obje_id`
- `newvital_2023.dta`: `obje_id` stored as numeric
- `gbcs_main.dta`: `obje_id` stored as character
- IDs are uniquely matched (1:1 after type harmonization to character/string)

## Exposure Definition (IC)

IC is not present in `newvital_2023.dta` and must be derived from `gbcs_main.dta` using the existing repo IC scoring logic (see `ic_wave_distribution_analysis.R`).

### Primary exposure (recommended)

- Baseline full IC composite (`ic5`, continuous)
- Baseline IC category (`ic3_classic`, categorical)

Rationale:

- Mortality follow-up is complete in the mortality file.
- Repeated IC at follow-up waves is sparse for some domains and may introduce selection bias.
- Baseline IC is a stable, interpretable primary exposure for prospective mortality analyses.

### Secondary exposure parameterizations

- Continuous IC (per 1-point lower/higher IC)
- Categorical IC (reference = highest IC category)
- Individual IC domain scores (domain-specific analyses; secondary)

## Outcomes

### Primary outcome: all-cause mortality

- Time: `futime`
- Event indicator:
  - `1` if `dead == 1`
  - `0` if `dead == 0`
- Exclude/handle separately:
  - `dead == 2` (Unknown)
  - missing `dead` or missing `futime`

### Secondary outcome: cause-specific mortality

- Use `und_icd` (ICD-10 underlying cause) among deaths
- Create prespecified ICD-10 cause groups (broad chapters / clinically meaningful groups), for example:
  - circulatory diseases
  - cancers
  - respiratory diseases
  - other causes
- Predefine rules for:
  - blank or missing `und_icd`
  - malformed/short codes
  - chapter grouping based on the leading ICD-10 letter and code range

## Analysis Cohort Construction (Step-by-Step)

1. Load `gbcs_main.dta` and `newvital_2023.dta`.
2. Harmonize `obje_id` type (convert both to character, trim spaces).
3. Merge 1:1 by `obje_id`.
4. Derive baseline IC exposure(s) from `gbcs_main.dta` using existing scoring logic.
5. Restrict to participants with valid baseline IC exposure.
6. Restrict to participants with valid mortality status (`dead` in `{0,1}`) and non-missing `futime`.
7. Create analysis flags and report exclusions at each step.

## Confounders (Prespecified)

The following choices are confirmed.

### Core confounders (primary adjustment set)

- `agec` (continuous; preferred age variable)
- `sex`
- `phase` (baseline phase)

### Socioeconomic / social confounders

- `edu` (education; categorical)
- `incomeh` (household annual income; categorical/ordinal as coded)
- `marital` (marital status)

### Lifestyle factors (choose lowest-missingness versions where overlap exists)

Select variables after a completeness comparison and prespecify before modeling.

- Smoking:
  - compare `sm_stat` vs `sm_stat2`
  - use the more complete variable as primary
  - if completeness is similar, prefer `sm_stat` (more detailed categories)
- Alcohol:
  - prefer `drk1` if completeness/coding is acceptable
- Physical activity:
  - choose an interpretable baseline activity/walking/exercise variable with minimal missingness (from available activity variables)

### Baseline comorbidity confounders (candidate pool)

Use one definition per construct to avoid redundancy/collinearity.

- Hypertension: `hxht` or `hyp`
- Diabetes: `hxdm` or `diab`
- CHD: `hxchd`
- Stroke: `hxstroke`
- Cancer history: `v12_1a`
- COPD / chronic respiratory disease: `cxcopd` (recode as needed)

## Pre-Model Data Checks (No Association Testing Yet)

1. Verify `futime` units and internal consistency using `reg_date` / `ddate` where available.
2. Check coding consistency between `dead` and `vstatus`.
3. Standardize `und_icd` formatting (uppercase, trim spaces, retain dot handling rule).
4. Summarize missingness for:
   - exposure(s)
   - outcome/time
   - all candidate confounders
5. Compare completeness of overlapping lifestyle variables (especially `sm_stat` vs `sm_stat2`) before finalizing covariates.

## Descriptive Analysis Plan (Before Modeling)

1. Create cohort flow summary (starting N, exclusions, final analytic N).
2. Produce baseline descriptive table by IC category (or IC quantiles/categories):
   - demographics
   - SES
   - lifestyle
   - comorbidities
3. Summarize deaths and person-time by IC category.
4. For cause-specific work, tabulate death counts by ICD cause group.

## Primary Statistical Models (All-Cause Mortality)

### Main approach

- Cox proportional hazards regression with `Surv(futime, dead == 1)`

### Sequential adjustment models

1. Model 1 (minimal):
   - IC + `agec` + `sex` + `phase`
2. Model 2 (SES/social):
   - Model 1 + `edu` + `incomeh` + `marital`
3. Model 3 (fully adjusted primary model):
   - Model 2 + selected smoking + alcohol + physical activity + baseline comorbidities

### Effect parameterization

- IC continuous: HR per 1-point lower (or higher) IC (specify direction consistently)
- IC categorical: HRs vs highest IC category + trend test

### Optional nonlinear assessment

- Restricted cubic splines for continuous IC (and for `agec` if needed)

## Cause-Specific Mortality Models

### Primary cause-specific approach (etiologic)

- Cause-specific Cox models for each prespecified ICD group
- Treat deaths from other causes as censored at death time

### Secondary competing-risk approach (if cumulative incidence is of interest)

- Fine-Gray subdistribution hazards (sensitivity/secondary)

### Event-count threshold rule

- Prespecify a minimum number of events before fitting a detailed cause-specific model (e.g., combine sparse groups if too few events)

## Model Diagnostics and Robustness Checks

1. Proportional hazards assumption:
   - Schoenfeld residual tests/plots
2. Functional form:
   - assess continuous IC and `agec`
3. Collinearity:
   - especially among SES and comorbidity variables
4. Influential observations / extreme follow-up values:
   - review diagnostics before interpreting estimates

## Missing Data Strategy

### Primary

- Complete-case analysis for the fully adjusted model

### Sensitivity

- Multiple imputation (if missingness is non-trivial and assumptions are reasonable)
- Compare estimates from complete-case vs imputed analyses

## Sensitivity Analyses

1. Exclude early deaths (e.g., first 1-2 years) to reduce reverse causation.
2. Alternate exposure coding:
   - continuous vs categorical IC
3. Alternate confounder definitions:
   - `sm_stat` vs `sm_stat2`
   - `hxht` vs `hyp`
   - `hxdm` vs `diab`
4. Cause-specific analysis restricted to deaths with valid `und_icd`.
5. Stratified analyses (sex and/or age groups) only if event counts support stable estimates.

## Outputs for the Analysis Stage

- Cohort flow diagram/table
- Baseline characteristics table by IC category
- Kaplan-Meier curves by IC category
- Cox model results table (Models 1-3)
- Cause-specific mortality results table with ICD grouping definitions
- Missingness summary and sensitivity analysis comparison table

## Notes / Decisions Locked In

- Use `futime` as follow-up duration.
- Use `dead` as vital status for all-cause mortality event definition.
- Use `und_icd` for cause-specific mortality classification among deaths.
- Prefer `agec` (continuous), `edu`, and `incomeh` as key confounders.
- For lifestyle covariates, choose the lowest-missingness variables after a formal completeness comparison (including `sm_stat` vs `sm_stat2`).
