# IC-Mortality Variable Dictionary and Model Formula Sheet (Concrete Spec)

## Scope

This document provides the concrete implementation specification for the planned analysis of intrinsic capacity (IC) and mortality:

- exact variable choices
- exact coding/recode rules
- derived analysis variables
- model formula sheet (primary and sensitivity models)

This is a data-prep and analysis-specification document only (no association analyses are run here).

## Source Datasets

- Mortality follow-up: `/Users/linxu/Documents/GBCS/newvital_2023.dta`
- Baseline IC + covariates: `/Users/linxu/Documents/GBCS/gbcs_main.dta`

## Linkage Rule

- Merge key: `obje_id`
- Harmonize type before merge:
  - `newvital_2023.dta::obje_id` is numeric
  - `gbcs_main.dta::obje_id` is character
- Exact rule:
  - convert both to character
  - `trimws()`
  - empty strings -> `NA`
- Merge expected: 1:1 by `obje_id`

## Final Variable Choices (Locked for Current Spec)

### Exposure

- Primary continuous: baseline `ic5` (derived from existing IC scoring logic)
- Primary categorical: baseline `ic3_classic` (derived from `ic5`)

### Outcomes

- All-cause mortality:
  - `futime`
  - `dead`
- Cause-specific mortality:
  - `und_icd`

### Confounders (primary model set)

- Core: `agec`, `sex`, `phase`
- SES/social: `edu`, `incomeh`, `marital`
- Lifestyle: `smk`, `drk1`, `v7_5a` (recode to binary walking)
- Adiposity: `bmi`
- Comorbidities: `hyp`, `diab`, `hxchd`, `hxstroke`, `v12_1a`

### Variables retained for sensitivity (not primary)

- Smoking alternative: `syn`, `sm_stat`, `sm_stat2`
- Hypertension alternative: `hxht`
- Diabetes alternative: `hxdm`
- COPD/CXR: `cxcopd` (very high missingness; not primary)
- Physical activity alternative coding: 3-level `v7_5a`

## Exact Coding / Recode Rules

## 1. Outcomes

### 1.1 All-cause mortality event (`event_allcause`)

Source: `newvital_2023.dta::dead`

Raw coding (labelled):

- `0` = Alive
- `1` = Dead
- `2` = Unknown

Recode:

- `event_allcause = 1` if `dead == 1`
- `event_allcause = 0` if `dead == 0`
- `event_allcause = NA` if `dead == 2` or missing

### 1.2 Survival time (`time_allcause`)

Source: `newvital_2023.dta::futime`

Recode:

- `time_allcause = futime`
- Set to `NA` if `futime` missing
- Do not transform in primary models

### 1.3 Cause-of-death code cleaning (`und_icd_clean`)

Source: `newvital_2023.dta::und_icd`

Recode:

- convert to upper case
- trim spaces
- internal spaces removed
- blank string -> `NA`

Derived helper:

- `icd3 = sub("\\\\..*$", "", und_icd_clean)` (ICD code stem before decimal)

## 2. Exposure (IC)

Exposure derivation should follow existing repo logic in `/Users/linxu/Documents/breast cancer/GBCS/ic_wave_distribution_analysis.R`.

### 2.1 Baseline full IC (`ic5`)

- Derived baseline 5-domain IC total
- Expected range: `0-10`
- Primary continuous exposure

### 2.2 Baseline IC category (`ic3_classic`)

Exact thresholds (from existing script):

- `0` if `ic5 >= 9` (highest IC; reference category)
- `1` if `6 <= ic5 < 9`
- `2` if `ic5 < 6` (lowest IC)
- `NA` if `ic5` missing

Derived reporting version (recommended):

- `ic3_classic_fct = factor(ic3_classic, levels = c(0,1,2), labels = c(\"high\", \"middle\", \"low\"))`

### 2.3 Continuous effect direction helper (`ic5_per1lower`)

For interpretable HRs:

- `ic5_per1lower = -1 * ic5`

Then HR > 1 means higher mortality risk per 1-point lower IC.

## 3. Confounders (Exact Rules)

## 3.1 Core confounders

### `agec` (continuous age; primary age variable)

Source: `gbcs_main.dta::agec`

Recode:

- Use as continuous numeric
- No categorization in primary model
- No winsorization in primary model

Derived scaling for interpretation (recommended):

- `agec_5 = agec / 5`

Use `agec_5` in model output if HR per 5-year increase is preferred.

### `sex`

Source: `gbcs_main.dta::sex`

Raw coding:

- `0` = female
- `1` = male

Recode:

- `sex_male = sex` (0/1)
- or `sex_fct = factor(sex, levels = c(0,1), labels = c(\"female\", \"male\"))`

Reference:

- `female`

### `phase`

Source: `gbcs_main.dta::phase`

Observed raw coding:

- `1`, `2`, `3` (no value labels)

Recode:

- `phase_fct = factor(phase, levels = c(1,2,3))`
- Set any non-`1/2/3` value to `NA` (if encountered)

Reference:

- `phase = 1`

## 3.2 SES / Social confounders

### `edu` (education)

Source: `gbcs_main.dta::edu`

Raw coding:

- `1` = `<primary`
- `2` = `primary`
- `3` = `junior middle`
- `4` = `senior middle`
- `5` = `junior college`
- `6` = `college`

Recode:

- Keep as categorical 6-level factor
- valid values: `1:6`
- any other code -> `NA`

Factor coding:

- `edu_fct = factor(edu, levels = c(6,5,4,3,2,1), labels = c(\"college\", \"junior_college\", \"senior_middle\", \"junior_middle\", \"primary\", \"lt_primary\"))`

Reference:

- `college` (highest education)

### `incomeh` (household income)

Source: `gbcs_main.dta::incomeh`

Raw coding (stored labels):

- `1` = `<5000 rmb`
- `2` = `5000-`
- `3` = `10000-`
- `4` = `20000-`
- `5` = `30000-`
- `6` = `>=50000`
- `7` = `notknow`

Recode:

- Keep code `7` (`notknow`) as a valid category (not missing)
- valid values: `1:7`
- any other code -> `NA`

Factor coding:

- `incomeh_fct = factor(incomeh, levels = c(6,5,4,3,2,1,7), labels = c(\"ge50000\", \"30000_band\", \"20000_band\", \"10000_band\", \"5000_band\", \"lt5000\", \"notknow\"))`

Reference:

- `>=50000` (code `6`)

Note:

- Middle-band labels are truncated in metadata export; preserve original raw codes and stored labels exactly.

### `marital` (high missingness; keep primary, add missing-category sensitivity)

Source: `gbcs_main.dta::marital`

Raw labelled coding:

- `0` = married
- `1` = separate
- `2` = widow
- `3` = never marry

Observed nonstandard unlabeled code:

- `56` (single record observed)

Primary recode:

- `marital_clean = marital`
- valid values: `0,1,2,3`
- set `56` and any other unlabeled code to `NA`

Factor coding:

- `marital_fct = factor(marital_clean, levels = c(0,1,2,3), labels = c(\"married\", \"separate\", \"widow\", \"never_married\"))`

Reference:

- `married`

Sensitivity recode (retain sample size):

- `marital5_fct = forcats::fct_explicit_na(marital_fct, na_level = \"missing_or_other\")`

## 3.3 Lifestyle confounders (minimal-missingness selection)

### Smoking (primary): `smk`

Source: `gbcs_main.dta::smk`

Why selected:

- Much lower missingness than `sm_stat`/`sm_stat2`
- Similar construct (smoking status)

Raw coding:

- `0` = never
- `1` = ex-occasional
- `2` = ex-daily
- `3` = current

Recode:

- valid values: `0,1,2,3`
- any other code -> `NA`

Factor coding:

- `smk_fct = factor(smk, levels = c(0,1,2,3), labels = c(\"never\", \"ex_occasional\", \"ex_daily\", \"current\"))`

Reference:

- `never`

Smoking sensitivity variables:

- `syn` (0=never, 1=ever) -> `syn_fct`
- `sm_stat` / `sm_stat2` for alternate definitions if needed

### Alcohol (primary): `drk1`

Source: `gbcs_main.dta::drk1`

Raw coding:

- `0` = never
- `1` = <1/month
- `2` = <1/week
- `3` = 1-4/week
- `4` = 5+/week
- `5` = ex-drinker

Recode:

- valid values: `0:5`
- any other code -> `NA`

Factor coding:

- `drk1_fct = factor(drk1, levels = c(0,1,2,3,4,5), labels = c(\"never\", \"lt1_month\", \"lt1_week\", \"1to4_week\", \"ge5_week\", \"ex_drinker\"))`

Reference:

- `never`

### Physical activity / walking (primary): `v7_5a` -> binary walking variable

Source: `gbcs_main.dta::v7_5a`

Raw coding (labelled `any walk in last week`):

- `1` = yes, >=10min/time
- `2` = yes, <10min/time
- `3` = no

Primary recode (binary, minimal sparse-category issue):

- `walk_any = 1` if `v7_5a %in% c(1,2)`
- `walk_any = 0` if `v7_5a == 3`
- otherwise `NA`

Factor coding:

- `walk_any_fct = factor(walk_any, levels = c(0,1), labels = c(\"no\", \"yes\"))`

Reference:

- `no`

Sensitivity recode (3-level):

- `walk3_fct = factor(v7_5a, levels = c(3,2,1), labels = c(\"no\", \"yes_lt10min\", \"yes_ge10min\"))`

## 3.4 Adiposity confounder

### `bmi`

Source: `gbcs_main.dta::bmi`

Recode:

- use raw continuous BMI (kg/m^2)
- `NA` stays `NA`
- no truncation/winsorization in primary model

Optional derived scaling:

- `bmi_5 = bmi / 5` (HR per 5 kg/m^2 increase)

Optional nonlinear sensitivity:

- spline term in Cox model

## 3.5 Baseline comorbidity confounders

### Hypertension (primary): `hyp`

Source: `gbcs_main.dta::hyp`

Raw coding:

- `0` = no
- `1` = yes

Recode:

- `hyp_bin = ifelse(hyp %in% c(0,1), hyp, NA)`
- factor version: `hyp_fct` with `no` as reference

Reference:

- `no`

Sensitivity alternative:

- `hxht` (self-reported hypertension)

### Diabetes (primary): `diab`

Source: `gbcs_main.dta::diab`

Raw coding:

- `0` = no
- `1` = yes

Recode:

- `diab_bin = ifelse(diab %in% c(0,1), diab, NA)`

Reference:

- `no`

Sensitivity alternative:

- `hxdm` (self-reported type II diabetes)

### CHD: `hxchd`

Source: `gbcs_main.dta::hxchd`

Raw coding:

- `0` = no
- `1` = yes

Recode:

- `hxchd_bin = ifelse(hxchd %in% c(0,1), hxchd, NA)`

Reference:

- `no`

### Stroke: `hxstroke`

Source: `gbcs_main.dta::hxstroke`

Raw coding:

- `0` = no
- `1` = yes

Recode:

- `hxstroke_bin = ifelse(hxstroke %in% c(0,1), hxstroke, NA)`

Reference:

- `no`

### Cancer history: `v12_1a`

Source: `gbcs_main.dta::v12_1a`

Raw coding:

- `1` = no
- `2` = yes

Recode:

- `cancerhx_bin = 0` if `v12_1a == 1`
- `cancerhx_bin = 1` if `v12_1a == 2`
- otherwise `NA`

Reference:

- `no`

### COPD / CXR abnormality (sensitivity only): `cxcopd`

Source: `gbcs_main.dta::cxcopd`

Raw coding:

- `30` = none
- `31` = mild
- `32` = moderate
- `33` = severe
- `34` = increased lung markings

Primary decision:

- Not included in primary fully adjusted model due very high missingness and sparse non-missing categories.

Sensitivity recode (binary):

- `cxcopd_any = 0` if `cxcopd == 30`
- `cxcopd_any = 1` if `cxcopd %in% c(31,32,33,34)`
- otherwise `NA`

## Exact Recode Implementation Template (R / dplyr)

```r
library(dplyr)
library(forcats)
library(stringr)

dat <- merged_dat %>%
  mutate(
    # Outcome
    event_allcause = case_when(
      dead == 1 ~ 1L,
      dead == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    time_allcause = futime,
    und_icd_clean = und_icd %>%
      as.character() %>%
      str_to_upper() %>%
      str_trim() %>%
      str_replace_all("\\\\s+", "") %>%
      na_if(""),
    icd3 = if_else(is.na(und_icd_clean), NA_character_, str_replace(und_icd_clean, "\\\\..*$", "")),

    # Exposure helpers
    ic5_per1lower = -1 * ic5,
    ic3_classic_fct = factor(ic3_classic, levels = c(0,1,2), labels = c("high", "middle", "low")),

    # Core confounders
    agec_5 = agec / 5,
    sex_fct = factor(if_else(sex %in% c(0,1), sex, NA_real_), levels = c(0,1), labels = c("female", "male")),
    phase_fct = factor(if_else(phase %in% c(1,2,3), phase, NA_real_), levels = c(1,2,3)),

    # SES/social
    edu_fct = factor(if_else(edu %in% 1:6, edu, NA_real_),
                     levels = c(6,5,4,3,2,1),
                     labels = c("college","junior_college","senior_middle","junior_middle","primary","lt_primary")),
    incomeh_fct = factor(if_else(incomeh %in% 1:7, incomeh, NA_real_),
                         levels = c(6,5,4,3,2,1,7),
                         labels = c("ge50000","30000_band","20000_band","10000_band","5000_band","lt5000","notknow")),
    marital_clean = if_else(marital %in% c(0,1,2,3), marital, NA_real_),
    marital_fct = factor(marital_clean, levels = c(0,1,2,3),
                         labels = c("married","separate","widow","never_married")),

    # Lifestyle
    smk_fct = factor(if_else(smk %in% c(0,1,2,3), smk, NA_real_),
                     levels = c(0,1,2,3),
                     labels = c("never","ex_occasional","ex_daily","current")),
    syn_fct = factor(if_else(syn %in% c(0,1), syn, NA_real_),
                     levels = c(0,1),
                     labels = c("never","ever")),
    drk1_fct = factor(if_else(drk1 %in% 0:5, drk1, NA_real_),
                      levels = c(0,1,2,3,4,5),
                      labels = c("never","lt1_month","lt1_week","1to4_week","ge5_week","ex_drinker")),
    walk_any = case_when(
      v7_5a %in% c(1,2) ~ 1L,
      v7_5a == 3 ~ 0L,
      TRUE ~ NA_integer_
    ),
    walk_any_fct = factor(walk_any, levels = c(0,1), labels = c("no","yes")),
    walk3_fct = factor(if_else(v7_5a %in% c(1,2,3), v7_5a, NA_real_),
                       levels = c(3,2,1),
                       labels = c("no","yes_lt10min","yes_ge10min")),

    # Comorbidities
    hyp_fct = factor(if_else(hyp %in% c(0,1), hyp, NA_real_), levels = c(0,1), labels = c("no","yes")),
    hxht_fct = factor(if_else(hxht %in% c(0,1), hxht, NA_real_), levels = c(0,1), labels = c("no","yes")),
    diab_fct = factor(if_else(diab %in% c(0,1), diab, NA_real_), levels = c(0,1), labels = c("no","yes")),
    hxdm_fct = factor(if_else(hxdm %in% c(0,1), hxdm, NA_real_), levels = c(0,1), labels = c("no","yes")),
    hxchd_fct = factor(if_else(hxchd %in% c(0,1), hxchd, NA_real_), levels = c(0,1), labels = c("no","yes")),
    hxstroke_fct = factor(if_else(hxstroke %in% c(0,1), hxstroke, NA_real_), levels = c(0,1), labels = c("no","yes")),
    cancerhx_bin = case_when(
      v12_1a == 1 ~ 0L,
      v12_1a == 2 ~ 1L,
      TRUE ~ NA_integer_
    ),
    cancerhx_fct = factor(cancerhx_bin, levels = c(0,1), labels = c("no","yes")),

    # Sensitivity-only COPD
    cxcopd_any = case_when(
      cxcopd == 30 ~ 0L,
      cxcopd %in% c(31,32,33,34) ~ 1L,
      TRUE ~ NA_integer_
    ),
    cxcopd_any_fct = factor(cxcopd_any, levels = c(0,1), labels = c("none","any_abnormal"))
  )
```

## Analytic Cohort Rules (Exact)

Apply in this order:

1. Successful linkage by non-missing `obje_id`
2. Non-missing IC exposure (`ic5` for continuous analyses; `ic3_classic` for categorical analyses)
3. Valid all-cause event coding (`event_allcause` non-missing)
4. Non-missing `time_allcause`
5. Positive follow-up time check (if any `futime <= 0`, flag and review; do not silently drop without logging)

## Model Formula Sheet (R Syntax + Statistical Meaning)

All primary time-to-event models use Cox proportional hazards regression.

Package:

- `survival::coxph`

Outcome:

- `Surv(time_allcause, event_allcause)`

## A. All-Cause Mortality (Primary)

### A1. Continuous IC (primary effect estimate)

Interpretation target: HR per 1-point lower baseline IC (`ic5_per1lower`)

### Model 1 (minimal)

```r
m1_ic5 <- survival::coxph(
  survival::Surv(time_allcause, event_allcause) ~
    ic5_per1lower + agec_5 + sex_fct + phase_fct,
  data = analysis_dat
)
```

### Model 2 (SES/social adjusted)

```r
m2_ic5 <- survival::coxph(
  survival::Surv(time_allcause, event_allcause) ~
    ic5_per1lower + agec_5 + sex_fct + phase_fct +
    edu_fct + incomeh_fct + marital_fct,
  data = analysis_dat
)
```

### Model 3 (fully adjusted primary)

```r
m3_ic5 <- survival::coxph(
  survival::Surv(time_allcause, event_allcause) ~
    ic5_per1lower + agec_5 + sex_fct + phase_fct +
    edu_fct + incomeh_fct + marital_fct +
    smk_fct + drk1_fct + walk_any_fct + bmi +
    hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct,
  data = analysis_dat
)
```

## B. All-Cause Mortality (Categorical IC)

Reference category: `ic3_classic_fct = "high"`

### Model 3 (fully adjusted; categorical IC)

```r
m3_ic3 <- survival::coxph(
  survival::Surv(time_allcause, event_allcause) ~
    ic3_classic_fct + agec_5 + sex_fct + phase_fct +
    edu_fct + incomeh_fct + marital_fct +
    smk_fct + drk1_fct + walk_any_fct + bmi +
    hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct,
  data = analysis_dat
)
```

Trend test option:

- fit ordinal numeric term `ic3_classic` (0/1/2) in separate model

## C. Nonlinearity Checks (Sensitivity / Extension)

### Continuous IC spline model

```r
m3_ic5_spline <- survival::coxph(
  survival::Surv(time_allcause, event_allcause) ~
    splines::ns(ic5, df = 4) + agec_5 + sex_fct + phase_fct +
    edu_fct + incomeh_fct + marital_fct +
    smk_fct + drk1_fct + walk_any_fct + splines::ns(bmi, df = 4) +
    hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct,
  data = analysis_dat
)
```

## D. Cause-Specific Mortality (Etiologic; Cause-Specific Cox)

Create broad cause groups from `und_icd_clean` (final grouping file can be documented separately).

Example broad groups:

- `cancer`: ICD-10 starts with `C`
- `circulatory`: starts with `I`
- `respiratory`: starts with `J`
- `other`: all other valid codes

Example event coding for circulatory mortality:

- `event_cvd = 1` if `event_allcause == 1` and `icd_group == "circulatory"`
- `event_cvd = 0` if alive/censored OR died from another cause (censored at death time)
- `NA` only if all-cause event/time invalid

Formula (fully adjusted):

```r
m3_cvd <- survival::coxph(
  survival::Surv(time_allcause, event_cvd) ~
    ic5_per1lower + agec_5 + sex_fct + phase_fct +
    edu_fct + incomeh_fct + marital_fct +
    smk_fct + drk1_fct + walk_any_fct + bmi +
    hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct,
  data = analysis_dat
)
```

Repeat for each prespecified cause group if event counts are adequate.

## E. Sensitivity Models (Pre-Specified)

### E1. Early-death exclusion (reverse-causation sensitivity)

```r
m3_ic5_lag1y <- survival::coxph(
  survival::Surv(time_allcause, event_allcause) ~
    ic5_per1lower + agec_5 + sex_fct + phase_fct +
    edu_fct + incomeh_fct + marital_fct +
    smk_fct + drk1_fct + walk_any_fct + bmi +
    hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct,
  data = dplyr::filter(analysis_dat, time_allcause > 1)
)
```

### E2. Alternate smoking definition (binary ever/never)

```r
m3_ic5_syn <- survival::coxph(
  survival::Surv(time_allcause, event_allcause) ~
    ic5_per1lower + agec_5 + sex_fct + phase_fct +
    edu_fct + incomeh_fct + marital_fct +
    syn_fct + drk1_fct + walk_any_fct + bmi +
    hyp_fct + diab_fct + hxchd_fct + hxstroke_fct + cancerhx_fct,
  data = analysis_dat
)
```

### E3. Alternate HTN/DM definitions (self-reported)

Replace:

- `hyp_fct` -> `hxht_fct`
- `diab_fct` -> `hxdm_fct`

### E4. Marital missing-category sensitivity

Replace `marital_fct` with `marital5_fct`.

### E5. COPD/CXR sensitivity (high missingness; exploratory)

Add `cxcopd_any_fct` to the fully adjusted model and report the reduced complete-case N.

## Missing-Data Handling Rules (Primary vs Sensitivity)

### Primary

- Complete-case per fitted model (`coxph` default `na.omit`)
- Report model-specific analytic N and number of events

### Required reporting (for every model)

- N included
- number of deaths/events
- number excluded due to missing covariates

### Sensitivity

- Multiple imputation (if implemented later)
- Marital missing-category model (`marital5_fct`)

## QC Checks Required Before Fitting

1. Confirm no duplicate `obje_id` after merge.
2. Confirm `event_allcause` only takes `0/1`.
3. Confirm `time_allcause > 0` for included records.
4. Confirm no invalid raw codes remain after recode for categorical confounders.
5. Report frequencies for each final factor level (including sparse levels).
6. Confirm reference levels are set as specified before model fitting.

## Notes on Variable Selection Changes vs Earlier Draft

- Smoking primary confounder is now `smk` (not `sm_stat`/`sm_stat2`) because `smk` has much lower missingness.
- Physical activity primary confounder is `v7_5a` (recode to `walk_any`) rather than `walk`, because `v7_5a` has much lower missingness.
- COPD/CXR (`cxcopd`) is not in the primary fully adjusted model because missingness is very high.
