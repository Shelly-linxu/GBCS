# IC-Mortality Analysis Results (Stepwise)

## Data and cohort
- Main rows: 30518
- Mortality rows: 30430
- Analytic base cohort (valid baseline IC + all-cause follow-up): 6550
- All-cause deaths in analytic base cohort: 1161

## IC category counts

  high middle    low 
   547   3904   2099 

## Primary Cox models (IC effect)
- M1 (continuous ic5_per1lower): 1.175 (1.127, 1.225), p=4.51e-14 [n=6550, events=1161]
- M2 (continuous ic5_per1lower): 1.165 (1.114, 1.218), p=1.66e-11 [n=6477, events=1145]
- M3 (continuous ic5_per1lower): 1.146 (1.089, 1.207), p=2.01e-07 [n=5083, events=902]
- M3 (IC category middle vs high): 1.104 (0.835, 1.460), p=0.486
- M3 (IC category low vs high): 1.497 (1.101, 2.035), p=0.01

## Diagnostics
- GLOBAL PH test p=0.0476

## Sensitivity models (continuous ic5_per1lower)
                      model_id       hr ci_lower ci_upper      p_value    n
1              S1_M3_ic5_lag1y 1.141750 1.083966 1.202613 5.652338e-07 5070
2                S2_M3_ic5_syn 1.148929 1.091223 1.209687 1.289269e-07 5083
3          S3_M3_ic5_hxht_hxdm 1.140962 1.083725 1.201222 5.116785e-07 5112
4 S4_M3_ic5_marital_missingcat 1.144170 1.087215 1.204109 2.344644e-07 5138
5        S5_M3_ic5_plus_cxcopd 1.159118 1.100300 1.221080 2.739055e-08 5022
  nevent
1    889
2    902
3    908
4    914
5    891

## Cause-specific Cox models (broad ICD groups; continuous ic5_per1lower)
  cause_group              model_id        hr ci_lower ci_upper      p_value
1 circulatory CS_circulatory_M3_ic5 1.2847118 1.180232 1.398440 7.081953e-09
2      cancer      CS_cancer_M3_ic5 0.9603277 0.880415 1.047494 3.611321e-01
3 respiratory CS_respiratory_M3_ic5 1.2117381 1.045612 1.404258 1.068526e-02
4       other       CS_other_M3_ic5 1.2987998 1.130244 1.492493 2.275901e-04
     n nevent
1 5083    340
2 5083    322
3 5083    113
4 5083    125

## Notes
- Smoking primary confounder used in fitted models: `smk` (lower missingness than `sm_stat`/`sm_stat2`).
- Physical activity primary confounder used in fitted models: binary `walk_any` from `v7_5a`.
- `marital` has high missingness; both primary complete-case and missing-category sensitivity models were run.
- `cxcopd` is sensitivity-only due high missingness.
