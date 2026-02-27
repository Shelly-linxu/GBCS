# Physical Activity Change and All-Cause Mortality: Results Summary

## 1) Cohort flow
- Main rows: 30518
- Matched to mortality: 30020
- Valid all-cause event coding: 30020
- Valid positive futime: 30097

## 2) Wave exposure completeness
- Baseline IPAQ non-missing: 30518
- F1 IPAQ non-missing: 18100
- F2 IPAQ non-missing: 10566
- F3 IPAQ non-missing: 5379
- Baseline MET non-missing: 30427
- F1 MET non-missing: 18100
- F2 MET non-missing: 10548
- F3 MET non-missing: 5369

## 3) Landmark interval samples
- BL->F1: N=18103, events=3944
- F1->F2: N=10547, events=1383
- F2->F3: N=5372, events=479

## 4) Primary model effects (M3)
Per-interval and pooled effects are in `cox_landmark_main_effects.csv` and `cox_landmark_pooled_main_effects.csv`.

### Pooled fixed-effect estimates (M3)
- ipaq_directionimproved: HR=0.991 (95% CI 0.915, 1.073), p=0.8210, intervals=3
- ipaq_directionworsened: HR=1.038 (95% CI 0.957, 1.127), p=0.3670, intervals=3
- met_change_sd: HR=1.023 (95% CI 0.959, 1.092), p=0.4883, intervals=3

## 5) Model failures
- No model fitting failures.
