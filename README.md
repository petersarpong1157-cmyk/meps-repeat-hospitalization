# MEPS 2024 Repeat Hospitalization Reproducibility Files

## Study
Selected Five-Condition Chronic Disease Count and Repeat Hospitalization Among U.S. Adults With Inpatient Hospital Use: Evidence from the 2024 Medical Expenditure Panel Survey

## Archived release
Version v3.0.3 of these reproducibility materials is archived on Zenodo.

**DOI:** 10.5281/zenodo.22928116

## Data
This analysis uses the 2024 Medical Expenditure Panel Survey (MEPS) Household Component Full-Year Consolidated public-use file (HC-256), maintained by the Agency for Healthcare Research and Quality (AHRQ).

The MEPS public-use dataset is not redistributed in this repository. Obtain HC-256 from AHRQ and place the Stata-format data file in the working directory. If its filename is not `h256.dta`, update the `data_file` value near the top of the R script.

## Main file
`MEPS2024_repeat_hospitalization_reproducibility.R`

The script reconstructs the final analytic workflow used for the manuscript, including:

- the adult hospitalized analytic domain and five selected diagnosed chronic conditions;
- the MEPS complex survey design;
- survey-weighted descriptive estimates;
- the primary continuous-count survey-weighted logistic regression;
- categorical and binary exposure specifications;
- the quadratic functional-form assessment;
- the survey-adjusted Rao-Scott likelihood-ratio comparison;
- the disease-combination count used in the limitations discussion; and
- Figures 1–3 and key output files.

## Required R packages
- haven
- dplyr
- survey
- ggplot2
- scales

The script intentionally does not install packages automatically.

## Expected verification checks
The script stops if it does not reproduce these core analytic-sample quantities:

- Analytic sample: N = 1,247
- Condition-count frequencies for 0–5 conditions: 226, 318, 337, 236, 112, 18
- Inpatient-discharge outcome counts: 963 with exactly one discharge and 284 with two or more
- Observed five-condition disease combinations: 31

These checks are intended to detect use of an incorrect file or unintended changes to the analytic sample.

## Reproducibility note
This is a cleaned version of the original working analysis. Exploratory commands, duplicated code, package-install commands, machine-specific paths, and superseded figure code were removed. The final analysis specifications were retained.

## License
The analysis code in this repository is released under the MIT License.
