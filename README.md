# High-BMI-Attributable Burden of Gout, Osteoarthritis, and Low Back Pain in China and G20 Countries

This repository contains the R analysis scripts for the manuscript
"High-BMI-Attributable Burden of Gout, Osteoarthritis, and Low Back Pain in
China and G20 Countries: A Multi-Disease Comparative Analysis of Intervention
Priorities and Health Inequities".

## Data source

Data were extracted from the Global Burden of Disease (GBD) 2023 study via the
GBD Results Tool (https://vizhub.healthdata.org/gbd-results/) on 26 January 2026.

Query parameters:

- GBD round: GBD 2023
- Risk: high body-mass index
- Causes: gout, osteoarthritis, low back pain
- Measure: DALYs
- Metric: rate and number
- Age groups: age-standardized, all ages, 5-year age groups
- Sex: male, female, both
- Locations: 19 sovereign G20 member states
- Years: 1990–2023

The 19 sovereign G20 member states included in the analysis are: Argentina,
Australia, Brazil, Canada, China, France, Germany, India, Indonesia, Italy,
Japan, Mexico, Russia, Saudi Arabia, South Africa, South Korea, Türkiye, the
United Kingdom, and the United States. The European Union and the African Union
were not included (the European Union is not a sovereign country and has no
separate SDI value in the GBD framework; the African Union joined the G20 only
in 2023, after the period covered by most analyses).

## Software environment

- R version 4.5.2
- Main packages: tidyverse, segmented, ggrepel, tidytext. The complete list of
  packages required by each script is given at the top of the corresponding R
  script.

## File structure

- `code/` — R scripts
  - `01_trend_analysis.R`
  - `02_age_sex_distribution.R`
  - `03_decomposition_analysis_Code_CHINA.R`
  - `04_decomposition_analysis_Code_G20.R`
  - `05_joinpoint_threshold_monte_carlo_and_crossing_year.R`
  - `06_inequality_computation.R`
  - `07_inequality_visualization.R`

## Description of scripts

| Script | Description |
|---|---|
| `01_trend_analysis.R` | Trend analysis of age-standardized DALY rates for China and the G20 median. |
| `02_age_sex_distribution.R` | Age- and sex-specific DALY rates attributable to high BMI for the five representative G20 countries. |
| `03_decomposition_analysis_Code_CHINA.R` | Das Gupta decomposition of the DALY increase in China into population growth, aging, and epidemiological change. |
| `04_decomposition_analysis_Code_G20.R` | Das Gupta decomposition for the 19 sovereign G20 member states, including SDI-stratified group summaries. |
| `05_joinpoint_threshold_monte_carlo_and_crossing_year.R` | Joinpoint regression, fixed absolute threshold, Monte Carlo projection for China (2024–2050), and threshold-crossing year estimation. |
| `06_inequality_computation.R` | Computation of the Slope Index of Inequality (SII) and Concentration Index (CI). |
| `07_inequality_visualization.R` | Visualization of the inequality results. |

## Data availability

The cleaned analysis dataset used in this study is derived from the GBD 2023
Results Tool and is subject to the IHME Free-of-Charge Non-Commercial User
Agreement. Because this agreement does not permit open redistribution of GBD
data, the cleaned dataset is **not included in this repository**.

The cleaned analysis dataset is available from the corresponding author upon
reasonable request. Access will be granted only to researchers from academic or
non-profit institutions who confirm in writing that the data will be used
exclusively for non-commercial research purposes and who agree not to
redistribute the data to third parties, in accordance with the IHME
Free-of-Charge Non-Commercial User Agreement.

## Contact

Song Yun Peng
The Second Affiliated Hospital of Guizhou University of Traditional Chinese Medicine
Guiyang, Guizhou, China
Email: 953648587@qq.com