# Infant-to-Adult Attachment Continuity Meta-Analysis

## Overview

This project conducts meta-analyses examining the continuity of attachment patterns from infancy/childhood to adulthood, using Strange Situation Procedure (SSP) or Attachment Q-Sort (AQS) measures in infancy and the Adult Attachment Interview (AAI) in adulthood.

## Six Meta-Analytic Syntheses

1. **2-way S/IS (Pearson's r)** - Secure/Insecure continuity, pooled SSP+AQS → AAI
2. **2-way O/D (Pearson's r)** - Organized/Disorganized continuity, pooled SSP+AQS → AAI
3. **3-way (Cohen's kappa)** - A/B/C classification continuity, pooled SSP+AQS → AAI
4. **4-way (Cohen's kappa)** - A/B/C/D classification continuity, pooled SSP+AQS → AAI
5. **SSP-only S/IS (Pearson's r)** - Secure/Insecure continuity, SSP → AAI only
6. **AQS-only S/IS (Pearson's r)** - Secure/Insecure continuity, AQS → AAI only

## Data Source

- **Original file**: `data/data_extraction_table.csv`

## Folder Structure

```
InfantAdultAttachment_MetaAnalysis/
├── Code/
│   ├── 01_PrepareData.R              # Parse Excel → clean CSV
│   ├── 02_CalculateEffectSizes.R     # Compute r and kappa
│   ├── 03_MetaAnalysis_2Way_SIS_r.R  # Analysis 1
│   ├── 04_MetaAnalysis_2Way_OD_r.R   # Analysis 2
│   ├── 05_MetaAnalysis_3Way_k.R      # Analysis 3
│   ├── 06_MetaAnalysis_4Way_k.R      # Analysis 4
│   ├── 07_MetaAnalysis_SSP_only.R    # Analysis 5
│   ├── 08_MetaAnalysis_AQS_only.R    # Analysis 6
│   └── Utilities/
│       ├── forest.col.r              # Forest plot functions (based on Opie2019)
│       ├── robu.custom.r             # Custom RVE (based on Opie2019)
│       ├── RVE.r                     # RVE wrappers (based on Opie2019)
│       ├── meta_analysis_helpers.r
│       └── effect_size_utils.r       # New: r and kappa calculations
├── Data/
├── Output/
│   ├── 2way_SIS/                     # Results from analysis 1
│   ├── 2way_OD/                      # Results from analysis 2
│   ├── 3way/                         # Results from analysis 3
│   ├── 4way/                         # Results from analysis 4
│   ├── SSP_only/                     # Results from analysis 5
│   └── AQS_only/                     # Results from analysis 6
├── Figures/                          # Forest plots, funnel plots
└── README.md

```

## Running the Analyses

### Run individual steps:
```r
source("Code/01_PrepareData.R")              # Step 1: Clean data
source("Code/02_CalculateEffectSizes.R")     # Step 2: Calculate effect sizes
source("Code/03_MetaAnalysis_2Way_SIS_r.R")  # Step 3: Run specific analysis
# ... etc
```

## Required R Packages

```r
install.packages(c(
  "readxl",       # Read Excel files
  "dplyr",        # Data manipulation
  "tidyr",        # Data reshaping
  "robumeta",     # RVE meta-analysis
  "metafor",      # Additional meta-analysis functions
  "psych"         # Kappa calculations
))
```

## Analysis Methods

### Effect Size Metrics
- **Pearson's r**: Calculated from 2×2 contingency tables using phi coefficient
- **Cohen's kappa**: Calculated from k×k contingency tables, accounting for chance agreement

### Meta-Analytic Method
- **Robust Variance Estimation (RVE)**: Accounts for dependencies among effect sizes from the same study
- **Random-effects models**: Assumes heterogeneity across studies
- **Small-sample corrections**

## Output Files

For each analysis:
- **Forest plot** (TIFF, 300 DPI)
- **Funnel plot** (TIFF, 300 DPI)
- **Summary table** (CSV and TXT)
- **Full results object** (RData)
- **Moderator results** (CSV)
