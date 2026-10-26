# Antarctic Penguin Population Trend Analysis

![CCAMLR](https://img.shields.io/badge/CCAMLR-State%20of%20the%20Antarctic%20Environment-teal) ![R](https://img.shields.io/badge/language-R-blue)

R script for monitoring and visualising population trends in Antarctic pygoscelid penguins using count data from [MAPPPD](https://www.penguinmap.com/mapppd/). Developed to support State of the Antarctic Environment (SOAE) reporting within the [CCAMLR](https://www.ccamlr.org) framework, where penguin trends serve as indicators of Southern Ocean ecosystem health.

## Species and regions

| Species | CCAMLR regions |
|---|---|
| Adélie, Chinstrap, Gentoo penguin | 48.1 (Antarctic Peninsula), 58.4.2 (East Antarctica), 88.1 (Ross Sea) |

## Statistical approach

Trend analysis compares a **recent period** (last 3 years) against a **historical period** (3–10 years prior to the latest survey) using a negative binomial mixed model (`glmer.nb`) with site as a random intercept. Marginal means are estimated via `emmeans`. Results significant at p < 0.05 are flagged with `**` in output figures.

## Requirements

```r
install.packages(c("MASS", "ggplot2", "tidyverse", "mgcv", "emmeans", "lme4"))
```

On first run the script downloads `MAPPPD_AllCounts.csv` automatically and reuses it on subsequent runs. Only sites with ≥5 years of observations are included.

## Usage and outputs

```bash
Rscript penguin_trend_analysis.R
```

Produces one SVG figure per CCAMLR region (e.g. `penguins_recent_change_48.1.svg`), faceted by species and count type, showing model-adjusted mean counts with 95% CIs and annotated percentage change.

## Data source

Humphries et al. (2017). MAPPPD: data and tools for dynamic management and decision support. *Polar Record*, 53(2), 160–166. https://doi.org/10.1017/S0032247417000055

---

*Data accessed via penguinmap.com. Maintained as part of the CCAMLR SAE reporting workflow.*
