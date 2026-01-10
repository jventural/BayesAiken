# BayesAiken <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/jventural/BayesAiken/workflows/R-CMD-check/badge.svg)](https://github.com/jventural/BayesAiken/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/BayesAiken)](https://CRAN.R-project.org/package=BayesAiken)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**BayesAiken** provides Bayesian estimation of Aiken's V coefficient for content validity studies, implementing both the standard Beta-Binomial model and a conservative Dirichlet-Multinomial alternative designed for small expert panels.

## Why Dirichlet-Multinomial?

The classical Beta-Binomial approach treats each rating point as an independent Bernoulli trial, effectively inflating the sample size. With a scale of 0-3 and 5 judges, this creates 15 pseudo-observations from just 5 judges.

The **Dirichlet-Multinomial** model counts each judge exactly once, producing:
- **Wider, more honest credible intervals** when panels are small (n < 10)
- **Conservative uncertainty quantification** that better reflects sparse data
- **Direct probability statements** about validity thresholds

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("jventural/BayesAiken")
```

## Quick Start

```r
library(BayesAiken)

# Example: 5 judges rate an item on a 0-3 scale
# Ratings: one judge chose 2, four chose 3
ratings <- c(2, 3, 3, 3, 3)

# Dirichlet-Multinomial model (recommended for small panels)
result <- V_aiken(ratings, l = 0, s = 3, model = "dirichlet")

# =======================================================
#          BAYESIAN AIKEN'S V - Dirichlet-Multinomial
# =======================================================
#
# DATA SUMMARY:
#   Number of judges:     5
#   Scale range:          [0, 3]
#   Category frequencies: 0, 0, 1, 4
#
# ESTIMATES:
#   Classical V:          0.9333
#   Posterior mean:       0.8519
#   Posterior median:     0.8704
#   Posterior SD:         0.1082
#
# CREDIBLE INTERVALS (95%):
#   HDI: [0.6179, 0.9951]  (width = 0.3773)
#   ETI: [0.5878, 0.9855]  (width = 0.3977)
#
# DECISION SUPPORT:
#   P(V >= 0.70 | data) = 0.8937
#   Interpretation:       Moderate evidence for adequate validity
```

## Compare Models

```r
# Compare Beta-Binomial vs Dirichlet-Multinomial
plot_V_comparison(ratings, l = 0, s = 3)
```

This produces a publication-ready figure showing how the Dirichlet model yields wider intervals that better reflect uncertainty with small panels.

## Key Functions

| Function | Description |
|----------|-------------|
| `V_aiken()` | Main function for Bayesian V estimation |
| `plot_V_comparison()` | Compare posterior densities from both models |
| `plot_V_forest()` | Forest plot for multiple items |

## Input Formats

```r
# Individual ratings
V_aiken(c(2, 3, 3, 3, 3), l = 0, s = 3)

# Frequency counts: c(n0, n1, n2, n3)
V_aiken(c(0, 0, 1, 4), l = 0, s = 3, input_type = "counts")
```

## Citation

If you use this package, please cite:

```
Ventura-León, J. (2026). *BayesAiken* [Software]. GitHub. https://github.com/jventural/BayesAiken
```

## References

- Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings. *Educational and Psychological Measurement*, 45(1), 131-142.
- Kruschke, J. K. (2018). Rejecting or accepting parameter values in Bayesian estimation. *Advances in Methods and Practices in Psychological Science*, 1(2), 270-280.

## Author

**Jose Ventura-León** - [ORCID](https://orcid.org/0000-0003-2996-4244)

## License

MIT License
