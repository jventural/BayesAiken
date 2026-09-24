# BayesAiken

Bayesian inference for the six rating coefficients of Aiken (1985, 1989):
**V** (content validity), **H** (homogeneity), **C** (congruence),
**R** (repeatability), **A** (rater agreement) and **I** (item coherence).
All six are functionals of the population distribution of rating
categories, so a single Dirichlet-Multinomial model gives their full
posterior. The package also provides the classical coefficients, the
Penfield-Giacobbi score confidence interval for V and the significance tests
of Aiken (1985).

## Installation

```r
# install.packages("remotes")
remotes::install_github("jventural/BayesAiken")
```

## Example

```r
library(BayesAiken)

# Ten judges rate one item on a 1-4 scale
ratings <- c(4, 4, 3, 4, 4, 3, 4, 4, 4, 3)

# Bayesian V with its credibility interval
aiken_bayes(ratings, coef = "V", l = 1, s = 4)

# Classical V, Penfield-Giacobbi interval and Aiken's (1985) test
k <- 3
aiken_V_classic(ratings - 1, k)
aiken_V_ic_pg(aiken_V_classic(ratings - 1, k), n = 10, k = k)
aiken_test(ratings, coef = "V", l = 1, s = 4)

# Interactive application
if (interactive()) run_aiken_app()
```

Version 0.1.0 (Bayesian V only, `V_aiken()`) is kept in the tag
[`v3-0.1.0`](https://github.com/jventural/BayesAiken/tree/v3-0.1.0).
