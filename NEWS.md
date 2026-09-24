# BayesAiken 2.0.0

First version submitted to CRAN.

## New features

* All six coefficients of Aiken (V, H, C, R, A, I) are now functionals of a
  single Dirichlet-Multinomial model (`aiken_bayes()`, `aiken_bayes_table()`),
  which treats each judge as one categorical observation and avoids the
  pseudo-replication of the Beta-Binomial formulation. R is modelled on the
  c-by-c joint table.
* `aiken_test()`: the significance tests of Aiken (1985) for V, R and H, with
  the exact right-tail probability (V and R, by convolution of the null
  distribution; H, by simulation) and the large-sample z test.
* Classical coefficients with the correct parity correction
  (`aiken_V_classic()`, `aiken_H_classic()`, ...) and the Penfield-Giacobbi
  score confidence interval for V (`aiken_V_ic_pg()`).
* APA-style reporting (`aiken_apa()`, `print_apa()`) and a local Shiny
  application (`run_aiken_app()`).

## Bug fixes (with respect to 1.0.0)

* The parity correction of H and C was fixed at 1, wrong with an even number
  of judges.
* The denominator of A omitted the number of criteria, and I used m = 2
  regardless of the design.
* `coef_C()` had no likelihood: it built pseudo-counts from its own point
  estimate.
* The five coefficients other than V used a Beta-Binomial model; for H this
  gave intervals about 3.4 times narrower than they should be.
