## Submission summary

New submission of BayesAiken (version 2.0.0): Bayesian inference for the six
rating coefficients of Aiken (1985, 1989) with a Dirichlet-Multinomial model,
the classical coefficients, the Penfield-Giacobbi score interval for V and
the significance tests of Aiken (1985).

## Test environments

* Local: Windows 11 x64, R 4.4.1 (R CMD check --as-cran --run-donttest,
  PDF manual built with pdflatex)
* win-builder: R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for the reviewer

* The Shiny application (`run_aiken_app()`) is in inst/shiny; shiny is in
  Suggests and the example only runs inside if (interactive()).
* Simulation-based functions do not set seeds; examples use small numbers of
  draws and the long simulation example is wrapped in \donttest{}.

## Downstream dependencies

There are currently no downstream dependencies.
