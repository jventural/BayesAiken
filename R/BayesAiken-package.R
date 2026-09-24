#' BayesAiken: Bayesian Content Validity Coefficients
#'
#' Bayesian inference for the six rating coefficients of Aiken (V, H, C, R,
#' A and I) with a Dirichlet-Multinomial model, together with the classical
#' coefficients, the Penfield-Giacobbi confidence interval for V and the
#' significance tests of Aiken (1985).
#'
#' @keywords internal
#' @importFrom grDevices adjustcolor rgb
#' @importFrom graphics abline axis barplot grid hist legend lines mtext par
#'   points polygon rect segments text
#' @importFrom stats aggregate approx dbeta density median optimize pbeta qbeta
#'   quantile rbeta rbinom rgamma sd var
#' @importFrom utils capture.output write.csv
"_PACKAGE"
