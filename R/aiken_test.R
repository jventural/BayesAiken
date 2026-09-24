#' Prueba de significacion de Aiken (1985) para V, R y H
#'
#' Contrasta si un coeficiente de Aiken es mayor que el esperado cuando los
#' jueces responden al azar, es decir, con calificaciones independientes y
#' uniformes sobre las \code{c} categorias de la escala. Es la prueba que
#' propone Aiken (1985): una probabilidad de cola derecha y, para muestras
#' grandes, una aproximacion normal. Es frecuentista y complementa, no
#' reemplaza, al intervalo de credibilidad de \code{\link{aiken_bayes}}.
#'
#' @param x Vector de calificaciones en la escala \code{l..s}. Para
#'   \code{coef = "R"}, calificaciones de la primera ocasion.
#' @param coef Coeficiente: \code{"V"} (validez de contenido), \code{"R"}
#'   (repetibilidad) o \code{"H"} (homogeneidad).
#' @param y Solo para \code{coef = "R"}: calificaciones de la segunda ocasion,
#'   en el mismo orden que \code{x}.
#' @param l,s Minimo y maximo de la escala (default 0 y 3).
#' @param B Numero de muestras simuladas para la probabilidad de H (default
#'   100000).
#'
#' @details
#' \describe{
#'   \item{V}{\eqn{V = S / [n(c - 1)]}, con \eqn{S} la suma de las
#'     calificaciones desplazadas a \code{0..k}. Bajo la hipotesis nula
#'     \eqn{S} es la suma de \eqn{n} uniformes discretas, cuya distribucion se
#'     obtiene exactamente por convolucion. La media nula es .5 y el desvio
#'     \eqn{\sqrt{(c + 1) / [12 n (c - 1)]}}.}
#'   \item{R}{\eqn{R = 1 - \sum |x_i - y_i| / [m(c - 1)]}. Bajo la hipotesis
#'     nula cada diferencia absoluta tiene
#'     \eqn{P(D = 0) = 1/c} y \eqn{P(D = d) = 2(c - d)/c^2}, y la distribucion
#'     de su suma tambien se obtiene exactamente. La media nula es
#'     \eqn{(2c - 1)/(3c)}.}
#'   \item{H}{Se calcula con \code{\link{aiken_H_classic}}. Su distribucion
#'     nula se estima por simulacion con \code{B} muestras, y la probabilidad
#'     lleva la correccion \eqn{(1 + b)/(B + 1)}.}
#' }
#' La columna \code{z} usa la media y el desvio de la distribucion nula; con
#' mas de 25 jueces o items coincide con la prueba de muestras grandes de
#' Aiken (1985). Con muestras pequenas conviene reportar la probabilidad
#' exacta (\code{p}).
#'
#' @return Un data frame de una fila con \code{coef}, el valor observado
#'   (\code{valor}), el numero de jueces o pares (\code{n}), la probabilidad de
#'   cola derecha (\code{p}), el metodo con que se obtuvo (\code{metodo}:
#'   \code{"exacta"} o \code{"simulacion"}), la media y el desvio nulos
#'   (\code{media_nula}, \code{de_nula}), el estadistico \code{z} y su
#'   probabilidad normal de cola derecha (\code{p_z}).
#'
#' @references
#' Aiken, L. R. (1985). Three coefficients for analyzing the reliability and
#' validity of ratings. \emph{Educational and Psychological Measurement,
#' 45}(1), 131-142. \doi{10.1177/0013164485451012}
#'
#' @seealso \code{\link{aiken_V_ic_pg}} para el intervalo de confianza de la
#'   V (Penfield y Giacobbi, 2004).
#'
#' @examples
#' # Cinco jueces califican un item en una escala de 0 a 3
#' aiken_test(c(3, 3, 2, 3, 3), coef = "V", l = 0, s = 3)
#'
#' # Repetibilidad: mismos 8 items calificados en dos ocasiones (escala 1-5)
#' aiken_test(c(5, 4, 4, 5, 3, 4, 5, 5), coef = "R",
#'            y = c(5, 4, 3, 5, 3, 4, 4, 5), l = 1, s = 5)
#'
#' # Homogeneidad (probabilidad por simulacion)
#' aiken_test(c(2, 3, 3, 3, 2, 3), coef = "H", l = 0, s = 3, B = 5000)
#'
#' @export
aiken_test <- function(x, coef = c("V", "R", "H"), y = NULL, l = 0, s = 3,
                       B = 100000) {
  coef <- match.arg(coef)
  sc <- .ab_scale(l, s)
  k <- sc$k
  c_cat <- sc$c

  conv_power <- function(pmf, times) {
    out <- 1
    for (i in seq_len(times)) out <- stats::convolve(out, rev(pmf), type = "open")
    pmax(out, 0)
  }

  if (coef == "V") {
    xi <- .ab_check_ratings(x, sc, "x")
    n <- length(xi)
    valor <- aiken_V_classic(xi, k)
    dist_S <- conv_power(rep(1 / c_cat, c_cat), n)   # S = 0..n*k
    S_obs <- sum(xi)
    p <- sum(dist_S[(S_obs + 1):length(dist_S)])
    media <- 0.5
    de <- sqrt((c_cat + 1) / (12 * n * (c_cat - 1)))
    metodo <- "exacta"
  } else if (coef == "R") {
    if (is.null(y)) stop("'y' es obligatorio para coef = \"R\".", call. = FALSE)
    if (length(x) != length(y))
      stop("'x' e 'y' deben tener la misma longitud.", call. = FALSE)
    ok <- !is.na(x) & !is.na(y)
    xi <- .ab_check_ratings(x[ok], sc, "x")
    yi <- .ab_check_ratings(y[ok], sc, "y")
    n <- length(xi)
    valor <- aiken_R_classic(xi, yi, k)
    pmf_D <- c(1 / c_cat, 2 * (c_cat - seq_len(k)) / c_cat^2)   # D = 0..k
    dist_S <- conv_power(pmf_D, n)
    S_obs <- sum(abs(xi - yi))
    p <- sum(dist_S[1:(S_obs + 1)])          # R alto <=> suma de diferencias baja
    ED <- sum(0:k * pmf_D)
    VD <- sum((0:k)^2 * pmf_D) - ED^2
    media <- 1 - ED / k
    de <- sqrt(VD / n) / k
    metodo <- "exacta"
  } else {
    xi <- .ab_check_ratings(x, sc, "x")
    n <- length(xi)
    valor <- aiken_H_classic(xi, k)
    sims <- vapply(seq_len(B), function(b) {
      aiken_H_classic(sample.int(c_cat, n, replace = TRUE) - 1L, k)
    }, numeric(1))
    p <- (1 + sum(sims >= valor - 1e-12)) / (B + 1)
    media <- mean(sims)
    de <- stats::sd(sims)
    metodo <- "simulacion"
  }

  z <- (valor - media) / de
  data.frame(coef = coef, valor = valor, n = n, p = min(1, p), metodo = metodo,
             media_nula = media, de_nula = de, z = z,
             p_z = stats::pnorm(z, lower.tail = FALSE),
             stringsAsFactors = FALSE)
}
