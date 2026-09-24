library(testthat)
library(BayesAiken)

# =============================================================================
# 1. CLASICOS: fidelidad a Aiken (1985, 1989)
# =============================================================================

test_that("V clasica reproduce Aiken (1985, p.133)", {
  # V = S/[n(c-1)], S = sum(r - lo)
  expect_equal(aiken_V_classic(c(3, 3, 3, 3), 3), 1)
  expect_equal(aiken_V_classic(c(0, 0, 0, 0), 3), 0)
  expect_equal(aiken_V_classic(c(2, 3, 3, 2, 3), 3), 13 / 15)
})

test_that("H clasica respeta la paridad j de Aiken (1985, p.140)", {
  # n par -> j = 0 ; n impar -> j = 1
  expect_equal(aiken_H_classic(c(3, 3, 3, 3, 3), 3), 1)      # acuerdo perfecto
  expect_equal(aiken_H_classic(c(0, 0, 3, 3), 3), 0)         # desacuerdo maximo, n par
  # n = 4, pares con el "2": 3 pares de diferencia 1 -> S = 3, j = 0
  expect_equal(aiken_H_classic(c(3, 3, 2, 3), 3), 1 - 4 * 3 / (3 * 16))
  # n = 5, S = 4: 1 - 16/(3*(25-1))
  expect_equal(aiken_H_classic(c(3, 3, 2, 3, 3), 3), 1 - 16 / (3 * 24))
})

test_that("C clasica respeta la paridad d de Aiken (1989, formula 3)", {
  expect_equal(aiken_C_classic(c(2, 2, 2, 2), 3), 1)         # varianza nula
  expect_equal(aiken_C_classic(c(0, 0, 3, 3), 3), 0)         # desacuerdo maximo, n par
  n <- 4; d <- 0
  x <- c(3, 3, 2, 3)
  expect_equal(aiken_C_classic(x, 3),
               1 - (4 * n * (n - 1) * var(x)) / ((n^2 - d) * 9))
})

test_that("R, A e I clasicos reproducen Aiken (1985)", {
  expect_equal(aiken_R_classic(c(3, 2, 3), c(3, 2, 3), 3), 1)
  expect_equal(aiken_R_classic(c(3, 3, 3), c(0, 0, 0), 3), 0)
  expect_equal(aiken_A_classic(3, c(3, 3, 3), 3), 1)
  expect_equal(aiken_A_classic(3, c(2, 3, 3), 3), 1 - 1 / 9)
  expect_equal(aiken_I_classic(3, c(2, 3, 3), 3), 1 - 1 / 9)
})

# =============================================================================
# 2. FUNCIONALES: coherencia entre coeficiente muestral y funcional poblacional
# =============================================================================

test_that("los funcionales alcanzan sus extremos teoricos", {
  k <- 3
  P_todo_max <- matrix(c(0, 0, 0, 1), nrow = 1)
  P_mitades  <- matrix(c(.5, 0, 0, .5), nrow = 1)
  expect_equal(BayesAiken:::.ab_T(P_todo_max, "V", k), 1)
  expect_equal(BayesAiken:::.ab_T(P_todo_max, "H", k), 1)
  expect_equal(BayesAiken:::.ab_T(P_todo_max, "C", k), 1)
  # desacuerdo maximo: mitad en cada extremo -> H = C = 0
  expect_equal(BayesAiken:::.ab_T(P_mitades, "H", k), 0)
  expect_equal(BayesAiken:::.ab_T(P_mitades, "C", k), 0)
})

test_that("el coeficiente muestral converge al funcional poblacional", {
  set.seed(99)
  k <- 3; pi <- c(.05, .10, .35, .50)
  r <- sample(0:k, 200000, TRUE, pi)
  P <- matrix(pi, nrow = 1)
  expect_equal(aiken_V_classic(r, k), BayesAiken:::.ab_T(P, "V", k), tolerance = 0.01)
  expect_equal(aiken_H_classic(r, k), BayesAiken:::.ab_T(P, "H", k), tolerance = 0.01)
  expect_equal(aiken_C_classic(r, k), BayesAiken:::.ab_T(P, "C", k), tolerance = 0.01)
})

test_that("R no es funcional de las marginales", {
  k <- 3
  PiA <- diag(rep(.25, 4))                                    # acuerdo perfecto
  PiB <- matrix(0, 4, 4); PiB[1,4] <- PiB[4,1] <- PiB[2,3] <- PiB[3,2] <- .25
  expect_equal(rowSums(PiA), rowSums(PiB))
  expect_equal(colSums(PiA), colSums(PiB))
  rA <- BayesAiken:::.ab_T(matrix(as.vector(PiA), nrow = 1), "R", k)
  rB <- BayesAiken:::.ab_T(matrix(as.vector(PiB), nrow = 1), "R", k)
  expect_equal(rA, 1)
  expect_true(abs(rA - rB) > 0.5)
})

# =============================================================================
# 3. API: aiken_bayes()
# =============================================================================

test_that("aiken_bayes devuelve una posterior coherente para los seis coeficientes", {
  r <- c(3, 3, 2, 3, 3)
  for (cf in c("V", "H", "C")) {
    z <- aiken_bayes(r, coef = cf, seed = 1, B = 4000, verbose = FALSE)
    expect_s3_class(z, "aiken_bayes")
    expect_true(z$media >= -1 && z$media <= 1)
    expect_true(z$ci[1] <= z$media && z$media <= z$ci[2])
    expect_length(z$draws, 4000)
  }
  zr <- aiken_bayes(r, y = c(3, 2, 2, 3, 3), coef = "R", seed = 1, B = 4000, verbose = FALSE)
  expect_equal(length(zr$conteos), 16)          # tabla 4x4
  za <- aiken_bayes(c(3, 3, 2, 3), focal = 3, coef = "A", seed = 1, B = 4000, verbose = FALSE)
  expect_true(za$media > 0.5)
  zi <- aiken_bayes(c(3, 3, 2), focal = 3, coef = "I", seed = 1, B = 4000, verbose = FALSE)
  expect_s3_class(zi, "aiken_bayes")
})

test_that("aiken_bayes reproduce el clasico correcto", {
  r <- c(3, 3, 2, 3, 3)
  expect_equal(aiken_bayes(r, coef = "V", seed = 1, B = 500, verbose = FALSE)$clasico,
               aiken_V_classic(r, 3))
  expect_equal(aiken_bayes(r, coef = "H", seed = 1, B = 500, verbose = FALSE)$clasico,
               aiken_H_classic(r, 3))
})

test_that("con mucha evidencia la posterior converge al valor poblacional", {
  set.seed(7)
  r <- sample(0:3, 4000, TRUE, c(.05, .10, .35, .50))
  z <- aiken_bayes(r, coef = "V", seed = 2, B = 4000, verbose = FALSE)
  expect_equal(z$media, aiken_V_classic(r, 3), tolerance = 0.01)
  expect_true(diff(z$ci) < 0.05)                # el intervalo se estrecha
})

test_that("el prior desplaza el resultado con paneles pequenos", {
  r <- c(3, 3, 2, 3, 3)
  vp <- aiken_bayes(r, coef = "V", prior = "perks",   seed = 1, B = 20000, verbose = FALSE)$media
  vu <- aiken_bayes(r, coef = "V", prior = "uniform", seed = 1, B = 20000, verbose = FALSE)$media
  expect_true(vp > vu)                          # el uniforme encoge hacia 0.5
  expect_true(vp - vu > 0.05)
})

test_that("los errores de escala son explicitos", {
  expect_error(aiken_bayes(c(1, 4, 4), l = 0, s = 3, verbose = FALSE), "fuera de la escala")
  expect_error(aiken_bayes(c(1, 2, 2), l = 3, s = 0, verbose = FALSE), "mayor que")
  expect_error(aiken_bayes(c(1, 2.5), l = 0, s = 3, verbose = FALSE), "no enteros")
  expect_error(aiken_bayes(c(1, 2), coef = "R", verbose = FALSE), "requiere 'y'")
  expect_error(aiken_bayes(c(1, 2), coef = "A", verbose = FALSE), "requieren 'focal'")
})

test_that("la escala 1-4 se maneja bien (el fallo clasico de AikenCalc)", {
  z <- aiken_bayes(c(4, 4, 4, 4, 4), l = 1, s = 4, coef = "V", seed = 1,
                   B = 2000, verbose = FALSE)
  expect_equal(z$clasico, 1)                    # y NO 1.33
  expect_true(z$media <= 1)
})

# =============================================================================
# 4. API: aiken_bayes_table()
# =============================================================================

test_that("aiken_bayes_table procesa un formato Item x Criterio x Jueces", {
  d <- data.frame(
    Item = paste0("It", rep(1:3, each = 2)),
    Criterio = rep(c("Relevancia", "Claridad"), 3),
    J1 = c(3,3,2,3,3,3), J2 = c(3,2,3,3,2,3),
    J3 = c(2,3,3,2,3,3), J4 = c(3,3,3,3,3,2)
  )
  res <- aiken_bayes_table(d, seed = 1, B = 2000, verbose = FALSE)
  expect_s3_class(res, "aiken_table")
  expect_equal(nrow(res$global), 6)
  # nombres cortos y estables: <coef>, <coef>_clasico, _li, _ls, _pNN
  expect_true(all(c("V_clasico", "V", "V_li", "V_ls", "V_p7", "V_p8",
                    "H", "C") %in% names(res$global)))
  expect_equal(res$jueces, c("J1", "J2", "J3", "J4"))
  expect_equal(nrow(res$acuerdo), 6 * 4)        # A por fila y juez
  expect_true(!is.null(res$coherencia))
  expect_true(all(res$global$V >= 0 & res$global$V <= 1))
})

test_that("fmt_apa sigue APA 7: sin cero a la izquierda", {
  expect_equal(fmt_apa(0.75), ".75")
  expect_equal(fmt_apa(1), "1.00")
  expect_equal(fmt_apa(-0.04), "-.04")
  expect_equal(fmt_apa(NA_real_), "—")
})

test_that("las tablas locales usan la MISMA nomenclatura que la global", {
  # Un renombrado a medias dejo la app pidiendo A_bayes/I_bayes y reventaba las
  # dos pestanas locales. Este test congela los nombres de las tres tablas.
  d <- data.frame(
    Item = paste0("It", rep(1:3, each = 2)),
    Criterio = rep(c("Rel", "Cla"), 3),
    J1 = c(3,3,2,3,3,3), J2 = c(3,2,3,3,2,3), J3 = c(2,3,3,2,3,3), J4 = c(3,3,3,3,3,2))
  res <- aiken_bayes_table(d, seed = 1, B = 800, verbose = FALSE)
  expect_equal(names(res$acuerdo),    c("Unidad", "Juez", "A_clasico", "A", "A_li", "A_ls"))
  expect_equal(names(res$coherencia), c("Item", "Juez", "I_clasico", "I", "I_li", "I_ls"))
  # el patron <coef>, <coef>_clasico, _li, _ls debe valer para las tres
  for (cf in c("V", "H", "C"))
    expect_true(all(c(cf, paste0(cf, c("_clasico", "_li", "_ls"))) %in% names(res$global)))
  # y la agregacion que hace la app debe funcionar sin errores
  expect_silent(aggregate(cbind(A_clasico, A, A_li, A_ls) ~ Juez,
                          data = res$acuerdo, FUN = mean))
  expect_silent(aggregate(cbind(I_clasico, I, I_li, I_ls) ~ Item,
                          data = res$coherencia, FUN = mean))
})

test_that("el IC de Penfield-Giacobbi reproduce el intervalo de puntuacion", {
  # Wilson sobre n*k ensayos: acotado en [0,1] y centrado cerca de V
  ic <- aiken_V_ic_pg(0.92, n = 5, k = 3)
  expect_true(ic[1] >= 0 && ic[2] <= 1)
  expect_true(ic[1] < 0.92 && 0.92 < ic[2])
  # con V = 1 el limite superior toca el techo
  expect_equal(aiken_V_ic_pg(1, n = 5, k = 3)[2], 1, tolerance = 1e-8)
  # mas jueces, intervalo mas estrecho
  expect_true(diff(aiken_V_ic_pg(0.9, 20, 3)) < diff(aiken_V_ic_pg(0.9, 5, 3)))
})

test_that("aiken_apa anade el IC de PG solo si se pide y solo a la V", {
  d <- data.frame(Item = c("I1", "I2"), J1 = c(3, 2), J2 = c(3, 3),
                  J3 = c(2, 3), J4 = c(3, 2))
  res <- aiken_bayes_table(d, seed = 1, B = 800, local = FALSE, verbose = FALSE)
  expect_true(all(c("V_pg_li", "V_pg_ls") %in% names(res$global)))
  sin <- aiken_apa(res)
  con <- aiken_apa(res, pg = TRUE)
  expect_false(any(grepl("IC PG", sin$V)))
  expect_true(all(grepl("IC PG", con$V)))
  expect_false(any(grepl("IC PG", con$H)))   # solo la V
})

test_that("aiken_apa colapsa a una columna por coeficiente", {
  d <- data.frame(Item = c("I1", "I2"), Criterio = c("Rel", "Rel"),
                  J1 = c(3, 2), J2 = c(3, 3), J3 = c(2, 3), J4 = c(3, 2))
  res <- aiken_bayes_table(d, seed = 1, B = 800, local = FALSE, verbose = FALSE)
  ap <- aiken_apa(res)
  expect_equal(names(ap), c("Item", "Criterio", "V", "H", "C"))
  expect_equal(nrow(ap), 2)
  expect_match(ap$V[1], "^\\.[0-9]{2} \\[")            # ".86 ["
  expect_match(ap$V[1], "P = ")
  expect_match(ap$V[1], "clásico")
  expect_equal(attr(ap, "umbral"), 0.70)
  # variante HTML para la app
  aph <- aiken_apa(res, html = TRUE)
  expect_match(aph$V[1], "ab-val")
  expect_match(aph$V[1], "<br>")
})

test_that("aiken_bayes_table detecta calificaciones fuera de escala", {
  d <- data.frame(Item = c("a", "b"), J1 = c(3, 9), J2 = c(2, 3))
  expect_error(aiken_bayes_table(d, l = 0, s = 3, verbose = FALSE), "fuera de la escala")
})
