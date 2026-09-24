test_that("aiken_test() V matches the exact enumeration", {
  r <- aiken_test(c(3, 3, 2, 3, 3), coef = "V", l = 0, s = 3)
  expect_equal(r$p, 0.005859375)
  expect_equal(r$valor, 14 / 15)
  expect_equal(aiken_test(c(4, 3, 4, 2), coef = "V", l = 1, s = 5)$p, 0.432)
})

test_that("aiken_test() V null SD equals Aiken (1985)", {
  r <- aiken_test(c(3, 2, 3, 3, 2, 3, 3, 3, 1, 3), coef = "V", l = 0, s = 3)
  expect_equal(r$de_nula, 0.5 * sqrt(5 / (3 * 10 * 3)))
})

test_that("aiken_test() R matches the exact enumeration and Aiken's mean", {
  r <- aiken_test(c(2, 1, 2, 0), coef = "R", y = c(2, 1, 1, 0), l = 0, s = 2)
  expect_equal(r$p, 0.0781893, tolerance = 1e-6)
  r2 <- aiken_test(c(1, 2, 3), coef = "R", y = c(1, 2, 2), l = 1, s = 5)
  expect_equal(r2$media_nula, (2 * 5 - 1) / (3 * 5))
})

test_that("aiken_test() H simulation is close to the exact probability", {
  set.seed(1)
  r <- aiken_test(c(2, 3, 3, 3, 2, 3), coef = "H", l = 0, s = 3, B = 20000)
  expect_equal(r$p, 0.03173828, tolerance = 0.005)
  expect_identical(r$metodo, "simulacion")
})

test_that("aiken_test() validates its input", {
  expect_error(aiken_test(c(1, 2), coef = "R", l = 0, s = 3), "'y'")
  expect_error(aiken_test(c(1, 5), coef = "V", l = 0, s = 3), "fuera de la escala")
})
