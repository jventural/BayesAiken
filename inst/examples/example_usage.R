# =============================================================================
# BayesAiken - Ejemplos de Uso
# =============================================================================

library(BayesAiken)

# =============================================================================
# EJEMPLO 1: Vector simple - Coeficiente V de Aiken
# =============================================================================

cat("\n=== EJEMPLO 1: Coeficiente V (vector) ===\n")

# 10 jueces calificaron un item en escala 0-3
ratings <- c(3, 3, 2, 3, 3, 2, 3, 3, 2, 3)

# Calcular con la funcion principal
result_v <- bayes_aiken(ratings, coefficient = "V", l = 0, s = 3)

# Ver resultado
print(result_v)

# Grafico
plot(result_v)

# =============================================================================
# EJEMPLO 2: Multiples coeficientes
# =============================================================================

cat("\n=== EJEMPLO 2: Multiples coeficientes ===\n")

ratings <- c(2, 2, 3, 2, 2, 3, 2, 2, 3, 2)

# Calcular V, H y C simultaneamente
result_multi <- bayes_aiken(ratings, coefficient = c("V", "H", "C"),
                            l = 0, s = 3, c = 4)

# Tabla resumen
summary_table(result_multi)

# =============================================================================
# EJEMPLO 3: Coeficiente R (Reproducibilidad test-retest)
# =============================================================================

cat("\n=== EJEMPLO 3: Coeficiente R (Reproducibilidad) ===\n")

# Calificaciones en dos momentos
tiempo_1 <- c(3, 3, 3, 3, 3, 2, 3, 3, 2, 3)
tiempo_2 <- c(3, 3, 3, 2, 3, 2, 3, 3, 3, 3)

# Usando funcion modular
result_r <- coef_R(tiempo_1, tiempo_2, c = 4)
print(result_r)

# =============================================================================
# EJEMPLO 4: Coeficiente A (Acuerdo del juez)
# =============================================================================

cat("\n=== EJEMPLO 4: Coeficiente A (Acuerdo del juez) ===\n")

# Calificaciones del juez focal
juez_focal <- c(3, 3, 2, 3)

# Calificaciones de otros jueces (matriz: filas=criterios, cols=otros jueces)
otros_jueces <- matrix(c(
  3, 3, 2, 3,  # Juez 2
  3, 2, 3, 3,  # Juez 3
  2, 3, 2, 3,  # Juez 4
  3, 3, 3, 2   # Juez 5
), nrow = 4, byrow = FALSE)

result_a <- coef_A(juez_focal, otros_jueces, c = 4)
print(result_a)

# =============================================================================
# EJEMPLO 5: Coeficiente I (Consistencia inter-item)
# =============================================================================

cat("\n=== EJEMPLO 5: Coeficiente I (Consistencia inter-item) ===\n")

# Dos items calificados por los mismos jueces
item_1 <- c(3, 3, 2, 3, 3)
item_2 <- c(3, 2, 3, 3, 3)

result_i <- coef_I(item_1, item_2, c = 4)
print(result_i)

# =============================================================================
# EJEMPLO 6: Data frame con multiples items
# =============================================================================

cat("\n=== EJEMPLO 6: Data frame con multiples items ===\n")

# Crear datos de ejemplo
data_items <- data.frame(
  item = rep(1:5, each = 3),
  criterio = rep(c("relevancia", "representatividad", "claridad"), 5),
  J1 = c(3,3,3, 2,3,3, 3,2,3, 3,3,2, 2,3,3),
  J2 = c(3,3,2, 3,3,3, 3,3,2, 3,2,3, 3,3,2),
  J3 = c(2,3,3, 3,2,3, 2,3,3, 2,3,3, 3,2,3),
  J4 = c(3,3,3, 3,3,2, 3,3,3, 3,3,3, 3,3,3),
  J5 = c(3,2,3, 2,3,3, 3,2,3, 3,3,2, 2,3,3)
)

# Calcular V para todos los items
result_df <- bayes_aiken(data_items, coefficient = "V",
                         item_col = "item", criterion_col = "criterio",
                         l = 0, s = 3)

# Ver resumen
print(result_df)

# Grafico de comparacion
plot(result_df)

# =============================================================================
# EJEMPLO 7: Comparar diferentes priors
# =============================================================================

cat("\n=== EJEMPLO 7: Comparacion de priors ===\n")

ratings <- c(3, 3, 2, 3, 3, 2, 3, 3, 2, 3)

prior_comparison <- compare_priors(ratings, l = 0, s = 3)
print(prior_comparison)

# =============================================================================
# EJEMPLO 8: Captura y exportacion
# =============================================================================

cat("\n=== EJEMPLO 8: Captura y exportacion ===\n")

# Capturar salida verbose
captured <- bayes_aiken_capture(ratings, coefficient = "V", l = 0, s = 3)
print(captured)

# Exportar a archivo (descomentar para usar)
# bayes_aiken_export(result_v, "resultados_V.txt", format = "both")

# =============================================================================
# EJEMPLO 9: Usando diferentes priors
# =============================================================================

cat("\n=== EJEMPLO 9: Priors informativos ===\n")

# Prior uniforme (default)
result_uniforme <- coef_V(ratings, l = 0, s = 3,
                          prior_alpha = 1, prior_beta = 1,
                          verbose = FALSE)

# Prior de Jeffreys
result_jeffreys <- coef_V(ratings, l = 0, s = 3,
                          prior_alpha = 0.5, prior_beta = 0.5,
                          verbose = FALSE)

# Prior esceptico (espera valores bajos)
result_esceptico <- coef_V(ratings, l = 0, s = 3,
                           prior_alpha = 1, prior_beta = 3,
                           verbose = FALSE)

cat("Con prior uniforme:  V =", result_uniforme$bayesiano_media, "\n")
cat("Con prior Jeffreys:  V =", result_jeffreys$bayesiano_media, "\n")
cat("Con prior esceptico: V =", result_esceptico$bayesiano_media, "\n")

# =============================================================================
# EJEMPLO 10: Funciones multi
# =============================================================================

cat("\n=== EJEMPLO 10: Funciones multi ===\n")

# Matriz de jueces por items
ratings_matrix <- matrix(c(
  3, 3, 2, 3, 3,  # Item 1
  2, 3, 3, 2, 3,  # Item 2
  3, 2, 3, 3, 2,  # Item 3
  3, 3, 3, 3, 3   # Item 4
), nrow = 4, byrow = TRUE)

# Acuerdo de cada juez vs los demas
result_a_multi <- coef_A_multi(t(ratings_matrix), c = 4)
print(result_a_multi)

# Consistencia inter-item
result_i_multi <- coef_I_multi(ratings_matrix, c = 4)
print(result_i_multi)

cat("\n=== Ejemplos completados ===\n")
