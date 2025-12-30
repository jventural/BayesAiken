# =============================================================================
# BayesAiken - Ejemplo Reproducible Completo
# =============================================================================
# Este script demuestra todas las funcionalidades de la libreria BayesAiken
# para calcular intervalos de credibilidad bayesianos para coeficientes de
# validez de contenido.
#
# Autor: Jose Ventura-Leon
# Fecha: 2024-12-30
# =============================================================================

# =============================================================================
# PASO 0: INSTALACION DE LA LIBRERIA
# =============================================================================
# Primero instalar devtools si no lo tienes
# install.packages("devtools")

# Instalar BayesAiken desde la ruta local
# devtools::install("D:/14. LIBRERIAS/BayesAiken")

# O si esta en GitHub:
# devtools::install_github("jventural/BayesAiken")

# =============================================================================
# PASO 1: CARGAR LA LIBRERIA
# =============================================================================
library(BayesAiken)

# =============================================================================
# SECCION 1: DATOS DE EJEMPLO
# =============================================================================
# Simulamos un estudio de validez de contenido donde 10 jueces expertos
# evaluan 5 items de un instrumento psicologico en una escala de 0-3
# (0=nada relevante, 1=poco relevante, 2=relevante, 3=muy relevante)

set.seed(2024)  # Para reproducibilidad

# Crear datos simulados realistas
datos_validez <- data.frame(
  item = rep(paste0("Item_", 1:5), each = 3),
  criterio = rep(c("Relevancia", "Representatividad", "Claridad"), 5),
  J1 = c(3,3,3, 3,3,2, 3,2,3, 2,3,3, 3,3,2),
  J2 = c(3,3,2, 3,2,3, 3,3,2, 3,3,3, 2,3,3),
  J3 = c(3,2,3, 2,3,3, 3,3,3, 3,2,3, 3,2,3),
  J4 = c(2,3,3, 3,3,3, 2,3,3, 3,3,2, 3,3,3),
  J5 = c(3,3,3, 3,3,2, 3,2,3, 3,3,3, 3,3,3),
  J6 = c(3,3,2, 3,2,3, 3,3,2, 2,3,3, 3,2,3),
  J7 = c(3,2,3, 2,3,3, 3,3,3, 3,3,3, 2,3,3),
  J8 = c(2,3,3, 3,3,3, 2,3,2, 3,2,3, 3,3,2),
  J9 = c(3,3,3, 3,3,2, 3,2,3, 3,3,3, 3,3,3),
  J10 = c(3,3,2, 3,2,3, 3,3,3, 3,3,2, 3,2,3)
)

cat("\n")
cat("================================================================\n")
cat("           EJEMPLO REPRODUCIBLE - BayesAiken v1.0.0\n")
cat("================================================================\n")
cat("\nDatos de ejemplo:\n")
print(head(datos_validez, 10))

# =============================================================================
# SECCION 2: COEFICIENTE V DE AIKEN (UN SOLO ITEM)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 2: COEFICIENTE V DE AIKEN (ITEM INDIVIDUAL)\n")
cat("================================================================\n")

# Extraer calificaciones del Item 1, criterio Relevancia
ratings_item1 <- as.numeric(datos_validez[datos_validez$item == "Item_1" &
                                           datos_validez$criterio == "Relevancia",
                                           paste0("J", 1:10)])

cat("\nCalificaciones Item 1 (Relevancia):", ratings_item1, "\n")

# Calcular V de Aiken con enfoque bayesiano
resultado_V <- coef_V(ratings_item1, l = 0, s = 3,
                      prior_alpha = 1, prior_beta = 1,
                      cred_level = 0.95, verbose = TRUE)

# =============================================================================
# SECCION 3: TODOS LOS COEFICIENTES PARA UN ITEM
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 3: MULTIPLES COEFICIENTES (V, H, C)\n")
cat("================================================================\n")

# Usar la funcion principal para calcular V, H y C simultaneamente
resultado_multiple <- bayes_aiken(ratings_item1,
                                   coefficient = c("V", "H", "C"),
                                   l = 0, s = 3, c = 4,
                                   verbose = TRUE)

# =============================================================================
# SECCION 4: COEFICIENTE R (REPRODUCIBILIDAD TEST-RETEST)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 4: COEFICIENTE R (REPRODUCIBILIDAD)\n")
cat("================================================================\n")

# Simular calificaciones en dos momentos temporales
tiempo_1 <- c(3, 3, 2, 3, 3, 2, 3, 3, 3, 3)
tiempo_2 <- c(3, 3, 3, 3, 2, 2, 3, 3, 3, 3)

cat("\nCalificaciones Tiempo 1:", tiempo_1)
cat("\nCalificaciones Tiempo 2:", tiempo_2, "\n")

resultado_R <- coef_R(tiempo_1, tiempo_2, c = 4, verbose = TRUE)

# =============================================================================
# SECCION 5: COEFICIENTE A (ACUERDO DEL JUEZ)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 5: COEFICIENTE A (ACUERDO DEL JUEZ)\n")
cat("================================================================\n")

# Evaluar el acuerdo del Juez 1 con los demas jueces
juez_focal <- c(3, 3, 2)  # Calificaciones del Juez 1 en 3 criterios

# Calificaciones de otros jueces (matriz: filas=criterios, columnas=jueces)
otros_jueces <- matrix(c(
  3, 3, 2, 3, 3, 2, 3, 2, 3,  # Criterio 1: J2-J10
  3, 2, 3, 3, 3, 2, 3, 3, 3,  # Criterio 2: J2-J10
  2, 3, 3, 3, 2, 3, 3, 3, 2   # Criterio 3: J2-J10
), nrow = 3, byrow = TRUE)

cat("\nCalificaciones Juez focal:", juez_focal)
cat("\nMatriz otros jueces:\n")
print(otros_jueces)

resultado_A <- coef_A(juez_focal, otros_jueces, c = 4, verbose = TRUE)

# =============================================================================
# SECCION 6: COEFICIENTE I (CONSISTENCIA INTER-ITEM)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 6: COEFICIENTE I (CONSISTENCIA INTER-ITEM)\n")
cat("================================================================\n")

# Comparar consistencia entre Item 1 e Item 2
item1_ratings <- c(3, 3, 3, 2, 3, 3, 3, 2, 3, 3)
item2_ratings <- c(3, 3, 2, 3, 3, 3, 2, 3, 3, 3)

cat("\nCalificaciones Item 1:", item1_ratings)
cat("\nCalificaciones Item 2:", item2_ratings, "\n")

resultado_I <- coef_I(item1_ratings, item2_ratings, c = 4, verbose = TRUE)

# =============================================================================
# SECCION 7: ANALISIS CON DATA.FRAME (MULTIPLES ITEMS)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 7: ANALISIS DE MULTIPLES ITEMS (DATA.FRAME)\n")
cat("================================================================\n")

# Analizar todos los items del data.frame
resultado_df <- bayes_aiken(datos_validez,
                             coefficient = "V",
                             item_col = "item",
                             criterion_col = "criterio",
                             l = 0, s = 3,
                             verbose = TRUE)

# =============================================================================
# SECCION 8: COMPARACION DE PRIORS
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 8: SENSIBILIDAD A DIFERENTES PRIORS\n")
cat("================================================================\n")

cat("\nComparando el efecto de diferentes distribuciones prior:\n")

comparacion_priors <- compare_priors(ratings_item1, l = 0, s = 3,
                                      priors = list(
                                        "Uniforme Beta(1,1)" = c(1, 1),
                                        "Jeffreys Beta(0.5,0.5)" = c(0.5, 0.5),
                                        "Informativo Beta(2,2)" = c(2, 2),
                                        "Esceptico Beta(1,3)" = c(1, 3),
                                        "Optimista Beta(3,1)" = c(3, 1)
                                      ))
print(comparacion_priors)

# =============================================================================
# SECCION 9: VISUALIZACION
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 9: VISUALIZACION DE RESULTADOS\n")
cat("================================================================\n")

# Grafico de la distribucion posterior para V (individual)
cat("\n9.1 Grafico de distribucion posterior individual:\n")
plot_posterior(resultado_V$post_alpha, resultado_V$post_beta,
               coef = "V",
               prior_alpha = 1, prior_beta = 1,
               cred_level = 0.95,
               show_prior = TRUE,
               show_hdi = TRUE,
               main = "Distribucion Posterior - V de Aiken (Item 1)")

# Grafico de multiples items en un panel
cat("\n9.2 Panel de multiples items (20 items simulados):\n")

# Crear datos para 20 items
items_data <- lapply(1:20, function(i) {
  list(
    post_alpha = 25 + sample(-5:10, 1),
    post_beta = 3 + sample(0:3, 1),
    item_name = paste0("Item ", i)
  )
})

# Panel de 4 filas x 5 columnas
plot_multiple_posteriors(items_data,
                         n_cols = 5,
                         main_title = "Posterior Distributions - Items 1-20",
                         coef = "V")

# =============================================================================
# SECCION 10: EXPORTACION DE RESULTADOS
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 10: EXPORTACION DE RESULTADOS\n")
cat("================================================================\n")

# Capturar salida verbose
cat("\nCapturando salida verbose...\n")
captura <- bayes_aiken_capture(ratings_item1, coefficient = "V", l = 0, s = 3)

# Mostrar texto capturado (primeras lineas)
cat("\nPrimeras 20 lineas del texto capturado:\n")
cat(paste(head(captura$text, 20), collapse = "\n"))
cat("\n...\n")

# Para exportar a archivo (descomentar para usar):
# bayes_aiken_export(resultado_multiple, "resultados_validez.txt", format = "both")

# =============================================================================
# SECCION 11: TABLA RESUMEN FINAL
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 11: TABLA RESUMEN FINAL\n")
cat("================================================================\n")

cat("\nResumen del analisis:\n")
tabla_resumen <- summary_table(resultado_multiple)
print(tabla_resumen)

# =============================================================================
# SECCION 12: INTERPRETACION DE RESULTADOS
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 12: GUIA DE INTERPRETACION\n")
cat("================================================================\n")

cat("
INTERPRETACION DE RESULTADOS BAYESIANOS:

1. ESTIMACION PUNTUAL:
   - V Clasico: Estimacion tradicional de Aiken
   - V Bayesiano: Media de la distribucion posterior

2. INTERVALOS DE CREDIBILIDAD:
   - ETI (Equal-Tailed Interval): Percentiles simetricos
   - HDI (Highest Density Interval): Region de mayor densidad
   - El HDI es preferible porque contiene los valores mas probables

3. PROBABILIDADES:
   - P(V > 0.70): Probabilidad de validez adecuada
   - P(V > 0.80): Probabilidad de validez excelente

4. CRITERIOS DE DECISION:
   - Si HDI_inferior >= 0.70: Evidencia fuerte de validez adecuada
   - Si P(V > 0.70) >= 0.95: Alta confianza en validez adecuada
   - Si P(V > 0.80) >= 0.95: Alta confianza en validez excelente

5. VENTAJAS DEL ENFOQUE BAYESIANO:
   - Proporciona probabilidades directas (no p-valores)
   - Incorpora incertidumbre de forma natural
   - Permite incorporar conocimiento previo via priors
   - Los intervalos tienen interpretacion probabilistica directa
")

cat("\n")
cat("================================================================\n")
cat("           FIN DEL EJEMPLO REPRODUCIBLE\n")
cat("================================================================\n")
