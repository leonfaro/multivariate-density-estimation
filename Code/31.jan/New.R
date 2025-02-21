# Beispielpakete laden
library(MASS)    # fitdistr() für parametrische Fits



# -------------------------------
# 1. Parameterfunktionen definieren
# -------------------------------
# f1: Mittelwert für Y2, abhängig von y1
f1 <- function(y) {
  y^2
}

# f2: Standardabweichung für Y2, abhängig von y1
f2 <- function(y) {
  exp(y)
}

# f3: α-Parameter für Y3, abhängig von y1 und y2
f3 <- function(y1, y2) {
  y1^2 + y2^3
}

# f4: β-Parameter für Y3, abhängig von y1 und y2
f4 <- function(y1, y2) {
  abs(y1 * y2) + 1
}

# -------------------------------
# 2. 3D-DGP Simulationsfunktion (d = 3)
# -------------------------------
simulate_3d_dgp <- function(N) {
  # Y1 ~ Chi²(df = 3)
  y1 <- rchisq(N, df = 3)
  
  # Initialisiere Vektoren für Y2 und Y3
  y2 <- numeric(N)
  y3 <- numeric(N)
  
  # Y2 | Y1 ~ Normal(mu = f1(y1), sd = f2(y1))
  for (i in 1:N) {
    mu    <- f1(y1[i])
    sigma <- f2(y1[i])
    y2[i] <- rnorm(1, mean = mu, sd = sigma)
  }
  
  # Y3 | (Y1, Y2) ~ Beta(alpha = f3(y1, y2), beta = f4(y1, y2))
  for (i in 1:N) {
    alpha <- f3(y1[i], y2[i])
    beta  <- f4(y1[i], y2[i])
    # Sicherstellen, dass Parameter positiv sind
    if (alpha <= 0) alpha <- 0.1
    if (beta  <= 0) beta  <- 0.1
    y3[i] <- rbeta(1, shape1 = alpha, shape2 = beta)
  }
  
  data.frame(y1 = y1, y2 = y2, y3 = y3)
}

# -------------------------------
# 3. Simulation ausführen
# -------------------------------
set.seed(123)
df <- simulate_3d_dgp(1000)

# -------------------------------
# 4. Erwartete Dichten und CDFs berechnen
# -------------------------------
# Y1: Gitter, erwartete Dichte und CDF (Chi²(df=3))
grid_y1 <- seq(min(df$y1), max(df$y1), length.out = 100)
dens_expected_y1 <- dchisq(grid_y1, df = 3)
cdf_expected_y1  <- pchisq(grid_y1, df = 3)

# Y2: Gitter, erwartete Dichte und CDF (als Mischung der bedingten Normalverteilungen)
grid_y2 <- seq(min(df$y2), max(df$y2), length.out = 100)
dens_expected_y2 <- sapply(grid_y2, function(x) {
  mean(dnorm(x, mean = f1(df$y1), sd = f2(df$y1)))
})
cdf_expected_y2 <- sapply(grid_y2, function(x) {
  mean(pnorm(x, mean = f1(df$y1), sd = f2(df$y1)))
})

# Y3: Gitter, erwartete Dichte und CDF (als Mischung der bedingten Beta-Verteilungen)
grid_y3 <- seq(0, 1, length.out = 100)
dens_expected_y3 <- sapply(grid_y3, function(x) {
  params_alpha <- f3(df$y1, df$y2)
  params_beta  <- f4(df$y1, df$y2)
  # Korrigiere negative oder zu kleine Parameter
  params_alpha[params_alpha <= 0] <- 0.1
  params_beta[params_beta <= 0]   <- 0.1
  mean(dbeta(x, shape1 = params_alpha, shape2 = params_beta))
})
cdf_expected_y3 <- sapply(grid_y3, function(x) {
  params_alpha <- f3(df$y1, df$y2)
  params_beta  <- f4(df$y1, df$y2)
  params_alpha[params_alpha <= 0] <- 0.1
  params_beta[params_beta <= 0]   <- 0.1
  mean(pbeta(x, shape1 = params_alpha, shape2 = params_beta))
})

# -------------------------------
# 5. Plots: 6 Plots in einem Fenster (2 Zeilen x 3 Spalten)
#    Oben: Histogramme mit erwarteter (rot) und empirischer Dichte (blau)
#    Unten: CDFs (empirisch in blau, erwartete in rot)
# -------------------------------
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

## Histogramm für Y1
hist(df$y1, probability = TRUE, main = "Histogram Y1", 
     xlab = "Y1", col = "lightgray")
lines(grid_y1, dens_expected_y1, col = "red", lwd = 2)
lines(density(df$y1), col = "blue", lwd = 2)
legend("topright", legend = c("Expected", "Empirical"), 
       col = c("red", "blue"), lty = 1, cex = 0.8)

## Histogramm für Y2
hist(df$y2, probability = TRUE, main = "Histogram Y2", 
     xlab = "Y2", col = "lightgray")
lines(grid_y2, dens_expected_y2, col = "red", lwd = 2)
lines(density(df$y2), col = "blue", lwd = 2)

## Histogramm für Y3
hist(df$y3, probability = TRUE, main = "Histogram Y3", 
     xlab = "Y3", col = "lightgray")
lines(grid_y3, dens_expected_y3, col = "red", lwd = 2)
lines(density(df$y3), col = "blue", lwd = 2)

## CDF für Y1
plot(ecdf(df$y1), main = "CDF Y1", xlab = "Y1", ylab = "CDF", 
     col = "blue", lwd = 2)
lines(grid_y1, cdf_expected_y1, col = "red", lwd = 2)

## CDF für Y2
plot(ecdf(df$y2), main = "CDF Y2", xlab = "Y2", ylab = "CDF", 
     col = "blue", lwd = 2)
lines(grid_y2, cdf_expected_y2, col = "red", lwd = 2)

## CDF für Y3
plot(ecdf(df$y3), main = "CDF Y3", xlab = "Y3", ylab = "CDF", 
     col = "blue", lwd = 2)
lines(grid_y3, cdf_expected_y3, col = "red", lwd = 2)






# =============================================================================
# DEBUGGING-SKRIPT FÜR 'data' must be of a vector type, was 'NULL' BEI tram::BoxCox
# =============================================================================

# 1) Kontrolle: Welche Pakete sind geladen, welche Funktionen sind maskiert?
cat("\n=== (1) Aktuelle Pakete und eventuell maskierte Funktionen ===\n")
print(search())             # Zeigt die Suchpfade
cat("Konflikte:\n")
print(conflicts(detail=TRUE))  # Zeigt Funktionsnamen, die mehrfach definiert sind

# 2) Welche 'BoxCox' und 'boxcox' -Funktionen existieren in deinem R-Setup?
cat("\n=== (2) Wo werden 'BoxCox' und 'boxcox' gefunden? ===\n")
cat("Ergebnisse von find('BoxCox'):\n")
print(find("BoxCox"))
cat("Ergebnisse von find('boxcox'):\n")
print(find("boxcox"))

# 3) Pakete laden (falls noch nicht):
#    Wir laden tram explizit nochmal und drucken dessen Version:
cat("\n=== (3) Lade explizit das 'tram'-Paket (nochmal) ===\n")
library(tram)
packageVersion("tram")

# 4) Struktur der Daten checken:
cat("\n=== (4) STR von df ===\n")
if (!exists("df")) {
  cat("WARNUNG: Objekt 'df' existiert nicht! Bitte lade/simuliere es zuerst.\n")
} else {
  str(df)
  
  # Gibt es y1 überhaupt?
  cat("\n=== (4a) Existiert df$y1 und ist numeric? ===\n")
  if (!("y1" %in% names(df))) {
    cat("FEHLER: Spalte 'y1' ist nicht in df vorhanden.\n")
  } else {
    cat("summary(df$y1):\n")
    print(summary(df$y1))
    cat("class(df$y1):\n")
    print(class(df$y1))
    cat("Min/Max df$y1:\n")
    print(range(df$y1))
    cat("NA-Werte:\n")
    print(any(is.na(df$y1)))
  }
}

# 5) Versuche den BoxCox-Fit explizit mit tram::BoxCox()
cat("\n=== (5) Versuch: tram::BoxCox(y1 ~ 1, data = df) ===\n")
if (!exists("df") || !("y1" %in% names(df))) {
  cat("ABBRUCH: 'df' existiert nicht oder df$y1 fehlt.\n")
} else {
  # Falls du sicher bist, dass y1 > 0 (oder hinreichend > 0):
  # Falls y1 0 oder negative Werte hat, verschiebe sie minimal.
  if (any(df$y1 <= 0)) {
    cat("WARNUNG: df$y1 enthält nicht-positive Werte.\n")
    cat("Verschiebe minimal...\n")
    df$y1[df$y1 <= 0] <- 1e-6
  }
  
  # Debug-Call mit tryCatch, damit wir Fehlermeldung abfangen und ausgeben
  res <- tryCatch({
    mod_y1_tram <- tram::BoxCox(y1 ~ 1, data = df)
    list(success=TRUE, model=mod_y1_tram)
  }, error=function(e) {
    list(success=FALSE, msg=e$message)
  })
  
  # Auswertung
  if (!res$success) {
    cat("FEHLER beim Fit:\n")
    cat(res$msg, "\n")
  } else {
    cat("\n+++ ERFOLG: mod_y1_tram passt, summary(mod_y1_tram): +++\n")
    print(summary(res$model))
  }
}

cat("\n=== ENDE: Debug-Skript ===\n")


# ------------------------------------------------------------
# 1) Objekt inspizieren ohne summary()
# ------------------------------------------------------------
cat("\n--- (A) str(mod_y1_tram) ---\n")
str(mod_y1_tram)

cat("\n--- (B) coef(mod_y1_tram) ---\n")
out_coef <- tryCatch(coef(mod_y1_tram), error=function(e) e)
print(out_coef)

cat("\n--- (C) vcov(mod_y1_tram) ---\n")
out_vcov <- tryCatch(vcov(mod_y1_tram), error=function(e) e)
print(out_vcov)

cat("\n--- (D) print(mod_y1_tram) ---\n")
out_print <- tryCatch({ print(mod_y1_tram) }, error=function(e) e)
print(out_print)

# ------------------------------------------------------------
# 2) Dann 'summary(mod_y1_tram)' in tryCatch
# ------------------------------------------------------------
cat("\n--- (E) summary(mod_y1_tram) ---\n")
res_sum <- tryCatch({
  summary(mod_y1_tram)
}, error=function(e) e)
print(res_sum)

cat("\n+++ Falls das wieder scheitert, dann bitte 'traceback()' +++\n")
traceback()
