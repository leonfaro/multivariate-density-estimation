
# 06_KLDiv_Viz.R
# Erweiterte Visualisierung der KL-Divergenz und Dichtevergleiche


# 1) Balkendiagramm mit optionalen Fehlerbalken
#    Barplot der KL-Divergenz
viz_kl_divergence_bar <- function(kl_list, kl_se=NULL) {
  methods <- names(kl_list)
  bp <- barplot(height = kl_list,
                main   = "KL-Divergenz (niedriger = besser)",
                xlab   = "Methoden",
                ylab   = "KL Divergence",
                col    = c("steelblue", "tomato", "forestgreen"),
                ylim   = c(0, 1.1 * max(kl_list + if(is.null(kl_se)) 0 else kl_se, na.rm=TRUE))
  )
  if(!is.null(kl_se)) {
    stopifnot(length(kl_se) == length(kl_list))
    arrows(x0 = bp, y0 = kl_list - kl_se,
           x1 = bp, y1 = kl_list + kl_se,
           angle = 90, code = 3, length = 0.05, col = "black")
  }
  axis(1, at = bp, labels = methods, tick = FALSE)
}

# 2) Linienplot der KL-Divergenz über verschiedene Stichprobengrößen
#    kl_values: Matrix, Zeilen=Methoden, Spalten=verschiedene N
viz_kl_vs_n <- function(sample_sizes, kl_values) {
  matplot(x = sample_sizes, y = t(kl_values), type = "b", pch = 19,
          lty = 1, col = c("steelblue", "tomato", "forestgreen"),
          xlab = "Stichprobengröße (N)", ylab = "KL Divergence",
          main = "KL-Divergenz vs. Sample-Size")
  legend("topright", legend = rownames(kl_values),
         col = c("steelblue", "tomato", "forestgreen"), pch = 19, lty = 1)
}

# 3) Overlay von true vs. estimated Dichte (1D)
#    Zeigt wahre vs. geschätzte Dichteverläufe
viz_density_overlay <- function(xvals, dens_true, dens_est, main_title = "True vs. Estimated") {
  plot(xvals, dens_true, type = "l", col = "red", lwd = 2,
       main = main_title, xlab = "x", ylab = "density")
  lines(xvals, dens_est, col = "blue", lty = 2, lwd = 2)
  legend("topright", legend = c("True", "Estimated"), col = c("red", "blue"),
         lty = c(1, 2), lwd = 2)
}

# 4) Heatmap/Contour der Differenz log(p)-log(q) in einem 2D-Gitter
#    Negative Werte => q > p, Positive => p > q
viz_logdiff_heatmap <- function(grid_df, logp, logq, nx, ny, title_str = "log(p)-log(q)") {
  logdiff <- logp - logq
  diffmat <- matrix(logdiff, nrow = nx, ncol = ny)
  xvals <- sort(unique(grid_df[, 1]))
  yvals <- sort(unique(grid_df[, 2]))
  
  image(xvals, yvals, diffmat, col = rev(heat.colors(50)),
        main = title_str, xlab = "x", ylab = "y")
  contour(xvals, yvals, diffmat, add = TRUE, col = "black")
  legend("topright",
         legend = c("log(p) > log(q)", "log(p) < log(q)"),
         fill = c(heat.colors(2)[2], heat.colors(2)[1]), bty = "n")
}

# 5) Einfache Wrapper-Funktion:
#    - Ruft barplot der KL-Werte auf
viz_kl_divergence <- function(kl_list) {
  viz_kl_divergence_bar(kl_list)
}

# --- Abschnitt zur Anzeige aller Plots ---

# Speichern der aktuellen Plot-Parameter
old_par <- par(no.readonly = TRUE)
# Layout für 4 Plots in einem Fenster
par(mfrow = c(2, 2))

# 1) Balkendiagramm der KL-Divergenz mit Fehlerbalken
kl_list <- c("TF" = 5.315, "KDE" = 7.851, "PM" = 2.814)
kl_se <- c("TF" = 0.3, "KDE" = 0.5, "PM" = 0.2)
viz_kl_divergence_bar(kl_list, kl_se)

# 2) Linienplot der KL-Divergenz über verschiedene Stichprobengrößen
sample_sizes <- c(100, 500, 1000, 5000)
kl_values <- matrix(c(5.5, 5.2, 5.3, 5.315,
                      8.0, 7.9, 7.85, 7.851,
                      3.0, 2.9, 2.8, 2.814),
                    nrow = 3, byrow = TRUE)
rownames(kl_values) <- c("TF", "KDE", "PM")
viz_kl_vs_n(sample_sizes, kl_values)

# 3) Overlay von true vs. estimated Dichte (1D)
xvals <- seq(0, 10, length.out = 100)
dens_true <- dnorm(xvals, mean = 5, sd = 1)
dens_est <- dnorm(xvals, mean = 5.2, sd = 1.1)
viz_density_overlay(xvals, dens_true, dens_est, main_title = "True vs. Estimated Density")

# 4) Heatmap/Contour der Differenz log(p)-log(q) in einem 2D-Gitter
nx <- 50; ny <- 50
xgrid <- seq(0, 1, length.out = nx)
ygrid <- seq(0, 1, length.out = ny)
grid_df <- expand.grid(x = xgrid, y = ygrid)
logp <- dnorm(grid_df$x, mean = 0.5, sd = 0.1, log = TRUE) + 
  dnorm(grid_df$y, mean = 0.5, sd = 0.1, log = TRUE)
logq <- dnorm(grid_df$x, mean = 0.55, sd = 0.1, log = TRUE) + 
  dnorm(grid_df$y, mean = 0.55, sd = 0.1, log = TRUE)
viz_logdiff_heatmap(grid_df, logp, logq, nx, ny, title_str = "log(p)-log(q) Heatmap")

# Zurücksetzen der ursprünglichen Plot-Parameter
par(old_par)

