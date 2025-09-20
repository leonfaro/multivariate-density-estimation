power_read_and_preprocess <- function(
    file_path,
    seed = 42L,
    include_time = FALSE,            # D=6: FALSE (MAF-Tabelle); D=7: TRUE (Zeit als Feature)
    add_noise = TRUE,                # Dequantisierung pro Feature
    test_frac = 0.10,                # 10% Test
    val_frac_of_remaining = 0.10,    # 10% von (Rest) als Validation
    write_csv = FALSE,               # optional: standardisierte Splits als CSV schreiben
    out_dir = file.path("experiments", "normalizing flow", "power")
) {
  stopifnot(file.exists(file_path))
  
  ## 1) Einlesen (UCI: Semikolon, Dezimalpunkt ".", fehlende Werte = "?" oder leer)
  raw <- utils::read.csv(
    file_path,
    sep = ";",
    na.strings = c("?", "NA", ""),
    stringsAsFactors = FALSE
  )
  
  req <- c(
    "Date","Time",
    "Global_active_power","Global_reactive_power","Voltage","Global_intensity",
    "Sub_metering_1","Sub_metering_2","Sub_metering_3"
  )
  miss <- setdiff(req, names(raw))
  if (length(miss)) stop("Fehlende Spalten: ", paste(miss, collapse=", "))
  
  ## 2) Numerische Spalten sicher numerisch parsen
  num_cols <- c("Global_active_power","Voltage","Global_intensity",
                "Sub_metering_1","Sub_metering_2","Sub_metering_3",
                "Global_reactive_power")
  for (nm in num_cols) raw[[nm]] <- as.numeric(raw[[nm]])
  
  ## 3) Optional: Zeit -> Minuten-im-Tag (Integer 0..1439). Rauschen kommt später zentral.
  if (include_time) {
    tt <- raw[["Time"]]
    hh <- suppressWarnings(as.integer(substr(tt, 1L, 2L)))
    mm <- suppressWarnings(as.integer(substr(tt, 4L, 5L)))
    minute_of_day <- hh * 60L + mm
  }
  
  ## 4) Verwerfen: Date, Global_reactive_power (Paper D.2)
  X <- raw[, c("Global_active_power","Voltage","Global_intensity",
               "Sub_metering_1","Sub_metering_2","Sub_metering_3"),
           drop = FALSE]
  if (include_time) X[["MinuteOfDay"]] <- minute_of_day
  
  ## 5) Zeilen mit NA in den verwendeten Features entfernen (robust ggü. "?" im Rohfile)
  X <- X[stats::complete.cases(X), , drop = FALSE]
  
  ## 6) Dequantisierung: pro Feature gleichverteiltes Rauschen U(0, ε_i)
  ##    ε_i = 0.5 * (minimale positive Differenz der eindeutigen Werte) 
  ##    Fallback: ε_i = max(1e-6, 1e-6 * sd), Time: ε=1 (Minutenraster).
  if (add_noise) {
    set.seed(seed + 20000L)
    for (j in seq_len(ncol(X))) {
      x <- X[[j]]
      eps <- NA_real_
      if (include_time && names(X)[j] == "MinuteOfDay") {
        eps <- 1.0
      } else {
        u <- sort(unique(x[is.finite(x)]))
        du <- diff(u)
        min_pos <- suppressWarnings(min(du[du > 0], na.rm = TRUE))
        if (!is.finite(min_pos) || min_pos <= 0) {
          s <- stats::sd(x, na.rm = TRUE)
          min_pos <- if (is.finite(s) && s > 0) s * 1e-6 else 1e-6
        }
        eps <- 0.5 * min_pos
      }
      # add uniform noise in [0, eps]
      X[[j]] <- x + stats::runif(length(x), min = 0, max = eps)
    }
  }
  
  ## 7) Splits: 10% Test; 10% des Rests Val; Rest Train (deterministisch)
  N <- nrow(X)
  if (N < 10) stop("Zu wenige Zeilen nach Bereinigung: N=", N)
  set.seed(seed)
  idx <- sample.int(N)
  n_test <- floor(test_frac * N)
  n_val  <- floor(val_frac_of_remaining * (N - n_test))
  
  idx_te <- idx[seq_len(n_test)]
  idx_va <- idx[seq.int(n_test + 1L, n_test + n_val)]
  idx_tr <- idx[-seq_len(n_test + n_val)]
  
  X_tr <- as.matrix(X[idx_tr, , drop = FALSE])
  X_va <- as.matrix(X[idx_va, , drop = FALSE])
  X_te <- as.matrix(X[idx_te, , drop = FALSE])
  
  ## 8) Train-only Standardisierung
  mu    <- colMeans(X_tr)
  sigma <- apply(X_tr, 2L, stats::sd)
  sigma[!is.finite(sigma) | sigma == 0] <- 1
  
  standardize <- function(M) {
    sweep(sweep(M, 2L, mu, FUN = "-"), 2L, sigma, FUN = "/")
  }
  X_tr_std <- standardize(X_tr)
  X_va_std <- standardize(X_va)
  X_te_std <- standardize(X_te)
  
  ## 9) Akzeptanztests (Formate, Numerik)
  stopifnot(
    ncol(X_tr_std) == (6L + as.integer(include_time)),
    all(is.finite(X_tr_std)), all(is.finite(X_va_std)), all(is.finite(X_te_std))
  )
  
  ## 10) Optional: persistieren
  if (write_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    utils::write.csv(X_tr_std, file = file.path(out_dir, "power_train.csv"), row.names = FALSE)
    utils::write.csv(X_va_std, file = file.path(out_dir, "power_val.csv"),   row.names = FALSE)
    utils::write.csv(X_te_std, file = file.path(out_dir, "power_test.csv"),  row.names = FALSE)
    saveRDS(list(mu = mu, sigma = sigma, features = colnames(X_tr_std)),
            file = file.path(out_dir, "power_standardization.rds"))
  }
  
  list(
    X_train  = X_tr_std,
    X_val    = X_va_std,
    X_test   = X_te_std,
    mu       = mu,
    sigma    = sigma,
    features = colnames(X_tr_std),
    index    = list(train = idx_tr, val = idx_va, test = idx_te)
  )
}



res <- power_read_and_preprocess(
  file_path   = "/Users/leonkiafaro/Downloads/power.txt",
  seed        = 42L,
  include_time = FALSE,   # D=6 (entspricht Tabelle 5)
  add_noise    = TRUE,
  write_csv    = TRUE,
  out_dir      = file.path("experiments", "normalizing flow", "power")
)

str(res$X_train)  # standardisierte Trainingsmatrix (N_tr x D)
res$features      # Spaltennamen
