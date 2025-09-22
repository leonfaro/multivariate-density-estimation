# ==== USER CONFIG (edit here) ==========================================
INPUT_POWER_TXT <- "household_power_consumption.txt"
N <- 1000L
SEED <- 42L
INCLUDE_TIME <- FALSE
ADD_NOISE <- TRUE
OUT_DIR <- "."
WRITE_VAL <- FALSE
# Checklist: set absolute input path, confirm row cap, toggle noise/time flags intentionally

suppressPackageStartupMessages({ library(utils); library(stats) })

power_make_splits_write_csv <- function(input_power_txt = INPUT_POWER_TXT,
                                        out_dir = OUT_DIR,
                                        seed = SEED,
                                        include_time = INCLUDE_TIME,
                                        add_noise = ADD_NOISE,
                                        n_rows_per_split = N,
                                        write_val = WRITE_VAL) {
  input_power_txt <- normalizePath(input_power_txt, winslash = "/", mustWork = TRUE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # 1) Einlesen
  raw <- utils::read.csv(input_power_txt, sep = ";",
                         na.strings = c("?", "NA", ""), stringsAsFactors = FALSE)
  req <- c("Date","Time","Global_active_power","Global_reactive_power","Voltage",
           "Global_intensity","Sub_metering_1","Sub_metering_2","Sub_metering_3")
  miss <- setdiff(req, names(raw))
  if (length(miss)) stop("Fehlende Spalten: ", paste(miss, collapse = ", "))

  num_cols <- c("Global_active_power","Voltage","Global_intensity",
                "Sub_metering_1","Sub_metering_2","Sub_metering_3","Global_reactive_power")
  for (nm in num_cols) raw[[nm]] <- as.numeric(raw[[nm]])

  # 2) Zeit -> Minute-im-Tag (Integer)
  if (include_time) {
    tt <- raw[["Time"]]
    hh <- suppressWarnings(as.integer(substr(tt, 1L, 2L)))
    mm <- suppressWarnings(as.integer(substr(tt, 4L, 5L)))
    minute_of_day <- hh * 60L + mm
  }

  # 3) Features wählen (D=6 oder 7)
  X <- raw[, c("Global_active_power","Voltage","Global_intensity",
               "Sub_metering_1","Sub_metering_2","Sub_metering_3"), drop = FALSE]
  if (include_time) X[["MinuteOfDay"]] <- minute_of_day

  # 4) NA-Zeilen raus
  X <- X[stats::complete.cases(X), , drop = FALSE]

  # 5) Dequantisierung: EINSEITIG [0, eps_j] (MAF)
  eps <- rep(0, ncol(X)); names(eps) <- colnames(X)
  if (add_noise) {
    set.seed(seed + 20000L)
    for (j in seq_len(ncol(X))) {
      x <- X[[j]]
      if (include_time && names(X)[j] == "MinuteOfDay") {
        eps[j] <- 1.0
        X[[j]] <- x + stats::runif(length(x), min = 0, max = eps[j])
      } else {
        u <- sort(unique(x[is.finite(x)]))
        du <- diff(u)
        min_pos <- suppressWarnings(min(du[du > 0], na.rm = TRUE))
        if (!is.finite(min_pos) || min_pos <= 0) {
          s <- stats::sd(x, na.rm = TRUE)
          min_pos <- if (is.finite(s) && s > 0) s * 1e-6 else 1e-6
        }
        eps[j] <- 0.5 * min_pos
        X[[j]] <- x + stats::runif(length(x), min = 0, max = eps[j])
      }
    }
  }

  # Checklist: ensure 10% test split, optional validation only when requested, reproducible shuffle
  # 6) Splits: 10% Test; Rest Train (kein Val, wenn write_val==FALSE)
  N_all <- nrow(X)
  set.seed(seed)
  idx <- sample.int(N_all)

  n_test <- floor(0.10 * N_all)
  n_val  <- if (write_val) floor(0.10 * (N_all - n_test)) else 0L

  idx_te <- idx[seq_len(n_test)]
  idx_va <- if (n_val > 0L) idx[seq.int(n_test + 1L, n_test + n_val)] else integer(0)
  idx_tr <- idx[-seq_len(n_test + n_val)]

  X_tr <- as.matrix(X[idx_tr, , drop = FALSE])
  X_va <- as.matrix(X[idx_va, , drop = FALSE])
  X_te <- as.matrix(X[idx_te, , drop = FALSE])

  # 7) Train-only Standardisierung
  mu    <- colMeans(X_tr)
  sigma <- apply(X_tr, 2L, stats::sd)
  sigma[!is.finite(sigma) | sigma == 0] <- 1
  standardize <- function(M) sweep(sweep(M, 2L, mu, "-"), 2L, sigma, "/")
  Z_tr <- standardize(X_tr)
  Z_va <- standardize(X_va)
  Z_te <- standardize(X_te)

  # 8) Limitierung pro Split (nach Standardisierung)
  limit <- function(M, n) if (!is.finite(n) || is.na(n)) M else utils::head(M, min(nrow(M), as.integer(n)))
  Z_tr_out <- limit(Z_tr, n_rows_per_split)
  Z_te_out <- limit(Z_te, n_rows_per_split)
  if (write_val) Z_va_out <- limit(Z_va, n_rows_per_split)

  # 9) Sanity-Ausgaben
  tie_ratio <- function(z) 1 - length(unique(round(z, 12))) / length(z)
  max_tie <- max(apply(Z_tr_out, 2, tie_ratio))
  message(sprintf("POWER prep | D=%d | train=%d | test=%d | max tie(train)≈%.2e",
                  ncol(Z_tr_out), nrow(Z_tr_out), nrow(Z_te_out), max_tie))

  # 10) Schreiben
  utils::write.csv(Z_tr_out, file = file.path(out_dir, "power_train.csv"), row.names = FALSE)
  utils::write.csv(Z_te_out, file = file.path(out_dir, "power_test.csv"),  row.names = FALSE)
  if (write_val) utils::write.csv(Z_va_out, file = file.path(out_dir, "power_val.csv"), row.names = FALSE)

  std_df <- data.frame(
    feature = colnames(Z_tr),
    mu = as.numeric(mu),
    sigma = as.numeric(sigma),
    eps = as.numeric(eps[colnames(Z_tr)]),
    stringsAsFactors = FALSE
  )
  utils::write.csv(std_df, file = file.path(out_dir, "power_standardization.csv"), row.names = FALSE)

  invisible(list(train = Z_tr_out, test = Z_te_out, standardization = std_df))
}

# 11) Main: nur wenn direkt ausgeführt
if (sys.nframe() == 0L) {
  res <- power_make_splits_write_csv()
  out_dir_norm <- normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE)
  message("Wrote: ", file.path(out_dir_norm, "power_train.csv"))
  message("Wrote: ", file.path(out_dir_norm, "power_test.csv"))
}
