# GAS (UCI) -> MAF-kompatible CSVs mit 10k Zeilen pro Split
# Erwartet: ethylene_CO.txt im Arbeitsordner (oder Pfad unten anpassen)

co_file <- "ethylene_CO.txt"
if (!file.exists(co_file)) stop("Datei nicht gefunden: ", co_file)

# --- Laden (schnell mit data.table::fread, sonst Base R) ---
use_fread <- requireNamespace("data.table", quietly = TRUE)
df <- if (use_fread) {
  data.table::fread(co_file, fill = TRUE)
} else {
  read.table(co_file, header = TRUE, sep = "", stringsAsFactors = FALSE)
}

if (ncol(df) < 19) stop("Unerwartete Spaltenanzahl. Erwartet 19, gefunden: ", ncol(df))

# Laut UCI: Spalten 1..3 = Zeit & Setpoints; 4..19 = 16 Sensoren
X <- as.matrix(df[, 4:19])

# --- Stark korrelierte Attribute entfernen (|Pearson| > 0.98), deterministisch von links nach rechts ---
thr <- 0.98
target_features <- 8L
CM <- cor(X)  # 16x16
CM[is.na(CM)] <- 0
keep <- rep(TRUE, ncol(X))
for (i in seq_len(ncol(X) - 1)) {
  if (!keep[i]) next
  for (j in (i + 1):ncol(X)) {
    if (keep[j] && abs(CM[i, j]) > thr) keep[j] <- FALSE
  }
}

kept_idx <- which(keep)
if (length(kept_idx) > target_features) {
  drop_order <- rev(kept_idx)
  for (idx_drop in drop_order) {
    if (length(kept_idx) <= target_features) break
    keep[idx_drop] <- FALSE
    kept_idx <- which(keep)
  }
}

if (length(kept_idx) < target_features) {
  combos <- utils::combn(seq_len(ncol(X)), target_features)
  found <- NULL
  for (col_idx in seq_len(ncol(combos))) {
    cand <- combos[, col_idx]
    sub_cor <- abs(CM[cand, cand])
    diag(sub_cor) <- 0
    if (all(sub_cor <= thr)) {
      found <- cand
      break
    }
  }
  if (is.null(found)) {
    stop("Keine Kombination mit ", target_features, " Features bei Schwelle ", thr, " gefunden.")
  }
  keep <- rep(FALSE, ncol(X))
  keep[found] <- TRUE
}

X <- X[, keep, drop = FALSE]
if (ncol(X) != target_features) stop("Es bleiben ", ncol(X), " Features übrig (erwartet 8).")

sel_names_raw <- colnames(df)[3 + which(keep)]  # 4:19 -> +3 Offset
sel_names <- make.names(sel_names_raw, unique = TRUE)
colnames(X) <- sel_names

# --- Split-Konfiguration (Standard 10k, override via GAS_ROWS_PER_SPLIT) ---
rows_env <- Sys.getenv("GAS_ROWS_PER_SPLIT", "")
n_per_split <- suppressWarnings(as.integer(rows_env))
if (is.na(n_per_split) || !is.finite(n_per_split) || n_per_split <= 0L) {
  n_per_split <- 10000L
}

# --- Subsampling & Splits nach Featureauswahl ---
set.seed(170507057)  # stabil reproduzierbar
n_total <- n_per_split * 3L
if (nrow(X) < n_total) stop("Zu wenig Zeilen: habe ", nrow(X), ", brauche ", n_total)

idx <- sample.int(nrow(X), n_total, replace = FALSE)
train_idx <- idx[seq_len(n_per_split)]
val_idx   <- idx[seq.int(n_per_split + 1L, 2L * n_per_split)]
test_idx  <- idx[seq.int(2L * n_per_split + 1L, n_total)]

train_raw <- X[train_idx, , drop = FALSE]
val_raw   <- X[val_idx, , drop = FALSE]
test_raw  <- X[test_idx, , drop = FALSE]

# --- Standardisierung pro Feature (fit auf Train, apply auf alle Splits) ---
mu <- colMeans(train_raw)
sdv <- apply(train_raw, 2, sd)
sdv[sdv == 0] <- 1

standardize <- function(mat, center, scale) {
  sweep(sweep(mat, 2, center, "-"), 2, scale, "/")
}

train <- standardize(train_raw, mu, sdv)
val   <- standardize(val_raw,   mu, sdv)
test  <- standardize(test_raw,  mu, sdv)

colnames(train) <- colnames(val) <- colnames(test) <- colnames(train_raw)

# --- Dateien schreiben (im gleichen Ordner wie Eingabedatei) ---
outdir <- dirname(normalizePath(co_file))
wpath  <- function(name) file.path(outdir, name)
write.csv(train, wpath("gas_train.csv"), row.names = FALSE)
write.csv(val,   wpath("gas_valid.csv"), row.names = FALSE)
write.csv(test,  wpath("gas_test.csv"),  row.names = FALSE)

# Dokumentation der gewählten Sensorkanäle
writeLines(sel_names, wpath("gas_selected_columns.txt"))

# Standardisierungsparameter dokumentieren
std_df <- data.frame(feature = sel_names,
                     mu = as.numeric(mu),
                     sigma = as.numeric(sdv),
                     stringsAsFactors = FALSE)
write.csv(std_df, wpath("gas_train_standardization.csv"), row.names = FALSE)

cat("Fertig.\n",
    "Geschrieben: gas_train.csv (", nrow(train), "), gas_valid.csv (", nrow(val), "), gas_test.csv (", nrow(test), ")\n",
    "Standardisierung: gas_train_standardization.csv\n",
    "Gewählte Spalten: ", paste(sel_names, collapse = ", "), "\n", sep = "")
