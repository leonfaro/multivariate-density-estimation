# prep_hepmass_10k_preprocess.R
# Zweck: 1) 1000_train/test laden  2) je TOP_N=10k Zeilen nehmen  3) Preprocessing  4) 2 CSV speichern
TOP_N <- 10000
OUT_TRAIN <- "hepmass_train_10k_preprocessed.csv"
OUT_TEST  <- "hepmass_test_10k_preprocessed.csv"

options(stringsAsFactors = FALSE)

install_if_missing <- function(pkgs){
  to_get <- pkgs[!pkgs %in% rownames(installed.packages())]
  if(length(to_get)) install.packages(to_get, repos = "https://cloud.r-project.org")
}
install_if_missing(c("data.table"))
library(data.table)

# bevorzugt lokale .csv, dann .csv.gz; sonst Download von UCI
resolve_path <- function(fname_csv, base="https://archive.ics.uci.edu/ml/machine-learning-databases/00347"){
  f_csv <- file.path(".", fname_csv)
  f_gz  <- paste0(f_csv, ".gz")
  if (file.exists(f_csv)) return(f_csv)
  if (file.exists(f_gz))  return(f_gz)
  url <- file.path(base, basename(f_gz))
  message("Lade: ", url)
  download.file(url, f_gz, mode = "wb", quiet = TRUE)
  return(f_gz)
}

# 1) Pfade ermitteln
path_train <- resolve_path("1000_train.csv")
path_test  <- resolve_path("1000_test.csv")

# 2) Jeweils nur TOP_N Zeilen einlesen (schnell, deterministisch)
tr_raw <- fread(path_train, nrows = TOP_N)
te_raw <- fread(path_test,  nrows = TOP_N)

# 3) Preprocessing
# 3.1 nur positive Klasse (Label==1 in erster Spalte), dann Label-Spalte entfernen
tr <- tr_raw[tr_raw[[1]] == 1][ , -1]
te <- te_raw[te_raw[[1]] == 1][ , -1]

# Trainings- und Testspalten auf konsistente Namen bringen
if (ncol(tr) != ncol(te)) {
  stop(sprintf("Spaltenanzahl inkonsistent: Train=%d, Test=%d", ncol(tr), ncol(te)))
}
feature_names <- paste0("f", seq_len(ncol(tr)) - 1L)
setnames(tr, feature_names)
setnames(te, feature_names)

# 3.2 diskrete Attribute eliminieren (ganzzahlig und <=10 Ausprägungen; basierend auf Train)
is_int_col <- function(v) all(is.finite(v)) && all(abs(v - round(v)) < 1e-8)
is_discr <- sapply(tr, function(v) is_int_col(v) && length(unique(v)) <= 10)
if (any(is_discr)) {
  tr <- tr[, !is_discr, with=FALSE]
  te <- te[, !is_discr, with=FALSE]
}

# 3.3 5 Features mit vielen Wiederholungen entfernen (Train-basiert)
max_prop <- function(v){ tb <- table(v); as.numeric(max(tb)) / length(v) }
props <- sapply(tr, max_prop)
drop5 <- if (length(props) > 0) order(props, decreasing=TRUE)[seq_len(min(5, length(props)))] else integer(0)
if (length(drop5)) {
  tr <- tr[, -drop5, with=FALSE]
  te <- te[, -drop5, with=FALSE]
}

# 3.4 hoch korrelierte Features entfernen (|ρ| > 0.98; Train-basiert)
drop_high_corr <- function(DT, thr=0.98){
  X <- as.matrix(DT)
  s <- apply(X, 2, sd)
  nz <- which(is.finite(s) & s > 0)
  if (length(nz) == 0) return(DT[, nz, with=FALSE])
  X <- X[, nz, drop=FALSE]
  if (ncol(X) <= 1) return(as.data.table(X))
  C <- cor(X); diag(C) <- 0
  keep <- rep(TRUE, ncol(C))
  while(TRUE){
    mx <- suppressWarnings(max(abs(C), na.rm=TRUE))
    if(!is.finite(mx) || mx <= thr) break
    ij <- which(abs(C) == mx, arr.ind=TRUE)[1,]
    mac <- colMeans(abs(C), na.rm=TRUE)
    drop <- if(mac[ij[1]] >= mac[ij[2]]) ij[1] else ij[2]
    keep[drop] <- FALSE
    C[drop,] <- 0; C[,drop] <- 0
  }
  kept_names <- colnames(X)[keep]
  as.data.table(X[, kept_names, drop=FALSE])
}
tr <- drop_high_corr(tr, 0.98)
kept_cols <- names(tr)
te <- te[, kept_cols, with=FALSE]

# 3.5 Standardisierung mit Train-Statistiken
mu <- sapply(tr, mean)
sdv <- sapply(tr, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
std <- function(DT, mu, sdv){
  as.data.table(mapply(function(x,m,s) (x - m)/s, DT, mu[names(DT)], sdv[names(DT)], SIMPLIFY=FALSE))
}
tr <- std(tr, mu, sdv)
te <- std(te, mu, sdv)

# Spaltennamen neutralisieren
setnames(tr, paste0("x", seq_len(ncol(tr))))
setnames(te, names(tr))

# 4) Speichern
fwrite(tr, OUT_TRAIN)
fwrite(te, OUT_TEST)

cat(sprintf("Gespeichert: %s (%d×%d), %s (%d×%d)\n",
            OUT_TRAIN, nrow(tr), ncol(tr),
            OUT_TEST,  nrow(te), ncol(te)))
