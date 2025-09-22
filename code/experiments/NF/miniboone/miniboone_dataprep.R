# MiniBooNE -> 2x10k CSVs (MAF-style preprocessing on first 20k positives)
# Expected input: MiniBooNE_PID.txt in working dir

suppressPackageStartupMessages({
  library(utils)
  library(stats)
})

txt <- "MiniBooNE_PID.txt"
if (!file.exists(txt)) stop("File not found: ", txt)

# --- Read header (supports 'n_signal n_background' on line 1 OR split across 2 lines) ---
hdr_line <- readLines(txt, n = 1L, warn = FALSE)
p1 <- scan(text = hdr_line, what = integer(), quiet = TRUE)
if (length(p1) == 2L) {
  n_signal <- p1[1]; n_background <- p1[2]
  skip <- 1L
} else {
  stop("Unrecognized header format: expecting two integers for counts")
}

message(sprintf("Header: n_signal=%d | n_background=%d", n_signal, n_background))

# --- Load only positive examples (signal block) ---
use_fread <- requireNamespace("data.table", quietly = TRUE)
raw_mat <- if (use_fread) {
  data.table::fread(txt, skip = skip, header = FALSE, data.table = FALSE,
                    fill = TRUE, showProgress = FALSE)
} else {
  utils::read.table(txt, skip = skip, header = FALSE, fill = TRUE)
}

if (nrow(raw_mat) < n_signal) {
  stop("File shorter than expected: have ", nrow(raw_mat), " rows, need ", n_signal)
}
if (ncol(raw_mat) < 50) stop("Expected >=50 columns, found: ", ncol(raw_mat))
if (ncol(raw_mat) > 50) raw_mat <- raw_mat[, seq_len(50)]
raw_mat[] <- lapply(raw_mat, function(x) as.numeric(as.character(x)))

pos <- raw_mat[seq_len(n_signal), , drop = FALSE]
ok <- stats::complete.cases(pos)
pos <- pos[ok, , drop = FALSE]

# --- Remove rows where all features equal -1000 (documented outliers) ---
is_out <- rowSums(pos == -1000) == ncol(pos)
pos <- pos[!is_out, , drop = FALSE]

# --- Drop rows containing sentinel placeholders (-1000 or -999) in any column ---
pos <- pos[apply(pos <= -999, 1L, function(r) !any(r)), , drop = FALSE]

# --- Trim to first 20k positives after cleanup (preserve order) ---
N_total <- 20000L
if (nrow(pos) < N_total) stop("Not enough rows after cleanup: have ", nrow(pos))
X20 <- pos[seq_len(N_total), , drop = FALSE]

# --- Remove 7 "spiky" features (highest mode frequency); expect D = 43 ---
mode_ratio <- vapply(X20, function(v) {
  tb <- table(v, useNA = "no")
  max(tb) / length(v)
}, numeric(1))
nuniq <- vapply(X20, function(v) length(unique(v)), integer(1))
drop_order <- order(-mode_ratio, nuniq)
drop_idx <- head(drop_order, 7L)
keep_idx <- setdiff(seq_len(ncol(X20)), drop_idx)
X20 <- as.matrix(X20[, keep_idx, drop = FALSE])
colnames(X20) <- paste0("X", sprintf("%02d", keep_idx))
if (ncol(X20) != 43L) warning("Resulting dims: ", ncol(X20), " (expected 43)")

# --- Split: first 10k -> train, next 10k -> test ---
train_raw <- X20[1:10000,  , drop = FALSE]
test_raw  <- X20[10001:20000, , drop = FALSE]

# --- Standardize using train stats only ---
mu  <- colMeans(train_raw)
sdv <- apply(train_raw, 2, sd); sdv[sdv == 0] <- 1
scale_fn <- function(M) sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
train <- scale_fn(train_raw)
test  <- scale_fn(test_raw)
colnames(train) <- colnames(test) <- colnames(train_raw)

# --- Write outputs ---
utils::write.csv(train, "miniboone_train.csv", row.names = FALSE)
utils::write.csv(test,  "miniboone_test.csv",  row.names = FALSE)
writeLines(paste0("X", sprintf("%02d", drop_idx)), "miniboone_dropped_columns.txt")
utils::write.csv(
  data.frame(feature = colnames(train), mu = as.numeric(mu), sigma = as.numeric(sdv)),
  "miniboone_train_standardization.csv", row.names = FALSE
)

cat("Done. Train/Test: 10k each. D=", ncol(train), "\n", sep = "")
