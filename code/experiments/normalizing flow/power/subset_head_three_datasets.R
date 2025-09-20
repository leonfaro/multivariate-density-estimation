# subset_head_three_datasets.R
# Zweck: Von drei Dateien (POWER, GAS, HEPMASS) die ersten p% Datenzeilen
#        (exkl. Header) extrahieren und als CSV speichern.

suppressWarnings(suppressMessages({
  have_dt <- requireNamespace("data.table", quietly = TRUE)
}))

# -- Hilfsfunktion: schnelle Zeilenzahl (Unix: wc -l; Fallback: binär lesen) --
count_lines_fast <- function(path) {
  stopifnot(file.exists(path))
  out <- tryCatch(
    system2("wc", c("-l", shQuote(path)), stdout = TRUE),
    error = function(e) NA_character_
  )
  if (!is.na(out[1])) {
    n <- suppressWarnings(as.integer(sub("^\\s*(\\d+).*$", "\\1", out[1])))
    if (is.finite(n)) return(n)
  }
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  n <- 0L
  repeat {
    buf <- readBin(con, what = "raw", n = 8L * 1024L * 1024L)
    if (!length(buf)) break
    n <- n + sum(buf == as.raw(10L)) # '\n'
  }
  n
}

# -- Hilfsfunktion: Prozentzahl in Bruch umwandeln (5 -> 0.05; 0.5 -> 0.5) --
as_fraction <- function(pct) {
  stopifnot(is.numeric(pct), length(pct) == 1L, is.finite(pct), pct > 0)
  if (pct > 1) pct/100 else pct
}

# -- Kernfunktion: p% Kopf speichern --
subset_head_percent_save <- function(in_file, pct, out_file,
                                     sep_in = NULL, header = TRUE,
                                     na = c("", "?", "NA"),
                                     out_sep = ",",
                                     use_fread = have_dt) {
  frac <- as_fraction(pct)
  total <- count_lines_fast(in_file)
  if (!is.finite(total) || total <= 0) stop("Konnte Zeilenzahl nicht bestimmen: ", in_file)
  header_rows <- if (isTRUE(header)) 1L else 0L
  data_rows_total <- max(0L, total - header_rows)
  n <- max(1L, floor(data_rows_total * frac))
  
  if (use_fread) {
    args <- list(input = in_file, nrows = n, header = header,
                 na.strings = na, showProgress = FALSE)
    if (!is.null(sep_in)) args$sep <- sep_in
    DT <- do.call(data.table::fread, args)
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(DT, out_file, sep = out_sep)
  } else {
    # Base-R Fallback
    read_fun <- utils::read.table
    if (!is.null(sep_in)) {
      DF <- read_fun(in_file, nrows = n, header = header, sep = sep_in,
                     na.strings = na, comment.char = "", quote = "\"",
                     stringsAsFactors = FALSE)
    } else {
      DF <- read_fun(in_file, nrows = n, header = header,
                     na.strings = na, comment.char = "", quote = "\"",
                     stringsAsFactors = FALSE)
    }
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(DF, out_file, row.names = FALSE)
  }
  invisible(list(total_lines = total, n_rows_written = n, out_file = out_file))
}

# -- Bequemer Wrapper für POWER/GAS/HEPMASS --
subset_three_datasets <- function(
    in_power, in_gas, in_hepmass,
    pct_power, pct_gas, pct_hepmass,
    out_dir = "subset_head",
    out_sep = ",",   # Ausgabe-CSV-Separator
    use_fread = have_dt
) {
  tag <- function(p) {
    x <- if (p <= 1) p*100 else p
    format(round(x, 6), trim = TRUE, scientific = FALSE)
  }
  
  # POWER: Eingabe mit ';' (Semikolon), Header vorhanden
  out_power <- file.path(out_dir, sprintf("power_head_%sperc.csv", tag(pct_power)))
  r1 <- subset_head_percent_save(in_power, pct_power, out_power,
                                 sep_in=";", header=TRUE, out_sep=out_sep,
                                 use_fread=use_fread)
  message(sprintf("[POWER] Gesamtzeilen ≈ %d; geschrieben: %d -> %s",
                  r1$total_lines, r1$n_rows_written, r1$out_file))
  
  # GAS: Whitespace (Auto-Erkennung), Header vorhanden
  out_gas <- file.path(out_dir, sprintf("gas_head_%sperc.csv", tag(pct_gas)))
  r2 <- subset_head_percent_save(in_gas, pct_gas, out_gas,
                                 sep_in=NULL, header=TRUE, out_sep=out_sep,
                                 use_fread=use_fread)
  message(sprintf("[GAS] Gesamtzeilen ≈ %d; geschrieben: %d -> %s",
                  r2$total_lines, r2$n_rows_written, r2$out_file))
  
  # HEPMASS: Eingabe mit ',', Header vorhanden
  out_hm <- file.path(out_dir, sprintf("hepmass_head_%sperc.csv", tag(pct_hepmass)))
  r3 <- subset_head_percent_save(in_hepmass, pct_hepmass, out_hm,
                                 sep_in=",", header=TRUE, out_sep=out_sep,
                                 use_fread=use_fread)
  message(sprintf("[HEPMASS] Gesamtzeilen ≈ %d; geschrieben: %d -> %s",
                  r3$total_lines, r3$n_rows_written, r3$out_file))
  
  invisible(list(power = r1, gas = r2, hepmass = r3))
}

# ----------------------------
# >>>>>>> KONFIGURATION <<<<<<
# (hier per Hand eintragen)
# ----------------------------
# Beispielpfade (bitte anpassen):
# in_power   <- "/Users/DEINNAME/Downloads/household_power_consumption.txt"
# in_gas     <- "/Users/DEINNAME/Downloads/ethylene_CO.txt"
# in_hepmass <- "/Users/DEINNAME/Downloads/1000_train.csv"
#
# Prozentangaben (entweder 5  == 5%  oder  0.05 == 5%):
# pct_power   <- 5
# pct_gas     <- 2.5
# pct_hepmass <- 0.5
#
# Ausgabeverzeichnis:
# out_dir <- "/Users/DEINNAME/Downloads/subsets"
#
# AUSFÜHRUNG:
# res <- subset_three_datasets(
#   in_power   = in_power,
#   in_gas     = in_gas,
#   in_hepmass = in_hepmass,
#   pct_power   = pct_power,
#   pct_gas     = pct_gas,
#   pct_hepmass = pct_hepmass,
#   out_dir     = out_dir,
#   out_sep     = ","   # oder ";" wenn gewünscht
# )
# str(res)
