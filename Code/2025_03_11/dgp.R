
## -- Zusätzliche Pakete
library("trtf")
library("VineCopula")

## -- (A) DGP + Transformation-Forest-Code (unverändert)
set.seed(29)

### config für 3 Dimensionen
config <- list(
  ### Y1
  list(distr = "norm", parm = function(data) {list(mean = 0, sd = 1)}),
  ### Y2
  list(distr = "exp", parm = function(data) {list(rate = 1 + 0.1 * data[["Y1"]])}),
  ### Y3
  list(distr = "gamma", parm = function(data) {list(shape = 1 + data[["Y2"]], rate = 1)}),
  ### Y4
  list(distr = "norm", parm = function(data) {list(mean = data[["Y3"]], sd = 1)}),
  ### Y5
  list(distr = "exp", parm = function(data) {list(rate = 1 + 0.1 * data[["Y4"]])}),
  ### Y6
  list(distr = "gamma", parm = function(data) {list(shape = 1 + data[["Y5"]], rate = 1)})
  # ### Y7
  # list(distr = "norm", parm = function(data) {list(mean = data[["Y6"]], sd = 1)}),
  # ### Y8
  # list(distr = "exp", parm = function(data) {list(rate = 1 + 0.1 * data[["Y7"]])}),
  # ### Y9
  # list(distr = "gamma", parm = function(data) {list(shape = 1 + data[["Y8"]], rate = 1)}),
  # ### Y10
  # list(distr = "norm", parm = function(data) {list(mean = data[["Y9"]], sd = 1)}),
  # ### Y11
  # list(distr = "exp", parm = function(data) {list(rate = 1 + 0.1 * data[["Y10"]])}),
  # ### Y12
  # list(distr = "gamma", parm = function(data) {list(shape = 1 + data[["Y11"]], rate = 1)}),
  # ### Y13
  # list(distr = "norm", parm = function(data) {list(mean = data[["Y12"]], sd = 1)}),
  # ### Y14
  # list(distr = "exp", parm = function(data) {list(rate = 1 + 0.1 * data[["Y13"]])}),
  # ### Y15
  # list(distr = "gamma", parm = function(data) {list(shape = 1 + data[["Y14"]], rate = 1)}),
  # ### Y16
  # list(distr = "norm", parm = function(data) {list(mean = data[["Y15"]], sd = 1)}),
  # ### Y17
  # list(distr = "exp", parm = function(data) {list(rate = 1 + 0.1 * data[["Y16"]])}),
  # ### Y18
  # list(distr = "gamma", parm = function(data) {list(shape = 1 + data[["Y17"]], rate = 1)}),
  # ### Y19
  # list(distr = "norm", parm = function(data) {list(mean = data[["Y18"]], sd = 1)}),
  # ### Y20
  # list(distr = "exp", parm = function(data) {list(rate = 1 + 0.1 * data[["Y19"]])}),
  # ### Y21
  # list(distr = "gamma", parm = function(data) {list(shape = 1 + data[["Y20"]], rate = 1)}),
  # ### Y22
  # list(distr = "norm", parm = function(data) {list(mean = data[["Y21"]], sd = 1)}),
  # ### Y23
  # list(distr = "exp", parm = function(data) {list(rate = 1 + 0.1 * data[["Y22"]])}),
  # ### Y24
  # list(distr = "gamma", parm = function(data) {list(shape = 1 + data[["Y23"]], rate = 1)}),
  # ### Y25
  # list(distr = "norm", parm = function(data) {list(mean = data[["Y24"]], sd = 1)})
)

## (1) dgp mit wahrer Log-Dichte
# dgp
# Erzeugt synthetische Daten
# Input: N, config
# Output: data.frame mit Y1..Yk + Attribut logd
# Algorithmus: ruft r* pro Dimension auf, summiert d*
dgp <- function(N = 100, config) {
  ret <- matrix(0, nrow = N, ncol = length(config))
  colnames(ret) <- paste0("Y", 1:length(config))
  ret <- as.data.frame(ret)
  logd <- 0
  for (j in seq_along(config)) {
    cfg  <- config[[j]]
    rfun <- paste0("r", cfg[["distr"]])  # Zufallszahlengenerator
    dfun <- paste0("d", cfg[["distr"]])  # Dichtefunktion
    
    # Parameter fürs Ziehen
    args <- cfg[["parm"]](ret)
    args$n <- N
    ret[[j]] <- do.call(rfun, args)
    
    # Parameter für Log-Dichte
    args$n   <- NULL
    args$x   <- ret[[j]]
    args$log <- TRUE
    logd <- logd + do.call(dfun, args)
  }
  attr(ret, "logd") <- logd
  ret
}

## (2) Transformation-Forest-Modell
# mytrtf
# Baut Transformation-Forest
# Input: data.frame
# Output: Liste mit ymod und forests
# Algorithmus: BoxCox(Y1) + traforest(Y2..Yk|Vorgänger)
mytrtf <- function(data) {
  # (2a) Marginales Modell für Y1
  ymod <- lapply(colnames(data), function(y) {
    fm <- as.formula(paste(y, "~ 1"))
    BoxCox(fm, data = data)
  })
  
  # (2b) Konditionale Transformation-Forests für Y2..Y25
  forests <- lapply(2:ncol(data), function(j) {
    fm <- as.formula(paste(colnames(data)[j], "~", colnames(data)[j - 1]))
    traforest(ymod[[j]], formula = fm, data = data, trace = FALSE)
  })
  
  ret <- list(ymod = ymod, forests = forests)
  class(ret) <- "mytrtf"
  ret
}

## (3) Vorhersage (bedingte Log-Dichten) für mytrtf
# predict.mytrtf
# Sagt Forest-Logdichten voraus
# Input: Modell-Objekt, neues data.frame
# Output: Vektor bedingter Log-Dichten
# Algorithmus: predict(Y1) + sum predict(Yj|Yj-1)
predict.mytrtf <- function(object, newdata, type = "logdensity") {
  # 3a) Marginaler Beitrag (Y1)
  ld1 <- predict(object$ymod[[1]], newdata = newdata, type = "logdensity")
  # 3b) Konditionelle Beiträge
  ld_cond <- lapply(seq_along(object$forests), function(idx) {
    frst  <- object$forests[[idx]]
    j <- idx + 1
    qvals <- newdata[[paste0("Y", j)]]
    pr <- predict(frst, newdata = newdata, type = "logdensity", q = qvals)
    diag(do.call("cbind", pr))
  })
  ld1 + Reduce("+", ld_cond)
}

## R-Vine fit

# Transformiert die Daten matrixweise per cdf, um U in [0,1] zu erhalten
# rvineCdfFromConfig
# Transformiert Daten in [0,1]
# Input: data.frame, config
# Output: U-Matrix
# Algorithmus: p* pro Spalte, shape/rate/sd=1 bei invalid

rvineCdfFromConfig <- function(data, config) { 
  U <- matrix(NA_real_, nrow = nrow(data), ncol = ncol(data)) # Allozierung der U-Matrix
  colnames(U) <- colnames(data) # Übernimmt Spaltennamen
  for (j in seq_len(ncol(data))) { # Schleife über Dimensionen
    cdfname <- paste0("p", config[[j]]$distr) # Bestimmt den cdf-Funktionsnamen
    parList <- config[[j]]$parm(data) # Erzeugt bedingte Parameter
    for (i in seq_len(nrow(data))) { # Schleife über Beobachtungen
      args_i <- lapply(parList, function(x) x[i]) # Extrahiert i-ten Wert pro Parameter
      args_i$q <- data[i, j] # cdf-Argument
      args_i$x <- NULL # Entfernt potenzielles x
      for (nm in names(args_i)) { val <- args_i[[nm]] # Überprüft shape/rate/sd
      if (is.na(val)) { args_i[[nm]] <- 1 } else if (nm %in% c("shape","rate","sd") && val <= 0) { args_i[[nm]] <- 1 } }
      U[i, j] <- tryCatch(do.call(cdfname, args_i), error = function(e) NA_real_) # Aufruf der cdf
    }
  }
  U # Gibt U zurück
}

# Berechnet zu jeder Beobachtung die Summe ihrer bedingten Log-PDFs
# rvineLogpdfFromConfig
# Summiert Log-PDF pro Beobachtung
# Input: data.frame, config
# Output: Vektor Summen-Log-PDF
# Algorithmus: d* pro Spalte (log=TRUE), shape/rate/sd=1
rvineLogpdfFromConfig <- function(data, config) { logpdf_vec <- numeric(nrow(data)) # Initialisiert Summenvektor
for (j in seq_len(ncol(data))) { # Schleife über Spalten
  pdfname <- paste0("d", config[[j]]$distr) # Bestimmt die pdf-Funktion
  parList <- config[[j]]$parm(data) # Parameter generieren
  vals <- numeric(nrow(data)) # Vektor pro Spalte
  for (i in seq_len(nrow(data))) {
    args_i <- lapply(parList, function(x) x[i]) # i-te Zeile extrahieren
    args_i$x <- data[i, j] # pdf-Argument
    args_i$log <- TRUE # Log-Argument
    for (nm in names(args_i)) { val <- args_i[[nm]]
    if (is.na(val)) { args_i[[nm]] <- 1 } else if (nm %in% c("shape","rate","sd") && val <= 0) { args_i[[nm]] <- 1 } }
    vals[i] <- tryCatch(do.call(pdfname, args_i), error = function(e) NA_real_) # PDF-Berechnung
  }
  logpdf_vec <- logpdf_vec + vals # Addiert Spaltenbeitrag
}
logpdf_vec # Gibt Summen zurück
}

# Passt ein R-Vine-Modell mit logLik-Kriterium an
# fitRvine
# Passt R-Vine-Modell an
# Input: data.frame, config
# Output: Liste mit rvine und config
# Algorithmus: cdf->U, filter, RVineStructureSelect(type=0,logLik)
fitRvine <- function(data, config) {
  U <- rvineCdfFromConfig(data, config) # Wandelt Daten in U-Scale um
  ok <- complete.cases(U) # Prüft Vollständigkeit
  if (!any(ok)) stop("No valid rows for R-Vine.") # Abbruch bei fehlenden Zeilen
  U_ok <- U[ok, , drop = FALSE] # Filtert valide Zeilen
  rvine <- RVineStructureSelect(
    data          = U_ok, # Übergibt U-Daten
    familyset     = NA,   # Alle Copula-Familien
    type          = 0,    # R-Vine
    selectioncrit = "logLik", # Auswahlkriterium
    indeptest     = FALSE
  )
  list(rvine = rvine, config = config) # Gibt Modellobjekt zurück
}

# Sagt per R-Vine die bedingte Dichte voraus
# predictRvine
# Sagt R-Vine-Logdichten voraus
# Input: R-Vine-Objekt, neues data.frame
# Output: Vektor bedingter Log-Dichten
# Algorithmus: summiere Log-PDF + log(RVinePDF(U))

predictRvine <- function(object, newdata) {
  
  logMarg <- rvineLogpdfFromConfig(newdata, object$config) # Berechnet Summen-Log-PDF
  Unew    <- rvineCdfFromConfig(newdata, object$config)    # Transformiert in [0,1]
  out     <- rep(NA_real_, nrow(Unew)) # Vektor für Ergebnisse
  ok      <- complete.cases(Unew)      # Prüft Komplettheit
  if (any(ok)) {
    out[ok] <- RVinePDF(Unew[ok, , drop=FALSE], object$rvine) # Berechnet Vine-PDF
  }
  logMarg + log(out) # Summiert bedingte Log-PDF + log(Vine-PDF)
}

trainRatio <- 0.7 # 70% für Training
Ntotal <- 50 # Gesamtdatenzahl
Ntrain <- floor(Ntotal * trainRatio) # Anzahl Trainingsdaten
Ntest  <- Ntotal - Ntrain # Anzahl Testdaten

trainData <- dgp(N = Ntrain, config = config) # Erzeugt Trainingsdaten
testData  <- dgp(N = Ntest,  config = config) # Erzeugt Testdaten
trueLog   <- attr(testData, "logd") # Wahre Logdichte

tf_model   <- mytrtf(trainData) # Transformation-Forest
tf_logdens <- predict(tf_model, newdata=testData) # Prognose Log-Dichte
tf_ll_marg  <- logLik(tf_model$ymod[[1]], newdata = testData) # LogLik marg Y1
tf_ll_conds <- sum(sapply(tf_model$forests, logLik, newdata = testData)) # Summe LogLik restl. Dimensionen
tf_ll_total <- tf_ll_marg + tf_ll_conds # Gesamte LogLik

cat("\n=== Transformation-Forest Out-of-Sample-Loglikelihood ===\n") # Ausgabe
cat("Sum:", tf_ll_total, "\n") # Zeigt Gesamt-LogLik

plot(tf_logdens, trueLog,
     xlab="Predicted TF Log-Density",
     ylab="True Log-Density",
     main="Transformation-Forest: Pred. vs. True") # Plot
abline(a=0, b=1, col="red", lwd=2) # Diagonale

rvine_model <- fitRvine(trainData, config) # R-Vine-Fit
rvine_logdens <- predictRvine(rvine_model, testData) # Prognose Log-Dichte R-Vine
rvine_sum_ll <- sum(rvine_logdens, na.rm=TRUE) # Summierte Log-Dichte

cat("\n=== R-Vine Out-of-Sample-Loglikelihood ===\n") # Ausgabe R-Vine
cat("Sum:", rvine_sum_ll, "\n")

plot(rvine_logdens, trueLog,
     xlab="Predicted R-Vine Log-Density",
     ylab="True Log-Density",
     main="R-Vine (logLik): Pred. vs. True") # Plot R-Vine
abline(a=0, b=1, col="red", lwd=2)

cat("\n=== Compare Summed Log-Density ===\n") # Vergleich
cat("Transformation-Forest Sum:", sum(tf_logdens, na.rm=TRUE), "\n")
cat("R-Vine Sum               :", rvine_sum_ll, "\n")
cat("True Summed Log-Density  :", sum(trueLog), "\n")