###############################################################################
#   1. Herunterladen & Einlesen                                              #
#   2. Nur positive Ereignisse (electron-neutrinos) behalten                 #
#   3. 11 Platzhalter-Zeilen (-1000/-999 in allen 50 Spalten) entfernen       #
#   4. 7 Spalten mit extrem häufig gleichem Wert (“Dirac-Peaks”) streichen    #
#   5. Standardisierung: (x-μ_train)/σ_train                                  #
#   6. Reproduzierbarer Split 80 %/10 %/10 % -> Train / Val / Test            #
#   7. Speichern als *.RData                                                 #
#                                                                            #
# Am Ende liegen drei Objekte vor: x_train, x_val, x_test                    #
# Dimensionen laut Paper:                                                    #
#   • Train: 29 556 × 43                                                     #
#   • Val:    3 284 × 43                                                     #
#   • Test:   3 648 × 43                                                     #
###############################################################################
# 
# # ------------------------- 1 | Download & Roh-Import -------------------------
# url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00199/MiniBooNE_PID.txt"
# tf  <- tempfile()
# download.file(url, tf, mode = "wb")          # ~87 MB
# 
# con   <- file(tf, "r")
# nums  <- scan(con, what = integer(), n = 2, quiet = TRUE)
# n_sig <- nums[1]                             # Zeile 1: #positiver Events
# # Zeile 2: #negativer Events (wird ignoriert)
# 
# raw   <- scan(con, what = numeric(), quiet = TRUE)  # restliche Zahlen
# close(con)
# 
# x <- matrix(raw, ncol = 50, byrow = TRUE)[1:n_sig, ]  # 50 Merkmale, nur positives
# 
# # --------------------- 2 | 11 vollständige Ausreißer löschen -----------------
# bad <- apply(x, 1, function(r) all(r == -1000) || all(r == -999))
# stopifnot(sum(bad) == 11)
# x <- x[!bad, ]                              # 36 488 Zeilen verbleiben
# 
# # ----------- 3 | 7 Features mit extrem häufiger Wiederholung streichen -------
# mode_ratio <- function(v) {                 # Anteil häufigster abgeschnittener Wert
# #   t <- table(signif(v, 6))
# #   max(t) / length(v)
# # }
# # ratios    <- apply(x, 2, mode_ratio)
# drop_cols <- order(ratios, decreasing = TRUE)[1:7]    # Top-7
# x         <- x[, -drop_cols]                          # 43 Dimensionen
# 
# # --------------------------- 4 | Standardisierung ----------------------------
# # Kennzahlen *nur* auf Basis des künftigen Trainings-Splits
# set.seed(42)
# N   <- nrow(x)
# idx <- sample(N)
# 
# train_idx <- idx[1:29556]
# val_idx   <- idx[29557:32840]
# test_idx  <- idx[32841:N]
# 
# mu  <- colMeans(x[train_idx, ])
# sig  <- apply(x[train_idx, ], 2, sd)
# 
# standardize <- function(a, mu, sd) sweep(sweep(a, 2, mu, `-`), 2, sd, `/`)
# x_train <- standardize(x[train_idx, ], mu, sig)
# x_val   <- standardize(x[val_idx,   ], mu, sig)
# x_test  <- standardize(x[test_idx,  ], mu, sig)
# 
# # --------------------------- 5 | Sanity-Checks -------------------------------
# stopifnot(dim(x_train) == c(29556, 43),
#           dim(x_val  ) == c( 3284, 43),
#           dim(x_test ) == c( 3648, 43))
# 
# # Keine Zeile mehr mit Platzhalter-Werten
# stopifnot(all(rowSums(x_train == -1000) == 0),
#           all(rowSums(x_val   == -1000) == 0),
#           all(rowSums(x_test  == -1000) == 0))
# 
# # Mittelwert ≈ 0 & SD ≈ 1 im Training
# stopifnot(max(abs(colMeans(x_train))) < 1e-10,
#           max(abs(apply(x_train, 2, sd) - 1)) < 1e-10)
# 
# # ------------------------------ 6 | Speichern --------------------------------
# # Wir packen die drei Splits in eine benannte Liste und sichern sie als .rds.
# saveRDS(
#   list(
#     x_train = x_train,
#     x_val   = x_val,
#     x_test  = x_test
#   ),
#   file = "miniboone.rds"
# )
# write.csv(x_train, "data/miniboone_train.csv", row.names = FALSE)
# write.csv(x_val,   "data/miniboone_val.csv",   row.names = FALSE)
# write.csv(x_test,  "data/miniboone_test.csv",  row.names = FALSE)


library("trtf")

mytrtf <- function(data, ntree = 50, mtry = floor(sqrt(ncol(data) - 1)),
                   minsplit = 25, minbucket = 20, maxdepth = 4, seed = 42){
  ymod <- lapply(colnames(data), function(y) BoxCox(as.formula(paste(y,"~1")), data=data))
  fm   <- lapply(2:ncol(data), function(j)
    as.formula(paste(colnames(data)[j], "~",
                     paste(colnames(data)[1:(j-1)], collapse="+"))))
  forests <- lapply(seq_along(fm), function(j)
    traforest(ymod[[j+1L]], formula=fm[[j]], data=data, trace=TRUE))
  structure(list(ymod=ymod, forests=forests), class="mytrtf")
}

predict.mytrtf <- function(object, newdata, type="logdensity"){
  ld1 <- predict(object$ymod[[1]], newdata=newdata, type=type)
  ld  <- lapply(object$forests, function(frst){
    q <- newdata[, variable.names(frst$model)[1]]
    diag(do.call(cbind, predict(frst, newdata=newdata, type=type, q=q))) })
  Reduce("+", ld) + ld1
}




n_train_keep <- 50   
d_keep       <- 2 
set.seed(42)                 

train <- read.csv("data/miniboone_train.csv")
val   <- read.csv("data/miniboone_val.csv")
test  <- read.csv("data/miniboone_test.csv")

keep_num <- sapply(train, is.numeric)
train <- train[ , keep_num]
val   <- val[ , keep_num]
test  <- test[ , keep_num]

if(d_keep < ncol(train)) {
  cols <- seq_len(d_keep)
  train <- train[ , cols]
  val   <- val  [ , cols]
  test  <- test [ , cols]
}

if (is.finite(n_train_keep) && !is.na(n_train_keep) && n_train_keep < nrow(train)) {
  train <- train[sample(nrow(train), n_train_keep), ]
}

(n_train <- nrow(train))
(n_val   <- min(floor(n_train * 0.1 / 0.8), nrow(val)))
(n_test  <- n_train)                                         # gleich groß wie val

val  <- val [sample(nrow(val),  n_val ), ]
test <- test[sample(nrow(test), n_test), ]

m  <- mytrtf(train)

ld     <- predict(m, newdata=test)        
avgLL  <- mean(ld)                        
err2SD <- 2 * sd(ld) / sqrt(length(ld))  # 2×Standard­fehler

cat(sprintf("Average test log-likelihood (nats): %.2f ± %.2f\n",
            round(avgLL, 2), round(err2SD, 2)))
