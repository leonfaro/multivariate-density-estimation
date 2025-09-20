#!/usr/bin/env Rscript

# Power dataset runner with no artifacts/logs. Wraps the embedded TRTF power script.
Sys.setenv(NO_ARTIFACTS = "1")
options(mde.no_artifacts = TRUE)

source(file.path("experiments", "tools", "trtf_power.R"))

invisible(trtf_power())
