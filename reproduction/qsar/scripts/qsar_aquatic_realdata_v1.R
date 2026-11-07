#!/usr/bin/env Rscript
# Local-only UCI QSAR nested-CV driver. Production is a user-terminal command.
driver <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]))
root <- dirname(dirname(driver))
source(file.path(root, "R/penalties.R"))
source(file.path(root, "R/lambert_family.R"))
source(file.path(root, "R/performance_models.R"))
source(file.path(root, "R/performance_run.R"))
source(file.path(root, "R/qsar_aquatic_realdata_v1.R"))
qsar_v1_main(root, commandArgs(trailingOnly = TRUE))
