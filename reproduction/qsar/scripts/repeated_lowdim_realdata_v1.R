#!/usr/bin/env Rscript
# Local repeated-CV runner. Production fits are a user-terminal command.
driver <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value = TRUE)[1L]))
root <- dirname(dirname(driver))
source(file.path(root, "R/penalties.R"))
source(file.path(root, "R/lambert_family.R"))
source(file.path(root, "R/performance_models.R"))
source(file.path(root, "R/performance_run.R"))
source(file.path(root, "R/qsar_aquatic_realdata_v1.R"))
source(file.path(root, "R/diabetes_realdata_v1.R"))
source(file.path(root, "R/hitters_realdata_v1.R"))
source(file.path(root, "R/repeated_lowdim_realdata_v1.R"))
rcv1_main(root, commandArgs(trailingOnly = TRUE))
