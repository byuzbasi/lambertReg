#!/usr/bin/env Rscript
# Production execution is reserved for the user's terminal.
driver <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]))
root <- dirname(dirname(driver))
for (path in c("R/penalties.R", "R/lambert_family.R", "R/performance_models.R",
    "R/performance_run.R", "R/qsar_aquatic_realdata_v1.R", "R/diabetes_realdata_v1.R",
    "R/hitters_realdata_v1.R", "R/lambert_fixed_c1_realdata_v1.R")) source(file.path(root, path))
runner <- lf1_environment(root)
runner$rcv1_main(root, commandArgs(trailingOnly = TRUE))
