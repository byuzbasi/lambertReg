#!/usr/bin/env Rscript
driver <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]))
root <- dirname(dirname(driver))
for (path in c("R/penalties.R", "R/lambert_family.R", "R/performance_models.R",
    "R/performance_run.R", "R/qsar_aquatic_realdata_v1.R", 
    "R/lambert_fixed_c1_realdata_v1.R",
    "R/lambert_fixed_c1_realdata_v2.R")) source(file.path(root, path))
lf2_environment(root)$rcv1_main(root, commandArgs(trailingOnly = TRUE))
