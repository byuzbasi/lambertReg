#!/usr/bin/env Rscript
# Verify and seal the original UCI archive and its unmodified CSV member.
local({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L || !args %in% c("seal", "validate"))
    stop("Use: Rscript --vanilla scripts/prepare_qsar_aquatic_v1.R seal|validate")
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]))
  root <- dirname(dirname(script))
  for (package in c("digest", "jsonlite"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
  raw_dir <- file.path(root, "data/qsar")
  archive <- file.path(raw_dir, "qsar_aquatic_toxicity.zip")
  csv <- file.path(raw_dir, "qsar_aquatic_toxicity.csv")
  receipt_path <- file.path(raw_dir, "ACQUISITION_RECEIPT.json")
  manifest_path <- file.path(raw_dir, "manifest.csv")
  expected <- c(qsar_aquatic_toxicity.zip =
    "09e45d1b2b66c5bc871a64b87ad40e452623e1c2dc375e517e4769379a9db9cd",
    qsar_aquatic_toxicity.csv =
    "cee55fadc4b021d3ff1567612a7ac32d7e212845a4ca9eb3ef983b78e0c17856")
  files <- file.path(raw_dir, names(expected))
  if (!all(file.exists(files)) || any(dir.exists(files))) stop("Raw files missing")
  hashes <- unname(vapply(files, digest::digest, character(1), file = TRUE, algo = "sha256"))
  if (!identical(hashes, unname(expected))) stop("Raw SHA-256 mismatch")
  listing <- utils::unzip(archive, list = TRUE)
  if (nrow(listing) != 1L || listing$Name[[1L]] != basename(csv) ||
      listing$Length[[1L]] != file.info(csv)$size) stop("Archive member mismatch")
  connection <- unz(archive, basename(csv), open = "rb")
  member <- readBin(connection, what = "raw", n = listing$Length[[1L]])
  close(connection)
  if (!identical(digest::digest(member, algo = "sha256", serialize = FALSE), hashes[[2L]]))
    stop("Extracted CSV differs from archive member")
  data <- utils::read.csv(csv, sep = ";", header = FALSE,
    colClasses = rep("numeric", 9L), check.names = FALSE)
  if (!identical(dim(data), c(546L, 9L)) || any(!is.finite(as.matrix(data))) ||
      any(apply(data[, 1:8, drop = FALSE], 2L, stats::sd) <= 0))
    stop("Raw data do not satisfy the frozen 546 x 9 numeric schema")
  manifest <- data.frame(file = names(expected), bytes = as.numeric(file.info(files)$size),
    sha256 = hashes, stringsAsFactors = FALSE)
  receipt <- list(schema = "uci_qsar_aquatic_acquisition_v1", dataset_id = 505L,
    source_url = "https://archive.ics.uci.edu/static/public/505/qsar%2Baquatic%2Btoxicity.zip",
    archive_member = basename(csv), archive_sha256 = hashes[[1L]],
    csv_sha256 = hashes[[2L]], observations = 546L, predictors = 8L,
    response_column = 9L, delimiter = ";", header = FALSE,
    response_transformation = "none")
  json <- paste(as.character(jsonlite::toJSON(receipt, auto_unbox = TRUE,
    pretty = TRUE, digits = 17)), collapse = "\n")
  if (args == "seal") {
    if (file.exists(receipt_path) || file.exists(manifest_path))
      stop("Provenance files already exist; refusing overwrite")
    writeLines(json, receipt_path)
    utils::write.csv(manifest, manifest_path, row.names = FALSE)
  } else {
    if (!file.exists(receipt_path) || !file.exists(manifest_path))
      stop("Provenance files missing")
    saved_manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    saved_receipt <- paste(readLines(receipt_path, warn = FALSE), collapse = "\n")
    if (!identical(saved_manifest$file, manifest$file) ||
        !identical(as.numeric(saved_manifest$bytes), manifest$bytes) ||
        !identical(saved_manifest$sha256, manifest$sha256) ||
        !identical(saved_receipt, json))
      stop("Provenance validation failed")
  }
  cat("QSAR_SOURCE_VALID: 546 observations, 8 predictors; raw SHA-256 verified; model fits: 0.\n")
})
