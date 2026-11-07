# TRUBA execution overlay for the frozen performance-v1 scientific workflow.
# This file changes scheduling, persistence and telemetry only.

perf_truba_release_id <- function() "performance_confirmation_truba_v1_a03"

perf_truba_source_paths <- function() unique(c(
  perf_source_paths(),
  "R/performance_truba_v1.R",
  "scripts/performance_truba_v1.R",
  "scripts/validate_performance_truba_v1.R",
  "tests/performance_truba_v1_checks.R",
  "protocol/performance_confirmation_v1.md",
  "protocol/performance_confirmation_truba_v1.md",
  "truba/run_performance_confirmation_v1.slurm",
  "truba/build_performance_confirmation_bundle_v1.R",
  "truba/README_TRUBA_PERFORMANCE_CONFIRMATION_V1.md"
))

perf_truba_path <- function(root, relative) {
  if (length(relative) != 1L || is.na(relative) || !nzchar(relative) ||
      startsWith(relative, "/") || grepl("(^|/)\\.\\.(/|$)", relative)) {
    stop("Unsafe relative path: ", relative)
  }
  file.path(root, relative)
}

perf_truba_inventory <- function(root, files) {
  files <- sort(unique(files), method = "radix")
  paths <- vapply(files, function(x) perf_truba_path(root, x), character(1))
  if (!length(files) || any(!file.exists(paths)) || any(dir.exists(paths))) {
    stop("Required inventory file is missing or is not a regular file")
  }
  data.frame(
    file = files,
    bytes = as.numeric(file.info(paths)$size),
    sha256 = unname(vapply(paths, perf_hash, character(1))),
    stringsAsFactors = FALSE
  )
}

perf_truba_sources <- function(root) {
  perf_truba_inventory(root, perf_truba_source_paths())
}

perf_truba_write_utf8 <- function(connection, value) {
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop("Canonical text value must be one non-missing character string")
  }
  bytes <- charToRaw(enc2utf8(value))
  writeBin(as.integer(length(bytes)), connection, size = 4L, endian = "little")
  if (length(bytes)) writeBin(bytes, connection)
  invisible(NULL)
}

perf_truba_write_atomic_vector <- function(connection, value) {
  type <- typeof(value)
  if (!type %in% c("integer", "double", "character") || anyNA(value) ||
      (type == "double" && any(!is.finite(value)))) {
    stop("Unsupported value in canonical portable identity")
  }
  perf_truba_write_utf8(connection, type)
  writeBin(as.integer(length(value)), connection, size = 4L, endian = "little")
  if (type == "integer") {
    writeBin(value, connection, size = 4L, endian = "little")
  } else if (type == "double") {
    writeBin(value, connection, size = 8L, endian = "little")
  } else {
    for (item in value) perf_truba_write_utf8(connection, item)
  }
  invisible(NULL)
}

perf_truba_tasks_hash <- function(tasks) {
  expected <- c("unit", "scenario", "repetition", "n", "p", "rho", "snr",
                "profile", "stream_index")
  if (!is.data.frame(tasks) || !identical(names(tasks), expected) ||
      nrow(tasks) < 1L || nrow(tasks) > .Machine$integer.max) {
    stop("Invalid task grid for canonical portable identity")
  }
  connection <- rawConnection(raw(), open = "wb")
  on.exit(close(connection), add = TRUE)
  perf_truba_write_utf8(connection, "performance_truba_tasks_canonical_v1")
  writeBin(as.integer(c(nrow(tasks), ncol(tasks))), connection,
           size = 4L, endian = "little")
  for (name in names(tasks)) {
    perf_truba_write_utf8(connection, name)
    perf_truba_write_atomic_vector(connection, tasks[[name]])
  }
  digest::digest(rawConnectionValue(connection), algo = "sha256", serialize = FALSE)
}

perf_truba_seeds_hash <- function(seeds) {
  if (!is.list(seeds) || !length(seeds) || length(seeds) > .Machine$integer.max) {
    stop("Invalid RNG stream list for canonical portable identity")
  }
  connection <- rawConnection(raw(), open = "wb")
  on.exit(close(connection), add = TRUE)
  perf_truba_write_utf8(connection, "performance_truba_seeds_canonical_v1")
  writeBin(as.integer(length(seeds)), connection, size = 4L, endian = "little")
  for (unit in seeds) {
    if (!is.list(unit) || length(unit) != 5L) {
      stop("Each work unit must contain five RNG substreams")
    }
    writeBin(as.integer(length(unit)), connection, size = 4L, endian = "little")
    for (state in unit) {
      if (!is.integer(state) || length(state) != 7L || anyNA(state)) {
        stop("Each RNG substream must be a seven-integer L'Ecuyer state")
      }
      perf_truba_write_atomic_vector(connection, state)
    }
  }
  digest::digest(rawConnectionValue(connection), algo = "sha256", serialize = FALSE)
}

perf_truba_portable_identity <- function(root, stage = "confirmation") {
  cfg <- perf_config(root, stage)
  tasks <- perf_tasks(cfg)
  seeds <- perf_seeds(cfg, tasks)
  list(
    schema = "performance_truba_portable_identity_v2",
    stage = stage,
    configuration_sha256 = digest::digest(perf_json(cfg), algo = "sha256",
                                           serialize = FALSE),
    tasks_sha256 = perf_truba_tasks_hash(tasks),
    seeds_sha256 = perf_truba_seeds_hash(seeds),
    scientific_sources = perf_sources(root),
    work_units = nrow(tasks),
    expected_method_rows = nrow(tasks) * length(cfg$methods),
    stream_range = range(tasks$stream_index),
    unique_streams = length(unique(tasks$stream_index))
  )
}

perf_truba_execution_spec <- function(root, stage) {
  list(
    schema = "performance_truba_execution_v1",
    release_id = perf_truba_release_id(),
    stage = stage,
    dispatcher = "parallel_mcparallel_dynamic_queue",
    worker_rng = "explicit_frozen_Lecuyer_states_mc_set_seed_false",
    queue_order = "repetition_then_scenario",
    shard_publication = "parent_validates_then_atomic_directory_rename",
    telemetry_writer = "single_parent_process",
    numerical_threads_per_worker = 1L,
    source_manifest = perf_truba_sources(root),
    portable_identity = perf_truba_portable_identity(root, stage)
  )
}

perf_truba_write_execution_metadata <- function(out, root, stage) {
  spec <- perf_truba_execution_spec(root, stage)
  rds <- file.path(out, "execution.rds")
  json <- file.path(out, "execution.json")
  signature <- file.path(out, "execution.sha256")
  perf_atomic(rds, function(path) saveRDS(spec, path, version = 3))
  perf_atomic(json, function(path) writeLines(perf_json(spec), path, useBytes = TRUE))
  perf_atomic(signature, function(path) writeLines(perf_hash(json), path, useBytes = TRUE))
  files <- c("execution.rds", "execution.json", "execution.sha256")
  manifest <- data.frame(
    file = files,
    sha256 = unname(vapply(file.path(out, files), perf_hash, character(1))),
    stringsAsFactors = FALSE
  )
  perf_atomic(file.path(out, "execution_manifest.csv"), function(path) {
    write.csv(manifest, path, row.names = FALSE)
  })
  invisible(spec)
}

perf_truba_verify_execution_metadata <- function(out, root, stage) {
  path <- file.path(out, "execution_manifest.csv")
  if (!file.exists(path)) stop("Missing execution metadata manifest")
  manifest <- read.csv(path, stringsAsFactors = FALSE)
  expected_files <- c("execution.rds", "execution.json", "execution.sha256")
  if (!identical(manifest$file, expected_files) || anyDuplicated(manifest$file) ||
      !identical(unname(vapply(file.path(out, manifest$file), perf_hash,
                               character(1))), manifest$sha256)) {
    stop("Execution metadata checksum mismatch")
  }
  stored <- readRDS(file.path(out, "execution.rds"))
  expected <- perf_truba_execution_spec(root, stage)
  if (!identical(perf_json(stored), perf_json(expected)) ||
      !identical(paste(readLines(file.path(out, "execution.json"), warn = FALSE),
                       collapse = "\n"), perf_json(stored))) {
    stop("Execution specification mismatch")
  }
  signature <- readLines(file.path(out, "execution.sha256"), warn = FALSE)
  if (length(signature) != 1L ||
      !identical(signature, perf_hash(file.path(out, "execution.json")))) {
    stop("Execution signature mismatch")
  }
  signature
}

perf_truba_verify_metadata_portable <- function(out, root) {
  manifest_path <- file.path(out, "metadata_manifest.csv")
  if (!file.exists(manifest_path)) stop("Missing scientific metadata manifest")
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  expected <- c("config.json", "config.rds", "sources.csv", "software.json",
                "tasks.csv", "seeds.rds", "signature.sha256", "session_info.txt")
  if (!identical(manifest$file, expected) || anyDuplicated(manifest$file) ||
      !identical(unname(vapply(file.path(out, manifest$file), perf_hash,
                               character(1))), manifest$sha256)) {
    stop("Scientific metadata checksum mismatch")
  }
  cfg <- readRDS(file.path(out, "config.rds"))
  if (!identical(perf_json(cfg), perf_json(perf_config(root, cfg$stage))) ||
      !identical(paste(readLines(file.path(out, "config.json"), warn = FALSE),
                       collapse = "\n"), perf_json(cfg))) {
    stop("Portable configuration mismatch")
  }
  sources <- read.csv(file.path(out, "sources.csv"), stringsAsFactors = FALSE)
  sources$bytes <- as.numeric(sources$bytes)
  if (!identical(sources, perf_sources(root))) stop("Portable source mismatch")
  software <- jsonlite::fromJSON(file.path(out, "software.json"), simplifyVector = FALSE)
  required <- c("glmnet", "ncvreg", "digest", "jsonlite", "knitr", "Rcpp")
  if (!all(required %in% names(software$packages)) ||
      !all(vapply(software$packages[required], function(x) length(x) == 1L &&
                    nzchar(as.character(x)), logical(1)))) {
    stop("Stored software record is incomplete")
  }
  thread_values <- unlist(software$threads, use.names = TRUE)
  if (length(thread_values) != 6L || any(thread_values != "1")) {
    stop("Stored numerical thread controls are invalid")
  }
  signature <- readLines(file.path(out, "signature.sha256"), warn = FALSE)
  if (length(signature) != 1L || !identical(signature, perf_signature(out))) {
    stop("Stored scientific signature mismatch")
  }
  tasks <- perf_tasks(cfg)
  stored_tasks <- read.csv(file.path(out, "tasks.csv"), stringsAsFactors = FALSE)
  if (!isTRUE(all.equal(tasks, stored_tasks, check.attributes = FALSE))) {
    stop("Portable task-grid mismatch")
  }
  seeds <- readRDS(file.path(out, "seeds.rds"))
  if (!identical(seeds, perf_seeds(cfg, tasks))) stop("Portable frozen RNG mismatch")
  list(cfg = cfg, tasks = tasks, seeds = seeds, signature = signature,
       stored_software = software)
}

perf_truba_receipt <- function(result, result_path) {
  list(
    schema = "performance_truba_shard_receipt_v1",
    unit = result$task$unit,
    scientific_signature = result$signature,
    execution_signature = result$execution_signature,
    result_bytes = as.numeric(file.info(result_path)$size),
    result_sha256 = perf_hash(result_path)
  )
}

perf_truba_write_unit <- function(out, result, state) {
  i <- result$task$unit
  perf_validate_unit(result, state$tasks[i, , drop = FALSE], state$cfg,
                     state$seeds[[i]], state$signature)
  if (!identical(result$execution_signature, state$execution_signature)) {
    stop("Result execution signature mismatch")
  }
  path <- perf_unit_path(out, i)
  if (file.exists(path)) stop("Refusing to overwrite a committed unit")
  partial <- paste0(path, ".partial-", Sys.getpid())
  if (!dir.create(partial)) stop("Cannot create atomic shard staging directory")
  committed <- FALSE
  on.exit(if (!committed && dir.exists(partial)) {
    preserved <- paste0(partial, ".failed-", format(Sys.time(), "%Y%m%d%H%M%S"))
    file.rename(partial, preserved)
  }, add = TRUE)
  result_path <- file.path(partial, "result.rds")
  saveRDS(result, result_path, version = 3, compress = "gzip")
  writeLines(perf_hash(result_path), file.path(partial, "result.sha256"),
             useBytes = TRUE)
  saveRDS(perf_truba_receipt(result, result_path),
          file.path(partial, "receipt.rds"), version = 3)
  if (!file.rename(partial, path)) stop("Cannot commit shard")
  committed <- TRUE
  invisible(path)
}

perf_truba_read_unit <- function(out, i, state) {
  path <- perf_unit_path(out, i)
  expected <- c("receipt.rds", "result.rds", "result.sha256")
  actual <- sort(list.files(path, all.files = TRUE, no.. = TRUE), method = "radix")
  if (!identical(actual, expected)) stop("Incomplete or unexpected shard files: ", path)
  result <- perf_read_unit(out, i, state)
  receipt <- readRDS(file.path(path, "receipt.rds"))
  wanted <- perf_truba_receipt(result, file.path(path, "result.rds"))
  if (!identical(receipt, wanted) ||
      !identical(result$execution_signature, state$execution_signature)) {
    stop("Shard receipt or execution signature mismatch: ", path)
  }
  result
}

perf_truba_read_attempts <- function(out, state) {
  directory <- file.path(out, "attempts")
  if (!dir.exists(directory)) stop("Missing attempt-record directory")
  files <- sort(list.files(directory, all.files = TRUE, no.. = TRUE), method = "radix")
  rds <- files[grepl("[.]rds$", files)]
  expected <- if (length(rds))
    sort(c(rds, paste0(rds, ".sha256")), method = "radix") else character()
  if ((length(files) || length(expected)) && !identical(files, expected)) {
    stop("Incomplete or unexpected attempt record")
  }
  attempts <- lapply(rds, function(name) {
    path <- file.path(directory, name)
    checksum <- readLines(paste0(path, ".sha256"), warn = FALSE)
    if (length(checksum) != 1L || !identical(checksum, perf_hash(path))) {
      stop("Attempt checksum mismatch: ", name)
    }
    item <- readRDS(path)
    if (!identical(item$schema, "performance_truba_attempt_v1") ||
        !identical(item$scientific_signature, state$signature) ||
        !identical(item$execution_signature, state$execution_signature) ||
        !identical(item$stage, state$cfg$stage) ||
        length(item$started_epoch) != 1L || !is.finite(item$started_epoch) ||
        length(item$finished_epoch) != 1L || !is.finite(item$finished_epoch) ||
        item$finished_epoch < item$started_epoch ||
        length(item$status) != 1L ||
        !item$status %in% c("paused", "failed", "complete") ||
        length(item$last_completed_unit) != 1L ||
        (!is.na(item$last_completed_unit) &&
         !item$last_completed_unit %in% state$tasks$unit)) {
      stop("Attempt identity mismatch: ", name)
    }
    item
  })
  attempts
}

perf_truba_write_attempt <- function(out, attempt) {
  path <- file.path(out, "attempts", paste0(attempt$invocation_id, ".rds"))
  perf_atomic(path, function(file) saveRDS(attempt, file, version = 3))
  perf_atomic(paste0(path, ".sha256"), function(file) {
    writeLines(perf_hash(path), file, useBytes = TRUE)
  })
  invisible(path)
}

perf_truba_scan <- function(out, root, complete = FALSE, portable = FALSE,
                            quiet = FALSE) {
  state <- if (portable) perf_truba_verify_metadata_portable(out, root) else
    perf_verify_metadata(out, root)
  state$execution_signature <- perf_truba_verify_execution_metadata(
    out, root, state$cfg$stage)
  state$attempts <- perf_truba_read_attempts(out, state)
  shard_root <- file.path(out, "shards")
  if (!dir.exists(shard_root)) stop("Missing shard directory")
  actual <- list.dirs(shard_root, recursive = FALSE, full.names = FALSE)
  partial <- actual[grepl("\\.partial-|\\.failed-", actual)]
  if (length(partial)) stop("Preserved partial shard requires inspection: ", partial[1L])
  present <- vapply(state$tasks$unit, function(i) dir.exists(perf_unit_path(out, i)),
                    logical(1))
  expected_dirs <- sprintf("unit_%05d", which(present))
  if (!setequal(actual, expected_dirs)) stop("Unexpected shard directories")
  results <- vector("list", sum(present))
  j <- 0L
  for (i in which(present)) {
    j <- j + 1L
    results[[j]] <- perf_truba_read_unit(out, i, state)
  }
  state$completed <- which(present)
  state$results <- results
  state$rows <- if (length(results)) do.call(rbind, lapply(results, `[[`, "rows")) else NULL
  if (complete) {
    if (!all(present) || !file.exists(file.path(out, "COMPLETED"))) stop("Run incomplete")
    if (length(list.files(file.path(out, "failures"), all.files = TRUE,
                          no.. = TRUE))) stop("Unresolved task failure evidence exists")
    perf_verify_manifest(out)
    aggregate <- readRDS(file.path(out, "results.rds"))
    if (!identical(aggregate, state$rows)) stop("Aggregate differs from validated shards")
    csv <- read.csv(file.path(out, "results.csv"), stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(csv, aggregate, tolerance = 1e-12,
                          check.attributes = FALSE))) stop("CSV differs from shard results")
    progress <- jsonlite::fromJSON(file.path(out, "progress.json"))
    if (!identical(progress$phase, "completed") ||
        progress$validated_completed != nrow(state$tasks) || progress$running != 0L ||
        progress$failed != 0L || progress$pending != 0L) {
      stop("Final parallel progress is inconsistent")
    }
    if (!identical(readLines(file.path(out, "COMPLETED"), warn = FALSE),
                   state$signature)) stop("Completion signature mismatch")
    if (!portable) perf_validate(out, root, complete = TRUE, quiet = TRUE)
    if (!quiet) cat(sprintf(
      paste0("COMPLETE: %d work units, %d method rows, parallel/CV/seed/source/",
             "manifest gates verified.\nFailed method fits retained: %d; partial searches: %d.\n"),
      nrow(state$tasks), nrow(aggregate), sum(!aggregate$ok),
      sum(aggregate$partial_search)))
  }
  invisible(state)
}

perf_truba_lock <- function(out) {
  path <- paste0(out, ".run_lock")
  if (!dir.create(path, showWarnings = FALSE)) {
    owner <- tryCatch(jsonlite::fromJSON(file.path(path, "owner.json")),
                      error = function(e) NULL)
    detail <- if (is.null(owner)) "unreadable owner" else paste(
      "host", owner$host, "pid", owner$pid, "job", owner$slurm_job_id)
    stop("Run lock exists (", detail, "). Verify the old job ended and preserve the lock before recovery.")
  }
  token <- paste(Sys.info()[["nodename"]], Sys.getpid(), Sys.time(), sep = ":")
  owner <- list(token = token, host = unname(Sys.info()[["nodename"]]),
                pid = Sys.getpid(), slurm_job_id = Sys.getenv("SLURM_JOB_ID", NA_character_),
                created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  writeLines(perf_json(owner), file.path(path, "owner.json"), useBytes = TRUE)
  list(path = path, token = token)
}

perf_truba_unlock <- function(lock) {
  if (is.null(lock) || !dir.exists(lock$path)) return(invisible(FALSE))
  owner <- tryCatch(jsonlite::fromJSON(file.path(lock$path, "owner.json")),
                    error = function(e) NULL)
  if (is.null(owner) || !identical(owner$token, lock$token)) {
    warning("Run-lock ownership changed; lock was preserved")
    return(invisible(FALSE))
  }
  unlink(lock$path, recursive = TRUE)
  invisible(!dir.exists(lock$path))
}

perf_truba_cost_class <- function(tasks) {
  with(tasks, paste(p, rho, snr, profile, sep = "/"))
}

perf_truba_progress_eta <- function(progress, now) {
  todo <- which(progress$status %in% c("running", "pending"))
  if (!length(todo)) return(0)
  observed <- which(is.finite(progress$durations) & progress$durations > 0)
  if (length(observed) < 2L) return(NA_real_)
  classes <- perf_truba_cost_class(progress$state$tasks)
  means <- tapply(progress$durations[observed], classes[observed], mean)
  costs <- as.numeric(means[classes[todo]])
  if (any(!is.finite(costs)) || any(costs <= 0)) return(NA_real_)
  running <- progress$status[todo] == "running"
  remaining <- costs[running] - (now - progress$task_started[todo[running]])
  if (any(!is.finite(remaining)) || any(remaining <= 0)) return(NA_real_)
  load <- c(remaining, rep(0, max(0L, progress$workers - length(remaining))))
  for (seconds in costs[!running]) {
    k <- which.min(load)
    load[k] <- load[k] + seconds
  }
  if (!length(load) || max(load) <= 0) NA_real_ else max(load)
}

perf_truba_progress_new <- function(out, state, workers, interval,
                                    now = as.numeric(Sys.time())) {
  progress <- new.env(parent = emptyenv())
  progress$out <- out
  progress$state <- state
  progress$workers <- workers
  progress$interval <- interval
  progress$started <- now
  previous_attempts <- state$attempts
  progress$previous_elapsed <- sum(vapply(previous_attempts, function(item) {
    if (is.null(item) || length(item$started_epoch) != 1L ||
        length(item$finished_epoch) != 1L || !is.finite(item$started_epoch) ||
        !is.finite(item$finished_epoch)) return(0)
    max(0, item$finished_epoch - item$started_epoch)
  }, numeric(1)))
  progress$last_emit <- -Inf
  progress$invocation_id <- paste0(
    format(as.POSIXct(now, origin = "1970-01-01", tz = "UTC"), "%Y%m%dT%H%M%S"),
    "_", Sys.getpid(), "_", Sys.getenv("SLURM_JOB_ID", "local"))
  progress$status <- rep("pending", nrow(state$tasks))
  progress$status[state$completed] <- "completed"
  progress$task_started <- rep(NA_real_, nrow(state$tasks))
  progress$durations <- rep(NA_real_, nrow(state$tasks))
  progress$new_completed <- integer()
  prior_latest <- if (length(previous_attempts)) {
    finished <- vapply(previous_attempts, function(item) item$finished_epoch, numeric(1))
    previous_attempts[[which.max(finished)]]$last_completed_unit
  } else NULL
  progress$latest <- if (length(prior_latest) == 1L && is.finite(prior_latest) &&
                         prior_latest %in% state$completed) as.integer(prior_latest) else
    if (length(state$completed)) {
      paths <- vapply(state$completed, function(i) file.path(
        perf_unit_path(out, i), "result.rds"), character(1))
      state$completed[which.max(as.numeric(file.info(paths)$mtime))]
    } else NA_integer_
  progress$phase <- "preflight"
  progress$detail <- ""
  progress
}

perf_truba_utc <- function(epoch) {
  if (length(epoch) != 1L || !is.finite(epoch)) return(NA_character_)
  format(as.POSIXct(epoch, origin = "1970-01-01", tz = "UTC"),
         "%Y-%m-%dT%H:%M:%SZ")
}

perf_truba_progress_emit <- function(progress, phase = progress$phase,
                                     detail = "", force = FALSE,
                                     now = as.numeric(Sys.time())) {
  if (!force && now - progress$last_emit < progress$interval) return(invisible(NULL))
  progress$last_emit <- now
  progress$phase <- phase
  progress$detail <- detail
  counts <- table(factor(progress$status,
                         levels = c("completed", "running", "failed", "pending")))
  total <- nrow(progress$state$tasks)
  if (sum(counts) != total) stop("Progress statuses do not partition work units")
  invocation_elapsed <- max(0, now - progress$started)
  elapsed <- progress$previous_elapsed + invocation_elapsed
  eta <- if (identical(phase, "computation")) perf_truba_progress_eta(progress, now) else
    if (identical(phase, "completed")) 0 else NA_real_
  deadline <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_END_TIME", NA_character_)))
  wall_left <- if (length(deadline) == 1L && is.finite(deadline))
    max(0, deadline - now) else NA_real_
  running_units <- which(progress$status == "running")
  snapshot <- list(
    schema = "performance_truba_progress_v1",
    run_id = basename(progress$out),
    scientific_signature = progress$state$signature,
    execution_signature = progress$state$execution_signature,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", NA_character_),
    invocation_id = progress$invocation_id,
    phase = phase,
    work_unit = "one scenario/replication including all methods, CV and validation",
    total = total,
    validated_completed = unname(counts[["completed"]]),
    running = unname(counts[["running"]]),
    failed = unname(counts[["failed"]]),
    pending = unname(counts[["pending"]]),
    percent_completed = 100 * unname(counts[["completed"]]) / total,
    elapsed_seconds = elapsed,
    invocation_elapsed_seconds = invocation_elapsed,
    previous_recorded_elapsed_seconds = progress$previous_elapsed,
    completed_this_invocation = length(progress$new_completed),
    throughput_units_per_second = if (invocation_elapsed > 0 && length(progress$new_completed))
      length(progress$new_completed) / invocation_elapsed else NA_real_,
    estimated_remaining_seconds = eta,
    estimated_completion_utc = perf_truba_utc(now + eta),
    eta_status = if (is.finite(eta)) "estimated" else "estimating_or_unknown",
    eta_scope = "computation_and_per_unit_validation_only; aggregation_and_final_verification_excluded",
    slurm_remaining_seconds = wall_left,
    last_heartbeat_utc = perf_truba_utc(now),
    current_work_units = if (length(running_units))
      paste(running_units, collapse = ",") else NA_character_,
    most_recently_completed_unit = progress$latest,
    detail = detail
  )
  perf_atomic_json(snapshot, file.path(progress$out, "progress.json"))
  history <- file.path(progress$out, "progress.tsv")
  write.table(as.data.frame(snapshot, stringsAsFactors = FALSE), history,
              sep = "\t", row.names = FALSE, col.names = !file.exists(history),
              append = file.exists(history), quote = TRUE, na = "NA")
  cat(sprintf(
    "[%s] %s | %d/%d (%.1f%%) | running %d failed %d | ETA %s | wall %s | %s\n",
    snapshot$last_heartbeat_utc, phase, snapshot$validated_completed, total,
    snapshot$percent_completed, snapshot$running, snapshot$failed,
    if (is.finite(eta)) paste0(round(eta), " s, computation only") else "estimating/unknown",
    if (is.finite(wall_left)) paste0(round(wall_left), " s remaining") else "unknown",
    detail))
  flush.console()
  invisible(snapshot)
}

perf_truba_failure_record <- function(out, progress, unit, item) {
  dir.create(file.path(out, "failures"), recursive = TRUE, showWarnings = FALSE)
  stem <- sprintf("%s_unit_%05d", progress$invocation_id, unit)
  record <- list(
    schema = "performance_truba_failure_v1",
    invocation_id = progress$invocation_id,
    unit = unit,
    task = progress$state$tasks[unit, , drop = FALSE],
    scientific_signature = progress$state$signature,
    execution_signature = progress$state$execution_signature,
    status = item$status,
    message = item$message,
    captured = item,
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
  perf_atomic(file.path(out, "failures", paste0(stem, ".rds")), function(path) {
    saveRDS(record, path, version = 3)
  })
  perf_atomic_json(record[c("schema", "invocation_id", "unit", "task",
                            "scientific_signature", "execution_signature",
                            "status", "message", "created_utc")],
                   file.path(out, "failures", paste0(stem, ".json")))
  invisible(record)
}

perf_truba_dispatch <- function(out, state, workers, max_seconds, task_limit,
                                interval) {
  progress <- perf_truba_progress_new(out, state, workers, interval)
  progress$phase <- "computation"
  perf_truba_progress_emit(progress, detail = "dispatcher initialized", force = TRUE)
  queue <- state$tasks$unit[order(state$tasks$repetition, state$tasks$scenario)]
  jobs <- list()
  outcomes <- list()
  launched <- 0L
  halt <- FALSE
  cleaned <- FALSE
  cleanup <- function() {
    if (cleaned) return(invisible(NULL))
    for (job in jobs) try(tools::pskill(job$pid), silent = TRUE)
    if (length(jobs)) try(parallel::mccollect(jobs, wait = FALSE), silent = TRUE)
    cleaned <<- TRUE
    invisible(NULL)
  }
  on.exit(cleanup(), add = TRUE)
  worker <- function(i) {
    tryCatch({
      result <- perf_run_task(state$tasks[i, , drop = FALSE], state$cfg,
                              state$seeds[[i]], function(...) invisible(NULL))
      list(status = "computed", result = result, message = "")
    }, error = function(e) list(status = "worker_error", result = NULL,
                                message = conditionMessage(e)))
  }
  accept <- function(i, item) {
    if (!identical(item$status, "computed") || is.null(item$result)) return(item)
    tryCatch({
      item$result$signature <- state$signature
      item$result$execution_signature <- state$execution_signature
      perf_truba_write_unit(out, item$result, state)
      readback <- perf_truba_read_unit(out, i, state)
      list(status = "passed", message = "", elapsed = readback$elapsed)
    }, error = function(e) list(status = "accept_error", result = item$result,
                                message = conditionMessage(e)))
  }
  repeat {
    now <- as.numeric(Sys.time())
    ready <- queue[progress$status[queue] == "pending"]
    within_budget <- now - progress$started < max_seconds
    while (!halt && within_budget && length(ready) && length(jobs) < workers &&
           launched < task_limit) {
      i <- ready[1L]
      ready <- ready[-1L]
      progress$status[i] <- "running"
      progress$task_started[i] <- now
      launched <- launched + 1L
      job <- parallel::mcparallel(local({ unit <- i; worker(unit) }),
                                  mc.set.seed = FALSE, silent = TRUE)
      attr(job, "unit") <- i
      jobs[[as.character(job$pid)]] <- job
      perf_truba_progress_emit(progress, detail = paste("starting work unit", i),
                               force = TRUE)
      now <- as.numeric(Sys.time())
      within_budget <- now - progress$started < max_seconds
    }
    if (!length(jobs)) break
    done <- parallel::mccollect(jobs, wait = FALSE)
    if (length(done)) for (pid in names(done)) {
      job <- jobs[[pid]]
      i <- attr(job, "unit")
      jobs[[pid]] <- NULL
      item <- done[[pid]]
      if (inherits(item, "try-error") || !is.list(item) || is.null(item$status)) {
        item <- list(status = "worker_error", result = NULL,
                     message = paste(as.character(item), collapse = " | "))
      }
      item <- accept(i, item)
      outcomes[[as.character(i)]] <- item[c("status", "message", "elapsed")]
      if (identical(item$status, "passed")) {
        progress$status[i] <- "completed"
        progress$durations[i] <- item$elapsed
        progress$new_completed <- c(progress$new_completed, i)
        progress$latest <- i
        perf_truba_progress_emit(progress, detail = paste("work unit validated", i),
                                 force = TRUE)
      } else {
        progress$status[i] <- "failed"
        perf_truba_failure_record(out, progress, i, item)
        halt <- TRUE
        perf_truba_progress_emit(progress, detail = paste("work unit failed", i,
                                                          item$message), force = TRUE)
      }
    }
    perf_truba_progress_emit(progress, detail = "workers active")
    if (length(jobs)) Sys.sleep(min(1, max(0.05, interval / 5)))
  }
  cleaned <- TRUE
  status <- if (any(progress$status == "failed")) "failed" else
    if (all(progress$status == "completed")) "complete" else "paused"
  list(status = status, progress = progress, outcomes = outcomes,
       launched = launched, elapsed = as.numeric(Sys.time()) - progress$started)
}

perf_truba_publish_aggregate <- function(out, rows) {
  rds <- file.path(out, "results.rds")
  csv <- file.path(out, "results.csv")
  if (file.exists(rds)) {
    if (!identical(readRDS(rds), rows)) stop("Existing aggregate RDS differs from shards")
  } else perf_atomic(rds, function(path) saveRDS(rows, path, version = 3))
  if (file.exists(csv)) {
    stored <- read.csv(csv, stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(stored, rows, tolerance = 1e-12,
                          check.attributes = FALSE))) stop("Existing aggregate CSV differs")
  } else perf_atomic(csv, function(path) write.csv(rows, path, row.names = FALSE))
  invisible(TRUE)
}

perf_truba_recover_manifest_completion_gap <- function(out, root, state) {
  if (!file.exists(file.path(out, "manifest.csv")) ||
      file.exists(file.path(out, "COMPLETED"))) return(FALSE)
  if (length(state$completed) != nrow(state$tasks)) {
    stop("Final manifest exists before all work units were committed")
  }
  perf_verify_manifest(out)
  if (!file.exists(file.path(out, "results.rds")) ||
      !file.exists(file.path(out, "results.csv"))) {
    stop("Final manifest exists without both aggregate results")
  }
  aggregate <- readRDS(file.path(out, "results.rds"))
  csv <- read.csv(file.path(out, "results.csv"), stringsAsFactors = FALSE)
  if (!identical(aggregate, state$rows) ||
      !isTRUE(all.equal(csv, aggregate, tolerance = 1e-12,
                        check.attributes = FALSE))) {
    stop("Manifest-gap aggregate differs from validated shards")
  }
  progress <- jsonlite::fromJSON(file.path(out, "progress.json"))
  if (!identical(progress$phase, "completed") ||
      progress$validated_completed != nrow(state$tasks) || progress$running != 0L ||
      progress$failed != 0L || progress$pending != 0L) {
    stop("Manifest-gap progress is not a completed state")
  }
  perf_atomic(file.path(out, "COMPLETED"), function(path) {
    writeLines(state$signature, path, useBytes = TRUE)
  })
  perf_truba_scan(out, root, complete = TRUE, quiet = TRUE)
  TRUE
}

perf_truba_operational_guard <- function(stage, workers, max_seconds, task_limit,
                                         interval, require_slurm = stage == "confirmation") {
  if (.Platform$OS.type != "unix") stop("The fork dispatcher requires Unix")
  if (length(workers) != 1L || !is.finite(workers) || workers != floor(workers) ||
      workers < 1L || workers > if (stage == "smoke") 2L else 56L) {
    stop("Invalid task-worker count")
  }
  if (stage == "confirmation" && workers != 56L) {
    stop("Confirmation requires exactly 56 task workers")
  }
  ceiling <- if (stage == "smoke") 300 else 21600
  if (length(max_seconds) != 1L || !is.finite(max_seconds) || max_seconds <= 0 ||
      max_seconds > ceiling) stop("Invalid soft dispatch budget")
  if (length(task_limit) != 1L || is.na(task_limit) || task_limit < 0 ||
      (!is.infinite(task_limit) && task_limit != floor(task_limit))) {
    stop("Invalid task limit")
  }
  if (stage == "confirmation" && !is.infinite(task_limit)) {
    stop("Confirmation task limit must be infinite")
  }
  if (length(interval) != 1L || !is.finite(interval) || interval <= 0 || interval > 60) {
    stop("Heartbeat interval must be positive and at most 60 seconds")
  }
  if (require_slurm) {
    if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Confirmation run requires a submitted SLURM job")
    if (!identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "56")) {
      stop("Confirmation run requires exactly 56 allocated CPUs")
    }
  }
  invisible(TRUE)
}

perf_truba_runtime_preflight <- function(root, stage, require_release = FALSE,
                                         require_truba_runtime = FALSE) {
  perf_preflight(compile = TRUE, root = root)
  if (require_truba_runtime) {
    if (!identical(as.character(getRversion()), "4.3.0") ||
        !grepl("linux", R.version$platform, ignore.case = TRUE)) {
      stop("TRUBA execution requires the Linux R 4.3.0 module")
    }
    wanted_library <- Sys.getenv("R_LIBS_USER")
    if (!identical(Sys.getenv("R_LIBS_USER"), wanted_library) ||
        !wanted_library %in% .libPaths()) stop("Unexpected TRUBA R user library")
  }
  source(file.path(root, "tests/performance_v1_checks.R"), local = .GlobalEnv)
  checks <- perf_checks(root)
  if (!all(checks$passed)) stop("Small scientific/API preflight failed")
  if (require_release) perf_truba_validate_release(root)
  software <- perf_software()
  cat("Runtime preflight accepted; model fits in confirmation workload: 0\n")
  cat(software$R, "\n", software$platform, "\n")
  print(unlist(software$packages))
  cat("Library paths:\n", paste(.libPaths(), collapse = "\n"), "\n")
  invisible(list(checks = checks, software = software))
}

perf_truba_validate_release <- function(root, quiet = FALSE) {
  path <- file.path(root, "release", perf_truba_release_id(), "LOCAL_VALIDATED.rds")
  if (!file.exists(path)) stop("Missing local validated TRUBA release receipt")
  receipt <- readRDS(path)
  if (!identical(receipt$schema, "performance_confirmation_truba_release_v2") ||
      !isTRUE(receipt$accepted) || !nrow(receipt$checks) ||
      !all(receipt$checks$passed %in% TRUE)) stop("Invalid TRUBA release receipt")
  if (!identical(receipt$sources, perf_truba_sources(root))) {
    stop("TRUBA release source hashes changed")
  }
  if (!identical(receipt$portable_identity,
                 perf_truba_portable_identity(root, "confirmation"))) {
    stop("Released confirmation design, sources or seeds changed")
  }
  if (!identical(receipt$evidence,
                 perf_truba_inventory(root, receipt$evidence$file))) {
    stop("Local validation evidence checksum mismatch")
  }
  if (!quiet) cat("Frozen local TRUBA release verified: TRUE; production fits: 0\n")
  invisible(receipt)
}

perf_truba_run <- function(root, out, stage, workers, max_seconds, task_limit = Inf,
                           interval = 30, resume = NULL) {
  stopifnot(stage %in% c("smoke", "confirmation"))
  if (is.null(resume)) resume <- dir.exists(out)
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) stop("Invalid resume mode")
  if (dir.exists(out) && !resume) stop("Output exists; refusing overwrite")
  if (!dir.exists(out) && resume) stop("Cannot resume a missing output")
  perf_truba_operational_guard(stage, workers, max_seconds, task_limit, interval)
  if (stage == "confirmation") perf_truba_validate_release(root)
  perf_truba_runtime_preflight(root, stage, require_release = stage == "confirmation",
                               require_truba_runtime = stage == "confirmation")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  lock <- perf_truba_lock(out)
  on.exit(perf_truba_unlock(lock), add = TRUE)
  if (!resume) {
    if (!dir.create(out) || !dir.create(file.path(out, "shards")) ||
        !dir.create(file.path(out, "attempts")) ||
        !dir.create(file.path(out, "failures"))) stop("Cannot create fresh output structure")
    cfg <- perf_config(root, stage)
    tasks <- perf_tasks(cfg)
    seeds <- perf_seeds(cfg, tasks)
    perf_metadata(out, root, cfg, tasks, seeds)
    perf_truba_write_execution_metadata(out, root, stage)
  }
  failure_files <- list.files(file.path(out, "failures"), all.files = TRUE, no.. = TRUE)
  if (length(failure_files)) stop("Unresolved failure evidence blocks automatic resume")
  state <- perf_truba_scan(out, root, complete = FALSE)
  if (file.exists(file.path(out, "COMPLETED"))) {
    perf_truba_scan(out, root, complete = TRUE)
    return(invisible(list(complete = TRUE, output = out)))
  }
  if (perf_truba_recover_manifest_completion_gap(out, root, state)) {
    perf_truba_scan(out, root, complete = TRUE)
    return(invisible(list(complete = TRUE, output = out,
                          recovered_manifest_completion_gap = TRUE)))
  }
  dispatched <- perf_truba_dispatch(out, state, workers, max_seconds, task_limit, interval)
  attempt <- list(
    schema = "performance_truba_attempt_v1",
    invocation_id = dispatched$progress$invocation_id,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", NA_character_),
    scientific_signature = state$signature,
    execution_signature = state$execution_signature,
    stage = stage,
    workers = workers,
    task_limit = task_limit,
    max_seconds_soft = max_seconds,
    interval_seconds = interval,
    launched = dispatched$launched,
    status = dispatched$status,
    last_completed_unit = dispatched$progress$latest,
    outcomes = dispatched$outcomes,
    started_epoch = dispatched$progress$started,
    finished_epoch = as.numeric(Sys.time()),
    runtime = perf_software(),
    session = capture.output(sessionInfo())
  )
  perf_truba_write_attempt(out, attempt)
  if (identical(dispatched$status, "failed")) {
    perf_truba_progress_emit(dispatched$progress, "failed",
                             "worker or shard validation failure; evidence preserved", TRUE)
    stop("Parallel confirmation stopped after a worker or shard validation failure")
  }
  if (!identical(dispatched$status, "complete")) {
    perf_truba_progress_emit(dispatched$progress, "checkpointed",
                             "soft stop or task limit reached; valid shards preserved", TRUE)
    return(invisible(list(complete = FALSE, output = out,
                          completed = sum(dispatched$progress$status == "completed"))))
  }
  perf_truba_progress_emit(dispatched$progress, "aggregation_and_validation",
                           "reading all validated shards", TRUE)
  checked <- perf_truba_scan(out, root, complete = FALSE)
  perf_truba_publish_aggregate(out, checked$rows)
  perf_truba_progress_emit(dispatched$progress, "completed",
                           "all shards aggregated and validated", TRUE)
  if (file.exists(file.path(out, "manifest.csv"))) {
    perf_verify_manifest(out)
  } else {
    perf_manifest(out)
    perf_verify_manifest(out)
  }
  if (!file.exists(file.path(out, "COMPLETED"))) {
    perf_atomic(file.path(out, "COMPLETED"), function(path) {
      writeLines(state$signature, path, useBytes = TRUE)
    })
  }
  perf_truba_scan(out, root, complete = TRUE)
  invisible(list(complete = TRUE, output = out))
}

perf_truba_main <- function(root) {
  args <- commandArgs(trailingOnly = TRUE)
  if (any(args %in% c("--help", "-h"))) {
    cat(paste0(
      "--stage=smoke|confirmation --action=design|preflight|run|verify|portable-verify\n",
      "--out=PATH --workers=N --max-seconds=N --task-limit=N|Inf --interval=N\n"))
    return(invisible(NULL))
  }
  if (length(args) && any(!grepl("^--[a-z-]+=.+$", args))) stop("Use --key=value arguments")
  keys <- sub("^--([^=]+)=.*$", "\\1", args)
  allowed <- c("stage", "action", "out", "workers", "max-seconds", "task-limit",
               "interval")
  if (anyDuplicated(keys) || any(!keys %in% allowed)) stop("Unknown or duplicate option")
  opts <- setNames(as.list(sub("^--[^=]+=", "", args)), keys)
  get <- function(name, fallback) if (is.null(opts[[name]])) fallback else opts[[name]]
  stage <- match.arg(get("stage", "confirmation"), c("smoke", "confirmation"))
  action <- match.arg(get("action", "design"),
                      c("design", "preflight", "run", "verify", "portable-verify"))
  if (action == "design") {
    cfg <- perf_config(root, stage)
    tasks <- perf_tasks(cfg)
    cat("DESIGN ONLY: zero model fits. Stage:", stage, "\n")
    cat("Tasks:", nrow(tasks), "Expected method rows:",
        nrow(tasks) * length(cfg$methods), "\n")
    cat("Streams:", paste(range(tasks$stream_index), collapse = "-"),
        "duplicates:", anyDuplicated(tasks$stream_index), "\n")
    return(invisible(NULL))
  }
  if (action == "preflight") {
    perf_truba_runtime_preflight(
      root, stage, require_release = stage == "confirmation",
      require_truba_runtime = identical(Sys.getenv("PERF_REQUIRE_TRUBA_RUNTIME"), "1"))
    return(invisible(NULL))
  }
  out_value <- get("out", if (stage == "confirmation")
    "outputs/performance/confirmation_truba_v1" else
      "outputs/performance/smoke_truba_v1")
  out <- if (startsWith(out_value, "/")) out_value else file.path(root, out_value)
  if (action == "verify") return(perf_truba_scan(out, root, complete = TRUE))
  if (action == "portable-verify") return(perf_truba_scan(
    out, root, complete = TRUE, portable = TRUE))
  workers <- as.numeric(get("workers", if (stage == "confirmation") "56" else "1"))
  max_seconds <- as.numeric(get("max-seconds", if (stage == "confirmation") "21600" else "300"))
  task_limit <- as.numeric(get("task-limit", "Inf"))
  interval <- as.numeric(get("interval", "30"))
  result <- perf_truba_run(root, out, stage, workers, max_seconds, task_limit, interval)
  if (identical(result$complete, FALSE)) {
    cat("Checkpointed incomplete run; resubmit the identical command to resume.\n")
    quit(status = 3L)
  }
  invisible(result)
}
