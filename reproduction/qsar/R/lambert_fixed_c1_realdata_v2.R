# Default-shape comparator overlay; no frozen v1 function is edited.
lf2_environment <- function(root) {
  env <- lf1_environment(root)
  old_engine <- env$rcv1_engine; old_config <- env$rcv1_config
  old_sources <- env$rcv1_source_paths
  env$rcv1_engine <- function(root, dataset) {
    engine <- old_engine(root, dataset); base_config <- engine$qsar_v1_config
    engine$qsar_v1_config <- function(root, mode) {
      cfg <- base_config(root, mode)
      cfg$schema <- paste0(cfg$schema, "_default_comparators_v2")
      cfg$mcp_gamma <- 3; cfg$scad_a <- 3.7
      cfg$comparator_shapes <- list(source = "ncvreg_3.16.0_defaults", mcp = 3, scad = 3.7)
      cfg
    }
    original_validate <- engine$qsar_v1_validate_unit
    engine$qsar_v1_validate_unit <- function(result, unit, state) {
      original_validate(result, unit, state)
      cv <- result$fit$cv
      stopifnot(identical(state$cfg$mcp_gamma, 3), identical(state$cfg$scad_a, 3.7),
        all(cv$shape[cv$family == "mcp"] == 3), all(cv$shape[cv$family == "scad"] == 3.7))
      for (m in c("mcp", "scad")) {
        row <- result$rows[result$rows$method == m, ]
        if (row$ok) stopifnot(row$shape == if (m == "mcp") 3 else 3.7)
      }
      invisible(TRUE)
    }
    engine
  }
  env$rcv1_config <- function(root, dataset, mode, engine = env$rcv1_engine(root, dataset)) {
    cfg <- old_config(root, dataset, mode, engine)
    cfg$schema <- "lambert_fixed_c1_realdata_v2_default_comparators"
    cfg
  }
  env$rcv1_source_paths <- function(engine) unique(c(old_sources(engine),
    "R/lambert_fixed_c1_realdata_v2.R", "scripts/lambert_fixed_c1_realdata_v2.R",
    "tests/lambert_fixed_c1_v2_checks.R", "protocol/lambert_fixed_c1_v2.md"))
  env
}
