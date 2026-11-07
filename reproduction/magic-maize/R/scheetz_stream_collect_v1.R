# Reuse the frozen collector's scientific expressions, with bounded record IO.
schm_legacy_collect <- sch_collect
schm_collect_parts <- local({
  b <- as.list(body(schm_legacy_collect))
  stopifnot(length(b)==14L,identical(b[[4]][[2]],quote(records)),
    identical(b[[8]],quote(rows <- sch_bind(rows))),
    identical(b[[9]],quote(predictions <- sch_bind(predictions))))
  group <- function(records,e)NULL
  body(group) <- as.call(c(list(as.name("{")),b[5:9],list(quote(list(
    results=rows,predictions=predictions,path_summary=sch_bind(diagnostics),
    candidate_diagnostics=sch_bind(candidate_rows))))))
  summary <- function(rows,predictions,diagnostics,candidate_rows,e)NULL
  body(summary) <- as.call(c(list(as.name("{")),b[8:14]))
  list(group=group,summary=summary)
})
sch_collect <- function(out,e,portable=FALSE,pulse=function(...)NULL) {
  state <- sch_scan(out,e,portable=portable,pulse=pulse)
  stopifnot(identical(state$completed,e$tasks$unit))
  chunks <- list();pending <- list();seen <- integer()
  for(batch in state$batches) {
    pack <- readRDS(file.path(out,cpk_pack_file(batch)))
    for(i in seq_along(pack$records)) {
      z <- pack$records[[i]];id <- z$task$unit;pending[[as.character(id)]] <- z
      ids <- e$tasks$unit[e$tasks$repeat_id==z$task$repeat_id & e$tasks$outer_fold==z$task$outer_fold]
      if(all(as.character(ids)%in%names(pending))) {
        part <- list2env(as.list(e),parent=parent.env(e))
        part$repeats <- z$task$repeat_id;part$active_outer <- z$task$outer_fold
        chunks[[length(chunks)+1L]] <- schm_collect_parts$group(pending[as.character(ids)],part)
        seen <- c(seen,ids);pending[as.character(ids)] <- NULL
      }
    }
    rm(pack,z);gc(verbose=FALSE)
  }
  stopifnot(!length(pending),identical(seen,e$tasks$unit))
  schm_collect_parts$summary(lapply(chunks,`[[`,"results"),lapply(chunks,`[[`,"predictions"),
    lapply(chunks,`[[`,"path_summary"),lapply(chunks,`[[`,"candidate_diagnostics"),e)
}
