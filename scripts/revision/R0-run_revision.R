## ---------------------------
##
## Script name: R0-run_revision.R
## Purpose of script: To run the R1 revision pipeline (unsupervised PCA and
##                    cross-fitted mixture indices) end to end
##
## Author: Yufan Gong
##
## Date Created: 2026-08-27
##
## Date Modified: 2026-08-27
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: Runs the revision analogues of scripts 3-7 in order:
##
##   R3-composites_revision.R    build the exposures
##   R4-mwas_revision.R          limma + PLS MWAS
##   R5-annotation_revision.R    attach compound annotations
##   R6-pathway_revision.R       mummichog + metapone
##   R7-visualization_revision.R figures
##
##        Each script runs in its OWN R process, because they all source
##        scripts/1-functions.R (which calls rm(list = ls())) and several hold
##        multi-GB metabolite matrices. One script failing does not abort the
##        run; the summary printed at the end says what failed and where the
##        log is.
##
##        From the shell:
##          Rscript scripts/revision/R0-run_revision.R
##          Rscript scripts/revision/R0-run_revision.R R4,R5
##
##        Smoke test first (minutes, NOT reportable):
##          SALSA_REVISION_QUICK=true \
##          SALSA_REVISION_ROOT=revision_output_smoke \
##          Rscript scripts/revision/R0-run_revision.R
##
##        From inside R:
##          STEPS <- c("R3", "R4")
##          source(here::here("scripts", "revision", "R0-run_revision.R"))
##
## ---------------------------

library(here)

REV_ROOT <- Sys.getenv("SALSA_REVISION_ROOT", unset = "revision_output")

revision_scripts <- c(
  R3 = "R3-composites_revision.R",
  R4 = "R4-mwas_revision.R",
  R5 = "R5-annotation_revision.R",
  R6 = "R6-pathway_revision.R",
  R7 = "R7-visualization_revision.R"
)

## Which steps to run: command line argument, the STEPS object, or all of them
if (!exists("STEPS")) {
  args <- commandArgs(trailingOnly = TRUE)
  STEPS <- if (length(args) > 0) {
    trimws(unlist(strsplit(args[1], ",")))
  } else {
    names(revision_scripts)
  }
}

unknown <- setdiff(STEPS, names(revision_scripts))
if (length(unknown) > 0) {
  stop("Unknown step(s): ", paste(unknown, collapse = ", "),
       ". Valid steps: ", paste(names(revision_scripts), collapse = ", "))
}

log_dir <- here::here(REV_ROOT, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

message(strrep("=", 78))
message("  SALSA air toxicants - R1 revision pipeline")
message("  output root : ", here::here(REV_ROOT))
message("  steps       : ", paste(STEPS, collapse = ", "))
message("  quick mode  : ", Sys.getenv("SALSA_REVISION_QUICK", unset = "false"))
message(strrep("=", 78))

run_summary <- STEPS |>
  lapply(function(step) {
    script <- revision_scripts[[step]]
    path   <- here::here("scripts", "revision", script)
    log    <- file.path(log_dir, paste0(stamp, "_", step, ".log"))

    message("\n>>> ", step, ": ", script)
    message("    log -> ", log)

    started <- Sys.time()
    status <- system2(
      file.path(R.home("bin"), "Rscript"),
      args = c("--vanilla", shQuote(path)),
      stdout = log, stderr = log
    )
    elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)

    message("    ", if (status == 0) "OK" else paste0("FAILED (status ",
                                                      status, ")"),
            " in ", elapsed, " min")

    data.frame(step = step, script = script, status = status,
               ok = status == 0, minutes = elapsed, log = log)
  }) |>
  do.call(what = rbind)

message("\n", strrep("=", 78))
message("  Run summary")
message(strrep("=", 78))
print(run_summary[, c("step", "ok", "minutes")], row.names = FALSE)

summary_path <- file.path(log_dir, paste0(stamp, "_run_summary.csv"))
utils::write.csv(run_summary, summary_path, row.names = FALSE)
message("\nSummary written to ", summary_path)

if (any(!run_summary$ok)) {
  message("\nFailed steps: ",
          paste(run_summary$step[!run_summary$ok], collapse = ", "))
  message("Check the logs listed above.")
}

#--------------------------------End of the code--------------------------------
