custom_format <- function(x) {
  ifelse(abs(x) < 0.01,
         format(x, scientific = TRUE, digits = 2),
         format(round(x, 2), nsmall = 2))
}

run_qgcomp_boot_parallel <- function(data,
                                     exposures_list,
                                     outcomes_list,
                                     covariates_list,
                                     q = 4, B = 200, seed = 42,
                                     family = NULL,
                                     rr = FALSE,         # used only for binomial
                                     workers = NULL,
                                     quiet = FALSE,
                                     progress_handler = c("txtprogressbar", "cli")) {
  
  # Dependencies (fail fast with clear message)
  pkgs <- c("dplyr", "purrr", "tibble", "qgcomp", "future", "furrr", "progressr")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Please install required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  
  progress_handler <- match.arg(progress_handler)
  
  stopifnot(is.data.frame(data))
  stopifnot(length(exposures_list) > 0, length(outcomes_list) > 0)
  
  # Validate columns
  all_needed <- unique(c(outcomes_list, exposures_list, covariates_list))
  missing_cols <- setdiff(all_needed, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in `data`: ", paste(missing_cols, collapse = ", "))
  }
  
  # ---- automatic family detection ----
  detect_family <- function(y) {
    y_nonmiss <- y[!is.na(y)]
    uniq <- unique(y_nonmiss)
    
    if (is.logical(y_nonmiss) ||
        is.factor(y_nonmiss) && length(uniq) == 2 ||
        (length(uniq) <= 2 && all(sort(uniq) %in% c(0,1)))) {
      return(stats::binomial())
    } else {
      return(stats::gaussian())
    }
  }

  # Reduce copying: keep only needed cols once
  data_small <- dplyr::select(data, dplyr::all_of(all_needed))
  
  # Shared RHS
  rhs <- paste(c(exposures_list, covariates_list), collapse = " + ")
  
  fit_one <- function(outcome, this_seed) {
    form <- stats::as.formula(paste(outcome, "~", rhs))
    
    dat_cc <- data_small |>
      dplyr::select(dplyr::all_of(c(outcome, exposures_list, covariates_list))) |>
      dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.x)))
    
    # Determine family
    fam <- if (is.null(family)) {
      detect_family(dat_cc[[outcome]])
    } else {
      family
    }
    
    if (!quiet) {
      message("Fitting qgcomp.boot model for: ", outcome, deparse(form), 
              " with family = ", fam$family)
    }
    
    
    model <- tryCatch({
      if (fam$family == "binomial") {
        qgcomp::qgcomp.glm.boot(
          f = form,
          data = dat_cc,
          expnms = exposures_list,
          q = q,
          B = B,
          seed   = this_seed,   # IMPORTANT: distinct per outcome
          family = fam,
          rr = rr
        )
      } else {
        qgcomp::qgcomp.glm.boot(
          f      = form,
          data   = dat_cc,
          expnms = exposures_list,
          q      = q,
          B      = B,
          seed   = this_seed,   # IMPORTANT: distinct per outcome
          family = fam
        )
      }
    }, error = function(e) e)
    
    if (inherits(model, "error")) {
      if (!quiet) message("Model failed for outcome ", outcome, ": ", model$message)
      return(list(
        result = tibble::tibble(
          Outcome     = outcome,
          N           = NA_integer_,
          Est_CI      = NA_character_,
          P_value_raw = NA_real_,
          Est         = NA_real_,
          Conf_low    = NA_real_,
          Conf_high   = NA_real_
        ),
        model = NULL
      ))
    }
    
    coefs <- summary(model)$coefficients
    if (!("psi1" %in% rownames(coefs))) {
      if (!quiet) message("Outcome ", outcome, ": 'psi1' not found in coefficients.")
      return(list(
        result = tibble::tibble(
          Outcome     = outcome,
          N           = nrow(dat_cc),
          Est_CI      = NA_character_,
          P_value_raw = NA_real_,
          Est         = NA_real_,
          Conf_low    = NA_real_,
          Conf_high   = NA_real_
        ),
        model = model
      ))
    }
    
    res <- coefs["psi1", , drop = TRUE]
    est       <- unname(res[["Estimate"]])
    conf_low  <- unname(res[["Lower CI"]])
    conf_high <- unname(res[["Upper CI"]])
    pval      <- unname(res[["Pr(>|z|)"]])
    
    est_ci <- sprintf("%.2f (%.2f, %.2f)", est, conf_low, conf_high)
    
    list(
      result = tibble::tibble(
        Outcome     = outcome,
        N           = nrow(dat_cc),
        Est_CI      = est_ci,
        P_value_raw = pval,
        Est         = est,
        Conf_low    = conf_low,
        Conf_high   = conf_high
      ),
      model = model
    )
  }
  
  # Decide workers (leave one core free by default)
  if (is.null(workers)) {
    avail <- future::availableCores()
    workers <- max(1, avail - 1)
  }
  
  # Set up future plan (restore after)
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  
  # Set progress handler (global on; restore afterwards)
  old_handlers <- progressr::handlers()
  on.exit(progressr::handlers(old_handlers), add = TRUE)
  progressr::handlers(global = TRUE)
  progressr::handlers(progress_handler)
  
  # Create per-outcome seeds (stable + distinct)
  # Use a simple deterministic scheme to avoid identical bootstrap streams.
  outcome_seeds <- seed + seq_along(outcomes_list) - 1L
  
  # Run in parallel with progress
  fits <- progressr::with_progress({
    p <- progressr::progressor(along = outcomes_list)
    
    furrr::future_map2(
      .x = outcomes_list,
      .y = outcome_seeds,
      .f = function(outcome, s) {
        p(sprintf("Outcome: %s", outcome))
        fit_one(outcome, this_seed = s)
      },
      .options = furrr::furrr_options(seed = TRUE) # reproducible parallel RNG
    )
  })
  
  names(fits) <- outcomes_list
  
  results <- purrr::map_dfr(fits, "result") |>
    dplyr::mutate(
      P_value_FDR_raw = stats::p.adjust(P_value_raw, method = "fdr"),
      P_value         = custom_format(P_value_raw),
      P_value_FDR     = custom_format(P_value_FDR_raw)
    ) |>
    dplyr::select(
      Outcome, N, Est_CI, P_value, P_value_FDR,
      P_value_raw, P_value_FDR_raw, Est, Conf_low, Conf_high
    )
  
  models <- purrr::map(fits, "model")
  
  list(results = results, models = models)
}

# qgcomp.glm.noboot for fixed weights
run_qgcomp_noboot_parallel <- function(data,
                                       exposures_list,
                                       outcomes_list,
                                       covariates_list,
                                       q = 4,
                                       seed = 42,
                                       family = NULL,
                                       id_cols = NULL,   # e.g., c("rand_id","blood_date") for merging
                                       workers = NULL,
                                       quiet = FALSE,
                                       progress_handler = c("txtprogressbar", "cli")) {
  
  # ---- deps ----
  pkgs <- c("dplyr", "purrr", "tibble", "qgcomp", "future", "furrr", "progressr", "rlang")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Please install required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  progress_handler <- match.arg(progress_handler)
  
  stopifnot(is.data.frame(data))
  stopifnot(length(exposures_list) > 0, length(outcomes_list) > 0)
  
  # ---- validate columns ----
  needed <- unique(c(exposures_list, outcomes_list, covariates_list, id_cols))
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in `data`: ", paste(missing_cols, collapse = ", "))
  }
  
  # ---- automatic family detection ----
  detect_family <- function(y) {
    y_nonmiss <- y[!is.na(y)]
    uniq <- unique(y_nonmiss)
    
    if (is.logical(y_nonmiss) ||
        is.factor(y_nonmiss) && length(uniq) == 2 ||
        (length(uniq) <= 2 && all(sort(uniq) %in% c(0,1)))) {
      return(stats::binomial())
    } else {
      return(stats::gaussian())
    }
  }
  
  # ID strategy for merging back
  if (is.null(id_cols) || length(id_cols) == 0) {
    data_id <- dplyr::mutate(data, row_id = dplyr::row_number())
    id_cols <- "row_id"
  } else {
    data_id <- data
  }
  
  # Keep only what we need for modeling + later composite construction
  data_small <- dplyr::select(data_id, dplyr::all_of(needed))
  
  rhs <- paste(c(exposures_list, covariates_list), collapse = " + ")
  
  # # ---- helper: compute q-category using cutpoints from training data ----
  # # Returns integer in 0..(q-1), NA if x is NA.
  # quantize_with_breaks <- function(x, breaks) {
  #   if (all(is.na(x))) return(rep(NA_integer_, length(x)))
  #   # cut returns factor; convert to integer 0..(q-1)
  #   out <- cut(x, breaks = breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)
  #   ifelse(is.na(out), NA_integer_, out - 1L)
  # }
  # 
  # # ---- helper: build per-exposure breaks from training data ----
  # # Uses quantiles; if ties collapse breaks too much, falls back to ntile for that exposure.
  # make_breaks <- function(x, q) {
  #   probs <- seq(0, 1, length.out = q + 1)
  #   br <- stats::quantile(x, probs = probs, na.rm = TRUE, type = 7)
  #   br <- unique(as.numeric(br))
  #   if (length(br) < 2) return(NULL)      # cannot cut
  #   if (length(br) < (q + 1)) return(br)  # fewer bins than q due to ties; still usable
  #   br
  # }
  
  # ---- helper: extract weights robustly across qgcomp versions ----
  get_weights <- function(mod, exposures_list) {
    # Common in qgcomp objects:
    # - mod$pos.weights and mod$neg.weights (named vectors)
    # Some versions also have mod$weights with sign; we handle both.
    
    if (!is.null(mod$pos.weights) && !is.null(mod$neg.weights)) {
      wpos <- mod$pos.weights
      wneg <- mod$neg.weights
      wpos <- wpos[exposures_list]; wpos[is.na(wpos)] <- 0
      wneg <- wneg[exposures_list]; wneg[is.na(wneg)] <- 0
      return(list(wpos = wpos, wneg = wneg))
    }
    
    if (!is.null(mod$weights)) {
      w <- mod$weights
      w <- w[exposures_list]
      w[is.na(w)] <- 0
      wpos <- pmax(w, 0)
      wneg <- pmax(-w, 0)
      return(list(wpos = wpos, wneg = wneg))
    }
    
    stop("Could not find weights in qgcomp model object. Inspect `names(model)` to locate weights.")
  }
  
  fit_one <- function(outcome) {
    form <- stats::as.formula(paste(outcome, "~", rhs))
    
    # complete cases for model fit
    dat_cc <- data_small |>
      dplyr::select(dplyr::all_of(c(id_cols, outcome, 
                                    exposures_list, covariates_list))) |>
      dplyr::filter(dplyr::if_all(dplyr::all_of(c(outcome, 
                                                  exposures_list, 
                                                  covariates_list)), 
                                  ~ !is.na(.x)))
    
    # Determine family
    fam <- if (is.null(family)) {
      detect_family(dat_cc[[outcome]])
    } else {
      family
    }
    
    if (!quiet) {
      message("Fitting qgcomp.noboot model for: ", outcome, deparse(form), 
              " with family = ", fam$family)
    }
    
    # Fit noboot model (fixed weights)
    model <- tryCatch(
      qgcomp::qgcomp.glm.noboot(
        f      = form,
        data   = dat_cc,
        expnms = exposures_list,
        q      = q,
        family = fam
      ),
      error = function(e) e
    )
    
    if (inherits(model, "error")) {
      if (!quiet) message("Model failed for outcome ", 
                          outcome, ": ", model$message)
      
      # composite for this outcome = all NA
      comp_tbl <- data_small |>
        dplyr::distinct(dplyr::across(dplyr::all_of(id_cols))) |>
        dplyr::mutate(!!paste0("comp_", outcome) := NA_real_)
      
      res_tbl <- tibble::tibble(
        Outcome     = outcome,
        N           = NA_integer_,
        Est_CI      = NA_character_,
        P_value_raw = NA_real_,
        Est         = NA_real_,
        Conf_low    = NA_real_,
        Conf_high   = NA_real_
      )
      
      return(list(result = res_tbl, model = NULL, composite = comp_tbl))
    }
    
    # Extract psi1 info (noboot often has SE/pvalue but may not have CI columns named the same)
    coefs <- summary(model)$coefficients
    if (!("psi1" %in% rownames(coefs))) {
      if (!quiet) message("Outcome ", outcome, ": 'psi1' not found in coefficients.")
      est <- NA_real_; pval <- NA_real_; conf_low <- NA_real_; conf_high <- NA_real_
      est_ci <- NA_character_
    } else {
      res <- coefs["psi1", , drop = TRUE]
      est  <- unname(res[["Estimate"]])
      pval <- unname(res[["Pr(>|z|)"]])
      
      # CI columns can differ; try common names, else NA
      conf_low  <- if ("Lower CI" %in% names(res)) unname(res[["Lower CI"]]) else NA_real_
      conf_high <- if ("Upper CI" %in% names(res)) unname(res[["Upper CI"]]) else NA_real_
      est_ci <- if (is.finite(conf_low) && is.finite(conf_high)) {
        sprintf("%.2f (%.2f, %.2f)", est, conf_low, conf_high)
      } else {
        sprintf("%.2f", est)
      }
    }
    
    res_tbl <- tibble::tibble(
      Outcome     = outcome,
      N           = nrow(dat_cc),
      Est_CI      = est_ci,
      P_value_raw = pval,
      Est         = est,
      Conf_low    = conf_low,
      Conf_high   = conf_high
    )
    
    # ---- Build composite exposure for ALL rows using training cutpoints + fixed weights ----
    w <- get_weights(model, exposures_list) # list(wpos, wneg)
    
    # # breaks per exposure computed from training complete-case data
    # breaks_list <- purrr::map(exposures_list, ~ make_breaks(dat_cc[[.x]], q))
    # names(breaks_list) <- exposures_list
    # 
    # # quantize each exposure in the full dataset using training breaks; fallback to ntile if breaks NULL
    # qmat <- purrr::map_dfc(exposures_list, function(xnm) {
    #   x_full <- data_small[[xnm]]
    #   br <- breaks_list[[xnm]]
    #   
    #   qx <- if (is.null(br)) {
    #     # extreme tie case: fall back to ntile on full data
    #     dplyr::ntile(x_full, q) - 1L
    #   } else {
    #     quantize_with_breaks(x_full, br)
    #   }
    #   
    #   tibble::tibble(!!xnm := as.integer(qx))
    # })
    
    # composite = sum(wpos * qx) - sum(wneg * qx)
    # qmat_num <- as.matrix(qmat)
    qmat_num <- as.matrix(model$qx)
    comp <- as.numeric(qmat_num %*% w$wpos) - as.numeric(qmat_num %*% w$wneg)
    
    comp_tbl <- data_small |>
      dplyr::distinct(dplyr::across(dplyr::all_of(id_cols))) |>
      dplyr::mutate(!!paste0("comp_", outcome) := comp)
    
    list(result = res_tbl, model = model, composite = comp_tbl)
  }
  
  # ---- workers + parallel plan ----
  if (is.null(workers)) {
    workers <- max(1, future::availableCores() - 1)
  }
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  
  old_handlers <- progressr::handlers()
  on.exit(progressr::handlers(old_handlers), add = TRUE)
  progressr::handlers(global = TRUE)
  progressr::handlers(progress_handler)
  
  # reproducibility: seed the parallel RNG streams
  set.seed(seed)
  
  fits <- progressr::with_progress({
    p <- progressr::progressor(along = outcomes_list)
    
    furrr::future_map(
      outcomes_list,
      ~ { p(sprintf("Outcome: %s", .x)); fit_one(.x) },
      .options = furrr::furrr_options(seed = TRUE)
    )
  })
  names(fits) <- outcomes_list
  
  results <- purrr::map_dfr(fits, "result") |>
    dplyr::mutate(
      P_value_FDR_raw = stats::p.adjust(P_value_raw, method = "fdr"),
      P_value         = custom_format(P_value_raw),
      P_value_FDR     = custom_format(P_value_FDR_raw)
    ) |>
    dplyr::select(
      Outcome, N, Est_CI, P_value, P_value_FDR,
      P_value_raw, P_value_FDR_raw, Est, Conf_low, Conf_high
    )
  
  models <- purrr::map(fits, "model")
  
  # Merge all per-outcome composite columns into one tibble keyed by id_cols
  composites <- purrr::map(fits, "composite") |>
    purrr::reduce(dplyr::left_join, by = id_cols)
  
  list(results = results, models = models, 
       composites = composites, id_cols = id_cols)
}