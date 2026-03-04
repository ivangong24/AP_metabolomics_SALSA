run_wqs <- function(data,
                    outcome,
                    mix_name,              # character vector of exposures
                    covariates,
                    id_cols,
                    q = 4,
                    validation = 0.6,
                    b = 200,
                    b1_pos = TRUE,
                    b_constr = FALSE,
                    rh = 5,
                    family = c("binomial", "gaussian"),
                    seed = 42,
                    ...) {
  
  stopifnot(is.data.frame(data), length(id_cols) > 0, length(mix_name) > 0)
  family <- match.arg(family)
  
  needed <- unique(c(id_cols, outcome, mix_name, covariates))
  miss <- setdiff(needed, names(data))
  if (length(miss) > 0) stop("Missing columns in `data`: ", paste(miss, collapse = ", "))
  
  # Complete cases for model variables
  dat_cc <- data |>
    dplyr::select(dplyr::all_of(needed)) |>
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.x)))
  
  # Formula: outcome ~ wqs + covariates
  rhs <- paste(c("wqs", covariates), collapse = " + ")
  form <- stats::as.formula(paste(outcome, "~", rhs))
  
  mod <- gWQS::gwqs(
    formula = form,
    mix_name = mix_name,      # exposures
    data = dat_cc,
    q = q,
    validation = validation,
    b = b,
    b1_pos = b1_pos,
    b1_constr = b1_constr,
    rh = rh,
    family = family,
    seed = seed,
    ...
  )
  
  # Extract WQS index (varies by gWQS version; these are common)
  wqs_vec <- NULL
  if (!is.null(mod$wqs)) wqs_vec <- mod$wqs
  if (is.null(wqs_vec) && !is.null(mod$data) && "wqs" %in% names(mod$data)) {
    wqs_vec <- mod$data$wqs
  }
  if (is.null(wqs_vec)) {
    stop("Couldn't find WQS index in model output. Try `names(mod)` or `str(mod, max.level=1)`.")
  }
  
  # Build y_wqs_df with IDs (merge-ready)
  wqs_df <- dat_cc |>
    dplyr::select(dplyr::all_of(id_cols)) |>
    dplyr::mutate(wqs = as.numeric(wqs_vec))
  
  list(model = mod, wqs_df = wqs_df)
}