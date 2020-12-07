# vim: foldmethod=marker

# Estimate -----------------------------------------------------------------{{{

#' Stratified clustered bootstrapping of (multivariate) linear models
#'
#' coef(fit, boot = TRUE) returns a 3d array where the first dimention is
#' the response variables, second dimension is the coefficients, and third
#' dimension is bootstrap replicates.
#'
#' @param observed A phyloseq object or phyloseq-otu_table object
#' @param data A data.frame or object coerceable to one
#' @param log_obs Set to TRUE if `observed` is already on the log scale
#'
#' @export
#'
# TODO: Allow the observation matrix to be in any orientation; and perhaps make
# it be passed separately rather than in the formula.
complm <- function(observed,
                   formula, 
                   data = sample_data(observed),
                   contrasts = NULL,
                   boot = FALSE,
                   times = 1000,
                   clusters = seq(nrow(data)),
                   strata = rep(1, nrow(data)),
                   log_obs = FALSE) {
  mf <- stats::model.frame(formula, 
    as(data, "data.frame") %>% tibble::remove_rownames()
  )
  mt <- stats::terms(mf)
  X <- stats::model.matrix(mt, mf, contrasts)
  ref <- stats::model.offset(mf)
  if (!is.null(ref)) stop("Offsets not currently supported")
  # Save the same names just in case
  sns <- sample_names(observed)

  # TODO: stopifnot observed is a phyloseq or otu_table
  ## Prepare the response matrix
  # Coerce phyloseq object to otu_table object
  observed <- otu_table(observed)
  # Convert observed to clr
  if (!log_obs) {
    stopifnot(all(observed > 0))
    observed <- transform_sample_counts(observed, log)
  }
  observed <- transform_sample_counts(observed, function(x) x - mean(x))
  # Get `observed` with taxa as rows, and `observed_t` with taxa as cols, and
  # coerce to plain matrices
  if (taxa_are_rows(observed)) {
    observed_t <- t(observed)
  } else {
    observed_t <- observed
    observed <- t(observed_t)
  }
  observed <- as(observed, "matrix")
  observed_t <- as(observed_t, "matrix")

  ## The main fit
  fit <- stats::lm.fit(X, observed_t, offset = ref)

  ## Bootstrap replicates
  if (boot) {
    wts <- draw_weights(data, times, clusters = clusters, strata = strata)
    bootreps <- wts %>%
      purrr::map(~stats::lm.wfit(X, observed_t, .x, offset = ref))
    extract_slot <- function(bootreps, slot, name = as.character(slot)) {
      bootreps %>%
        # Extract item in specified slot
        purrr::map(slot) %>%
        # Combine matrices into a 3-d array
        simplify2array %>%
        # Move responses/taxa from cols to rows
        aperm(c(2,1,3)) %>%
        # Set names of dimensions
        {names(dimnames(.)) <- c("response", name, ".draw"); .}
    }
    # Replace `bootreps` with a list of slot arrays
    bootreps <- 
      tibble::tibble(
        slot = c("coefficients", "residuals", "fitted.values"),
        name = c("term", "sample_idx", "sample_idx"),
      ) %>%
      purrr::pmap(extract_slot, bootreps = bootreps) %>%
      rlang::set_names(c("coefficients", "residuals", "fitted"))
  } else 
    bootreps <- NULL

  structure(class = "mc_complm_fit",
    list(
      fit = fit,
      offset = ref,
      contrasts = attr(X, "contrasts"),
      # model frame
      model = mf,
      # model terms
      terms = mt,
      # factor levels
      xlevels = stats::.getXlevels(mt, mf),
      #> data = data,
      sample_names = sns,
      observed = observed,
      coefficients = coef(fit) %>% t,
      residuals = residuals(fit) %>% t,
      fitted = fitted(fit) %>% t,
      bootreps = bootreps
    )
  )
}

# Function for getting the observed values in a desired form; not sure if
# actually needed for anything

#' @export
complm_observed <- function(object, units = "clr") {
  stopifnot(units %in% c("clr", "proportion"))
  if (units == "clr") 
    object$observed
  else if (units == "proportion") 
    object$observed %>% exp %>% apply(2, metacal::close_elts)
}

#' @export
fitted.mc_complm_fit <- function(object, boot = FALSE, units = "clr") {
  if (boot) stop("boot option not yet implemented")
  stopifnot(units %in% c("clr", "proportion"))
  if (units == "clr") 
    object$fitted
  else if (units == "proportion") 
    object$fitted %>% exp %>% apply(2, metacal::close_elts)
}

#' @export
residuals.mc_complm_fit <- function(object, boot = FALSE, units = "clr") {
  if (boot) stop("boot option not yet implemented")
  stopifnot(units %in% c("clr", "proportion"))
  if (units == "clr") 
    object$residuals
  else if (units == "proportion") {
    # Compute observed and fitted proportions, then subtract
    list(observed = object$observed, fitted = object$fitted) %>%
      purrr::map(exp) %>%
      purrr::map(apply, 2, metacal::close_elts) %>%
      {.$observed - .$fitted}
  }
}

#' @export
coef.mc_complm_fit <- function(object, boot = FALSE) {
  if (!boot) {
    object$coefficients
  } else {
    if (is.null(object$bootreps))
      stop("Bootstrap replicates are not included")
    object$bootreps$coefficients
  }
}

#' @export
terms.mc_complm_fit <- function(x) {
  x$terms
}

#' @export
predict.mc_complm_fit <- function(object, newdata, boot = FALSE) {
  # terms <- delete.response(terms(object))
  mf <- stats::model.frame(stats::terms(object), as(newdata, "data.frame"), 
    xlev = object$xlevels)
  X_t <- stats::model.matrix(stats::terms(object), mf, 
    contrasts.arg = object$contrasts) %>% t
  colnames(X_t) <- NULL
  if (!boot) {
    (object$coefficients %*% X_t) %>%
    {names(dimnames(.)) <- c(".response", ".prediction"); .}
  } else {
    if (is.null(object$bootreps))
      stop("Bootstrap replicates are not included")
      brc <- object$bootreps$coefficients
      purrr::map(seq(dim(brc)[3]), ~brc[,,.x] %*% X_t) %>%
        simplify2array %>%
        # Reset names of dimensions
        {names(dimnames(.)) <- c(".response", ".prediction", ".draw"); .}
  }
}

#' @export
print.mc_complm_fit <- function(x, coef = FALSE, digits = 2) {
  cat('A metacal complm fit.', fill = TRUE)
  if (coef) {
    cat(fill = TRUE)
    cat('Estimated coefficients:', fill = TRUE)
    print(x$coefficients, digits = digits)
  }
  if (!is.null(x$bootreps)) {
    cat(fill = TRUE)
    cat('Contains', dim(x$bootreps$coefficients)[3], 
      'bootstrap replicates.', fill = TRUE)
  }
  invisible(x)
}

#> summary.mc_complm_fit <- function(object) {
#>   s <- structure(list(), class = "mc_complm_fit_summary")
#> }

# Broom generics

# Consider computing the bootrep summary stats in a summary method that is
# called upon here; or perhaps precomputing in the original call to complm so
# that this is faster. Or: compute with calls to apply on the bootrep arrays,
# so that don't need to melt the giant 3d array.

# Which summary stats to report? Perhaps median and median sad are better?

#' @export
tidy.mc_complm_fit <- function(x) {
  # point estimate
  tb <- x %>%
    coef %>%
    tibble::as_tibble(rownames = "response") %>%
    tidyr::pivot_longer(-response, names_to = "term", values_to = "estimate")
  # mean and standard error from bootreps
  if (!is.null(x$bootreps)) {
    bootrep <- x %>% 
      coef(boot = TRUE) %>% 
      data.table::as.data.table() %>%
      tibble::as_tibble() %>%
      dplyr::with_groups(c(response, term), dplyr::summarize, 
        dplyr::across(value, list(mean = mean, sd = stats::sd))
      ) %>%
      dplyr::rename(std.error = value_sd)
    tb <- tb %>% 
      dplyr::left_join(bootrep, by = c("response", "term")) %>%
      dplyr::mutate(bias = estimate - value_mean) %>%
      dplyr::select(-value_mean) %>%
      dplyr::relocate(bias, .before = std.error)
  }
  tb
}

#' @export
augment.mc_complm_fit <- function(x, units = "clr") {
  stopifnot(units %in% c("clr", "proportion"))
  #> .observed <- x %>%
  #>   complm_observed(units = units) %>%
  #>   as_tibble(rownames = ".response") %>%
  #>   pivot_longer(-.response, names_to = ".sample", values_to = ".observed")
  .fitted <- x %>%
    fitted(units = units) %>%
    tibble::as_tibble(rownames = ".response") %>%
    tidyr::pivot_longer(-.response, names_to = ".sample", values_to = ".fitted")
  .resid <- x %>% 
    resid(units = units) %>%
    tibble::as_tibble(rownames = ".response") %>%
    tidyr::pivot_longer(-.response, names_to = ".sample", values_to = ".resid")
  x$model %>% 
    tibble::as_tibble() %>%
    dplyr::mutate(.sample = colnames(fitted(x))) %>%
    dplyr::left_join(.fitted, by = c(".sample")) %>%
    dplyr::left_join(.resid, by = c(".sample", ".response"))
}


# # Tidybayes functions
#
# # Problem: Supposed to have one column per variable; Not sure how to deal with
# # the multivariate outcome
# tidy_draws.complm_fit <- function(model) {
#   if (!requireNamespace("tidybayes", quietly = TRUE)) {
#     stop("The tidybayes package must be installed to use `tidy_draws()`.", 
#       call. = FALSE)
#   }
#
#   if (is.null(model$bootreps))
#     stop("Bootstrap replicates are not included")
#   coef(model, boot = TRUE) %>% 
#     data.table::as.data.table(key = ".draw") %>%
#     as_tibble %>%
#     mutate(.chain = 1L, .iteration = .draw)
# }

# }}}


# Helpers ------------------------------------------------------------------{{{

#' Set the reference taxon/taxa
#'
#' @examples
#' smaller array for testing
#' fit0 <- fit[1:5,,1:2]
#' set_denom(fit0, 1)
#' set_denom(fit0, 1:2)
#' set_denom(fit0, c("T1", "T2"))
#'
#' @export
set_ref_taxa <- function(x, ref) {
  apply(x, c(2,3), function(x) {x - mean(x[ref])})
}

#' Add a reference coefficient to the array
#'
#' fit1 <- fit0 %>% add_ref_coef("protocolP1")
#'
#' @export
add_ref_coef <- function(fit, name) {
  ref_coef <- array(0, 
    dim = c(dim(fit)[1], 1, dim(fit)[3]),
    dimnames = c(dimnames(fit)[1], list(coefficient = name), dimnames(fit)[3])
  )
  a <- abind::abind(ref_coef, fit, along = 2)
  names(dimnames(a)) <- names(dimnames(fit))
  a
}

#' Set the reference taxa or coefficient
#'
#' Note, I'm not sure what conditions its actually safe to do this
#'
#' @examples
#' set_ref(fit1, 1, 1)
#' set_ref(fit1, 1, c("T1", "T2"))
#' set_ref(fit1, 2, c("protocolP2"))
#' set_ref(fit1, 2, str_c("protocol", protocols))
#' set_ref(fit1, 3, 1) 
#'
#' @export
set_ref <- function(x, margin, ref) {
  stopifnot(identical(length(dim(x)), 3L))
  stopifnot(margin %in% 1:2)
  if (is.character(ref))
    stopifnot(all(ref %in% dimnames(x)[[margin]]))
  x <- apply(x, setdiff(seq(3), margin), function(y) {y - mean(y[ref])})
  if (margin == 2)
    x <- aperm(x, c(2,1,3))
  x
}

#' Convert response (taxa) to pairwise logratios
#'
#' `upper` and `diag` arguments patterned after `dist`, and refer to the
#' corresponding entries in the pairwise difference matrix, even though we're
#' giving the pairwise differences stacked in a vector.
#'
#' @examples
#' fit[1:3,,1:2] %>% set_pairwise
#' fit[1:3,,1:2] %>% set_pairwise(upper = TRUE)
#'
#' @export
set_pairwise <- function(x, diag = FALSE, upper = FALSE) {
  idx <- seq(dim(x)[1])
  prs <- expand.grid(i = idx, j = idx)
  if (!diag) prs <- filter(prs, i != j)
  if (!upper) prs <- filter(prs, i <= j)

  if (is.matrix(x)) 
    y <- x[prs$i,] - x[prs$j,]
  else
    y <- x[prs$i,,] - x[prs$j,,]
  nms <- dimnames(x)[[1]]
  dimnames(y)[[1]] <- stringr::str_c(nms[prs$i], nms[prs$j], sep = ":")
  y
}

# TODO: Make s.t. clusters and strata can optionally be the names of variables
# in `data`
draw_weights <- function(data, 
                         times, 
                         clusters = seq(nrow(data)),
                         strata = rep(1, nrow(data))) {
  clusters <- clusters %>% as.factor %>% as.integer
  n_clusters <- clusters %>% unique %>% length
  cluster_weights <- stats::rexp(n_clusters * times) %>%
    matrix(nrow = n_clusters)
  # Column names are needed before converting matrix to tibble
  colnames(cluster_weights) <- stringr::str_c("...", seq(ncol(cluster_weights)))
  # Expand the cluster weights into sample weights, then normalize within each
  # stratum to have the appropriate stratum weights
  cluster_weights[clusters,] %>%
    tibble::as_tibble() %>%
    tibble::add_column(stratum = strata) %>%
    dplyr::group_by(stratum) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
        ~metacal::close_elts(.) * dplyr::n())) %>%
    dplyr::ungroup() %>%
    dplyr::select(-stratum) %>%
    # Return as a list of weight vectors
    as.list %>%
    unname
}

# }}}

