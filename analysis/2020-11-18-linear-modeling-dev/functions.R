# vim: set foldmethod=marker
# Simulate -----------------------------------------------------------------{{{

# 
# Note, assumes that specimen is nested within condition, like the vivo-vitro
# experiment. Should think of condition as a treatment-var of the specimens.
# The first condition is always named "Ref". Tries to even splits specimen
# among the conditions. (Perhaps I should just call the conditions C1, ...)
#
# The returned `coefficients` are the coefficents relative to the Ref-Ref
# protocol, condition combination.
#
# TODO: Give option of returning a phyloseq otu-table object on the
# counts/abundance
simulate_bias_dataset <- function(n, sd) {
  n_taxa <- n[[1]]
  n_specimens <- n[[2]]
  n_protocols <- n[[3]]
  n_conditions <- n[[4]] 
  n_replicates <- n[[5]]

  taxa <- str_c("T", seq(n_taxa)) %>% fct_inorder
  specimens <- str_c("S", seq(n_specimens)) %>% fct_inorder
  protocols <- c("Ref", str_c("P", seq(n_protocols))) %>% fct_inorder
  conditions <- c("Ref", str_c("C", seq(n_conditions))) %>% fct_inorder

  condition_of_specimens <- sample(
    rep_len(conditions, n_specimens),
    n_specimens, replace = FALSE
  )
  
  sam <- crossing(
    nesting(specimen = specimens, condition = condition_of_specimens),
    replicate = seq(n_replicates),
    protocol = protocols,
  ) %>%
    mutate(
      sample = str_glue("{specimen}_{protocol}_{condition}_{replicate}"),
    )

  X <- model.matrix(
    ~ 0 + specimen + protocol + condition + protocol:condition, 
    data = sam
  )
  rownames(X) <- sam$sample

  # The observed dataset does not have the reference protocol measurements;
  # let's also create corresponding data frames and model matrices restricted
  # to just those.
  sam.obs <- sam %>%
    filter(protocol != "Ref") %>%
    mutate(across(where(is.factor), fct_drop))
  X.obs <- X[sam.obs$sample,]

  # Generate coefficients
  cff_tb <- crossing(coefficient = colnames(X) %>% fct_inorder, taxon = taxa) %>%
    mutate(
      mean = 0,
      sd = case_when(
        # Actual abundances
        str_detect(coefficient, "^specimen[^:]+$") ~ sd[[1]],
        # Actual bias in condition A
        str_detect(coefficient, "^protocol[^:]+$") ~ sd[[2]],
        # Bias associated with the conditions
        str_detect(coefficient, "^condition[^:]+$") ~ sd[[3]],
        # Change in protocol bias associated with the conditions
        str_detect(coefficient, "^protocol[^:]+:condition[^:]+$") ~ sd[[4]]
      ),
      value = rnorm(n(), mean, sd)
    ) %>%
    with_groups(coefficient, mutate, across(value, ~ . - mean(.)))
  B <- cff_tb %>%
    build_matrix(taxon, coefficient, value)
  stopifnot(identical(colnames(X), colnames(B)))

  # Generate the predicted and observed data on just the non-ref-protocol
  # samples
  Y.pred <- B %*% t(X.obs)
  noise <- rnorm(Y.pred, sd = sd[[5]]) %>%
    matrix(nrow = nrow(Y.pred)) %>%
    apply(2, function(x) x - mean(x))
    # scale(scale = FALSE)
  Y <- Y.pred + noise
  # Y should have one column per sample
  stopifnot(identical(ncol(Y), nrow(sam.obs)))
  # The samples should already be centered (CLR values)
  stopifnot(all.equal(rep(0, ncol(Y)) , colSums(Y), check.attributes = FALSE))

  list(n = n, sd = sd, data = sam.obs, coefficients = B,
    observed = Y, predicted = Y.pred)
}

# }}}

# Helpers ------------------------------------------------------------------{{{

#' Set the reference taxon/taxa
#'
#' @example
#' smaller array for testing
#' fit0 <- fit[1:5,,1:2]
#' set_denom(fit0, 1)
#' set_denom(fit0, 1:2)
#' set_denom(fit0, c("T1", "T2"))
set_ref_taxa <- function(x, ref) {
  apply(x, c(2,3), function(x) {x - mean(x[ref])})
}

#' Add a reference coefficient to the array
#'
#' fit1 <- fit0 %>% add_ref_coef("protocolP1")
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
  dimnames(y)[[1]] <- str_c(nms[prs$i], nms[prs$j], sep = ":")
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
  cluster_weights <- rexp(n_clusters * times) %>%
    matrix(nrow = n_clusters)
  # Column names are needed before converting matrix to tibble
  colnames(cluster_weights) <- str_c("...", seq(ncol(cluster_weights)))
  # Expand the cluster weights into sample weights, then normalize within each
  # stratum to have the appropriate stratum weights
  cluster_weights[clusters,] %>%
    as_tibble %>%
    add_column(stratum = strata) %>%
    group_by(stratum) %>%
    mutate(across(everything(), ~close_elts(.) * n()),) %>%
    ungroup %>%
    select(-stratum) %>%
    # Return as a list of weight vectors
    as.list %>%
    unname
}

# }}}

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
  mf <- model.frame(formula, as(data, "data.frame") %>% tibble::remove_rownames())
  mt <- terms(mf)
  X <- model.matrix(mt, mf, contrasts)
  ref <- model.offset(mf)
  if (!is.null(ref)) stop("Offsets not currently supported")
  # Save the same names just in case
  sns <- sample_names(observed)

  # TODO: stopifnot observed is a phyloseq or otu_table
  ## Prepare the response matrix
  # Coerce phyloseq object to otu_table object
  observed <- otu_table(observed)
  # Convert observed to clr
  if (log_obs) {
    observed <- transform_sample_counts(observed, ~. - mean(.))
  } else {
    stopifnot(all(observed > 0))
    observed <- transform_sample_counts(observed, clr)
  }
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

  if (!is.matrix(observed)) stop("The response must be a matrix")
  # NOTE: Currently model.response will give a 1-d vector if obs is a column
  # matrix, which makes the next line moot
  if (!dim(observed_t)[[1]] > 1) stop("There must be more than 1 taxon")

  ## The main fit
  fit <- lm.fit(X, observed_t, offset = ref)

  ## 
  if (boot) {
    wts <- draw_weights(data, times, clusters = clusters, strata = strata)
    bootreps <- wts %>%
      map(~lm.wfit(X, observed_t, .x, offset = ref))
    extract_slot <- function(bootreps, slot, name = as.character(slot)) {
      bootreps %>%
        # Extract item in specified slot
        map(slot) %>%
        # Combine matrices into a 3-d array
        simplify2array %>%
        # Move responses/taxa from cols to rows
        aperm(c(2,1,3)) %>%
        # Set names of dimensions
        {names(dimnames(.)) <- c("response", name, ".draw"); .}
    }
    # Replace `bootreps` with a list of slot arrays
    bootreps <- 
      tibble(
        slot = c("coefficients", "residuals", "fitted.values"),
        name = c("coefficient", "residual", "fitted"),
      ) %>%
      pmap(extract_slot, bootreps = bootreps) %>%
      set_names(c("coefficients", "residuals", "fitted"))
  } else 
    bootreps <- NULL

  structure(class = "mc_complm_fit",
    list(
      fit = fit,
      offset <- ref,
      contrasts <- attr(X, "contrasts"),
      # model frame
      model = mf,
      # model terms
      terms = mt,
      # factor levels
      xlevels = .getXlevels(mt, mf),
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
complm_observed <- function(object, units = "clr") {
  stopifnot(units %in% c("clr", "proportion"))
  if (units == "clr") 
    object$observed
  else if (units == "proportion") 
    object$observed %>% exp %>% apply(2, close_elts)
}

# TODO: Implement the bootrep option for fitted and residuals

#' @param boot If TRUE, use the bootstrap replicates; otherwise use the point
#' estimate
#' @export
fitted.mc_complm_fit <- function(object, boot = FALSE, units = "clr") {
  if (boot) stop("boot option not yet implemented")
  stopifnot(units %in% c("clr", "proportion"))
  if (units == "clr") 
    object$fitted
  else if (units == "proportion") 
    object$fitted %>% exp %>% apply(2, close_elts)
}

#' @param boot If TRUE, use the bootstrap replicates; otherwise use the point
#' estimate
#' @export
residuals.mc_complm_fit <- function(object, boot = FALSE, units = "clr") {
  if (boot) stop("boot option not yet implemented")
  stopifnot(units %in% c("clr", "proportion"))
  if (units == "clr") 
    object$residuals
  else if (units == "proportion") {
    # Compute observed and fitted proportions, then subtract
    list(observed = object$observed, fitted = object$fitted) %>%
      map(exp) %>%
      map(apply, 2, close_elts) %>%
      {.$observed - .$fitted}
  }
}

#' @param boot If TRUE, use the bootstrap replicates; otherwise use the point
#' estimate
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

terms.mc_complm_fit <- function(object) {
  object$terms
}

predict.mc_complm_fit <- function(object, newdata, boot = FALSE) {
  # terms <- delete.response(terms(object))
  mf <- model.frame(terms(object), as(newdata, "data.frame"), 
    xlev = object$xlevels)
  X_t <- model.matrix(terms(object), mf, contrasts.arg = object$contrasts) %>% t
  colnames(X_t) <- NULL
  if (!boot) {
    (object$coefficients %*% X_t) %>%
    {names(dimnames(.)) <- c(".response", ".prediction"); .}
  } else {
    if (is.null(object$bootreps))
      stop("Bootstrap replicates are not included")
      brc <- object$bootreps$coefficients
      map(seq(dim(brc)[3]), ~brc[,,.x] %*% X_t) %>%
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
tidy.mc_complm_fit <- function(x) {
  # point estimate
  tb <- x %>%
    coef %>%
    as_tibble(rownames = "response") %>%
    pivot_longer(-response, names_to = "term", values_to = "estimate")
  # mean and standard error from bootreps
  if (!is.null(x$bootreps)) {
    bootrep <- x %>% 
      coef(boot = TRUE) %>% 
      data.table::as.data.table() %>%
      as_tibble %>%
      with_groups(c(response, coefficient), summarize, 
        across(value, list(mean = mean, sd = sd))
      ) %>%
      rename(term = coefficient, std.error = value_sd)
    tb <- tb %>% 
      left_join(bootrep, by = c("response", "term")) %>%
      mutate(bias = estimate - value_mean) %>%
      select(-value_mean) %>%
      relocate(bias, .before = std.error)
  }
  tb
}

augment.mc_complm_fit <- function(x, units = "clr") {
  stopifnot(units %in% c("clr", "proportion"))
  #> .observed <- x %>%
  #>   complm_observed(units = units) %>%
  #>   as_tibble(rownames = ".response") %>%
  #>   pivot_longer(-.response, names_to = ".sample", values_to = ".observed")
  .fitted <- x %>%
    fitted(units = units) %>%
    as_tibble(rownames = ".response") %>%
    pivot_longer(-.response, names_to = ".sample", values_to = ".fitted")
  .resid <- x %>% 
    resid(units = units) %>%
    as_tibble(rownames = ".response") %>%
    pivot_longer(-.response, names_to = ".sample", values_to = ".resid")
  x$data %>% 
    mutate(.sample = colnames(fitted(x))) %>%
    # left_join(.observed, by = ".sample") %>%
    left_join(.fitted, by = c(".sample")) %>%
    left_join(.resid, by = c(".sample", ".response"))
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

