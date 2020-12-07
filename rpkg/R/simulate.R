# vim: foldmethod=marker

#' Simulate a differential-bias dataset
#'
#' Note, assumes that specimen is nested within condition, like the vivo-vitro
#' experiment. Should think of condition as a treatment-var of the specimens.
#' The first condition is always named "Ref". Tries to even splits specimen
#' among the conditions. (Perhaps I should just call the conditions C1, ...)
#'
#' The returned `coefficients` are the coefficents relative to the Ref-Ref
#' protocol, condition combination.
#'
#' @param n 5-vector specifying the numbers of taxa, specimens, protocols,
#' conditions, and replicates
#' @param sd 5-vector specifying the standard deviation for the draws of the
#' specimen, protocol, condition, and protocol:condition effects and finally
#' the residuals
#
# TODO: Give option of returning a phyloseq otu-table object on the
# counts/abundance
simulate_dataset <- function(n, sd, log = TRUE, phyloseq = FALSE) {
  n_taxa <- n[[1]]
  n_specimens <- n[[2]]
  n_protocols <- n[[3]]
  n_conditions <- n[[4]] 
  n_replicates <- n[[5]]

  taxa <- stringr::str_c("T", seq(n_taxa)) %>% forcats::fct_inorder()
  specimens <- stringr::str_c("S", seq(n_specimens)) %>% forcats::fct_inorder()
  protocols <- c("Ref", stringr::str_c("P", seq(n_protocols))) %>%
    forcats::fct_inorder()
  conditions <- c("Ref", stringr::str_c("C", seq(n_conditions))) %>%
    forcats::fct_inorder()

  condition_of_specimens <- sample(
    rep_len(conditions, n_specimens),
    n_specimens, replace = FALSE
  )
  
  sam <- tidyr::crossing(
    tidyr::nesting(specimen = specimens, condition = condition_of_specimens),
    replicate = seq(n_replicates),
    protocol = protocols,
  ) %>%
    dplyr::mutate(
      sample = stringr::str_glue("{specimen}_{protocol}_{condition}_{replicate}"),
    )

  X <- stats::model.matrix(
    ~ 0 + specimen + protocol + condition + protocol:condition, 
    data = sam
  )
  rownames(X) <- sam$sample

  # The observed dataset does not have the reference protocol measurements;
  # let's also create corresponding data frames and model matrices restricted
  # to just those.
  sam.obs <- sam %>%
    dplyr::filter(protocol != "Ref") %>%
    dplyr::mutate(dplyr::across(where(is.factor), forcats::fct_drop))
  X.obs <- X[sam.obs$sample,]

  # Generate coefficients
  cff_tb <- tidyr::crossing(
    coefficient = colnames(X) %>% forcats::fct_inorder(),
    taxon = taxa
  ) %>%
    dplyr::mutate(
      mean = 0,
      sd = dplyr::case_when(
        # Actual abundances
        stringr::str_detect(coefficient, "^specimen[^:]+$") ~ sd[[1]],
        # Actual bias in condition A
        stringr::str_detect(coefficient, "^protocol[^:]+$") ~ sd[[2]],
        # Bias associated with the conditions
        stringr::str_detect(coefficient, "^condition[^:]+$") ~ sd[[3]],
        # Change in protocol bias associated with the conditions
        stringr::str_detect(coefficient, "^protocol[^:]+:condition[^:]+$") ~ sd[[4]]
      ),
      value = stats::rnorm(dplyr::n(), mean, sd)
    ) %>%
    dplyr::with_groups(coefficient, dplyr::mutate, dplyr::across(value, ~ . - mean(.)))
  B <- cff_tb %>%
    metacal::build_matrix(taxon, coefficient, value)
  stopifnot(identical(colnames(X), colnames(B)))

  # Generate the predicted and observed data on just the non-ref-protocol
  # samples
  Y.pred <- B %*% t(X.obs)
  noise <- stats::rnorm(Y.pred, sd = sd[[5]]) %>%
    matrix(nrow = nrow(Y.pred)) %>%
    apply(2, function(x) x - mean(x))
    # scale(scale = FALSE)
  Y <- Y.pred + noise
  # Y should have one column per sample
  stopifnot(identical(ncol(Y), nrow(sam.obs)))
  # The samples should already be centered (CLR values)
  stopifnot(all.equal(rep(0, ncol(Y)) , colSums(Y), check.attributes = FALSE))

  if (!log) {
    Y <- exp(Y)
    Y.pred <- exp(Y.pred)
  }

  if (phyloseq) {
    ps_sam <- sam.obs %>%
      tibble::column_to_rownames("sample") %>%
      phyloseq::sample_data()
    ps_otu <- phyloseq::otu_table(Y, taxa_are_rows = TRUE)
    ps <- phyloseq::phyloseq(ps_sam, ps_otu)
  } else
    ps <- NULL

  list(n = n, sd = sd, log = log, data = sam.obs, coefficients = B,
    observed = Y, predicted = Y.pred, phyloseq = ps)
}
