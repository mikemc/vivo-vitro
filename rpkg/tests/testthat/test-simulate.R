test_that("simulation and fitting work", {
  set.seed(42)
  sim <- simulate_dataset(
    n = c(5, 5, 3, 2, 2),
    sd = c(1, 1, 0.2, 0.2, 0.3),
    phyloseq = TRUE
  )
  expect_s4_class(sim$phyloseq, "phyloseq")

  ## Check that fit on deterministic sim gives exact coefficients
  obs <- sim$predicted %>% exp %>% otu_table(taxa_are_rows = TRUE)
  fit <- complm(obs,
    ~ 0 + specimen + protocol + condition + protocol:condition, 
    data = sim$data, 
    boot = TRUE, times = 10, 
    clusters = sim$data$specimen,
    strata = sim$data$condition
  )
  expect_equal(fitted(fit), sim$predicted)
  B <- sim$coefficients
  # Differential bias of P2:P1 in the Ref condition
  expect_equal(
    coef(fit)[,"protocolP2"],
    B[, "protocolP2"] - B[, "protocolP1"]
  )
  # Change in P2/P1 differential bias moving from Ref condition to C1
  expect_equal(
    coef(fit)[,"protocolP2:conditionC1"],
    B[, "protocolP2:conditionC1"] - B[, "protocolP1:conditionC1"]
  )

  expect_s3_class(tidy(fit), "tbl_df")
  expect_s3_class(augment(fit), "tbl_df")
})
