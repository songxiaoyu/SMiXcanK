# Helper functions shared by the simulation drivers.
#
# This version uses the current SMiXcan package API:
#   - MiXcan_train_K() for prediction model training
#   - run_gwas() for SNP-level summary statistics
#   - SMiXcan_assoc_test_K() for summary-statistics S-MiXcan


# Main function: run one simulation analysis function - including
#  prediction to association for multiple methods.
run_one_simulation <- function(sim_data,
                               family = c("gaussian", "binomial"),
                               regularization = "fixed",
                               reg_scale = 0.1,
                               condition_max = 100,
                               reg_min_scale = 0.001,
                               reg_max_scale = 0.1,
                               alpha = 0.5,
                               nfolds = 10,
                               smixcan_ld_reference = c("train", "test")) {
  family <- match.arg(family)
  smixcan_ld_reference <- match.arg(smixcan_ld_reference)
  if (family == "gaussian") {
    sim_data <- sim_center(sim_data)
    family_obj <- stats::gaussian()
    outcome <- as.numeric(sim_data$D.test[, 1])
  } else {
    sim_data <- sim_center_binom(sim_data)
    family_obj <- stats::binomial()
    outcome <- as.numeric(sim_data$D.test[, 1])
  }

  p_results <- train_prediction_model_latest(
    y.train = sim_data$y.train,
    x.train = sim_data$x.train,
    pi.train = sim_data$pi.train,
    alpha = alpha,
    nfolds = nfolds
  )

  W <- p_results$W
  selected_snp <- p_results$selected_snp

  gwas_results <- if (family == "binomial") {
    SMiXcan::run_gwas(
      X = sim_data$x.test[, selected_snp, drop = FALSE],
      D = outcome,
      family0 = family_obj,
      method_binomial = "glm"
    )
  } else {
    SMiXcan::run_gwas(
      X = sim_data$x.test[, selected_snp, drop = FALSE],
      D = outcome,
      family0 = family_obj
    )
  }

  tilde_y_1 <- as.numeric(sim_data$x.test[, selected_snp, drop = FALSE] %*% W[, 1])
  tilde_y_2 <- as.numeric(sim_data$x.test[, selected_snp, drop = FALSE] %*% W[, 2])

  predixcan <- run_predixcan_latest(
    outcome = outcome,
    x.test = sim_data$x.test,
    tissue_weights = p_results$tissue_weights,
    tissue_intercept = p_results$tissue_intercept,
    family = family_obj
  )
  mixcan_sep <- run_mixcan_separate(outcome, tilde_y_1, tilde_y_2, family_obj)
  mixcan_join <- run_mixcan_latest(outcome, tilde_y_1, tilde_y_2, family)

  n0 <- n1 <- NULL
  if (family == "binomial") {
    n0 <- sum(outcome == 0, na.rm = TRUE)
    n1 <- sum(outcome == 1, na.rm = TRUE)
  }

  smixcan <- run_smixcan_latest(
    W = W,
    selected_snp = selected_snp,
    sim_data = sim_data,
    gwas_results = gwas_results,
    family = family,
    n0 = n0,
    n1 = n1,
    regularization = regularization,
    reg_scale = reg_scale,
    condition_max = condition_max,
    reg_min_scale = reg_min_scale,
    reg_max_scale = reg_max_scale,
    ld_reference = smixcan_ld_reference
  )

  list(
    p_s_sep_1 = smixcan$p_1_sep,
    p_s_sep_2 = smixcan$p_2_sep,
    p_s_sep = smixcan$p_sep,
    p_m_sep_1 = mixcan_sep$cell1_p,
    p_m_sep_2 = mixcan_sep$cell2_p,
    p_m_sep = mixcan_sep$p_combined,
    p_s_join_1 = smixcan$p_1_join,
    p_s_join_2 = smixcan$p_2_join,
    p_s_join = smixcan$p_join,
    p_m_join_1 = mixcan_join$cell1_p,
    p_m_join_2 = mixcan_join$cell2_p,
    p_m_join = mixcan_join$p_combined,
    p_predixcan = predixcan$p,
    m_mode = mixcan_join$mode,
    s_reg_scale = smixcan$reg_scale_selected,
    s_reg_condition = smixcan$reg_condition,
    s_mode = smixcan$mode
  )
}

# Main function: run simulation studies for large data.
run_simulation <- function(Data_batch,
                           sim_result_10,
                           i,
                           family = "gaussian",
                           regularization = "fixed",
                           reg_scale = 0.1,
                           condition_max = 100,
                           reg_min_scale = 0.001,
                           reg_max_scale = 0.1,
                           alpha = 0.5,
                           nfolds = 10,
                           workers = 1,
                           smixcan_ld_reference = c("train", "test")) {
  smixcan_ld_reference <- match.arg(smixcan_ld_reference)
  run_one_safe <- function(j) {
    tryCatch(
      run_one_simulation(
        sim_data = Data_batch[[j]],
        family = family,
        regularization = regularization,
        reg_scale = reg_scale,
        condition_max = condition_max,
        reg_min_scale = reg_min_scale,
        reg_max_scale = reg_max_scale,
        alpha = alpha,
        nfolds = nfolds,
        smixcan_ld_reference = smixcan_ld_reference
      ),
      error = function(e) {
        message("Simulation batch item ", j, " failed: ", e$message)
        NULL
      }
    )
  }

  if (workers > 1) {
    sim_values_list <- parallel::mclapply(seq_len(B), run_one_safe, mc.cores = workers)
  } else {
    sim_values_list <- lapply(seq_len(B), run_one_safe)
  }

  for (j in seq_len(B)) {
    row_id <- (i - 1) * B + j
    sim_values <- sim_values_list[[j]]
    if (!is.null(sim_values)) {
      sim_result_10 <- write_if_column_exists(sim_result_10, row_id, sim_values)
    }
  }
  sim_result_10
}


# Prediction Function
train_prediction_model_latest <- function(y.train, x.train, pi.train,
                                          alpha = 1,
                                          nfolds = 10,
                                          weight_eps = 1e-12) {
  set.seed(1L)
  foldid <- sample(seq_len(nfolds), nrow(as.matrix(x.train)), replace = TRUE)
  fit <- SMiXcan::MiXcan_train_K(
    y = y.train,
    x = x.train,
    pi_k = make_pi_k(pi.train),
    foldid = foldid,
    alpha = alpha
  )

  W_all <- as.matrix(fit$W)
  keep <- rowSums(abs(W_all), na.rm = TRUE) > weight_eps
  if (!any(keep)) {
    keep <- rep(TRUE, nrow(W_all))
  }

  list(
    W1 = W_all[keep, 1],
    W2 = W_all[keep, 2],
    W = W_all[keep, , drop = FALSE],
    selected_snp = which(keep),
    tissue_weights = as.numeric(fit$beta.all.models[-1, "Tissue"]),
    tissue_intercept = as.numeric(fit$beta.all.models[1, "Tissue"])
  )
}

# Association Analysis Function - MiXcan Separate
run_mixcan_separate <- function(outcome, cell1, cell2, family) {
  dat <- stats::na.omit(data.frame(
    outcome = as.numeric(outcome),
    cell1 = as.numeric(cell1),
    cell2 = as.numeric(cell2)
  ))
  if (nrow(dat) == 0) {
    return(list(cell1_p = NA_real_, cell2_p = NA_real_, p_combined = NA_real_))
  }

  fit_one <- function(x) {
    fit <- tryCatch(
      stats::glm(dat$outcome ~ x, family = family),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NA_real_)
    }
    cs <- summary(fit)$coefficients
    if (nrow(cs) < 2 || !"Std. Error" %in% colnames(cs)) {
      return(NA_real_)
    }
    z <- cs[2, "Estimate"] / cs[2, "Std. Error"]
    2 * stats::pnorm(abs(z), lower.tail = FALSE)
  }

  p1 <- fit_one(scale(dat$cell1))
  p2 <- fit_one(scale(dat$cell2))
  list(
    cell1_p = p1,
    cell2_p = p2,
    p_combined = SMiXcan::safe_ACAT(c(p1, p2))
  )
}

# Association Analysis Function - MiXcan
run_mixcan_latest <- function(outcome, cell1, cell2, family) {
  fit <- tryCatch(
    SMiXcan::MiXcan_assoc_test(
      outcome = outcome,
      cell1 = cell1,
      cell2 = cell2,
      family = family
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(cell1_p = NA_real_, cell2_p = NA_real_, p_combined = NA_real_, mode = "error"))
  }
  fit
}

# Association Analysis Function - PrediXcan
run_predixcan_latest <- function(outcome, x.test, tissue_weights,
                                 tissue_intercept = 0, family) {
  tissue_weights <- as.numeric(tissue_weights)
  keep <- is.finite(tissue_weights) & abs(tissue_weights) > 1e-12
  if (!any(keep)) {
    return(list(p = NA_real_))
  }

  pred_expr <- as.numeric(tissue_intercept + as.matrix(x.test[, keep, drop = FALSE]) %*% tissue_weights[keep])
  dat <- stats::na.omit(data.frame(
    outcome = as.numeric(outcome),
    pred_expr = as.numeric(scale(pred_expr))
  ))
  if (nrow(dat) == 0 || stats::var(dat$pred_expr) == 0) {
    return(list(p = NA_real_))
  }

  fit <- tryCatch(
    stats::glm(dat$outcome ~ pred_expr, family = family),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(p = NA_real_))
  }

  cs <- summary(fit)$coefficients
  if (!"pred_expr" %in% rownames(cs)) {
    return(list(p = NA_real_))
  }
  z <- cs["pred_expr", "Estimate"] / cs["pred_expr", "Std. Error"]
  list(p = 2 * stats::pnorm(abs(z), lower.tail = FALSE))
}

# Association Analysis Function - S-MiXcan
run_smixcan_latest <- function(W,
                               selected_snp,
                               sim_data,
                               gwas_results,
                               family,
                               n0 = NULL,
                               n1 = NULL,
                               regularization = "fixed",
                               reg_scale = 0.1,
                               condition_max = 100,
                               reg_min_scale = 0.001,
                               reg_max_scale = 0.1,
                               ld_reference = c("train", "test")) {
  ld_reference <- match.arg(ld_reference)
  x_ref <- if (ld_reference == "test") sim_data$x.test else sim_data$x.train
  x_g <- as.matrix(x_ref[, selected_snp, drop = FALSE])
  W <- as.matrix(W)
  colnames(W) <- paste0("Cell", seq_len(ncol(W)))

  fit <- SMiXcan::SMiXcan_assoc_test_K(
    W = W,
    gwas_results = gwas_results,
    x_g = x_g,
    n0 = n0,
    n1 = n1,
    family = family,
    regularization = regularization,
    reg_scale = reg_scale,
    condition_max = condition_max,
    reg_min_scale = reg_min_scale,
    reg_max_scale = reg_max_scale
  )

  list(
    p_1_sep = fit$p_sep[1],
    p_2_sep = fit$p_sep[2],
    p_sep = SMiXcan::safe_ACAT(fit$p_sep),
    p_1_join = fit$p_join_vec[1],
    p_2_join = fit$p_join_vec[2],
    p_join = fit$p_join,
    reg_scale_selected = fit$reg_scale_selected,
    reg_condition = fit$reg_condition,
    mode = fit$mode
  )
}


# Helper function
sim_center_binom <- function(sim_data){
  sim_data$y.train = sim_data$y.train - mean(sim_data$y.train)
  sim_data$y1.train = sim_data$y1.train - mean(sim_data$y1.train)
  sim_data$y2.train = sim_data$y2.train - mean(sim_data$y2.train)

  sim_data$y.test = sim_data$y.test - mean(sim_data$y.test)
  sim_data$y1.test = sim_data$y1.test - mean(sim_data$y1.test)
  sim_data$y2.test = sim_data$y2.test - mean(sim_data$y2.test)

  return(sim_data)
}

# Helper function
make_pi_k <- function(pi.train) {
  if (is.null(dim(pi.train))) {
    pi_vec <- as.numeric(pi.train)
    return(cbind(Cell1 = pi_vec, Cell2 = 1 - pi_vec))
  }
  as.matrix(pi.train)
}

# Helper function
write_if_column_exists <- function(result, row_id, values) {
  for (nm in names(values)) {
    if (nm %in% colnames(result)) {
      if (is.null(values[[nm]]) || length(values[[nm]]) == 0L) {
        next
      }
      result[row_id, nm] <- values[[nm]]
    }
  }
  result
}

