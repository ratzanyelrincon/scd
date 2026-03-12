#' Synthetic Control with Differencing Estimator
#'
#' @param data Data frame containing the panel data.
#' @param yname Name of the outcome variable.
#' @param dname Name of the treatment indicator variable.
#' @param tname Name of the time variable.
#' @param gname Name of the group variable.
#' @param idname Name of the individual identifier variable (required for panel data).
#' @param gtreated Identifier for the treated group.
#' @param guntreated Vector of identifiers for untreated (donor) groups.
#' @param t_labels Labels for time periods.
#' @param guntreated_labels Labels for donor groups.
#' @param wname Name of the weight variable (optional).
#' @param lambda Lambda specification ("DID", "Unif", "SC", or custom vector).
#' @param data_type Data structure type ("auto", "panel", or "RC"). Default "auto" detects automatically.
#' @param L Number of post-treatment periods (0 = all available).
#' @param wdecimals Decimal places for rounding weights.
#' @param radius Radius for cloud point generation around hat_w.
#' @param seed Random seed.
#' @param alpha Significance level.
#' @param kappa Kappa level for Bonferroni method.
#' @param length_grid_w Number of weight grid points.
#' @param n_wpoints_per_edge Points per simplex edge.
#' @param tolerance Tolerance for binding constraint detection.
#' @param parallel Whether to use parallel computation.
#' @param post_treatment_only Only compute CIs for post-treatment periods (default TRUE).
#'
#' @return List containing estimation results.
#' @export
att_scd <- function(data, yname, dname, tname, gname, gtreated, guntreated,
                    t_labels, guntreated_labels, idname = NULL, wname = NULL,
                    lambda = "DID",
                    data_type = c("auto", "panel", "RC"),
                    L = 0, wdecimals = 10,
                    radius = 0.1, seed = 123456, alpha = 0.05, kappa = 0.005,
                    length_grid_w = 500, n_wpoints_per_edge = 20,
                    tolerance = 1e-04, parallel = TRUE,
                    post_treatment_only = TRUE) {

  data_type <- match.arg(data_type)

  # Input validation
  if (!is.data.frame(data)) stop("Error: 'data' must be a data.frame.")
  required_vars <- c(yname, dname, tname, gname)
  if (!all(required_vars %in% names(data))) stop("Error: Missing required variables.")
  if (!is.null(wname) && !wname %in% names(data)) stop("Error: Weight variable not found.")
  if (!is.null(idname) && !idname %in% names(data)) stop("Error: ID variable not found.")
  if (length(gtreated) != 1) stop("Error: 'gtreated' must be a single identifier.")

  all_groups <- unique(data[[gname]])
  if (!gtreated %in% all_groups) stop("Error: Treated group not found.")
  if (anyDuplicated(guntreated) > 0) stop("Error: 'guntreated' must be unique.")
  if (any(guntreated %in% gtreated)) stop("Error: 'guntreated' cannot overlap with 'gtreated'.")
  if (!all(guntreated %in% all_groups)) stop("Error: Donor groups not found.")

  max_t <- max(data[[tname]], na.rm = TRUE)
  if (length(t_labels) != max_t) stop("Error: 't_labels' length mismatch.")
  if (length(guntreated_labels) != length(guntreated)) stop("Error: 'guntreated_labels' length mismatch.")

  start_time <- Sys.time()
  set.seed(seed)

  # Rename variables
  data <- data %>%
    dplyr::rename(Y = !!dplyr::sym(yname), D = !!dplyr::sym(dname),
                  G_org = !!dplyr::sym(gname), t = !!dplyr::sym(tname))

  if (!is.null(idname)) {
    data <- data %>% dplyr::rename(id = !!dplyr::sym(idname))
  }

  if (is.null(wname)) { data$weight <- 1 } else { data <- data %>% dplyr::rename(weight = !!dplyr::sym(wname)) }

  if (!is.null(idname)) {
    data <- data %>% dplyr::select(id, Y, D, G_org, t, weight)
  } else {
    data <- data %>% dplyr::select(Y, D, G_org, t, weight)
  }

  if (all(data$D == 0)) stop("Error: No treatment observed.")

  # Filter data to only include treated and specified donor groups
  data <- data %>% dplyr::filter(G_org %in% c(gtreated, guntreated))
  if (nrow(data) == 0) stop("Error: No data remaining after filtering to specified groups.")

  # Remap groups: treated -> 0, donors -> 1..K (user-supplied order)
  group_order <- c(gtreated, guntreated)

  observed_groups <- unique(data$G_org)
  missing_groups <- setdiff(group_order, observed_groups)
  if (length(missing_groups) > 0) {
    stop(paste0("Error: The following supplied groups have no observations after filtering: ",
                paste(missing_groups, collapse = ", ")))
  }

  data$G <- match(data$G_org, group_order) - 1L
  if (any(is.na(data$G))) stop("Error: Group remapping failed.")
  data <- data %>% dplyr::arrange(t, G)

  if (nrow(data) == 0) stop("Error: Empty data.")

  T_total <- max(data$t, na.rm = TRUE)
  treated_times <- unique(data$t[data$D == 1 & data$G == 0])
  if (length(treated_times) == 0) stop("Error: No treatment for treated group.")
  T_star <- min(treated_times, na.rm = TRUE)
  T_0 <- T_star - 1

  full_post_len <- T_total - T_star + 1
  if (L == 0) L <- full_post_len - 1
  T_1 <- L + 1

  # =========== DETECT DATA TYPE ===========
  if (data_type == "panel" && is.null(idname)) {
    stop("Error: data_type='panel' requires 'idname' to be specified. Please provide the name of the individual ID variable.")
  }

  if (data_type == "auto") {
    is_panel <- detect_panel_data(data, idname, T_total)
  } else if (data_type == "panel") {
    # User forced panel mode - verify it's actually panel data
    is_panel <- detect_panel_data(data, idname, T_total)
    if (!is_panel) {
      stop("Error: data_type='panel' specified but data is not a balanced panel. Each individual must appear exactly once in each time period.")
    }
  } else {
    # data_type == "RC"
    is_panel <- FALSE
  }

  if (is_panel && is.null(idname)) {
    stop("Error: Panel data detected but 'idname' not provided. Please specify the individual ID variable.")
  }

  message(paste0("Data type: ", ifelse(is_panel, "Balanced Panel", "Repeated Cross-Sections")))

  # =========== CRITICAL: CORRECT DEFINITION OF n ===========
  # Paper (line 866): "We let N = ∪_t N_t and n = |N|"
  # - For panel: N_t = N for all t, so n = |N| = number of UNIQUE INDIVIDUALS
  # - For RC: N_t disjoint, so n = Σ|N_t| = total observations
  #
  # This affects the influence function formula: psi* = (n/n_jt) * (Y-m)
  # and the variance formula: V = (1/n) * Σ_i z_i z_i'
  # ============================================================
  if (is_panel) {
    # Panel: n = number of unique individuals
    n <- length(unique(data$id))
    message(paste0("Panel n = ", n, " unique individuals"))
  } else {
    # RC: n = total observations (sum of weights for survey data)
    n <- sum(data$weight)
    message(paste0("RC n = ", n, " total observations"))
  }

  if (is.character(lambda)) {
    lambda_vec <- if (lambda == "DID") { c(rep(0, T_0 - 1), 1) }
    else if (lambda == "Unif") { rep(1 / T_0, T_0) }
    else if (lambda == "SC") { rep(0, T_0) }
    else stop("Unknown lambda")
  } else { lambda_vec <- lambda }

  K <- length(guntreated)
  if (K > T_0) stop("Error: K > T_0.")

  means_ordered <- collapse::collap(data, Y ~ t + G, w = ~ weight) %>% dplyr::arrange(G, t) %>% dplyr::pull(Y)
  if (length(means_ordered) != T_total * (K + 1)) stop("Error: Incomplete means.")

  group_sample_averages <- cbind(
    matrix(means_ordered[1:T_total], ncol = 1),
    matrix(means_ordered[(T_total + 1):((K + 1) * T_total)], nrow = T_total, ncol = K)
  )

  pre <- group_sample_averages[1:T_0, , drop = FALSE]
  post <- group_sample_averages[T_star:(T_star + T_1 - 1), , drop = FALSE]

  within_transf_pars <- as.vector(lambda_vec %*% pre)
  hat_gamma_pre <- pre[, 1] - within_transf_pars[1]
  hat_gamma_post <- post[, 1] - within_transf_pars[1]

  rep_untr_pre <- matrix(within_transf_pars[-1], nrow = T_0, ncol = K, byrow = TRUE)
  rep_untr_post <- matrix(within_transf_pars[-1], nrow = T_1, ncol = K, byrow = TRUE)
  hat_Gamma_pre <- pre[, -1, drop = FALSE] - rep_untr_pre
  hat_Gamma_post <- post[, -1, drop = FALSE] - rep_untr_post

  if (Matrix::rankMatrix(hat_Gamma_pre)[1] != K) stop("Error: hat_Gamma_pre not full rank.")

  obj_func <- function(w) {
    diff_vec <- hat_gamma_pre - hat_Gamma_pre %*% w
    list(objective = crossprod(diff_vec), gradient = -2 * crossprod(hat_Gamma_pre, diff_vec))
  }
  eq_constr <- function(w) list(constraints = sum(w) - 1, jacobian = rep(1, length(w)))

  opt_res <- nloptr::nloptr(
    x0 = rep(1/K, K), eval_f = obj_func, eval_g_eq = eq_constr,
    lb = rep(0, K), ub = rep(1, K),
    opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-10, ftol_abs = 1e-14, maxeval = 100000, print_level = 0)
  )

  hat_w <- round(opt_res$solution, wdecimals)
  if (abs(sum(hat_w) - 1) > 1e-5) {
    hat_w[which.max(hat_w)] <- hat_w[which.max(hat_w)] + (1 - sum(hat_w))
  }

  hat_theta_scd_pre <- hat_gamma_pre - hat_Gamma_pre %*% hat_w
  hat_theta_scd <- hat_gamma_post - hat_Gamma_post %*% hat_w
  hat_theta_full <- c(hat_theta_scd_pre, hat_theta_scd)

  hat_H <- (1 / T_0) * crossprod(hat_Gamma_pre)
  hat_h <- (1 / T_0) * crossprod(hat_Gamma_pre, hat_gamma_pre)

  n_jt_df <- data %>% dplyr::group_by(G, t) %>% dplyr::summarise(n_jt_val = sum(weight, na.rm = TRUE), .groups = "drop") %>% dplyr::arrange(t, G)
  if(any(n_jt_df$n_jt_val == 0)) stop("Error: Zero weight cells.")
  n_jt <- matrix(n_jt_df %>% dplyr::arrange(t, G) %>% dplyr::pull(n_jt_val), nrow = K + 1, ncol = T_total)
  group_means <- matrix(collapse::collap(data, Y ~ t + G, w = ~weight) %>% dplyr::arrange(t, G) %>% dplyr::pull(Y), nrow = K + 1, ncol = T_total)

  I_K <- diag(K)
  M_mat <- I_K - (1 / K) * matrix(1, K, K)
  eig <- eigen(M_mat, symmetric = TRUE)
  B2 <- eig$vectors[, -which(abs(eig$values) < 1e-10), drop = FALSE]


  # =========== COMPUTE VARIANCE ESTIMATORS ===========
  panel_precomp <- NULL
  if (is_panel) {
    panel_precomp <- precompute_panel_all(hat_w, data, hat_Gamma_pre, lambda_vec, B2,
                                          group_means, n_jt, K, n, T_star, T_total)
    hat_V <- gen_hat_V_panel_fast(panel_precomp)
    hat_sigma_squared <- gen_hat_sigma_squared_panel_fast(panel_precomp)
  } else {
    hat_V <- gen_hat_V_RC(hat_w, hat_Gamma_pre, data %>% dplyr::filter(t <= T_star - 1),
                          lambda_vec, B2, group_means, n_jt, K, n, T_star)
    hat_sigma_squared <- gen_hat_sigma_squared_RC(hat_w, data, lambda_vec, group_means,
                                                  n_jt, K, n, T_star, L)
  }

  # No regularization needed - OSQP handles PSD matrices
  if (det(hat_V) <= 0 || any(eigen(hat_V)$values <= 0)) stop("Error: hat_V not positive definite.")
  precomp_bf <- B2 %*% solve(hat_V) %*% t(B2)

  cl <- NULL
  if (parallel) {
    n_cores <- parallel::detectCores(logical = TRUE)
    cl <- parallel::makeCluster(n_cores)
    on.exit({if (!is.null(cl)) parallel::stopCluster(cl)}, add = TRUE)
  }

  # Grid generation - put hat_w first for better early acceptance
  if (K > 1) {
    ij_pairs <- expand.grid(i = 1:(K - 1), j = 2:K)
    ij_pairs <- ij_pairs[ij_pairs$j > ij_pairs$i, ]
    edge_pts <- do.call(rbind, lapply(1:nrow(ij_pairs), function(k) {
      ij <- ij_pairs[k, ]; i <- ij$i; j <- ij$j
      v_i <- rep(0, K); v_i[i] <- 1; v_j <- rep(0, K); v_j[j] <- 1
      t(sapply(seq(0, 1, length.out = n_wpoints_per_edge), function(tt) tt * v_i + (1 - tt) * v_j))
    }))
  } else { edge_pts <- matrix(numeric(0), nrow = 0, ncol = K) }

  cloud_pts <- t(replicate(length_grid_w, {
    candidate <- hat_w + rnorm(K, 0, radius)
    u <- sort(candidate, decreasing = TRUE)
    cssv <- cumsum(u)
    rho <- max(which(u > (cssv - 1) / seq_along(u)))
    theta <- (cssv[rho] - 1) / rho
    pmax(candidate - theta, 0)
  }))

  # Put hat_w first in grid
  grid_w <- unique(rbind(hat_w, cloud_pts, gtools::rdirichlet(length_grid_w, rep(1, K)), edge_pts))

  ci_matrix <- matrix(NA_real_, nrow = T_0 + T_1, ncol = 2)
  conf_set_size <- NA

  # Use OSQP-based in_C function
  in_C_func <- function(w) in_C_osqp(w, kappa, hat_H, hat_h, precomp_bf, n, K, tolerance)

  if (parallel) {
    parallel::clusterExport(cl, varlist = c("grid_w", "precomp_bf", "kappa", "hat_H", "hat_h",
                                            "n", "K", "tolerance",
                                            "in_C_osqp"), envir = environment())
    parallel::clusterEvalQ(cl, library(osqp))
    parallel::clusterEvalQ(cl, library(Matrix))

    in_set <- unlist(parallel::parLapply(cl, seq_len(nrow(grid_w)), function(i) {
      in_C_osqp(as.numeric(grid_w[i, ]), kappa, hat_H, hat_h, precomp_bf, n, K, tolerance)
    }))
  } else {
    in_set <- unlist(lapply(seq_len(nrow(grid_w)), function(i) in_C_func(as.numeric(grid_w[i, ]))))
  }

  conf_set <- grid_w[in_set, , drop = FALSE]
  conf_set_size <- nrow(conf_set)
  z_val <- qnorm(1 - (alpha - kappa) / 2)

  # Determine which periods to compute
  if (post_treatment_only) {
    periods_to_compute <- (T_0 + 1):(T_0 + T_1)
  } else {
    periods_to_compute <- seq_len(T_0 + T_1)
  }

  for (tt in periods_to_compute) {
    # Degenerate case: variance = 0 (reference period with DID)
    if (hat_sigma_squared[tt] < 1e-12) {
      ci_matrix[tt, ] <- c(0, 0)
      next
    }

    if (tt <= T_0) {
      hat_mu_t <- hat_Gamma_pre[tt, ]; hat_mu_0t <- hat_gamma_pre[tt]
    } else {
      hat_mu_t <- hat_Gamma_post[tt - T_0, ]; hat_mu_0t <- hat_gamma_post[tt - T_0]
    }
    if (nrow(conf_set) == 0L) {
      ci_matrix[tt, ] <- NA_real_
    } else {
      synth_vals <- conf_set %*% hat_mu_t
      diffs <- hat_mu_0t - synth_vals
      se_tt <- sqrt(hat_sigma_squared[tt] / n)
      ci_matrix[tt, 1] <- min(diffs) - z_val * se_tt
      ci_matrix[tt, 2] <- max(diffs) + z_val * se_tt
    }
  }

  output <- list(
    donorpool_weights = data.frame(donors = guntreated,
                                   label = guntreated_labels,
                                   hat_w = hat_w),
    target_parameter = data.frame(period = seq_len(T_0 + T_1), hat_theta_scd = hat_theta_full,
                                  lci_hat_theta_scd = ci_matrix[, 1], uci_hat_theta_scd = ci_matrix[, 2]),
    synthetic_outcome = data.frame(period = seq_len(T_0 + T_1), Y = group_means[1, seq_len(T_0 + T_1)],
                                   Y_synth = group_means[1, seq_len(T_0 + T_1)] - hat_theta_full),
    hat_sigma = sqrt(hat_sigma_squared), n = n, K = K, T_0 = T_0, T_1 = T_1, T_star = T_star, L = L,
    lambda = lambda, data_type = ifelse(is_panel, "panel", "RC"),
    post_treatment_only = post_treatment_only,
    conf_set_size = conf_set_size,
    time_elapsed = Sys.time() - start_time
  )
  return(output)
}


#' Detect if data is balanced panel
#'
#' Checks whether data is a balanced panel where each individual
#' appears exactly once in each time period.
#'
#' @param data Data frame with id, t columns.
#' @param idname Name of individual ID variable (before renaming).
#' @param T_total Total number of time periods.
#'
#' @return TRUE if balanced panel, FALSE otherwise.
#' @keywords internal
detect_panel_data <- function(data, idname, T_total) {
  # If no idname provided, cannot be panel

  if (is.null(idname) || !"id" %in% names(data)) {
    return(FALSE)
  }

  # Check if each individual appears exactly T_total times with distinct periods
  id_summary <- data %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      n_obs = dplyr::n(),
      n_periods = dplyr::n_distinct(t),
      .groups = "drop"
    )

  # Balanced panel requirements:
  # 1. Each individual has exactly T_total observations
  # 2. Each individual has exactly T_total distinct periods (no duplicates within period)
  all_have_all_periods <- all(id_summary$n_periods == T_total)
  no_duplicates <- all(id_summary$n_obs == id_summary$n_periods)

  return(all_have_all_periods && no_duplicates)
}
