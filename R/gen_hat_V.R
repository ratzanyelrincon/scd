#' Precompute all panel quantities
#'
#' @param w Vector of weights (K x 1).
#' @param data Full data frame sorted by (id, t).
#' @param hat_Gamma_pre Pre-treatment DIFFERENCED means matrix (T_0 x K).
#' @param lambda Vector of lambdas (T_0).
#' @param B2 Orthonormal basis matrix (K x (K-1)).
#' @param group_means Matrix of group means ((K+1) x T_total).
#' @param n_jt Matrix of weighted sums ((K+1) x T_total).
#' @param K Number of donors.
#' @param n Total sample size.
#' @param T_star First post-treatment period.
#' @param T_total Total periods.
#'
#' @return List with all precomputed quantities
#' @export
precompute_panel_all <- function(w, data, hat_Gamma_pre, lambda, B2,
                                 group_means, n_jt, K, n, T_star, T_total) {

  T_0 <- T_star - 1

  # ============================================================
  # STEP 1: Sort and validate
  # ============================================================
  data <- data[order(data$id, data$t), ]
  ids <- unique(data$id)
  N <- length(ids)

  if (nrow(data) != N * T_total) {
    stop("Unbalanced panel: expected ", N * T_total, " rows, got ", nrow(data))
  }

  # ============================================================
  # STEP 2: Compute psi* for ALL observations - O(n)
  #
  # TIME-VARYING WEIGHTS: absorb pi_{i,t} into the influence function:
  #   psi*_{i,t} = pi_{i,t} * (n / n_{G_i,t}) * (Y_{i,t} - m_{G_i,t})
  #
  # This generalizes the constant-weight case (where pi_i factors
  # out) and allows weights to vary across periods.
  # ============================================================
  G_vec <- data$G
  t_vec <- data$t
  Y_vec <- data$Y
  wt_vec <- data$weight
  G_idx <- G_vec + 1L

  lin_idx <- G_idx + (t_vec - 1L) * (K + 1L)

  psi_star <- wt_vec * (n / n_jt[lin_idx]) * (Y_vec - group_means[lin_idx])

  # ============================================================
  # STEP 3: Reshape to N x T_total matrices
  # ============================================================
  Psi_star_mat <- matrix(psi_star, nrow = N, ncol = T_total, byrow = TRUE)

  # Get group for each individual (first observation per individual)
  # Data is sorted by (id, t), so individual i's first obs is at row (i-1)*T_total + 1
  G_by_id <- G_vec[seq(1L, N * T_total, by = T_total)]

  # ============================================================
  # STEP 4: Compute DIFFERENCED psi for Z computation
  # psi_{i,t} = psi*_{i,t} - sum_s lambda_s * psi*_{i,s}
  # ============================================================
  Psi_pre <- Psi_star_mat[, 1:T_0, drop = FALSE]  # N x T_0

  # diff_term_z = sum_s lambda_s * psi*_{i,s} for each individual
  diff_term_z <- as.vector(Psi_pre %*% lambda)  # N x 1

  # Psi_diff[i,t] = psi*_{i,t} - diff_term_z[i] for pre-treatment periods
  Psi_diff_pre <- Psi_pre - diff_term_z  # Broadcasting: N x T_0

  # ============================================================
  # STEP 5: Compute Z_i for hat_V using DIFFERENCED psi
  # z_i = alpha_i * (1/T_0) sum_t mu_t * psi_{i,G_i,t}
  # ============================================================
  alpha_lookup <- c(-1, w)
  alpha_vec <- alpha_lookup[G_by_id + 1L]

  # Weight by alpha and compute Z_raw
  Psi_diff_weighted <- alpha_vec * Psi_diff_pre  # N x T_0

  # Z_raw = (1/T_0) * Psi_diff_weighted %*% hat_Gamma_pre
  # [N x T_0] * [T_0 x K] = [N x K]
  Z_raw <- (1 / T_0) * (Psi_diff_weighted %*% hat_Gamma_pre)

  # Transform to (K-1) space
  Z_B2 <- Z_raw %*% B2  # N x (K-1)

  # ============================================================
  # STEP 6: Compute components for psi_{i,t,theta} (sigma_squared)
  # psi_{i,t,theta} = beta_i * psi_{i,G_i,t}
  # where psi_{i,G_i,t} = psi*_{i,t} - sum_s lambda_s * psi*_{i,s}
  # ============================================================
  beta_lookup <- c(1, -w)
  beta_vec <- beta_lookup[G_by_id + 1L]

  # For sigma_squared, we need psi_{i,t} for ALL periods (including post)
  # For post-treatment: psi_{i,t} = psi*_{i,t} - sum_s lambda_s * psi*_{i,s}
  # (same formula, but t > T_0)

  # Psi_theta_raw stores beta_i * psi*_{i,t} before differencing
  Psi_theta_raw <- beta_vec * Psi_star_mat  # N x T_total

  # diff_term for theta = beta_i * sum_s lambda_s * psi*_{i,s}
  diff_term_theta <- beta_vec * diff_term_z  # N x 1

  # ============================================================
  # RETURN
  # ============================================================
  return(list(
    Z_B2 = Z_B2,                        # N x (K-1) - using DIFFERENCED psi
    Psi_theta_raw = Psi_theta_raw,       # N x T_total - beta_i * psi*_{i,t}
    diff_term_theta = diff_term_theta,   # N x 1 - beta_i * sum_s lambda_s psi*_{i,s}
    N = N, K = K, n = n,
    T_0 = T_0, T_total = T_total,
    B2 = B2,
    is_panel = TRUE
  ))
}


#' Generate hat_V for Panel Data
#'
#' Since pi_{i,t} is absorbed into psi*, the variance is simply:
#'   hat_V = (1/n) * B2' (sum_i z_i z_i') B2 = (1/n) * crossprod(Z_B2)
#'
#' @param precomp List from precompute_panel_all
#' @return hat_V matrix ((K-1) x (K-1))
#' @export
gen_hat_V_panel_fast <- function(precomp) {
  # Z_B2 already contains pi_{i,t} inside each psi*_{i,t}
  return((1 / precomp$n) * crossprod(precomp$Z_B2))
}


#' Generate hat_sigma_squared for Panel Data
#'
#' Since pi_{i,t} is absorbed into psi*, the variance is simply:
#'   hat_sigma^2_t = (1/n) sum_i psi_{i,t,theta}^2
#'
#' @param precomp List from precompute_panel_all
#' @return hat_sigma_squared (T_total x 1 matrix)
#' @export
gen_hat_sigma_squared_panel_fast <- function(precomp) {
  Psi_theta <- precomp$Psi_theta_raw - precomp$diff_term_theta  # N x T_total
  return(matrix(
    colSums(Psi_theta^2) / precomp$n,
    ncol = 1
  ))
}


# ============================================================================
# RC (REPEATED CROSS-SECTIONS) ESTIMATORS
# ============================================================================

#' Generate hat_V for Repeated Cross-Sections
#'
#' Updated to absorb pi_{i,t} into the influence function for consistency.
#'
#' @export
gen_hat_V_RC <- function(w, hat_Gamma_pre, data, lambda, B2, group_means, n_jt, K, n, T_star) {
  T_0 <- T_star - 1
  bar_mu <- colMeans(hat_Gamma_pre)
  weight_vec <- c(1, -w)  # treated = 1, donor k = -w[k]
  hat_V_sum <- matrix(0, nrow = K-1, ncol = K-1)

  for (pre_t in 1:T_0) {
    idx_t <- which(data$t == pre_t)
    if (length(idx_t) == 0) next

    G_idx <- data$G[idx_t] + 1L
    lin_idx <- G_idx + (pre_t - 1L) * (K + 1L)

    # psi*_{i,t} = pi_{i,t} * (n / n_{j,t}) * (Y - m_{j,t})
    psi_star <- data$weight[idx_t] * (n / n_jt[lin_idx]) * (data$Y[idx_t] - group_means[lin_idx])

    # sum_j w_j psi*_{ij,t} - psi*_{i0,t}
    psi_wt <- psi_star * weight_vec[G_idx]

    # (1/n) sum_i psi_wt^2  (no separate pi^2 needed)
    hat_V_it <- sum(psi_wt^2) / n

    # M_t matrix
    mu_t <- hat_Gamma_pre[pre_t, ]
    M_t <- (1/T_0) * tcrossprod(mu_t) -
      lambda[pre_t] * (tcrossprod(bar_mu, mu_t) + tcrossprod(mu_t, bar_mu)) +
      T_0 * lambda[pre_t]^2 * tcrossprod(bar_mu)

    hat_V_sum <- hat_V_sum + hat_V_it * crossprod(B2, M_t) %*% B2
  }
  return((1/T_0) * hat_V_sum)
}


#' Generate hat_sigma_squared for Repeated Cross-Sections
#'
#' Updated to absorb pi_{i,t} into the influence function for consistency.
#'
#' @export
gen_hat_sigma_squared_RC <- function(w, data, lambda, group_means, n_jt, K, n, T_star, L) {
  T_0 <- T_star - 1
  T_total <- T_star + L

  # delta[t] = lambda[t] for t <= T_0, 0 otherwise
  delta <- c(lambda, rep(0, T_total - T_0))

  valid <- data$t >= 1 & data$t <= T_total
  G_idx <- data$G[valid] + 1L
  t_vals <- data$t[valid]
  lin_idx <- G_idx + (t_vals - 1L) * (K + 1L)

  # psi*_{i,t} = pi_{i,t} * (n / n_jt) * (Y - m)
  psi_star <- data$weight[valid] * (n / n_jt[lin_idx]) * (data$Y[valid] - group_means[lin_idx])

  # beta_i = 1 for treated (G=0), -w_k for donor k
  wt_mult <- ifelse(G_idx == 1L, 1, -w[G_idx - 1L])

  # (beta_i * psi*_{i,t})^2
  base2 <- (psi_star * wt_mult)^2

  # sum_vec[t] = sum_{i in N_t} beta_i^2 * psi*^2
  sum_by_t <- tapply(base2, t_vals, sum, default = 0)
  sum_vec <- numeric(T_total)
  sum_vec[as.integer(names(sum_by_t))] <- sum_by_t

  # total_const = sum_s delta_s^2 * sum_vec[s]
  total_const <- sum(base2 * delta[t_vals]^2)

  # sigma_t^2 = (1/n) * [total_const + (1 - 2*delta_t) * sum_vec[t]]
  return(matrix((1/n) * (total_const + sum_vec * (1 - 2*delta)), ncol = 1))
}
