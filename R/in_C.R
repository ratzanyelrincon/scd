#' Check Membership in Confidence Set C_{1-kappa}
#'
#' The optimization problem:
#'   min_{r} (phi - r)' P (phi - r)  subject to w'r = 0, r >= 0
#'
#' Rewritten in standard OSQP form:
#'   min_{r} (1/2) r' (2P) r - (2P*phi)' r
#'   subject to: [w'; -I] r <= [0; 0]  and  [w'; I] r >= [0; 0]
#'   i.e., w'r = 0, r >= 0
#'
#' @param w Vector of weights (K x 1).
#' @param kappa Kappa level.
#' @param hat_H hat_H matrix.
#' @param hat_h hat_h vector.
#' @param precomp Precomputed matrix B2 %*% solve(hat_V) %*% t(B2).
#' @param n Total sample size.
#' @param K Number of untreated groups.
#' @param tolerance Numerical tolerance.
#'
#' @return Logical: TRUE if in C, FALSE otherwise.
#' @export
in_C_osqp <- function(w, kappa, hat_H, hat_h, precomp, n, K, tolerance) {

  # Compute deviation
  phi_diff <- hat_H %*% w - hat_h

  if (all(w > 1e-8)) {
    # Interior point: r = 0 is optimal
    hat_r <- rep(0, K)
  } else {
    # Boundary point: solve with OSQP
    #
    # OSQP solves: min (1/2) x' P x + q' x
    #              subject to l <= A x <= u
    #
    # Our problem: min (phi - r)' precomp (phi - r)
    #            = r' precomp r - 2 phi' precomp r + const
    #
    # So: P_osqp = 2 * precomp (PSD, rank K-1)
    #     q_osqp = -2 * precomp %*% phi
    #
    # Constraints: w'r = 0, r >= 0
    # In OSQP form: A = [w', I_K]'
    #               l = [0, 0, ..., 0]'
    #               u = [0, Inf, ..., Inf]'

    P_osqp <- 2 * precomp
    q_osqp <- -2 * as.vector(precomp %*% phi_diff)

    # Constraint matrix: first row is w' (equality), next K rows are I (r >= 0)
    A_osqp <- rbind(w, diag(K))
    l_osqp <- c(0, rep(0, K))
    u_osqp <- c(0, rep(Inf, K))

    # Convert to sparse matrix format for osqp
    P_sparse <- Matrix::Matrix(P_osqp, sparse = TRUE)
    A_sparse <- Matrix::Matrix(A_osqp, sparse = TRUE)

    # OSQP settings
    settings <- osqp::osqpSettings(
      verbose = FALSE,
      eps_abs = 1e-8,
      eps_rel = 1e-8,
      max_iter = 10000,
      polish = TRUE  # Improves solution accuracy
    )

    # Solve
    result <- tryCatch({
      model <- osqp::osqp(P = P_sparse, q = q_osqp, A = A_sparse,
                          l = l_osqp, u = u_osqp, pars = settings)
      model$Solve()
    }, error = function(e) {
      warning(paste("OSQP failed in in_C_osqp:", e$message))
      return(NULL)
    })

    if (is.null(result) || result$info$status_val != 1) {
      # OSQP failed or didn't converge - use fallback
      # Conservative: return FALSE (exclude from confidence set)
      return(FALSE)
    }

    hat_r <- result$x
    hat_r <- pmax(hat_r, 0)  # Numerical cleanup
  }

  # Compute test statistic
  diff_vec <- phi_diff - hat_r
  T_stat <- as.numeric(n * (t(diff_vec) %*% precomp %*% diff_vec))

  # Compute degrees of freedom using the A42 definition of gamma
  gamma_vec <- precomp %*% diff_vec


  # Find binding constraints: w is effectively 0 AND gamma is effectively 0
  hat_d <- sum(abs(gamma_vec) < tolerance & w ==0)

  # Degrees of freedom (Note: K - 1 for the Bonferroni method per A42)
  hat_k <- max(K - 1 - hat_d, 1)

  # Critical value
  hat_c_1_kappa <- qchisq(1 - kappa, df = hat_k)

  return(T_stat <= hat_c_1_kappa)
}
