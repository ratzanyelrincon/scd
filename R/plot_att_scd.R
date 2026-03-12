#' Plot SCD ATT Estimates with Confidence Intervals
#'
#' @param results List output from \code{\link{att_scd}}.
#' @param date_start Optional \code{\link[Date]{Date}} for x-axis (default: "2003-01-01" for CPS-like data).
#' @param breaks_every Integer for x-breaks (default: 24, biennial for monthly).
#' @param save_path Optional file path to save (default: NULL; prints to device). Supports PDF (default) or TikZ (".tex").
#' @param width,height Numeric dimensions (inches; default: 4.5 x 3 for AER panels).
#'
#' @return A \code{ggplot} object (invisibly); saves if \code{save_path} provided.
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline scale_x_continuous scale_color_manual labs theme_bw theme element_text element_blank margin
#' @importFrom lubridate seq mdy year
#' @importFrom grDevices pdf dev.off
#' @importFrom tikzDevice tikz
#' @export
plot_att_scd <- function(results,
                         date_start = as.Date("2003-01-01"),
                         breaks_every = 24,
                         save_path = NULL,
                         width = 4.5,
                         height = 3) {

  # Input validation (rigorous for replicability)
  if (!inherits(results, "list") || is.null(results$target_parameter) || is.null(results$T_star)) {
    stop("Error: 'results' must be valid output from att_scd() with target_parameter and T_star.")
  }
  if (!inherits(date_start, "Date")) {
    warning("Warning: 'date_start' not a Date; falling back to numeric x-axis.")
    date_start <- NULL
  }
  if (!is.numeric(breaks_every) || breaks_every <= 0) {
    stop("Error: 'breaks_every' must be a positive integer.")
  }

  # Construct plotting data: Pre/post factor (stable matching diagnostic)
  df <- results$target_parameter %>%
    mutate(
      treatment_period = factor(
        ifelse(period <= results$T_star - 1, "Pre", "Post"),
        levels = c("Pre", "Post")  # Enforce legend order
      )
    )

  # X-axis setup: Date-mapped if provided, else numeric
  if (!is.null(date_start)) {
    period_dates <- seq(date_start, by = "month", length.out = nrow(df))
    df$x_axis <- period_dates
    period_breaks <- seq(date_start, max(period_dates), by = paste(breaks_every, "months"))
    period_labels <- format(period_breaks, "%Y")
    x_scale <- scale_x_date(breaks = period_breaks, labels = period_labels)
  } else {
    df$x_axis <- df$period
    period_breaks <- seq(1, max(df$period), by = breaks_every)
    period_labels <- as.character(period_breaks)  # Numeric fallback
    x_scale <- scale_x_continuous(breaks = period_breaks, labels = period_labels)
  }

  # Core plot: ATT points + CIs, zero line, pre/post colors
  p <- ggplot(df, aes(x = x_axis, y = hat_theta_scd, color = treatment_period)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lci_hat_theta_scd, ymax = uci_hat_theta_scd), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    x_scale +
    scale_color_manual(values = c("Pre" = "steelblue", "Post" = "firebrick")) +
    labs(
      x = "",
      y = "",
      color = NULL
    ) +
    theme_bw(base_size = 7.5) +
    theme(
      legend.position = "bottom",
      legend.margin = margin(t = -5),
      legend.justification = "center",
      legend.text = element_text(size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  # Save handling: PDF or TikZ (for LaTeX integration)
  if (!is.null(save_path)) {
    if (grepl("\\.tex$", save_path)) {
      tikz(save_path, width = width, height = height)
      print(p)
      dev.off()
    } else {
      pdf(save_path, width = width, height = height)
      print(p)
      dev.off()
    }
    message("Plot saved to: ", save_path)
  } else {
    print(p)  # Interactive display
  }

  invisible(p)  # Return for further customization
}
