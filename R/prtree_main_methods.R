# Last revision: July/2026

#' @title Print method for PRTree objects
#'
#' @description
#' Prints a brief summary of a fitted Probabilistic Regression Tree model.
#'
#' @param x An object of class `prtree`, as returned by [pr_tree].
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the original object.
#'
#' @export
print.prtree <- function(x, ...) {
  cat("\n========================================================\n")
  cat("Probabilistic Regression Tree (PRTree)\n")
  cat("========================================================\n")

  # Model dimensions
  n_obs <- length(x$yhat)
  n_tn <- length(x$gamma)

  cat(sprintf("Number of observations: %d\n", n_obs))
  cat(sprintf("Number of features: %d\n", x$n_feat))
  cat(sprintf("Number of terminal nodes: %d\n", n_tn))
  cat(sprintf("Distribution: %s\n", x$dist$dist_name))

  # Show distribution parameter if applicable
  if (!is.null(x$dist$par_name) && x$dist$par_name != "par_value") {
    cat(sprintf(
      "  with %s = %.4f\n",
      x$dist$par_name,
      x$dist[[x$dist$par_name]]
    ))
  }

  # Check for test/validation sets using mse names (more robust)
  has_test <- "test" %in% names(x$mse) && !is.na(x$mse["test"])
  has_validation <- "validation" %in% names(x$mse) && !is.na(x$mse["validation"])

  # mse
  cat("\nMean Squared Error:\n")
  cat(sprintf("  Train: %.6f\n", x$mse["train"]))

  if (has_validation) {
    cat(sprintf("  Validation: %.6f\n", x$mse["validation"]))
    cat(sprintf("  Global: %.6f\n", x$mse["global"]))
  } else if (has_test) {
    cat(sprintf("  Test: %.6f\n", x$mse["test"]))
    cat(sprintf("  Global: %.6f\n", x$mse["global"]))
  }

  # Sigma
  cat("\nSelected sigma:\n")
  if (length(x$sigma) == 1) {
    cat(sprintf("  %.4f\n", x$sigma))
  } else {
    cat("  ", .format_vector(x$sigma, max_print = 5, digits = 4), "\n")
  }

  cat("========================================================\n")
  return(invisible(x))
}

#' @title Summary method for PRTree objects
#'
#' @description
#' Provides a detailed summary of a fitted Probabilistic Regression Tree model.
#'
#' @param object An object of class `prtree`, as returned by [pr_tree].
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list of class `summary.prtree` containing summary statistics.
#'
#' @importFrom stats median quantile sd
#'
#' @export
summary.prtree <- function(object, ...) {
  # Basic information
  n_obs <- object$n_obs
  n_feat <- object$n_feat
  n_tn <- length(object$gamma)

  # Check for test/validation sets using mse names (more robust)
  has_test <- "test" %in% names(object$mse) && !is.na(object$mse["test"])
  has_validation <- "validation" %in% names(object$mse) && !is.na(object$mse["validation"])

  # mse statistics
  mse <- object$mse

  # Sigma statistics
  sigma <- object$sigma

  # Terminal node statistics
  gamma <- object$gamma[1:n_tn]

  # Probability matrix statistics
  # remove names to perform calculations
  P_matrix <- object$P
  colnames_saved <- colnames(P_matrix)
  colnames(P_matrix) <- NULL

  # Probability matrix statistics - quartis em ordem
  P_min <- apply(P_matrix, 2, min, na.rm = TRUE)
  P_q25 <- apply(P_matrix, 2, quantile, probs = 0.25, na.rm = TRUE)
  P_median <- apply(P_matrix, 2, median, na.rm = TRUE)
  P_q75 <- apply(P_matrix, 2, quantile, probs = 0.75, na.rm = TRUE)
  P_max <- apply(P_matrix, 2, max, na.rm = TRUE)
  P_mean <- colMeans(P_matrix, na.rm = TRUE)
  P_sd <- apply(P_matrix, 2, sd, na.rm = TRUE)

  # Distribution info
  dist_name <- object$dist$dist_name
  dist_par <- if (object$dist$par_name != "par_value") {
    object$dist[[object$dist$par_name]]
  } else {
    NULL
  }

  # Create summary list
  result <- list(
    n_obs = n_obs,
    n_feat = n_feat,
    n_terminal_nodes = n_tn,
    has_test = has_test,
    has_validation = has_validation,
    distribution = list(
      name = dist_name,
      parameter = dist_par,
      parameter_name = object$dist$par_name
    ),
    fill_type = object$fill_type,
    mse = mse,
    sigma = sigma,
    terminal_nodes = data.frame(
      node = 1:n_tn,
      gamma = gamma,
      P_mean = P_mean,
      P_sd = P_sd,
      P_min = P_min,
      P_q25 = P_q25,
      P_median = P_median,
      P_q75 = P_q75,
      P_max = P_max
    )
  )

  class(result) <- "summary.prtree"
  return(result)
}


#' @title Print method for summary of PRTree objects
#'
#' @description
#' Prints a summary of a fitted Probabilistic Regression Tree model.
#'
#' @param x An object of class `summary.prtree`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the original object.
#'
#' @export
print.summary.prtree <- function(x, ...) {
  cat("\n========================================================\n")
  cat("Summary of Probabilistic Regression Tree (PRTree)\n")
  cat("========================================================\n")

  # Model dimensions
  cat(sprintf("Number of observations: %d\n", x$n_obs))
  cat(sprintf("Number of features: %d\n", x$n_feat))
  cat(sprintf("Number of terminal nodes: %d\n", x$n_terminal_nodes))
  cat(sprintf("Distribution: %s\n", x$distribution$name))
  if (!is.null(x$distribution$parameter)) {
    cat(sprintf(
      "Distribution parameter (%s): %.4f\n",
      x$distribution$parameter_name, x$distribution$parameter
    ))
  }
  cat(sprintf("Fill type: %d\n", x$fill_type))

  # mse
  cat("\nMean Squared Error:\n")
  cat(sprintf("  Train: %.6f\n", x$mse["train"]))

  if (x$has_validation) {
    cat(sprintf("  Validation: %.6f\n", x$mse["validation"]))
    cat(sprintf("  Global: %.6f\n", x$mse["global"]))
  } else if (x$has_test) {
    cat(sprintf("  Test: %.6f\n", x$mse["test"]))
    cat(sprintf("  Global: %.6f\n", x$mse["global"]))
  }

  # Sigma
  cat("\nSelected sigma:\n")
  if (length(x$sigma) == 1) {
    cat(sprintf("  %.4f\n", x$sigma))
  } else {
    cat("  ", .format_vector(x$sigma, max_print = 5, digits = 4), "\n")
  }

  # Terminal nodes information
  cat("\nTerminal node statistics:\n")
  print(x$terminal_nodes, digits = 4, row.names = FALSE)

  cat("========================================================\n")
  return(invisible(x))
}

#' @title Plot method for PRTree objects
#'
#' @description
#' Provides visualization of a fitted Probabilistic Regression Tree model.
#' Note that the original response vector `y` is only required for plots 1-4.
#' Plot 5 (terminal node probabilities) can be generated without `y`.
#'
#' @param x An object of class `prtree`, as returned by [pr_tree].
#'
#' @param y Original response vector used to fit the model. Required only for
#'   plots 1-4 (Fitted vs Observed, Residuals, Histogram, Q-Q plot).
#'   Can be omitted if only plot 5 is requested.
#'
#' @param which Which plots to produce. Options are:
#'   \itemize{
#'     \item `1`: \strong{Ovserved vs Fitted} - Scatter plot comparing
#'           observe responses against fitted values. Requires `y`.
#'     \item `2`: \strong{Residuals (index plot)} - Residuals plotted
#'           against observation index. Useful for detecting patterns or outliers.
#'           Requires `y`.
#'     \item `3`: \strong{Histogram of residuals} - Distribution of
#'           residuals with Freedman-Diaconis bin selection. Requires `y`.
#'     \item `4`: \strong{Q-Q plot of residuals} - Quantile-quantile plot
#'           to assess normality of residuals. Requires `y`.
#'     \item `5`: \strong{Terminal node probabilities} - Boxplot of
#'           \eqn{P(R|X)} for each terminal node. Does not require `y`.
#'   }
#'   If `NULL` (default), shows all plots (1-5).
#'
#' @param ncol Integer. Number of columns in the plot layout. If `NULL` (default),
#'   an appropriate number of columns is chosen automatically based on the number
#'   of plots.
#'
#' @param ... Further arguments passed to ggplot2 (e.g., `color`, `size`, `alpha`).
#'
#' @return Invisibly returns the original object.
#'
#' @importFrom rlang .data
#' @export
plot.prtree <- function(x, y = NULL, which = NULL, ncol = NULL, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for plotting.")
}

  # Determine which plots to show
  if (is.null(which)) which <- 1:5

  # Check if y is needed
  plots_need_y <- intersect(which, 1:4)
  if (length(plots_need_y) > 0) {
    if (is.null(y)) {
      .stop.message(
        paste0("Argument 'y' is required for plots ",
               paste(plots_need_y, collapse = ", "))
      )
    }
    if (length(y) != length(x$yhat)) {
      .stop.message("Length of y does not match number of observations")
    }
    residuals <- y - x$yhat
  }

  # Get sample types
  types <- names(x$yhat)
  if (is.null(types)) types <- rep("train", length(x$yhat))

  # Create data frame
  if (length(plots_need_y) > 0) {
    df_points <- data.frame(
      fitted = x$yhat,
      observed = y,
      residuals = residuals,
      index = seq_along(y),
      type = types
    )
  }

  # Create plot list
  plot_list <- list()
  plot_counter <- 1

  # Plot 1: Fitted vs Observed
  if (1 %in% which) {
    p <- ggplot2::ggplot(df_points, 
                         ggplot2::aes(x = .data$observed, 
                                      y = .data$fitted, 
                                      color = .data$type)) +
      ggplot2::geom_point(...) +
      ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      ggplot2::labs(
        x = "Observed values", y = "Fitted values", 
        title = "Observed vs Fitted", color = "Data set"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      ) +
      ggplot2::scale_color_manual(
        values = c(train = "black", test = "blue", validation = "darkgreen")
      )
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }

  # Plot 2: Residuals vs Index
  if (2 %in% which) {
    p <- ggplot2::ggplot(df_points, 
                         ggplot2::aes(x = .data$index, 
                                      y = .data$residuals, 
                                      color = .data$type)) +
      ggplot2::geom_point(...) +
      ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      ggplot2::labs(
        x = "Index", y = "Residuals",
        title = "Residuals", color = "Data set"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      ) +
      ggplot2::scale_color_manual(
        values = c(train = "black", test = "blue", validation = "darkgreen")
      )
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }

  # Plot 3: Histogram of residuals
  if (3 %in% which) {
    p <- ggplot2::ggplot(df_points, ggplot2::aes(x = .data$residuals)) +
      ggplot2::geom_histogram(
        bins = grDevices::nclass.FD(residuals),
        fill = "lightgray", color = "gray", ...
      ) +
      ggplot2::labs(x = "Residuals", y = "Count", title = "Histogram of Residuals") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1

  }

  # Plot 4: Q-Q plot of residuals
  if (4 %in% which) {
    p <- ggplot2::ggplot(df_points, ggplot2::aes(sample = .data$residuals)) +
      ggplot2::geom_qq(...) +
      ggplot2::geom_qq_line(color = "red", linetype = "dashed") +
      ggplot2::labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "Q-Q Plot") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }

  # Plot 5: Terminal node probabilities
  if (5 %in% which) {
    P_df <- as.data.frame(x$P)
    colnames(P_df) <- paste0("R", 1:ncol(P_df))
    P_long <- tidyr::pivot_longer(
      P_df,
      cols = tidyselect::everything(),
      names_to = "node",
      values_to = "probability"
    )
    P_long$node <- factor(P_long$node,
                          levels = paste0("R", 1:ncol(P_df)),
                          ordered = TRUE)

    p <- ggplot2::ggplot(P_long, ggplot2::aes(x = .data$node, y = .data$probability)) +
      ggplot2::geom_boxplot(fill = "lightblue", color = "blue", ...) +
      ggplot2::geom_hline(yintercept = c(0.25, 0.5, 0.75),
                          color = "gray", linetype = "dotted") +
      ggplot2::labs(
        x = "Terminal Node",
        y = expression(P(R ~ "|" ~ X)),
        title = "Terminal Node Probabilities"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }

  # Arrange plots
  n_plots <- length(plot_list)
  if (n_plots == 0) {
    message("No plots to display.")
    return(invisible(x))
  }

  if (is.null(ncol)) ncol <- min(n_plots, 3)
  nrow <- ceiling(n_plots / ncol)

  do.call(gridExtra::grid.arrange, c(plot_list, ncol = ncol, nrow = nrow))

  invisible(x)
}

#' @title Plot a PRTree as a dendrogram with gamma and P distributions
#'
#' @description
#' Draws a tree representation of a fitted PRTree model with root at top.
#' Terminal nodes show gamma coefficients and a mini heatmap of the
#' distribution of the corresponding column of the probability matrix P.
#'
#' @param x An object of class `prtree`.
#' @param heatmap_width Numeric. Width of heatmaps as fraction of node spacing.
#'   Default = 0.8 (80% of spacing).
#' @param heatmap_height Numeric. Height of heatmaps. Default = 0.08.
#' @param max_bins Integer. Maximum number of bins in the heatmap. Default = 6.
#' @param ... Further arguments passed to plot.
#'
#' @return Invisibly returns the original object.
#'
#' @export
#'
plot_tree <- function(x, heatmap_width = 0.8, heatmap_height = 0.08,
                      max_bins = 6, ...) {
  # === BASIC VALIDATION ===
  if (!inherits(x, "prtree")) stop("Object must be of class 'prtree'")
  if (!is.numeric(max_bins) || length(max_bins) != 1 || max_bins < 2 ||
      max_bins != round(max_bins)) {
    stop("'max_bins' must be a single integer greater than or equal to 2.")
  }
  n_bins <- max_bins

  required <- c("nodes_matrix_info", "features", "gamma", "P")
  missing <- required[!required %in% names(x)]
  if (length(missing) > 0) stop("Missing components: ", paste(missing, collapse = ", "))

  nodes <- x$nodes_matrix_info
  if (!all(c("node", "isTerminal", "fatherNode", "depth", "feature") %in% colnames(nodes))) {
    stop("nodes_matrix_info missing required columns")
  }

  # === IDENTIFY NODES ===
  terminal_idx <- which(nodes$isTerminal == 1)
  internal_idx <- which(nodes$isTerminal == 0)
  n_terminal <- length(terminal_idx)
  n_nodes <- nrow(nodes)
  max_depth <- max(nodes$depth)

  if (n_terminal == 0) stop("Tree has no terminal nodes")
  if (n_terminal != length(x$gamma)) {
    warning("Dimension mismatch: adjusting")
    n_terminal <- min(n_terminal, length(x$gamma))
    terminal_idx <- terminal_idx[1:n_terminal]
    x$gamma <- x$gamma[1:n_terminal]
    x$P <- x$P[, 1:n_terminal, drop = FALSE]
  }

  # === NODE POSITIONING ===
  # Y-coordinate: depth (children below parents)
  # Greater depth corresponds to smaller y values
  y_pos <- max_depth - nodes$depth

  # X-coordinate: exact placement of terminal nodes
  # Formula: (k - 0.5) / n_terminal
  x_pos_terminal <- (1:n_terminal - 0.5) / n_terminal

  # Initialize x positions
  x_pos <- numeric(n_nodes)

  # Recursive function for node positioning
  position_node <- function(node_idx, left_bound, right_bound) {
    if (node_idx %in% terminal_idx) {
      # Terminal node: position already assigned
      return(x_pos[node_idx])
    }

    # Find children
    children <- which(nodes$fatherNode == node_idx)

    if (length(children) == 2) {
      # Position left and right children
      x_left <- position_node(children[1], left_bound, x_pos[node_idx])
      x_right <- position_node(children[2], x_pos[node_idx], right_bound)

      # Current node is placed at the midpoint of its children
      x_pos[node_idx] <<- (x_left + x_right) / 2

    } else if (length(children) == 1) {
      # Only one child
      x_child <- position_node(children[1], left_bound, right_bound)
      x_pos[node_idx] <<- x_child

    } else {
      # No children (should not occur for internal nodes)
      x_pos[node_idx] <<- (left_bound + right_bound) / 2
    }

    return(x_pos[node_idx])
  }

  # First assign positions to terminal nodes
  # We first determine their left-to-right order in the tree
  # using an in-order traversal
  terminal_positions <- list()

  # Function to collect terminal nodes in left-to-right order
  collect_terminals <- function(node_idx) {
    if (node_idx %in% terminal_idx) {
      return(list(node_idx))
    }
    children <- which(nodes$fatherNode == node_idx)
    if (length(children) == 2) {
      left <- collect_terminals(children[1])
      right <- collect_terminals(children[2])
      return(c(left, right))
    } else if (length(children) == 1) {
      return(collect_terminals(children[1]))
    } else {
      return(list())
    }
  }

  # Collect terminal nodes in left-to-right order
  terminal_order <- unlist(collect_terminals(1))

  # Assign x positions according to the collected order
  for (k in 1:length(terminal_order)) {
    x_pos[terminal_order[k]] <- x_pos_terminal[k]
  }

  # Recursively position internal nodes starting from the root
  x_pos[1] <- position_node(1, 0, 1)

  # === COMPUTE NODE SPACING ===
  terminal_x_sorted <- sort(x_pos[terminal_idx])
  if (n_terminal > 1) {
    node_spacing <- min(diff(terminal_x_sorted))
  } else {
    node_spacing <- 0.5
  }

  # Heatmap width proportional to terminal-node spacing
  actual_heatmap_width <- node_spacing * heatmap_width

  # === DYNAMIC HISTOGRAM BINS ===
  breaks <- seq(0, 1, length.out = n_bins + 1)

  # === PLOT ===
  graphics::par(mfrow = c(1,1))
  graphics::plot.new()
  bottom_margin <- 0.4 + heatmap_height * 1.5
  graphics::plot.window(xlim = c(0, 1), ylim = c(-bottom_margin, max_depth + 0.5))

  # Draw edges
  for (i in 1:n_nodes) {
    if (nodes$fatherNode[i] > 0) {
      parent <- nodes$fatherNode[i]
      graphics::lines(c(x_pos[i], x_pos[parent]),
                      c(y_pos[i], y_pos[parent]),
                      col = "gray", lwd = 1)
    }
  }

  # Internal nodes
  node_cex <- 5
  for (i in internal_idx) {
    if (!is.na(x_pos[i])) {
      graphics::points(x_pos[i], y_pos[i], pch = 22, cex = node_cex,
                       col = "black", bg = "lightgray", lwd = 1)

      feature_name <- x$features[nodes$feature[i]]
      thr <- nodes$threshold[i]
      thr_txt <- if (abs(thr) > 1e4) sprintf("%.2e", thr) else sprintf("%.2f", thr)
      graphics::text(x_pos[i], y_pos[i], paste0(feature_name, "\n<", thr_txt),
                     cex = 0.5, font = 2)
    }
  }

  # Terminal nodes and heatmaps
  for (j in 1:n_terminal) {
    i <- terminal_idx[j]
    if (is.na(x_pos[i])) next

    # Terminal nodes
    graphics::points(x_pos[i], y_pos[i], pch = 21, cex = node_cex * 1.1,
                     col = "black", bg = "lightblue", lwd = 1)

    # Terminad-node weight
    gamma_txt <- if (abs(x$gamma[j]) > 1e4) sprintf("%.2e", x$gamma[j]) else sprintf("%.3f", x$gamma[j])
    graphics::text(x_pos[i], y_pos[i], gamma_txt, cex = 0.6, font = 2)

    # Heatmap
    p_vals <- x$P[, j]
    hist_counts <- graphics::hist(p_vals, breaks = breaks, plot = FALSE)$counts
    total <- sum(hist_counts)
    if (total > 0) hist_counts <- hist_counts / total

    # Vertical position below the terminal node
    heat_y <- y_pos[i] - 0.4 - heatmap_height/2

    # Center heatmap beneath the node
    x_left <- x_pos[i] - actual_heatmap_width/2

    # Draw histogram bins
    bin_w <- actual_heatmap_width / n_bins
    for (b in 1:n_bins) {
      x1 <- x_left + (b-1) * bin_w
      x2 <- x_left + b * bin_w
      gray_level <- 0.9 - hist_counts[b] * 0.7
      graphics::rect(x1, heat_y,
                     x2, heat_y + heatmap_height,
                     col = grDevices::rgb(gray_level, gray_level, gray_level),
                     border = NA)
    }

    # Draw heatmap border
    graphics::rect(x_left, heat_y,
                   x_left + actual_heatmap_width, heat_y + heatmap_height,
                   col = NA, border = "gray50", lwd = 0.5)
  }

  # Title and axes
  graphics::title(main = "PRTree Structure",
                  xlab = paste("Terminal nodes:", n_terminal),
                  ylab = "Depth")
  graphics::axis(2, at = max_depth:0, labels = 0:max_depth, las = 1)

  invisible(x)
}