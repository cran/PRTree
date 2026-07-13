# Last revision: July/2026

#' @title Print method for PRTree Cross-Validation objects
#'
#' @description
#' Prints a brief summary of cross-validation results for Probabilistic
#' Regression Trees.
#'
#' @param x An object of class `prtree.cv`, as returned by [pr_tree_cv].
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the original object.
#'
#' @export
#'
print.prtree.cv <- function(x, ...) {
  method <- x$config$method
  config <- x$config
  n_rows <- nrow(x$sigma_matrix)

  cat("\n========================================================\n")
  cat("PRTree Cross-Validation Results\n")
  cat("========================================================\n")

  # --- Method Summary ---
  if (method == "kfold") {
    cat(sprintf("Method: %d-fold CV\n", config$n_rep))
  } else {
    cat(sprintf("Method: Monte Carlo (%d iterations)\n", config$n_rep))
  }

  # --- Test set information ---
  if (config$only_sigma) {
    cat("Test set: none (sigma selection only).\n")
  } else if (method == "kfold") {
    cat(sprintf("Test set: 1 fold per iteration (~%.1f%% of data).\n", 100 / config$n_rep))
  } else if (!is.null(config$prop_test) && config$prop_test > 0) {
    cat(sprintf("Test proportion: %.2f\n", config$prop_test))
  } else {
    cat("Test set: none.\n")
  }

  # --- Validation info ---
  has_validation <- !is.null(x$rmse_by_rep$validation)

  if (method == "kfold" && config$only_sigma) {
    cat(sprintf("Validation (for sigma): 1 fold per iteration (~%.1f%% of data).\n", 100 / config$n_rep))
  } else if (!is.null(config$prop_valid) && config$prop_valid > 0) {
    cat(sprintf("Validation (for sigma): %.2f of the training data/folds.\n", config$prop_valid))
  } else if (method == "kfold" && !has_validation) {
    cat("Validation (for sigma): none (sigma selected on training folds).\n")
  } else if (method == "montecarlo" && has_validation) {
    cat("  Validation: performed (internal split)\n")
  } else if (method == "montecarlo" && !has_validation && nrow(x$sigma_matrix) > 1) {
    cat("  Validation: none (sigma selected on training data)\n")
  } else {
    cat("  Validation: none\n")
  }

  # --- General Info ---
  cat(sprintf("Features: %d\n", ncol(x$sigma_matrix)))

  # --- Sigma Selection Info ---
  if (has_validation) {
    cat("\n--- Sigma Selection ---\n")
    avg_val_rmse <- mean(x$rmse_by_rep$validation)
    sd_val_rmse <- stats::sd(x$rmse_by_rep$validation)
    cat(sprintf("Validation RMSE (mean +/- sd): %.6f +/- %.6f\n", avg_val_rmse, sd_val_rmse))
  } else {
    if (nrow(x$sigma_matrix) == 1) {
      cat("\n--- Sigma fixed ---\n")
      cat("Single sigma candidate used.\n")
    } else {
      cat("\n--- Sigma Selection ---\n")
      cat("No validation set used (sigma selected on training data).\n")
    }
  }

  # --- Test Error Info ---
  if (!is.null(x$rmse_by_rep$test)) {
    cat("\n--- Test Error ---\n")
    avg_test_rmse <- mean(x$rmse_by_rep$test)
    sd_test_rmse <- stats::sd(x$rmse_by_rep$test)
    cat(sprintf("Test RMSE (mean +/- sd): %.6f +/- %.6f\n", avg_test_rmse, sd_test_rmse))
  } else {
    cat("\n--- Test Error ---\n")
    cat("No test set used.\n")
  }

  # --- First 5 results ---
  cat("\nFirst 5 iterations/folds (sigma values and RMSE):\n")
  n_show <- min(5, n_rows)
  if (n_show > 0) {
    df_show <- data.frame(
      iter = 1:n_show,
      x$sigma_matrix[1:n_show, , drop = FALSE]
    )

    if (!is.null(x$rmse_by_rep$validation)) {
      df_show$RMSE_val <- x$rmse_by_rep$validation[1:n_show]
    }

    if (!is.null(x$rmse_by_rep$test)) {
      df_show$RMSE_test <- x$rmse_by_rep$test[1:n_show]
    }

    print(df_show, row.names = FALSE, ...)
  }

  cat("========================================================\n")
  return(invisible(x))
}


#' @title Summary method for PRTree Cross-Validation objects
#'
#' @description
#' Provides a detailed summary of cross-validation results for Probabilistic
#' Regression Trees, including statistics for RMSE and selected sigma values.
#'
#' @param object An object of class `prtree.cv`, as returned by `[pr_tree_cv]`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list of class `summary.prtree.cv` containing:
#'   \item{method}{Cross-validation method used.}
#'   \item{config}{Complete configuration list.}
#'   \item{n_features}{Number of features.}
#'   \item{sigma_summary}{Matrix with summary statistics for selected sigma values
#'     (Min, Mean, Median, Max, SD) for each feature.}
#'   \item{rmse_summary}{List with summary statistics for validation and/or test
#'     RMSE (Min, Q1, Median, Mean, Q3, Max, SD).}
#'
#' @export
summary.prtree.cv <- function(object, ...) {
  method <- object$config$method
  config <- object$config
  n_features <- ncol(object$sigma_matrix)

  # --- Create summary list ---
  summ <- list(
    method = method,
    config = config,
    n_features = n_features,
    sigma_summary = NULL,
    sigma_fixed = FALSE,
    sigma_value = NULL,
    rmse_summary = list()
  )

  # --- Summarize RMSE Validation ---
  if (!is.null(object$rmse_by_rep$validation)) {
    rmse_val <- object$rmse_by_rep$validation
    summ$rmse_summary$validation <- data.frame(
      Min = min(rmse_val),
      Q1 = unname(stats::quantile(rmse_val, 0.25)),
      Median = stats::median(rmse_val),
      Mean = mean(rmse_val),
      Q3 = unname(stats::quantile(rmse_val, 0.75)),
      Max = max(rmse_val),
      SD = stats::sd(rmse_val),
      N = length(rmse_val),
      check.names = FALSE
    )
  }

  # --- Summarize RMSE Test ---
  if (!is.null(object$rmse_by_rep$test)) {
    rmse_test <- object$rmse_by_rep$test
    summ$rmse_summary$test <- data.frame(
      Min = min(rmse_test),
      Q1 = unname(stats::quantile(rmse_test, 0.25)),
      Median = stats::median(rmse_test),
      Mean = mean(rmse_test),
      Q3 = unname(stats::quantile(rmse_test, 0.75)),
      Max = max(rmse_test),
      SD = stats::sd(rmse_test),
      N = length(rmse_test),
      check.names = FALSE
    )
  }

  # --- Summarize Sigma ---
  sigma_mat <- object$sigma_matrix
  if (!is.null(sigma_mat) && nrow(sigma_mat) > 0) {
    # Verificar se sigma foi selecionado (múltiplas linhas) ou é fixo
    # Para ser selecionado, precisa ter mais de uma linha OU track_sigma = TRUE
    sigma_selected <- nrow(sigma_mat) > 1

    if (sigma_selected) {
      # Sigma foi selecionado - mostrar estatísticas
      sigma_summary_list <- list()
      for (j in 1:n_features) {
        feature_name <- if (!is.null(colnames(sigma_mat))) colnames(sigma_mat)[j] else paste0("V", j)
        sigma_j <- sigma_mat[, j]
        sigma_summary_list[[feature_name]] <- c(
          Min = min(sigma_j),
          Mean = mean(sigma_j),
          Median = stats::median(sigma_j),
          Max = max(sigma_j),
          SD = stats::sd(sigma_j)
        )
      }
      summ$sigma_summary <- do.call(rbind, sigma_summary_list)
      summ$sigma_fixed <- FALSE
    } else {
      # Sigma é fixo - guardar o valor
      summ$sigma_fixed <- TRUE
      sigma_value <- sigma_mat[1, ]
      if (!is.null(colnames(sigma_mat))) {
        names(sigma_value) <- colnames(sigma_mat)
      }  else {
        names(sigma_value) <- paste0("V", 1:n_features)
      }
      summ$sigma_value <- sigma_value
    }
  }

  class(summ) <- "summary.prtree.cv"
  return(summ)
}


#' @title Print method for summary of PRTree Cross-Validation objects
#'
#' @description
#' Prints a detailed summary of cross-validation results.
#'
#' @param x An object of class `summary.prtree.cv`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the original object.
#'
#' @export
print.summary.prtree.cv <- function(x, ...) {
  cat("\n========================================================\n")
  cat("Summary of PRTree Cross-Validation\n")
  cat("========================================================\n")

  # --- Configuration ---
  cat("Cross-Validation Configuration:\n")
  if (x$config$method == "kfold") {
    cat(sprintf("  Method: %d-fold CV\n", x$config$n_rep))
    if (x$config$only_sigma) {
      cat("  Test set: none (sigma selection only).\n")
      cat(sprintf("  Validation (for sigma): 1 fold per iteration (~%.1f%% of data).\n", 100 / x$config$n_rep))
    } else {
      cat(sprintf("  Test set: 1 fold per iteration (~%.1f%% of data).\n", 100 / x$config$n_rep))
      if (!is.null(x$config$prop_valid) && x$config$prop_valid > 0) {
        cat(sprintf("  Validation (for sigma): %.2f of the training folds.\n", x$config$prop_valid))
      } else {
        cat("  Validation: none (sigma selected on training fold)\n")
      }
    }
  } else {
    cat(sprintf("  Method: Monte Carlo (%d iterations)\n", x$config$n_rep))
    if (x$config$only_sigma) {
      cat("  Test set: none (sigma selection only)\n")
    } else if (!is.null(x$config$prop_test) && x$config$prop_test > 0) {
      cat(sprintf("  Test proportion: %.2f\n", x$config$prop_test))
    } else {
      cat("  Test set: none\n")
    }

    if (!is.null(x$config$prop_valid) && x$config$prop_valid > 0) {
      cat(sprintf("  Validation proportion: %.2f (internal split)\n", x$config$prop_valid))
    } else if (!is.null(x$rmse_summary$validation)) {
      cat("  Validation: performed (internal split)\n")
    } else {
      cat("  Validation: none\n")
    }
  }
  cat(sprintf("  Features: %d\n", x$n_features))
  cat(sprintf("  Stratification: %s\n", toupper(as.character(x$config$stratify))))

  # --- RMSE Summary ---
  cat("\n--- RMSE Summary ---\n")
  if (!is.null(x$rmse_summary$validation)) {
    cat("Validation Set (used for sigma selection):\n")
    print(x$rmse_summary$validation, row.names = FALSE,...)
  }
  if (!is.null(x$rmse_summary$test)) {
    cat("\nTest Set (unseen data):\n")
    print(x$rmse_summary$test, row.names = FALSE,...)
  }
  if (is.null(x$rmse_summary$validation) && is.null(x$rmse_summary$test)) {
    cat("  No RMSE data available.\n")
  }

  # --- Sigma Info ---
  if (!is.null(x$sigma_fixed) && x$sigma_fixed) {
    cat("\n--- Sigma fixed ---\n")
    cat("Single sigma candidate used:\n")
    print(x$sigma_value, ...)
  } else if (!is.null(x$sigma_summary)) {
    cat("\n--- Selected Sigma Summary (by feature) ---\n")
    print(x$sigma_summary, ...)
  } else {
    cat("\n--- No sigma information available ---\n")
  }

  cat("========================================================\n")
  return(invisible(x))
}

#' @title Plot method for PRTree Cross-Validation objects
#'
#' @description
#' Provides visualizations for cross-validation results of Probabilistic
#' Regression Trees.
#'
#' @param x An object of class `prtree.cv`, as returned by [pr_tree_cv].
#'
#' @param which Which plots to produce. Options are:
#'   \itemize{
#'     \item `1`: \strong{RMSE across iterations/folds} - Shows validation and test
#'           RMSE for each iteration/fold.
#'
#'     \item `2`: \strong{Selected sigma values} - Values of selected sigma
#'           across iterations/folds. Colors distinguish features.
#'
#'     \item `3`: \strong{Distribution of sigma} - For single feature: histogram;
#'           for multiple features: boxplots colored by feature.
#'
#'     \item `4`: \strong{Sigma vs Validation RMSE} - Scatter plot with LOESS smooth
#'           to show relationship between selected sigma and validation RMSE.
#'   }
#'   If `NULL` (default), shows all applicable plots.
#'
#' @param ncol Integer. Number of columns in the plot layout. If `NULL` (default),
#'   an appropriate number of columns is chosen automatically.
#'
#' @param ... Further arguments passed to ggplot2 (e.g., `size`, `alpha`).
#'
#' @return Invisibly returns the original object.
#'
#' @importFrom rlang .data
#' @export
#'
plot.prtree.cv <- function(x, which = NULL, ncol = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  # Determine available data
  has_validation <- !is.null(x$rmse_by_rep$validation)
  has_test <- !is.null(x$rmse_by_rep$test)
  n_iter <- nrow(x$sigma_matrix)
  n_features <- ncol(x$sigma_matrix)

  # Feature names
  feature_names <- colnames(x$sigma_matrix)
  if (is.null(feature_names)) {
    feature_names <- paste0("Feature_", 1:n_features)
  }

  # Color palette for features
  feature_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                      "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  if (n_features > length(feature_colors)) {
    feature_colors <- grDevices::rainbow(n_features)
  }
  names(feature_colors) <- feature_names[1:n_features]

  # Determine which plots are available
  available_plots <- 1:2
  if (n_iter > 1 && n_features > 0) available_plots <- c(available_plots, 3)
  if (has_validation) available_plots <- c(available_plots, 4)

  # Determine which plots to show
  if (is.null(which)) {
    which <- available_plots
  } else {
    invalid <- which[!which %in% 1:4]
    if (length(invalid) > 0) {
      stop(paste("Invalid plot number(s):", paste(invalid, collapse = ", "),
                 "\nValid options are 1, 2, 3, or 4."))
    }
    unavailable <- which[!which %in% available_plots]
    if (length(unavailable) > 0) {
      warning(paste("Plot(s)", paste(unavailable, collapse = ", "),
                    "not available. Skipping."))
      which <- which[which %in% available_plots]
    }
  }

  n_plots <- length(which)
  if (n_plots == 0) {
    message("No plots to display.")
    return(invisible(x))
  }

  # Prepare data for plots
  xlab_method <- ifelse(x$config$method == "kfold", "Fold", "Iteration")
  iter_names <- 1:n_iter

  # Data for plot 1 (RMSE)
  if (1 %in% which && (has_validation || has_test)) {
    df_rmse <- data.frame(iteration = iter_names)
    if (has_validation) {
      df_rmse$validation <- x$rmse_by_rep$validation
    }
    if (has_test) {
      df_rmse$test <- x$rmse_by_rep$test
    }
    df_rmse_long <- tidyr::pivot_longer(
      df_rmse,
      cols = -.data$iteration,
      names_to = "type",
      values_to = "rmse"
    )
  }

  # Data for plot 2 (Sigma values)
  if (2 %in% which && n_features >= 1) {
    df_sigma <- as.data.frame(x$sigma_matrix)
    colnames(df_sigma) <- feature_names
    df_sigma$iteration <- iter_names
    df_sigma_long <- tidyr::pivot_longer(
      df_sigma,
      cols = -.data$iteration,
      names_to = "feature",
      values_to = "sigma"
    )
  }

  # Data for plot 3 (Distribution)
  if (3 %in% which && n_iter > 1) {
    df_dist <- as.data.frame(x$sigma_matrix)
    colnames(df_dist) <- feature_names
  }

  # Data for plot 4 (Sigma vs RMSE)
  if (4 %in% which && has_validation) {
    df_relation <- data.frame()
    for (j in 1:n_features) {
      df_temp <- data.frame(
        rmse = x$rmse_by_rep$validation,
        sigma = x$sigma_matrix[, j],
        feature = feature_names[j]
      )
      df_relation <- rbind(df_relation, df_temp)
    }
  }

  # Create plot list
  plot_list <- list()
  plot_counter <- 1

  # Plot 1: RMSE across iterations/folds
  if (1 %in% which && (has_validation || has_test)) {
    p <- ggplot2::ggplot(df_rmse_long, ggplot2::aes(x = .data$iteration, 
                                                    y = .data$rmse, 
                                                    color = .data$type)) +
      ggplot2::geom_point(size = 2, ...) +
      ggplot2::geom_hline(
        data = stats::aggregate(rmse ~ type, df_rmse_long, mean),
        ggplot2::aes(yintercept = .data$rmse, 
                     color = .data$type),
        linetype = "dashed", linewidth = 0.5
      ) +
      ggplot2::labs(
        x = xlab_method,
        y = "RMSE",
        title = "RMSE by Iteration/Fold",
        color = "Type"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      ) +
      ggplot2::scale_color_manual(
        values = c(validation = "blue", test = "red")
      )
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }

  # Plot 2: Selected sigma values
  if (2 %in% which && n_features >= 1) {
    if (n_features == 1) {
      p <- ggplot2::ggplot(df_sigma_long, ggplot2::aes(x = .data$iteration, 
                                                       y = .data$sigma)) +
        ggplot2::geom_point(size = 2, color = "blue", ...) +
        ggplot2::geom_hline(
          yintercept = mean(df_sigma_long$sigma),
          color = "red", linetype = "dashed", linewidth = 0.5
        ) +
        ggplot2::labs(
          x = xlab_method,
          y = expression(sigma),
          title = "Selected Sigma Values"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    } else {
      p <- ggplot2::ggplot(df_sigma_long, ggplot2::aes(x = .data$iteration, 
                                                       y = .data$sigma, 
                                                       color = .data$feature)) +
        ggplot2::geom_point(size = 2, ...) +
        ggplot2::geom_hline(
          data = stats::aggregate(sigma ~ feature, df_sigma_long, mean),
          ggplot2::aes(yintercept = .data$sigma, color = .data$feature),
          linetype = "dashed", linewidth = 0.5, alpha = 0.7
        ) +
        ggplot2::labs(
          x = xlab_method,
          y = expression(sigma),
          title = "Selected Sigma Values by Feature",
          color = "Feature"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom"
        ) +
        ggplot2::scale_color_manual(values = feature_colors)
    }
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }

  # Plot 3: Distribution of selected sigma values
  if (3 %in% which && n_iter > 1) {
    if (n_features == 1) {
      df_hist <- df_dist
      colnames(df_hist)[1] <- "sigma_value"

      p <- ggplot2::ggplot(df_hist, ggplot2::aes(x = .data$sigma_value)) +
        ggplot2::geom_histogram(
          bins = grDevices::nclass.FD(df_hist$sigma_value),
          fill = "lightgray", color = "gray",
          ...
        ) +
        ggplot2::geom_vline(
          xintercept = mean(df_hist$sigma_value),
          color = "red", linetype = "dashed", linewidth = 0.5
        ) +
        ggplot2::labs(
          x = expression(sigma),
          y = "Count",
          title = "Distribution of Selected Sigma"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    } else {
      df_box <- tidyr::pivot_longer(
        df_dist,
        cols = tidyselect::everything(),
        names_to = "feature",
        values_to = "sigma"
      )
      p <- ggplot2::ggplot(df_box, ggplot2::aes(x = .data$feature, 
                                                y = .data$sigma, 
                                                fill = .data$feature)) +
        ggplot2::geom_boxplot(alpha = 0.7, linewidth = 0.3, ...) +
        ggplot2::geom_hline(yintercept = 0, color = "gray", linetype = "dotted", linewidth = 0.3) +
        ggplot2::labs(
          x = "Feature",
          y = expression(sigma),
          title = "Distribution of Selected Sigma",
          fill = "Feature"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "none"
        ) +
        ggplot2::scale_fill_manual(values = feature_colors)
    }
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }

  # Plot 4: Sigma vs Validation RMSE - Com LOESS ajustado para poucos pontos
  if (4 %in% which && has_validation) {
    n_points_per_feature <- nrow(df_relation) / n_features

    p <- ggplot2::ggplot(df_relation,
                         ggplot2::aes(x = .data$rmse, 
                                      y = .data$sigma, 
                                      color = .data$feature)) +
      ggplot2::geom_point(size = 2.5, alpha = 0.8, ...) +
      ggplot2::geom_smooth(
        method = "loess",
        formula = y ~ x,
        se = TRUE,
        linewidth = 0.8,
        alpha = 0.15,
        span = min(1.5, 2.0 * n_points_per_feature^(-0.3)),  # Ajuste dinâmico do span
        method.args = list(degree = 1, family = "symmetric")
      ) +
      ggplot2::geom_hline(
        data = stats::aggregate(sigma ~ feature, df_relation, mean),
        ggplot2::aes(yintercept = .data$sigma, color = .data$feature),
        linetype = "dotted", linewidth = 0.4, alpha = 0.5
      ) +
      ggplot2::labs(
        x = "Validation RMSE",
        y = expression(sigma),
        title = "Sigma vs Validation RMSE",
        subtitle = paste0("LOESS smooth (", n_points_per_feature, " points per feature)"),
        color = "Feature"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, face = "italic"),
        legend.position = "bottom"
      ) +
      ggplot2::scale_color_manual(values = feature_colors)

    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }

  # Arrange plots
  n_plots <- length(plot_list)

  if (n_plots == 0) {
    message("No plots to display.")
    return(invisible(x))
  }

  if (is.null(ncol)) {
    ncol <- min(n_plots, 2)
  }
  nrow <- ceiling(n_plots / ncol)

  # Arrange all plots in a single grid
  do.call(gridExtra::grid.arrange, c(plot_list, ncol = ncol, nrow = nrow))

  invisible(x)
}