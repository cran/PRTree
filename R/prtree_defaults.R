# Last revision: July/2026

#' @include prtree_rules.R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal. mapping names and code for distributions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# optionsion values for arguments #----
.pr.mapping <- list(
  dist =  list(
    norm  = list(dist_code = 1L, par_name = NULL),
    lnorm = list(dist_code = 2L, par_name = "sdlog"),
    t     = list(dist_code = 3L, par_name = "df"),
    gamma = list(dist_code = 4L, par_name = "shape")
  ),
  proxy_crit =  list(mean  = 1L, var = 2L, both = 3L)
)
attr(.pr.mapping, "dist_names") <- names(.pr.mapping$dist)
attr(.pr.mapping, "dist_par_names") <- unname(
  unlist(
    sapply(
      names(.pr.mapping$dist),
      function(n) .pr.mapping$dist[[n]]$par_name
    )
  )
)
attr(.pr.mapping, "proxy_crit_names") <- names(.pr.mapping$proxy_crit)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Parameter metadata (list format for easy maintenance) #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.pr.meta_list <- c(
  # distribution parameters: df, shape, sdlog
  sapply(attr(.pr.mapping, "dist_par_names"),
         function(n) {
           list(
             type = "numeric",
             rule = ".rule.param",
             dots_remove = TRUE,
             depends = "dist"
           )
         },
         simplify = FALSE),
  # current
  list(
    # debuggin option (debuggin purposes during code construction)
    debug = list(type = "logical", rule = ".rule.debug", default = FALSE),

    # Printing
    iprint = list(
      type = "integer", rule = ".rule.iprint",
      default = -1L, control_pr = TRUE, control_cv = TRUE
    ),
    verbose = list(
      type = "logical", rule = ".rule.verbose",
      default = TRUE, control_pr = TRUE, control_cv = TRUE
    ),

    # deprecated
    perc_x = list(type = "numeric", rule = ".rule.perc_x",
                  dots_remove = TRUE, deprecated = TRUE),
    perc_test = list(type = "numeric", rule = ".rule.perc_test",
                     dots_remove = TRUE, deprecated = TRUE),

    # Data
    y = list(type = "numeric", rule = ".rule.y"),
    X = list(type = "numeric", rule = ".rule.X", depends = "y"),
    n_obs = list(type = "integer", rule = ".rule.n_obs", depends = "y"),
    n_feat = list(type = "integer", rule = ".rule.n_feat", depends = "X"),

    # Grid parameters
    sigma_grid = list(
      type = "numeric", rule = ".rule.sigma_grid",
      control_pr = TRUE, control_cv = TRUE,
      depends = "n_feat"
    ),
    tiny_sigma = list(
      type = "numeric", rule = ".rule.tiny_sigma",
      control_pr = TRUE, control_cv = TRUE,
      depends = "sigma_grid"
    ),
    grid_size = list(
      type = "integer", rule = ".rule.grid_size",
      default = 8L, control_pr = TRUE, control_cv = TRUE,
      depends = "tiny_sigma"
    ),
    min_mult = list(
      type = "numeric", rule = ".rule.min_mult",
      default = 0.0, control_pr = TRUE, control_cv = TRUE,
      depends = "sigma_grid"
    ),
    max_mult = list(
      type = "numeric", rule = ".rule.max_mult",
      default = 2.0, control_pr = TRUE, control_cv = TRUE,
      depends = "min_mult"
    ),

    # Stopping criteria
    max_terminal_nodes = list(
      type = "integer", rule = ".rule.max_terminal_nodes",
      default = 15L, stopping = TRUE, control_pr = TRUE, control_cv = TRUE
    ),
    max_depth = list(
      type = "integer", rule = ".rule.max_depth",
      stopping = TRUE, control_pr = TRUE, control_cv = TRUE,
      depends = "max_terminal_nodes"
    ),
    cp = list(
      type = "numeric", rule = ".rule.cp",
      default = 0.01, stopping = TRUE, control_pr = TRUE, control_cv = TRUE
    ),
    n_min = list(
      type = "integer", rule = ".rule.n_min",
      default = 5L, stopping = TRUE, control_pr = TRUE, control_cv = TRUE
    ),
    prop_x = list(
      type = "numeric", rule = ".rule.prop_x",
      default = 0.1, stopping = TRUE, control_pr = TRUE, control_cv = TRUE,
      depends = "perc_x"
    ),
    p_min = list(
      type = "numeric", rule = ".rule.p_min",
      default = 0.05, stopping = TRUE, control_pr = TRUE, control_cv = TRUE
    ),

    # Missing data
    fill_type = list(
      type = "integer", rule = ".rule.fill_type",
      options = c(0L, 1L, 2L), default = 2L, control_pr = TRUE, control_cv = TRUE
    ),
    proxy_crit = list(
      type = "character", rule = ".rule.proxy_crit",
      options = attr(.pr.mapping, "proxy_crit_names"), default = "both",
      control_pr = TRUE, control_cv = TRUE
    ),

    # Split search
    n_candidates = list(
      type = "integer", rule = ".rule.n_candidates",
      default = 3L, control_pr = TRUE, control_cv = TRUE
    ),
    by_node = list(
      type = "logical", rule = ".rule.by_node",
      default = FALSE, control_pr = TRUE, control_cv = TRUE
    ),

    # Distribution
    dist = list(
      type = "character", rule = ".rule.dist",
      options = attr(.pr.mapping, "dist_names"), default = "norm",
      control_pr = TRUE, control_cv = TRUE
    ),
    dist_pars = list(
      type = "list", rule = ".rule.dist_pars",
      control_pr = TRUE, control_cv = TRUE,
      depends = attr(.pr.mapping, "dist_par_names")
    ),

    # Cross-validation
    method = list(
      type = "character", rule = ".rule.method",
      options = c("montecarlo", "kfold"), default = "montecarlo", control_cv = TRUE
    ),
    n_rep = list(
      type = "integer", rule = ".rule.n_rep",
      default = 10L, control_cv = TRUE, depends = "method"
    ),
    fold_idx = list(
      type = "integer", rule = ".rule.fold_idx",
      control_cv = TRUE, depends = c("method", "n_rep", "n_obs")
    ),
    only_sigma = list(
      type = "logical", rule = ".rule.only_sigma",
      default = FALSE, control_cv = TRUE, depends = "grid_size"
    ),
    stratify = list(
      type = "logical", rule = ".rule.stratify",
      default = FALSE, control_cv = TRUE
    ),
    is_NA = list(type = "logical", rule = ".rule.is_NA", depends = c("stratify", "n_obs")),
    update_final = list(
      type = "logical", rule = ".rule.update_final",
      default = TRUE, control_cv = TRUE
    ),

    # Data splitting
    prop_hold = list(
      type = "numeric", rule = ".rule.prop_hold",
      default = 0.2, control_pr = TRUE, depends = c("perc_test")
    ),
    prop_test = list(
      type = "numeric", rule = ".rule.prop_test",
      default = 0.2, control_cv = TRUE,
      depends = c("only_sigma", "n_rep")
    ),
    prop_valid = list(
      type = "numeric", rule = ".rule.prop_valid",
      default = 0.2, control_cv = TRUE, depends = "only_sigma"
    ),
    idx_train = list(
      type = "integer", rule = ".rule.idx_train", control_pr = TRUE,
      depends = "n_obs"
    )
  ))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Calculate validation order based on dependencies #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.sort_by_depends <- function(meta_list) {
  # Get all parameter names
  params <- names(meta_list)

  # Build dependency graph
  deps <- list()
  for (p in params) {
    if (!is.null(meta_list[[p]]$depends)) {
      deps[[p]] <- meta_list[[p]]$depends
    }
  }

  # Topological sort
  result <- character()
  remaining <- params

  while (length(remaining) > 0) {
    available <- remaining[sapply(remaining, function(p) {
      if (is.null(deps[[p]])) {
        return(TRUE)
      }
      all(deps[[p]] %in% result)
    })]

    if (length(available) == 0) {
      warning("Circular dependency detected. Using original order.")
      return(params)
    }

    result <- c(result, available[1])
    remaining <- setdiff(remaining, available[1])
  }

  result
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: helper functions #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.name <- function(x) rownames(.pr.meta)[x]
.rule <- function(name) attr(.pr.meta, "rules")[[name]]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Convert parameter list to data.frame for easy filtering #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.pr.meta <- do.call(rbind, lapply(names(.pr.meta_list), function(n) {
  data.frame(
    type = .null.default(.pr.meta_list[[n]]$type, "none"),
    options = I(list(.null.default(.pr.meta_list[[n]]$options, NULL))),
    default = I(list(.null.default(.pr.meta_list[[n]]$default, NULL))),
    control_pr = .null.default(.pr.meta_list[[n]]$control_pr, FALSE),
    control_cv = .null.default(.pr.meta_list[[n]]$control_cv, FALSE),
    stopping = .null.default(.pr.meta_list[[n]]$stopping, FALSE),
    dots_remove = .null.default(.pr.meta_list[[n]]$dots_remove, FALSE),
    stringsAsFactors = FALSE,
    row.names = n
  )
}))
attr(.pr.meta, "check_order") <- .sort_by_depends(.pr.meta_list)
attr(.pr.meta, "stopping") <- .get.name(.pr.meta$stopping)
attr(.pr.meta, "dots_remove") <- .get.name(.pr.meta$dots_remove)
attr(.pr.meta, "control_pr") <- .get.name(.pr.meta$control_pr)
attr(.pr.meta, "control_cv") <- .get.name(.pr.meta$control_cv)
attr(.pr.meta, "pr_int") <- .get.name(.pr.meta$type == "integer" & .pr.meta$control_pr)
attr(.pr.meta, "pr_dble") <- .get.name(.pr.meta$type == "numeric" & .pr.meta$control_pr)
attr(.pr.meta, "rules") <- sapply(rownames(.pr.meta),
                                  function(name) {get(.pr.meta_list[[name]]$rule)},
                                  simplify = FALSE)

