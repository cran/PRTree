# Last revision: July/2026

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal. prints a given message and stops #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.stop.message <- function(msg) {
  stop(paste0(
    "\n-----------------------------------------------------------\n",
    msg,
    "\n-----------------------------------------------------------\n"
  ), call. = FALSE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal. prints a given message (warning) #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.warning.message <- function(msg) {
  warning(paste0(
    "\n-----------------------------------------------------------\n",
    msg,
    "\n-----------------------------------------------------------\n"
  ), call. = FALSE, immediate. = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Format vector for printing with truncation #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.format_vector <- function(x, max_print = 5, digits = 5, collapse = ", ") {
  if (is.null(x) || length(x) == 0) {
    return("NULL")
  }

  if (length(x) == 1) {
    return(sprintf(paste0("%.", digits, "f"), x))
  }

  formatted <- sprintf(paste0("%.", digits, "f"), x)

  if (length(x) <= max_print) {
    return(paste(formatted, collapse = collapse))
  } else {
    first <- paste(formatted[1:max_print], collapse = collapse)
    return(paste0(first, collapse, "... (", length(x), " total)"))
  }
}
