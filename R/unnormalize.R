#' execute un-normalize fit results
#'
#' @noRd
execute_unnormalize <- function(data.normalized,
                                tangency.point,
                                shift,
                                fit) {
  # shorthands for calculations
  y.tangent <- tangency.point$y.tangent
  x.tangent <- tangency.point$x.tangent
  x.shift <- shift$x.shift
  y.shift <- shift$y.shift
  y.min <- data.normalized[data.normalized$otr.idx == fit$otr.idx.start,
                           "y.normalized"]
  y.max <- data.normalized[data.normalized$otr.idx == fit$otr.idx.end,
                           "y.normalized"]
  x.min <- data.normalized[data.normalized$otr.idx == fit$otr.idx.start,
                           "x.normalized"]
  x.max <- data.normalized[data.normalized$otr.idx == fit$otr.idx.end,
                           "x.normalized"]
  m <- fit$m
  b <- fit$b

  # Un-normalize fit
  result <- list("finalSlope" = m * y.tangent / x.tangent,
       "finalIntercept" = b * y.tangent,
       "trueIntercept" = b * y.tangent - x.shift * m * y.tangent / x.tangent + y.shift,
       "y.lowerBound" = y.min * y.tangent + y.shift,
       "y.upperBound" = y.max * y.tangent + y.shift,
       "x.lowerBound" = x.min * x.tangent + x.shift,
       "x.upperBound" = x.max * x.tangent + x.shift)

  # return result
  return(result)
}

#' execute un-normalize fit results
#'
#' @noRd
execute_unnormalize.data <- function(data.normalized,
                                     tangency.point,
                                     shift) {
  # shorthands for calculations
  y.tangent <- tangency.point$y.tangent
  x.tangent <- tangency.point$x.tangent
  x.shift <- shift$x.shift
  y.shift <- shift$y.shift

  # Un-normalize data
  data.normalized %>%
    dplyr::mutate(
      x.unnormalized = .data$x.normalized * x.tangent + x.shift) %>%
    dplyr::mutate(
      y.unnormalized = .data$y.normalized * y.tangent + y.shift) %>%
    dplyr::select(-c("x.normalized", "y.normalized")) %>%
    return()
}

#' un-normalize fit results
#'
#' @note
#' Refer to chapter 6.9 in ASTM E3076-18
#'
#' @description
#' un-normalize slope of optimum fit to engineering units and assemble report
#'
#' @noRd
unnormalize_results <- function(normalized_data,
                                otr.info = NULL,
                                dataQualityMetrics,
                                fit,
                                fitQualityMetrics,
                                verbose = FALSE,
                                warn = TRUE) {
  # get units if available
  if (!is.null(otr.info)) {
    x.unit <- otr.info$x$unit
    y.unit <- otr.info$y$unit
    if (!is.null(x.unit) && !is.null(y.unit)) {
      slope.unit <- paste0(y.unit, "/", x.unit)
    } else if (is.null(x.unit) && !is.null(y.unit)) {
      slope.unit <- y.unit
    } else {
      slope.unit <- ""
    }
    x.unit <- ifelse(!is.null(x.unit), x.unit, "")
    y.unit <- ifelse(!is.null(y.unit), y.unit, "")
  } else {
    x.unit <- ""
    y.unit <- ""
    slope.unit <- ""
  }

  # Un-normalize fit
  fit.unnormalized <- execute_unnormalize(
    data.normalized = normalized_data$data.normalized,
    tangency.point = normalized_data$tangency.point,
    shift = normalized_data$shift,
    fit = fit)

  finalSlope <- fit.unnormalized[["finalSlope"]]
  finalIntercept <- fit.unnormalized[["finalIntercept"]]
  trueIntercept <- fit.unnormalized[["trueIntercept"]]
  y.lowerBound <- fit.unnormalized[["y.lowerBound"]]
  y.upperBound <- fit.unnormalized[["y.upperBound"]]
  x.lowerBound <- fit.unnormalized[["x.lowerBound"]]
  x.upperBound <- fit.unnormalized[["x.upperBound"]]

  # assemble results
  all.checks.passed <- dataQualityMetrics[["passed.check"]] && fitQualityMetrics[["passed.check"]]

  # print messages
  if (verbose) {
    message("  Un-normalized fit")
    message(paste0("      Final Slope:             ", finalSlope, " ", slope.unit))
    message(paste0("      True Intercept:          ", trueIntercept, " ", y.unit))
    if (length(y.unit) > 0) {
      message(paste0(
        "      y-Range:                 ", y.lowerBound, " ", y.unit,
        " - ", y.upperBound, " ", y.unit, "\n"
      ))
    } else {
      message(paste0(
        "      y-Range:                 ", y.lowerBound,
        " - ", y.upperBound, "\n"
      ))
    }
  }

  # give a warning when quality checks have failed
  if (warn) {
    if (!all.checks.passed) {
      warning(
        "Quality checks failed: result might be inaccurate!\n",
        call. = FALSE
      )
    }
  }

  # return results
  list(
    "finalSlope" = finalSlope,
    "finalSlope.unit" = slope.unit,
    "trueIntercept" = trueIntercept,
    "trueIntercept.unit" = y.unit,
    "y.lowerBound" = y.lowerBound,
    "y.upperBound" = y.upperBound,
    "y.unit" = y.unit,
    "x.lowerBound" = x.lowerBound,
    "x.upperBound" = x.upperBound,
    "x.unit" = x.unit,
    "passed.check" = all.checks.passed
  ) %>%
    return()
}
