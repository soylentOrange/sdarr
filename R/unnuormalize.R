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
                                verbose = F,
                                warn = T) {

  # shorthands for calculations
  y.tangent <- normalized_data[["tangency.point"]][["y.tangent"]]
  x.tangent <- normalized_data[["tangency.point"]][["x.tangent"]]
  x.shift <- normalized_data[["shift"]][["x.shift"]]
  y.shift <- normalized_data[["shift"]][["y.shift"]]
  y.min <- normalized_data[["data.normalized"]][fit[["otr.idx.start"]], "y.normalized"]
  y.max <- normalized_data[["data.normalized"]][fit[["otr.idx.end"]], "y.normalized"]
  x.min <- normalized_data[["data.normalized"]][fit[["otr.idx.start"]], "x.normalized"]
  x.max <- normalized_data[["data.normalized"]][fit[["otr.idx.end"]], "x.normalized"]
  m <- fit[["m"]]
  b <- fit[["b"]]

  # get units if available
  if(!is.null(otr.info)) {
    x.unit <- otr.info$x$unit
    y.unit <- otr.info$y$unit
    if(!is.null(x.unit) && !is.null(y.unit)) {
      slope.unit <- paste0(y.unit , "/", x.unit)
    } else if(is.null(x.unit) && !is.null(y.unit)) {
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
  finalSlope <- m * y.tangent / x.tangent
  finalIntercept <- b * y.tangent
  trueIntercept <- finalIntercept - x.shift * finalSlope + y.shift
  y.lowerBound <- y.min * y.tangent + y.shift
  y.upperBound <- y.max * y.tangent + y.shift
  x.lowerBound <- x.min * x.tangent + x.shift
  x.upperBound <- x.max * x.tangent + x.shift

  # assemble results
  all.checks.passed <- dataQualityMetrics[["passed.check"]] && fitQualityMetrics[["passed.check"]]

  # print messages
  if(verbose) {
    message("  Un-normalized fit")
    message(paste0("      Final Slope:             ", finalSlope, " ", slope.unit))
    message(paste0("      True Intercept:          ", trueIntercept, " ", y.unit))
    if(length(y.unit) > 0) {
      message(paste0("      y-Range:                 ", y.lowerBound, " ", y.unit,
                     " - ", y.upperBound, " ", y.unit, "\n"))
    } else {
      message(paste0("      y-Range:                 ", y.lowerBound,
                     " - ", y.upperBound, "\n"))
    }
  }

  # give a warning when quality checks have failed
  if(warn) {
    if(!all.checks.passed) {
      warning("Some quality checks failed: result might be inaccurate or completely wrong!\n")
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
    "passed.check" = all.checks.passed %>%
      return()
  )
}
