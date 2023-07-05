#' Run the SDAR algorithm on prepared data
#' @noRd
sdarr_execute <- function(prepared_data, otr.info,
                          verbose.all = FALSE,
                          verbose.report = TRUE,
                          showPlots.all = FALSE,
                          showPlots.report = TRUE,
                          savePlots = FALSE) {
  # maybe the offset for step 1 needs to be raised later
  raise_offset_times <- 0
  optimum_fit_is_found <- FALSE

  # Give a welcome message
  if (verbose.report) {
    message("SDAR-algorithm\n")
  }

  # do step 1, check data quality, and get optimum fit
  # repeat until the upper index is not the last point in the optimum region
  while (!optimum_fit_is_found) {
    # do step 1:
    normalized_data <- normalize_data(
      data = prepared_data,
      otr.info = otr.info,
      raise_offset_times = raise_offset_times,
      denoise.x = FALSE,
      denoise.y = FALSE,
      verbose = verbose.all,
      showPlots = showPlots.all,
      savePlots = savePlots
    )

    # check data quality
    data_quality_metrics <- check_data_quality(
      normalized_data$data.normalized,
      verbose.all,
      showPlots.all,
      savePlots
    )

    # find best fit
    optimum_fit <- find_linear_regressions(
      normalized_data$data.normalized,
      verbose.all
    )

    # check if the offset needs to be raised
    optimum_fit_is_found <- !optimum_fit$offset_raise_required

    if (!optimum_fit_is_found) {
      raise_offset_times <- raise_offset_times + 1
    }

    # next offset is out of data range
    if (raise_offset_times > 17) {
      stop("Data is unfit for processing: upper index of the optimum fit region is beyond the last point in the truncated test record.")
    }
  }

  # check fit quality
  fit_quality_metrics <- check_fit_quality(
    normalized_data$data.normalized,
    optimum_fit,
    verbose.all,
    showPlots.all,
    savePlots
  )

  # un-normalize data and summarize
  assembled_results <- assemble_report(
    normalized_data,
    otr.info,
    data_quality_metrics,
    optimum_fit,
    fit_quality_metrics,
    verbose.report,
    showPlots.report,
    savePlots
  )

  # re-do the final fit to find confidence intervals for slope and intercepts
  final.fit.confint <- prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end, ] %>%
    {
      prepared_data.fit.lm <- stats::lm(y.data ~ x.data, data = .)
    } %>%
    stats::confint()

  results.labels <- list(
    "finalSlope" = assembled_results$Slope_Determination_Results$finalSlope.unit,
    "finalSlope.conf.low" = assembled_results$Slope_Determination_Results$finalSlope.unit,
    "finalSlope.conf.high" = assembled_results$Slope_Determination_Results$finalSlope.unit,
    "trueIntercept" = assembled_results$Slope_Determination_Results$trueIntercept.unit,
    "trueIntercept.conf.low" = assembled_results$Slope_Determination_Results$trueIntercept.unit,
    "trueIntercept.conf.high" = assembled_results$Slope_Determination_Results$trueIntercept.unit,
    "y.lowerBound" = assembled_results$Slope_Determination_Results$y.unit,
    "y.upperBound" = assembled_results$Slope_Determination_Results$y.unit,
    "x.lowerBound" = assembled_results$Slope_Determination_Results$x.unit,
    "x.upperBound" = assembled_results$Slope_Determination_Results$x.unit,
    "Num.Obs.fit" = "(points in final fit)",
    "Num.Obs.normalized" = "(points in normalized data)",
    "passed.check" = "(all checks)",
    "passed.check.data" = "(data quality checks)",
    "passed.check.fit" = "(fit quality checks)"
  )

  sdarr.results <- data.frame(
    "finalSlope" = assembled_results$Slope_Determination_Results$finalSlope,
    "finalSlope.conf.low" = final.fit.confint["x.data", 1],
    "finalSlope.conf.high" = final.fit.confint["x.data", 2],
    "trueIntercept" = assembled_results$Slope_Determination_Results$trueIntercept,
    "trueIntercept.conf.low" = final.fit.confint["(Intercept)", 1],
    "trueIntercept.conf.high" = final.fit.confint["(Intercept)", 2],
    "y.lowerBound" = assembled_results$Slope_Determination_Results$y.lowerBound,
    "y.upperBound" = assembled_results$Slope_Determination_Results$y.upperBound,
    "x.lowerBound" = assembled_results$Slope_Determination_Results$x.lowerBound,
    "x.upperBound" = assembled_results$Slope_Determination_Results$x.upperBound,
    "Num.Obs.fit" = nrow(prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end, ]),
    "Num.Obs.normalized" = assembled_results$Data_Quality_Metrics$Num.Obs.normalized,
    "passed.check" = assembled_results$Slope_Determination_Results$passed.check,
    "passed.check.data" = assembled_results$Data_Quality_Metrics$passed.check,
    "passed.check.fit" = assembled_results$Fit_Quality_Metrics$passed.check
  ) %>%
    labelled::set_variable_labels(.labels = results.labels)

  # return full results
  results <- list(
    "sdar" = sdarr.results,
    "Data_Quality_Metrics" = assembled_results$Data_Quality_Metrics,
    "Fit_Quality_Metrics" = assembled_results$Fit_Quality_Metrics
  )

  # append/sort plots to the results
  if (savePlots) {
    results <- results %>% append(list("plots" = list(
      "final.fit" = assembled_results$plots$plot.fit,
      "otr" = assembled_results$plots$plot.otr,
      "normalized" = assembled_results$plots$plot.normalized,
      "hist.x" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.x,
      "hist.y" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.y,
      "residuals" = assembled_results$Fit_Quality_Metrics$plots$plot.residuals
    )))

    results$Data_Quality_Metrics$digitalResolution$plots <- NULL
    results$Fit_Quality_Metrics$plots <- NULL
  }

  # return results
  return(results)
}


#' @title SDAR-algorithm
#'
#' @description Run the SDAR algorithm as standardized in ASTM E3076-18. Will
#'   use numerous linear regressions (`.lm.fit()` from the stats-package) and
#'   can be painfully slow for test data with high resolution.
#'
#' @note The function can use parallel processing via the
#'   [`furrr-package`][furrr::furrr-package]. To use this feature, set up a plan
#'   other than the default sequential strategy beforehand.
#'
#' @references Lucon, E. (2019). Use and validation of the slope determination
#'   by the analysis of residuals (SDAR) algorithm (NIST TN 2050; p. NIST TN
#'   2050). National Institute of Standards and Technology.
#'   https://doi.org/10.6028/NIST.TN.2050
#'
#' @references Standard Practice for Determination of the Slope in the Linear
#'   Region of a Test Record (ASTM E3076-18). (2018).
#'   https://doi.org/10.1520/E3076-18
#'
#' @references Graham, S., & Adler, M. (2011). Determining the Slope and Quality
#'   of Fit for the Linear Part of a Test Record. Journal of Testing and
#'   Evaluation - J TEST EVAL, 39. https://doi.org/10.1520/JTE103038
#'
#' @param data Data record to analyze. Labels of the data columns will be used
#'   as units.
#'
#' @param x,y Use <[`tidy-selections`][dplyr::dplyr_tidy_select]> to specify x
#'   and y within the data record.
#'
#' @param verbose,showPlots Give informational messages and plots during
#'   computation. Defaults to `"report"` to only show a report and a plot of the
#'   final fit. Set to `"all"` to also give messages from the individual steps.
#'   Set to `"none"` to be quiet. Can be abbreviated.
#'
#' @param savePlots Set to `TRUE` to get plot-functions with the result for
#'   later use.
#'
#' @seealso [sdarr.lazy()] for the random sub-sampling modification of the
#'   SDAR-algorithm.
#'
#' @examples
#' # Synthesize a test record resembling Al 6060 T66
#' # (Values according to Metallic Material Properties
#' # Development and Standardization (MMPDS) Handbook).
#' # Explicitly set names to "strain" and "stress",
#' # set effective number of bits in the x-data to 12
#' # to limit the number of data points.
#'
#' Al_6060_T66 <- synthesize_test_data(
#'   slope = 68000,
#'   yield.y = 160,
#'   ultimate.y = 215,
#'   ultimate.x = 0.091,
#'   x.name = "strain",
#'   y.name = "stress",
#'   enob.x = 12
#' )
#'
#'
#' # use sdarr() to analyze the synthetic test record
#' # will print a report and give a plot of the final fit
#' \donttest{
#' result <- sdarr(Al_6060_T66, strain, stress)
#' }
#'
#' @returns A list containing a data-frame with the results of the final fit,
#'   lists with the quality- and fit-metrics, and a list containing the
#'   <[`crated`][carrier::crate]> plot-functions (if `savePlots = TRUE`).
#'
#' @export
sdarr <- function(data,
                  x,
                  y,
                  verbose = "report",
                  showPlots = "report",
                  savePlots = FALSE) {
  # to be furrr-safe, enquote the tidy arguments here
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)

  # determine verbosity level
  verbose.all <- ifelse(verbose == "a" || verbose == "all", TRUE, FALSE)
  verbose.report <- ifelse(verbose == "r" || verbose == "report", TRUE, FALSE)
  verbose.report <- verbose.report || verbose.all

  # determine plot level
  showPlots.all <- ifelse(showPlots == "a" || showPlots == "all", TRUE, FALSE)
  showPlots.report <- ifelse(showPlots == "r" || showPlots == "report", TRUE, FALSE)
  showPlots.report <- showPlots.report || showPlots.all

  # get units for data
  x.label.unit <- data %>%
    dplyr::select(!!x) %>%
    labelled::var_label() %>%
    {
      .[[1]]
    }
  y.label.unit <- data %>%
    dplyr::select(!!y) %>%
    labelled::var_label() %>%
    {
      .[[1]]
    }

  # get names of data
  x.name <- data %>%
    dplyr::select(!!x) %>%
    names() %>%
    {
      .[[1]]
    }
  y.name <- data %>%
    dplyr::select(!!y) %>%
    names() %>%
    {
      .[[1]]
    }

  # prepare data, add an index and remove NA from data
  prepared_data <- data.frame(
    x.data = data %>% dplyr::select(!!x) %>% unlist(TRUE, FALSE),
    y.data = data %>% dplyr::select(!!y) %>% unlist(TRUE, FALSE),
    otr.idx = seq.int(nrow(data))
  ) %>%
    tidyr::drop_na()

  # assemble information of the original test record
  otr.info <- list(
    "x" = list(
      "name" = x.name,
      "unit" = x.label.unit
    ),
    "y" = list(
      "name" = y.name,
      "unit" = y.label.unit
    )
  )

  # Give a welcome message
  if (verbose.report) {
    message("Determination of Slope in the Linear Region of a Test Record:")
  }

  # execute the SDAR-algorithm
  sdarr_execute(
    prepared_data,
    otr.info,
    verbose.all,
    verbose.report,
    showPlots.all,
    showPlots.report,
    savePlots
  ) %>%
    return()
}
