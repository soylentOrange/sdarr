#' Run the SDAR algorithm on prepared data
#' @noRd
sdar_execute <- function(prepared_data,
                         otr.info,
                         verbose.all = FALSE,
                         verbose = TRUE,
                         plot.all = FALSE,
                         plot = TRUE,
                         plotFun.all = FALSE,
                         plotFun = FALSE,
                         Nmin_factor = 0.2) {
  # maybe the offset for step 1 needs to be raised later
  raise_offset_times <- 0
  optimum_fit_is_found <- FALSE

  # Give a welcome message
  if (verbose) {
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
      plot = plot.all,
      plotFun = plotFun.all
    )

    # check data quality
    data_quality_metrics <- check_data_quality(
      data.normalized = normalized_data$data.normalized,
      verbose = verbose.all,
      plot = plot.all,
      plotFun = plotFun.all
    )

    # find best fit
    optimum_fit <- find_linear_regressions(
      normalized_data = normalized_data,
      verbose = verbose.all,
      Nmin_factor = Nmin_factor
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
    data.normalized = normalized_data$data.normalized,
    fit = optimum_fit,
    verbose = verbose.all,
    plot = plot.all,
    plotFun = plotFun.all
  )

  # un-normalize data and summarize
  assembled_results <- assemble_report(
    normalized_data = normalized_data,
    otr.info = otr.info,
    dataQualityMetrics = data_quality_metrics,
    fit = optimum_fit,
    fitQualityMetrics = fit_quality_metrics,
    verbose = verbose,
    plot = plot,
    plotFun.all = plotFun.all,
    plotFun = plotFun
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

  sdar.results <- data.frame(
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
    "sdar" = sdar.results,
    "Data_Quality_Metrics" = assembled_results$Data_Quality_Metrics,
    "Fit_Quality_Metrics" = assembled_results$Fit_Quality_Metrics
  )

  # append plot
  if (plotFun || plotFun.all) {
    results <- results %>% append(list("plots" = list(
      "final.fit" = assembled_results$plots$plot.fit
    )))
  }

  # append more plots
  if (plotFun.all) {
    results$plots <- results$plots %>% append(list(
      "otr" = assembled_results$plots$plot.otr,
      "normalized" = assembled_results$plots$plot.normalized,
      "hist.x" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.x,
      "hist.y" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.y,
      "residuals" = assembled_results$Fit_Quality_Metrics$plots$plot.residuals
    ))

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
#'   can be painfully slow for test data with high resolution. See the article
#'   on [validation](
#'   https://soylentorange.github.io/sdarr/articles/sdarr_validation.html) on
#'   the [package-website](https://soylentorange.github.io/sdarr/) for further
#'   information.
#'
#' @note The function can use parallel processing via the
#'   [furrr-package](https://furrr.futureverse.org/). To use this feature, set
#'   up a plan other than the default sequential strategy beforehand.
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
#' @param x,y <[`tidy-select`][dplyr::dplyr_tidy_select]> Columns with x and y
#'   within data.
#'
#' @param verbose,plot Give a summarizing report / show a plot of the final fit.
#'
#' @param plotFun Set to `TRUE` to get a plot-function for the final fit with
#'   the results for later use.
#'
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> Pass parameters to downstream
#'  functions: set `verbose.all`, `plot.all` and `plotFun.all` to `TRUE` to get
#'  additional diagnostic information during processing data.
#'
#' @seealso [sdar.lazy()] for the random sub-sampling modification of the
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
#'   toe.start.y = 3, toe.end.y = 10,
#'   toe.start.slope = 13600,
#'   enob.x = 12
#' )
#'
#'
#' # use sdar() to analyze the synthetic test record
#' # will print a report and give a plot of the final fit
#' \donttest{
#' result <- sdar(Al_6060_T66, strain, stress)
#' }
#'
#' @returns A list containing a data.frame with the results of the final fit,
#'   lists with the quality- and fit-metrics, and a list containing the crated
#'   plot-function(s) (if `plotFun = TRUE`).
#'
#' @export
sdar <- function(data, x, y,
                 verbose = TRUE,
                 plot = TRUE,
                 plotFun = FALSE,
                 ...) {
  # to be furrr-safe, enquote the tidy arguments here
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)

  # take care of dynamic dots
  additional_parameters <- rlang::list2(...)
  Nmin_factor <- try(
    {
      additional_parameters$Nmin_factor
    },
    silent = TRUE
  )
  if (!is.numeric(Nmin_factor)) Nmin_factor <- 0.2
  verbose.all <- try(
    {
      additional_parameters$verbose.all
    },
    silent = TRUE
  )
  if (!is.logical(verbose.all)) verbose.all <- FALSE
  plot.all <- try(
    {
      additional_parameters$plot.all
    },
    silent = TRUE
  )
  if (!is.logical(plot.all)) plot.all <- FALSE
  plotFun.all <- try(
    {
      additional_parameters$plotFun.all
    },
    silent = TRUE
  )
  if (!is.logical(plotFun.all)) plotFun.all <- FALSE

  # save final fit plot, when plotFun.all is set
  plotFun <- plotFun || plotFun.all

  # give messages for report, when verbose.all is set
  verbose <- verbose || verbose.all

  # show final fit plot, when plot.all is set
  plot <- plot || plot.all

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
    y.data = data %>% dplyr::select(!!y) %>% unlist(TRUE, FALSE)
  ) %>%
    dplyr::mutate("otr.idx" = as.numeric(rownames(.))) %>%
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
  if (verbose) {
    message("Determination of Slope in the Linear Region of a Test Record:")
  }

  # execute the SDAR-algorithm
  result <- sdar_execute(
    prepared_data = prepared_data,
    otr.info = otr.info,
    verbose.all = verbose.all,
    verbose = verbose,
    plot.all = plot.all,
    plot = plot,
    plotFun.all = plotFun.all,
    plotFun = plotFun,
    Nmin_factor = Nmin_factor
  ) %>%
    return()
}
