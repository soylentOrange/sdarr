#' Run the SDAR algorithm on prepared data
#'
#' @noRd
sdarr_execute <- function(prepared_data, otr.info,
                          verbose.all = F,
                          verbose.report = T,
                          showPlots.all = F,
                          showPlots.report = T,
                          savePlots = F) {

  # maybe the offset for step 1 needs to be raised later
  raise_offset_times <- 0
  optimum_fit_is_found <- FALSE

  # Give a welcome message
  if(verbose.report) {
    message("SDAR-algorithm\n")
  }

  # do step 1, check data quality, and get optimum fit
  # repeat until the upper index is not the last point in the optimum region
  while(!optimum_fit_is_found) {
    # do step 1:
    normalized_data <- normalize_data(prepared_data,
                                      otr.info,
                                      raise_offset_times,
                                      verbose.all,
                                      showPlots.all,
                                      savePlots)

    # check data quality
    data_quality_metrics <- check_data_quality(normalized_data$data.normalized,
                                               verbose.all,
                                               showPlots.all,
                                               savePlots)

    # find best fit
    optimum_fit <- find_linear_regressions(normalized_data$data.normalized,
                                           verbose.all)

    # check if the offset needs to be raised
    optimum_fit_is_found <- !optimum_fit$offset_raise_required

    if(!optimum_fit_is_found) {
      raise_offset_times <- raise_offset_times + 1
    }

    # next offset is out of data range
    if(raise_offset_times > 17) {
      stop("Data is unfit for processing: upper index of the optimum fit region is beyond the last point in the truncated test record.")
    }
  }

  # check fit quality
  fit_quality_metrics <- check_fit_quality(normalized_data$data.normalized,
                                           optimum_fit,
                                           verbose.all,
                                           showPlots.all,
                                           savePlots)

  # un-normalize data and summarize
  assembled_results <- assemble_report(normalized_data,
                                       otr.info,
                                       data_quality_metrics,
                                       optimum_fit,
                                       fit_quality_metrics,
                                       verbose.report,
                                       showPlots.report,
                                       savePlots)

  # re-do the final fit to find confidence intervals for slope and intercepts
  final.fit.confint <- prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end,] %>% {
    prepared_data.fit.lm <- stats::lm(y.data ~ x.data, data = .)
  } %>% stats::confint()

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
    "Num.Obs.fit" = nrow(prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end,]),
    "Num.Obs.normalized" = assembled_results$Data_Quality_Metrics$Num.Obs.normalized,
    "passed.check" = assembled_results$Slope_Determination_Results$passed.check,
    "passed.check.data" = assembled_results$Data_Quality_Metrics$passed.check,
    "passed.check.fit" = assembled_results$Fit_Quality_Metrics$passed.check
  ) %>%
    labelled::set_variable_labels(.labels = results.labels)

  # return full results
  results <- list("sdar" = sdarr.results,
                  "Data_Quality_Metrics" = assembled_results$Data_Quality_Metrics,
                  "Fit_Quality_Metrics" = assembled_results$Fit_Quality_Metrics)

  # append/sort plots to the results
  if(savePlots) {
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
